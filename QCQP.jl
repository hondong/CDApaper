include("CDAutil.jl")
import CDA
using MAT
#using JLD

#matname = "/Users/hdong/workspace/Manuscripts/cda/data/qcqp/unitbox_c_10_10_1_50.mat";
#matname = "/Users/hdong/workspace/Manuscripts/cda/data/qcqp/unitbox_c_50_100_1_50.mat";
#matname = "/Users/hdong/workspace/Manuscripts/cda/data/qcqp/unitbox_c_40_40_3_50.mat";
#matname = "/Users/hdong/workspace/Manuscripts/cda/data/qcqp/unitbox_c_30_45_2_50.mat";
#matname = "/Users/hdong/workspace/Manuscripts/cda/data/qcqp/unitbox_c_30_60_3_50.mat";
matname = "/Users/hdong/workspace/Manuscripts/cda/data/qcqp/unitbox_c_40_60_1_100.mat";

function ReadQCQPMat(matname::String)
    MatStr = matread(matname);
    EPS = 1E-14;
    m   = Int(MatStr["m"]);
    n   = Int(MatStr["n"])
    P = CDA.MIQCP(Int(MatStr["n"]))
    P.Qobj = MatStr["Qobj"];

    P.cobj = squeeze(MatStr["Cobj"],2);
    P.lb   = squeeze(MatStr["lb"],  2);
    P.ub   = squeeze(MatStr["ub"],  2);
    P.Ql   = Array{AbstractMatrix{Float64},1}(m);
    P.A    = zeros(m, n);
    P.LHS  =Array{Float64,1}(m)
    P.RHS  =Array{Float64,1}(m)
    if (m>0)
        for k=1:m
            P.Ql[k] = sparse(MatStr["Ql"][:, :, k]);
            P.A[k,:] = (MatStr["Cl"][:,k])';
        end
        P.A     = sparse(P.A);
        P.LHS   = squeeze(MatStr["cl"], 2);
        P.RHS   = squeeze(MatStr["cu"], 2);
    end
    P.VarType = Array{CDA.VARTYPE,1}(n);
    for i=1:n
        if (MatStr["xtype"][i]=='C')
            P.VarType[i] = CDA.CONTINUOUS;
        elseif (MatStr["xtype"][i]=='B')
            P.VarType[i] = CDA.BINARY;
        elseif (MatStr["xtype"][i]=='I')
            P.VarType[i] = CDA.INTEGER;
        else
            error("Something wrong with the vartype?");
        end
    end

    # Clearing small entries
    P.Qobj[abs.(P.Qobj).< EPS] = 0.0;
    P.cobj[abs.(P.cobj).< EPS] = 0.0;
    P.lb[abs.(P.lb).< EPS]     = 0.0;
    P.ub[abs.(P.ub).< EPS]     = 0.0;
    if (m>0)
        P.A[abs.(P.A).< EPS]       = 0.0;
        P.LHS[abs.(P.LHS).< EPS]   = 0.0;
        P.RHS[abs.(P.RHS).< EPS]   = 0.0;
        for k=1:m
            (P.Ql[k])[abs.(P.Ql[k]).< EPS]   = 0.0;
        end
    end
    CDA.setSetIndices!(P)
    return P ;
end

# ###
# ########## Code block: Run QCQP with CDA ######################
function RandQCQPTest(FolderList, summary_file)
    fsummary = open(summary_file, "w");
    write(fsummary, "InstName             NumVars   NumCons     BestObj       BestBound      NumIter     Time \n")
    ops = CDA.CDAoption();
    ops.TimeLimit = 1000;
    ops.display = 1;
    ops.GapTol = 1E-4;
    ops.NumRefine = 20;
    AllRet1 = Dict{String, Any}()
    for folder in FolderList
        for matname in readdir(folder)
            (instname,extname) = splitext(matname);
            if (extname==".mat")
                try
                    P = ReadQCQPMat(joinpath(folder, matname))
                    open(joinpath(folder, string(instname,".cdaout")), "w") do outf
                        redirect_stdout(outf) do
                            println("==============Instane Name: ", instname, "==================")
                            CDA.displayMIQCP(P)
                            println("==============End of Instane ", instname, "=================\n")
                            RET = CDA_MIQCP_Driver(P, "Gurobi", ops);
                            push!(AllRet1, instname=>RET)
                            @printf(fsummary, "%20s", instname);
                            @printf(fsummary, "%6d ", P.n);
                            @printf(fsummary, "%6d ", length(P.Ql));
                            @printf(fsummary, "%16.5f", (RET["BestObjVal"]));
                            @printf(fsummary, "%16.5f", (RET["BestBound"]));
                            @printf(fsummary, "%6d", (count(RET["Iters"][i].solvername=="Gurobi" for i in 1:length(RET["Iters"]))) );
                            @printf(fsummary, "%12.2f\n", (RET["TotalTime"]));
                        end
                    end
                catch
                    println("Error processing instance ", instname, ".")
                    write(fsummary, "Error processing instance \n")
                    continue;
                end
                flush(fsummary)
            end
        end
    end

    close(fsummary)
end

# stubname = "newstub2.nl";
# CDA.GenNL(P,stubname)
# run(`couenne $stubname`)

# ops = CDA.CDAoption();
# ops.TimeLimit = 120;
# #ops.initialNu = 1;
# ops.display = 0;
# ops.GapTol = 1E-4;
# ops.NumRefine = 20;
# CDA_MIQCP_Driver(P, "Gurobi", ops);


# ops.TimeLimit = 1000;
# ops.UseCQ = 1;
# ops.NumRefine = 20;
# ops.initialNu = 2;
# ops.display = 0;
# CDA_MIQCP_Driver(P, "Gurobi", ops)

#(obj, dx) = Couenne_Driver(P);


"Code block: Generating nl files from mat files for Couenne runs"
function QCQPMat2NL(FolderList)
    for folder in FolderList
         for matname in readdir(folder)
             (instname,extname) = splitext(matname);
             if (extname==".mat")
                try
                     P = ReadQCQPMat(joinpath(folder, matname))
                     println("==============Instane Name: ", instname, "==================")
                     CDA.displayMIQCP(P)
                     println("==============End of Instane ", instname, "=================\n")
                     stubname = string(instname, ".nl");
                     CDA.GenNL(P,joinpath(folder, stubname))
                catch
                     println("Warning: instance ", instname, "readmat failed.")
                     continue;
                end
            end
        end
    end
end
