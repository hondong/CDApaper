include("CDAutil.jl")
using CPLEX

"Parse BoxQP input files, return (n,c,Q)
optimization problem is:
    max 0.5*x^T Q x + c'x       s.t., 0<=x<=1
 "
function ParseBoxQP(f)
    n = parse(Int, readline(f))
    println("n=", n)
    c = parse.(split(readline(f)))
    c = Array{Float64,1}(c);
    if (length(c) != n)
        println("Warning: length(c) error. c=",c)
    end
    Q = zeros(n,n);
    for i in 1:n
        Q[:,i] = parse.(split(readline(f)))
    end
    return (n,c,Q)
end

"Solve BoxQP instances with CPLEX, return an object of type CDA.SolverInfo.
"
function BoxQP_CplexDriver(n,c,Q)
    m = Model(solver=CplexSolver(CPX_PARAM_OPTIMALITYTARGET=3, CPX_PARAM_TILIM=1200))
    #m = Model(solver=GurobiSolver())
    lb = 0;
    ub = 1;
    @variable(m, 0<=x[1:n]<=1);
    @objective(m, Min, -0.5*x'*Q*x - c'*x)
    #varlist = SetUpApprox(m, x, y, lb, ub, 4);
    #println(m)
    println("Start solving full model: ")
    status = solve(m)
    println("Solution status: ", status)
    println("Objective value: ", getobjectivevalue(m))
    ret = CDA.SolverInfo("CPLEX", string(status), getsolvetime(m), getnodecount(m),
                            getobjgap(m), getobjbound(m), getobjectivevalue(m), -1);
    return ret;
end

function CheckBoxQPViolation(dx::AbstractVector{Float64},dy::AbstractVector{Float64}, S, aug::AbstractVector{Float64})
    viol = [];
    allviol = 0.0;
    for i in S
        #if (inst.VarType[i]==CDA.CONTINUOUS)
            push!( viol,(i,aug[i]*(dy[i]-dx[i]^2)) )
            allviol = allviol + aug[i]*(dy[i]-dx[i]^2);
        #end
    end
    sort!(viol, by=(x->x[2]), rev=true)
    return (allviol, viol)
end

"Set up MIQCP structure for BoxQP instances"
function SetUpBoxQP_CDA(n::Int, c::AbstractVector{Float64}, Q::AbstractMatrix{Float64})
    inst  = CDA.MIQCP(n) # last one is the dummy variable representing obj
    CDA.setVarType!(inst, fill(CDA.CONTINUOUS,n))
    CDA.setBounds!(inst, zeros(n), ones(n))
    CDA.setObj!(inst, -c, -0.5*Q)
    CDA.setSetIndices!(inst)
    return inst;
end

###
function BoxQPTest(FolderList, summary_file)
    fsummary = open(summary_file, "w");
    write(fsummary, "==================================== CPLEX Results ======================================= \n")
    @printf(fsummary, "%15s %7s %10s %10s %9s %9s\n", "InstName", "Time",
        "NodeCount", "Gap", "Bound", "objval")
    for folder in FolderList
        for fname in readdir(folder)
            (instname,extname) = splitext(fname);
            if (extname==".in")
                try
                    open(joinpath(folder,fname)) do f
                        (n,c,Q) = ParseBoxQP(f)
                        info = BoxQP_CplexDriver(n,c,Q)
                        @printf(fsummary, "%15s %7.3f %10d %10.5f %9.3e %9.3e\n", instname, info.time,
                        info.nodecount, info.gap, info.bound, info.objval)
                        flush(fsummary)
                    end
                catch
                    @printf(fsummary, "%15s Error processing instance \n", instname)
                end
            end
        end
    end
    write(fsummary, "================================= End of CPLEX Results ===================================== \n")

    write(fsummary, "==================================== CDA Results ======================================= \n")
    @printf(fsummary, "%15s %8s %8s %12s %12s %6s %8s\n", "InstName", "NumVars", "NumCons", "BestObj", "BestBound", "NumIter", "Time");

    ops = CDA.CDAoption();
    ops.TimeLimit = 1200;
    ops.display = 1;
    ops.GapTol = 1E-4;
    ops.NumRefine = 20;
    AllRet1 = Dict{String, Any}()
    for folder in FolderList
        for fname in readdir(folder)
            (instname,extname) = splitext(fname);
            if (extname==".in")
                try
                    f = open(joinpath(folder, fname));
                    (n,c,Q) = ParseBoxQP(f)
                    close(f);
                    P = SetUpBoxQP_CDA(n,c,Q);

                    open(joinpath(folder, string(instname,".cdaout")), "w") do outf
                        redirect_stdout(outf) do
                            println("==============Instane Name: ", instname, "==================")
                            CDA.displayMIQCP(P)
                            println("==============End of Instane ", instname, "=================\n")
                            RET = CDA_MIQCP_Driver(P, "Gurobi", ops);
                            push!(AllRet1, instname=>RET)
                            @printf(fsummary, "%15s", instname);
                            @printf(fsummary, "%8d ", P.n);
                            @printf(fsummary, "%8d ", length(P.Ql));
                            @printf(fsummary, "%12.3f", (RET["BestObjVal"]));
                            @printf(fsummary, "%12.3f", (RET["BestBound"]));
                            @printf(fsummary, "%6d", (count(RET["Iters"][i].solvername=="Gurobi" for i in 1:length(RET["Iters"]))) );
                            @printf(fsummary, "%8.2f\n", (RET["TotalTime"]));
                        end
                    end
                catch
                    println("Error processing instance ", instname, ".")
                    write(fsummary, "Error processing instance \n")
                end
                flush(fsummary)
            end
        end
    end

    @printf(fsummary, "================================ End of CDA Results ==================================== \n")

    close(fsummary)
end

###
########## Code block: Run BoxQP instances individually with CDA ######################
#fname = "/Users/hdong/workspace/Manuscripts/cda/data/boxqp/basic/spar020-100-1.in"
#fname = "/Users/hdong/workspace/Manuscripts/cda/data/boxqp/basic/spar060-020-3.in"
#fname = "/Users/hdong/workspace/Manuscripts/cda/data/boxqp/basic/spar050-050-1.in"
#fname = "/Users/hdong/workspace/Manuscripts/cda/data/boxqp/extended/spar070-050-1.in"
#fname = "/Users/hdong/workspace/Manuscripts/cda/data/boxqp/extended/spar070-075-1.in"
#fname = "/Users/hdong/workspace/Manuscripts/cda/data/boxqp/extended/spar100-025-3.in"
#fname = "/Users/hdong/workspace/Manuscripts/cda/data/boxqp/extended/spar100-075-3.in"
#fname = "/Users/hdong/workspace/Manuscripts/cda/data/boxqp/extended2/spar125-075-3.in"
# f = open(fname);
# (n,c,Q) = ParseBoxQP(f)
# close(f);
# inst = SetUpBoxQP_CDA(n,c,Q);
# ops = CDA.CDAoption();
# ops.TimeLimit = 1200;
# ops.initialNu = 1;
# ops.display = 0;
# ops.GapTol = 1E-4;
# CDA_MIQCP_Driver(inst, "Gurobi", ops);
#CDA_MIQCP_Driver(inst, "Gurobi", ops);
#BoxQP_CDADriver(n, Array{Float64,1}(c), Q, 300.0, "Gurobi")
#BoxQP_CDA_FixedPrecision(n, Array{Float64,1}(c), Q, 300.0, "Gurobi",6)

#inst = SetUpBoxQP_CDA(n,c,Q);
#CDA_MIQCP_Driver(inst, "Gurobi", 4, 120.0)
########## End Code block: Run BoxQP instances with CDA ######################
