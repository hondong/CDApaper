#"This file contains utility functions"
using JuMP
using Gurobi
using Ipopt
#using CPLEX
#using Mosek

include("CDA.jl")
import CDA

include("DiagPerturbSDP.jl")

" Set up the compact disjunctive approximation for pair (x,y) "
function SetUpApprox(m::JuMP.Model, x, y, l, u, nu::Int)
    #println("Setting up approximations. nu = ", nu)
    Tmin = (l==0)?(-pi/2):((l>0)? atan((l^2-1)/(2*l)) : (atan((l^2-1)/(2*l))-pi));
    Tmax = (u==0)?(-pi/2):((u>0)? atan((u^2-1)/2,u) : (atan((u^2-1)/2,u)-pi));
    Tmid  = (Tmin+Tmax)/2;
    Tdiff = Tmax - Tmin;
    #println("SetUpApprox: Tmin = ", Tmin, ", Tmax = ", Tmax);
    #println("SetUpApprox: Tmid = ", Tmid, ", Tdiff = ", Tdiff);
    C = max(l^2,u^2)/2+0.5;
    xi = @variable(m, [1:nu])
    eta = @variable(m, [1:nu])
    lam1 = @variable(m, [1:nu], lowerbound = 0)
    lam2 = @variable(m, [1:nu], lowerbound = 0)
    z    = @variable(m, [1:nu], category=:Bin)
    @constraint(m, xi[1] == cos(Tmid)*x + sin(Tmid)*(0.5*y-0.5))
    @constraint(m, (lam2[1]-lam1[1])*C == -sin(Tmid)*x + cos(Tmid)*(0.5*y-0.5))
    @constraint(m, eta[1] == (lam1[1]+lam2[1])*C)
    @constraint(m, lam1[1] <= 1-z[1])
    @constraint(m, lam2[1] <= z[1])
    for j in 1:(nu-1)
        tt = Tdiff/(2^(j+1))
        Cj = C*sin(tt)
        @constraint(m, xi[j+1] == cos(tt)*xi[j] + sin(tt)*eta[j])
        @constraint(m, (lam2[j+1]-lam1[j+1])*Cj == -sin(tt)*xi[j] + cos(tt)*eta[j])
        @constraint(m, eta[j+1] == (lam1[j+1]+lam2[j+1])*Cj)
        @constraint(m, lam1[j+1] <= 1-z[j+1])
        @constraint(m, lam2[j+1] <= z[j+1])
    end
    Ang = Tdiff/(2^nu);
    HalfAng = Ang/2;
    @constraint(m, xi[nu]*cos(HalfAng)+eta[nu]*sin(HalfAng) >= (0.5*y+0.5)*cos(HalfAng))
    @constraint(m, xi[nu]*cos(Ang)+eta[nu]*sin(Ang) <= (0.5*y+0.5))
    @constraint(m, eta[nu] <= 0.5*y+0.5)
    dalist = CDA.DAlist(nu, l, u, xi, eta, lam1, lam2, z);
    return dalist;
end

function RefineApprox!(m::JuMP.Model, dalist::CDA.DAlist, x, y)
    l = dalist.l
    u = dalist.u
    Tmin = (l==0)?(-pi/2):((l>0)? atan((l^2-1)/(2*l)) : (atan((l^2-1)/(2*l))-pi));
    Tmax = (u==0)?(-pi/2):((u>0)? atan((u^2-1)/2,u) : (atan((u^2-1)/2,u)-pi));
    Tmid  = (Tmin+Tmax)/2;
    Tdiff = Tmax - Tmin;
    C = max(l^2,u^2)/2+0.5;

    xi = @variable(m)
    eta = @variable(m)
    lam1 = @variable(m, lowerbound = 0)
    lam2 = @variable(m, lowerbound = 0)
    z    = @variable(m, category=:Bin)
    if (dalist.nu==0)
        dalist.nu = 1;
        @constraint(m, xi == cos(Tmid)*x + sin(Tmid)*(0.5*y-0.5))
        @constraint(m, (lam2-lam1)*C == -sin(Tmid)*x + cos(Tmid)*(0.5*y-0.5))
        @constraint(m, eta == (lam1+lam2)*C)
        @constraint(m, lam1 <= 1-z)
        @constraint(m, lam2 <= z)

        Ang = Tdiff/(2^1);
        HalfAng = Ang/2;
        #push!(dalist.DynCons, @constraint(m, xi*cos(HalfAng)+eta*sin(HalfAng) >= (0.5*y+0.5)*cos(HalfAng)))
        #push!(dalist.DynCons, @constraint(m, xi*cos(Ang)+eta*sin(Ang) <= (0.5*y+0.5)))
        #push!(dalist.DynCons, @constraint(m, eta <= 0.5*y+0.5))
        @constraint(m, xi*cos(HalfAng)+eta*sin(HalfAng) >= (0.5*y+0.5)*cos(HalfAng))
        @constraint(m, xi*cos(Ang)+eta*sin(Ang) <= (0.5*y+0.5))
        @constraint(m, eta <= 0.5*y+0.5)
    else
        #j = dalist.nu;
        dalist.nu = dalist.nu+1;
        tt = Tdiff/(2^(dalist.nu))
        Cj = C*sin(tt)
        @constraint(m, xi == cos(tt)*dalist.xi[end] + sin(tt)*dalist.eta[end])
        @constraint(m, (lam2-lam1)*Cj == -sin(tt)*dalist.xi[end] + cos(tt)*dalist.eta[end])
        @constraint(m, eta == (lam1+lam2)*Cj)
        @constraint(m, lam1 <= 1-z)
        @constraint(m, lam2 <= z)

        Ang = Tdiff/(2^(dalist.nu));
        HalfAng = Ang/2;
        #push!(dalist.DynCons, @constraint(m, xi*cos(HalfAng)+eta*sin(HalfAng) >= (0.5*y+0.5)*cos(HalfAng)))
        #push!(dalist.DynCons, @constraint(m, xi*cos(Ang)+eta*sin(Ang) <= (0.5*y+0.5)))
        #push!(dalist.DynCons, @constraint(m, eta <= 0.5*y+0.5))
        @constraint(m, xi*cos(HalfAng)+eta*sin(HalfAng) >= (0.5*y+0.5)*cos(HalfAng))
        @constraint(m, xi*cos(Ang)+eta*sin(Ang) <= (0.5*y+0.5))
        @constraint(m, eta <= 0.5*y+0.5)
    end

    push!(dalist.xi, xi)
    push!(dalist.eta, eta)
    push!(dalist.lam1, lam1)
    push!(dalist.lam2, lam2)
    push!(dalist.z, z)
end

function DiagPerturb(Q::AbstractMatrix, Si, isLeftBounded::Bool, isRightBounded::Bool)
    n      = size(Q, 1);
    #println("Diagonal Perturbation via SDP: size(Si)=", length(Si),", leftbounded=",isLeftBounded, ", rightbounded=", isRightBounded)
    (isLeftBounded||isRightBounded)?():(error("DiagPerturb: unbounded both sides."));
    pvec_L = []; pvec_R = [];
    lambda = eigvals(full(Q[Si, Si]));
    if (minimum(lambda) > 1E-5)
        Convexity = 1;
    elseif (maximum(lambda) < -1E-5)
        Convexity = -1;
    else
        Convexity = 0;
    end
    T = 0.0;
    if (isLeftBounded && Convexity!=-1)
        (d,Trun) = ComputeDiagPerturbByDSDP(Array{Float64,1}(ones(Si)), -Q[Si,Si]);
        pvec_L = zeros(n);
        pvec_L[Si] = d + (1E-5);
        T = T + Trun;
    end
    if (isRightBounded && Convexity!=1)
        (d,Trun) = ComputeDiagPerturbByDSDP(Array{Float64,1}(ones(Si)), Q[Si,Si]);
        pvec_R = zeros(n);
        pvec_R[Si] = d + (1E-5);
        T = T + Trun;
    end
    (isempty(pvec_L))?(sumL=0):(sumL = sum(pvec_L))
    (isempty(pvec_R))?(sumR=0):(sumR = sum(pvec_R))
    #@printf("sum of perturbation: (%10.3f, %10.3f) \n", sumL, sumR )
    return (Convexity, pvec_L, pvec_R, T);
end

# function CheckViolation(dx::AbstractVector{Float64},dy::AbstractVector{Float64},inst::CDA.MIQCP, S)
#     var_viol = [];
#     for i in S
#         if (inst.VarType[i]==CDA.CONTINUOUS)
#             push!(var_viol,(i,dy[i]-dx[i]^2))
#         end
#     end
#     sort!(var_viol, by=(x->abs(x[2])), rev=true)
#     #sort!(var_viol, by=(x->abs(x[2])), rev=true)
#     #println(var_viol)
#     return var_viol
# end

function CheckViolation(dx::AbstractVector{Float64},dy::AbstractVector{Float64},P::CDA.MIQCP)
    var_viol = Array{Tuple{Int, Float64}}(0);
    for i in intersect(P.S, P.S_cont)
        push!(var_viol, (i, dy[i]-dx[i]^2) )
    end
    sort!(var_viol, by=(x->abs(x[2])), rev=true)

    obj_viol = 0.0;
    if (!iszero(P.Qobj))
        for j=1:length(var_viol)
            i = var_viol[j][1];
            obj_viol = obj_viol + P.AugList[0].right[i] * var_viol[j][2];
        end
    end

    cons_viol = 0.0;
    for i in P.S_quad # actual quadratic constraints
        if (!isempty(P.AugList[i].left) && !iszero(P.AugList[i].left))
            for j=1:length(var_viol) # violated variables
                k = var_viol[j][1]; # variable index
                cons_viol = cons_viol + P.AugList[i].left[k] * var_viol[j][2];
            end
        end
        if (!isempty(P.AugList[i].right) && !iszero(P.AugList[i].right))
            for j=1:length(var_viol) # violated variables
                k = var_viol[j][1]; # variable index
                cons_viol = cons_viol + P.AugList[i].right[k] * var_viol[j][2];
            end
        end
    end
    #println(var_viol)
    return (var_viol, obj_viol, cons_viol)
end

" Bound Tightening for unbounded continuous variables by feasiblity"
function BoundTightByFeas(P::CDA.MIQCP, m::Model, x, y)
    S_cont = find(P.VarType[i]==CDA.CONTINUOUS for i=1:P.n);
    MaxPass = 3;
    println("Start Bound Tightening ... ")
    for pass = 1:MaxPass
        NumChange = 0;
        for i in S_cont
            old_lb = P.lb[i]; old_ub = P.ub[i];
            @objective(m, Min, x[i])
            status = solve(m, relaxation=true);
            #println(status)
            (status==:Optimal)?(min_xi = getobjectivevalue(m)):(min_xi=-Inf);
            (min_xi > P.lb[i]+1E-5)?
            (P.lb[i]=min_xi; setlowerbound(x[i], min_xi); NumChange = NumChange+1;):();

            @objective(m, Max, x[i])
            status = solve(m, relaxation=true);
            #println(status)
            (status==:Optimal)?(max_xi = getobjectivevalue(m)):(max_xi=Inf);
            (max_xi < P.ub[i]-1E-5)?
            (P.ub[i]=max_xi; setupperbound(x[i], max_xi); NumChange = NumChange+1;):();

            @printf("Pass: %4d, i=%4d, old = (%f, %f), new = (%f, %f), NumChange=%4d\n",
                pass, i, old_lb, old_ub, min_xi, max_xi, NumChange)
        end
        (NumChange>0)?(@printf("Pass %3d: %d bounds tightened.\n", pass, NumChange)):
            (@printf("Pass %3d: no bounds tightened.\n", pass); break;);
    end
end

function GetSolverInfo(m::Model, solvername::String, solverstatus, ExpectBnB::Bool=true)
    ret = CDA.SolverInfo()
    ret.solvername = solvername;
    ret.status     = string(solverstatus);
    try
        ret.time   = getsolvetime(m);
    catch
        if (solvername!="Ipopt") # Ipopt has no getsolvetime ...
            println("Warning: got error when getsolvetime.")
        end
        ret.time   = 0.0;
    end

    try
        ret.objval = getobjectivevalue(m);
    catch
        println("Warning: got error when getobjectivevalue.")
        ret.objval = +Inf;
    end

    if (ExpectBnB)
        try
            ret.nodecount = getnodecount(m);
        catch
            println("Warning: expect bound but got error when getnodecount.")
            ret.nodecount = 0;
        end
        try
            ret.bound = getobjectivebound(m);
        catch
            println("Warning: expect bound but got error when getobjectivebound.")
            ret.bound = -Inf;
        end
    else
        ret.nodecount = 0;
        ret.bound = ret.objval;
    end

    ret.gap        = (ret.bound - ret.objval)/((1E-6)+abs(ret.objval))
    return ret;
end

function TerminationCheck(BestObjVal::Float64, BestBound::Float64, ops::CDA.CDAoption)
    if abs(BestObjVal-BestBound)/(1E-6 + abs(BestObjVal) ) < ops.GapTol
        println("Desired accuracy achieved.")
        return true;
    else
        return false;
    end
end

function HeuristicRun(m_heu::Model, RET::Dict{String, Any}, xH, dx::AbstractVector{Float64}=Vector{Float64}(0))
    if (!isempty(dx))
        for i in 1:length(dx)
            setvalue(xH[i], dx[i])
        end
    end
    status = solve(m_heu; suppress_warnings=true)
    info = GetSolverInfo(m_heu, "Ipopt", status, false)
    if (info.objval<RET["BestObjVal"])
        RET["BestObjVal"] = info.objval;
        RET["BestX"] = getvalue(xH)
    end
    push!(RET["Iters"], info)
    return info
end

function RelaxationRun(P::CDA.MIQCP, m_rlx::Model, x, y, RET::Dict{String, Any}, env_rlx, TimeLeft::Float64, ops::CDA.CDAoption, ExpectBnB::Bool=true)
    setparam!(env_rlx, "TimeLimit", TimeLeft);
    setparam!(env_rlx, "Cutoff", RET["BestObjVal"]+(2E-2)*abs(RET["BestObjVal"]) );
    setparam!(env_rlx, "BestBdStop", RET["BestObjVal"]-((ops.GapTol)*abs(RET["BestObjVal"])) );
    setsolver(m_rlx, GurobiSolver(env_rlx));

    status = solve(m_rlx; suppress_warnings=true)
    dx = getvalue(x)
    dy = fill(NaN, P.n);
    for i in P.S
        dy[i] = getvalue(y[i])
    end
    (var_viol, obj_viol, cons_viol) = CheckViolation(dx, dy, P)

    info = GetSolverInfo(m_rlx, "Gurobi", status, ExpectBnB)
    RET["BestBound"] = max(RET["BestBound"], info.bound);
    #println("info.bound =  ", info.bound)
    push!(RET["Iters"], info)
    return (var_viol, obj_viol, cons_viol, info)
end

"Currently iterative refinement. To do: add rounding and heuristics."
function CDA_MIQCP_Driver(P::CDA.MIQCP, solvername::String, ops::CDA.CDAoption = CDA.CDAoption())
    tstart = time_ns()
    assert(solvername=="Gurobi")
    assert(length(P.S_bin) + length(P.S_cont) == P.n)

    RET = Dict{String,Any}();
    push!(RET,"Iters"=>Array{CDA.SolverInfo}(0));
    RET["BestObjVal"] = Float64(+Inf);
    RET["BestBound"]  = Float64(-Inf);

    NumQC = length(P.Ql) # number of constraints
    ### Start Building relaxation model
    env_rlx = Gurobi.Env();
    Gurobi.setparam!(env_rlx, "OutputFlag", ops.display);
    Gurobi.setparam!(env_rlx, "MIQCPMethod", 0);
    Gurobi.setparam!(env_rlx, "PreQLinearize", 2);
    Gurobi.setparam!(env_rlx, "MIPFocus", 2);

    m_rlx = Model(solver=GurobiSolver(env_rlx));

    @variable(m_rlx, x[i=1:P.n], lowerbound = P.lb[i], upperbound = P.ub[i]);
    @variable(m_rlx, y[i in P.S], lowerbound = (P.lb[i]<=0&&P.ub[i]>=0)?(0):(min(P.lb[i]^2, P.ub[i]^2)),
        upperbound = max(P.lb[i]^2, P.ub[i]^2));

    for i in P.S_bin
        if (i in P.S)
            @constraint(m_rlx, y[i]==x[i])
        end
        setcategory(x[i], :Bin)
    end

    # add diagonal RLT and convex part
    for i in intersect(P.S_cont, P.S)
        if ( P.lb[i]>-Inf && P.ub[i]<Inf )
            @constraint(m_rlx, y[i]<=(P.ub[i]+P.lb[i])*x[i]-P.ub[i]*P.lb[i])
            if (ops.UseCQ==1)
                @constraint(m_rlx, x[i]^2 - y[i] <= 0)
            end
        end
    end

    # Add linear Constraints
    for i in P.S_lin
        if (P.LHS[i]>-Inf)
            @constraint( m_rlx, P.A[i,:]'*x >= P.LHS[i] );
        end
        if (P.RHS[i]<Inf)
            @constraint( m_rlx, P.A[i,:]'*x <= P.RHS[i] );
        end
    end
    # Add diagonal constraints
    for i in P.S_diag
        if (P.LHS[i]>-Inf)
            @constraint( m_rlx, sum(P.Ql[i][j,j]*y[j] for j in P.QuadIdx[i]) + P.A[i,:]'*x >= P.LHS[i] );
        end
        if (P.RHS[i]<Inf)
            @constraint( m_rlx, sum(P.Ql[i][j,j]*y[j] for j in P.QuadIdx[i]) + P.A[i,:]'*x <= P.RHS[i] );
        end
    end

    NumUnbdCont = count((P.lb[i]==-Inf || P.ub[i]==Inf) for i in 1:P.n)
    if (NumUnbdCont>0) # Bound inference if necssary
        @printf("Warning: %4d Unbounded vars exist. \n", NumUnbdCont)
        BoundTightByFeas(P, m_rlx, x, y);
        NumUnbdCont = count((P.lb[i]==-Inf || P.ub[i]==Inf) for i in P.S_cont);
        if (NumUnbdCont>0)
            error("Error: fail to infer bounds for $NumUnbdCont cont. vars. \n")
        end
        for i in P.S # TODO: this should checked inside BoundTightByFeas
            setlowerbound( y[i], (P.lb[i]<0&&P.ub[i]>0)?(0):(min(P.lb[i]^2, P.ub[i]^2)));
            setupperbound( y[i], max(P.lb[i]^2, P.ub[i]^2) );
        end
    end
    T_sdp = 0.0;
    for i in P.S_quad # processing QCs
        Si = P.QuadIdx[i];
        isLeftBounded  = (P.LHS[i]!=-Inf);
        isRightBounded = (P.RHS[i]!=Inf);
        (Convexity, pvec_L, pvec_R, Trun) = DiagPerturb(full(P.Ql[i]), Si,
                (P.LHS[i]!=-Inf), (P.RHS[i]!=Inf));
        T_sdp = T_sdp + Trun;
        if (isLeftBounded)
            if (Convexity==-1)
                @constraint( m_rlx, x'*(P.Ql[i])*x + P.A[i,:]'*x >= P.LHS[i] );
            else
                @constraint( m_rlx, x'*(P.Ql[i]-diagm(pvec_L))*x + P.A[i,:]'*x
                    >= P.LHS[i] - sum(y[j]*pvec_L[j] for j in Si) );
            end
        end
        if (isRightBounded)
            if (Convexity==1)
                @constraint( m_rlx, x'*P.Ql[i]*x + P.A[i,:]'*x <= P.RHS[i] );
            else
                @constraint( m_rlx, x'*(P.Ql[i]+diagm(pvec_R))*x + P.A[i,:]'*x
                    <= P.RHS[i] + sum(y[j]*pvec_R[j] for j in Si) );
            end
        end
        push!(P.AugList, i=>CDA.VectorPair(pvec_L, pvec_R))
    end
    # set objective!
    if (iszero(P.Qobj))
        @objective(m_rlx, Min, P.cobj'*x + P.constant)
    else
        (Val, _) = findmax(abs.(P.Qobj),2)
        Si = find(Val)
        (Convexity, _, pvec_R,Trun) = DiagPerturb(full(P.Qobj), Si, false, true);
        T_sdp = T_sdp + Trun;
        if (Convexity==1)
            @objective(m_rlx, Min, x'*P.Qobj*x + P.cobj'*x + P.constant)
        else
            @objective(m_rlx, Min, x'*(P.Qobj+diagm(pvec_R))*x + P.cobj'*x -
                        sum( y[j]*pvec_R[j] for j in Si) );
            push!(P.AugList, 0=>CDA.VectorPair([], pvec_R));
        end
    end
    @printf("Perturbing quad diagonals... DSDP time: %6.3f\n", T_sdp)

    _dalist = Dict{Int, CDA.DAlist}();
    for j in intersect(P.S_cont, P.S)
        #push!(_dalist, j=>SetUpApprox(m_rlx, x[j], y[j], inst.lb[j], inst.ub[j],
        # ops.initialNu));
        push!(_dalist, j=>CDA.DAlist(0,P.lb[j], P.ub[j],[],[],[],[],[]))
    end

    # Set up the heuristics model
    (ops.display>0)?(ipopt_print_level=2):(ipopt_print_level=0)
    m_heu = Model(solver=IpoptSolver(print_level=ipopt_print_level))
    @variable(m_heu, xH[i=1:P.n], lowerbound = P.lb[i], upperbound = P.ub[i]);
    for i in 1:length(P.Ql) # processing QCs
        if (P.LHS[i]!=-Inf)
            @constraint( m_heu, xH'*P.Ql[i]*xH + P.A[i,:]'*xH >= P.LHS[i] );
        end
        if (P.RHS[i]!=Inf)
            @constraint( m_heu, -xH'*P.Ql[i]*xH - P.A[i,:]'*xH >= -P.RHS[i] );
        end
    end
    @objective(m_heu, Min, xH'*P.Qobj*xH + P.cobj'*xH + P.constant)

    # First heuristics run, no initial point
    info = HeuristicRun(m_heu, RET, xH);
    Gap  = (RET["BestObjVal"]-RET["BestBound"])/(1E-6+abs(RET["BestObjVal"]));
    @printf("Iter     STATUS     Time    NodeCount      OBJVAL      BOUND      Gap       ObjViol     ConsViol  TotalTime\n" )
    @printf("%4c  %12s  %7.3f  %7d    %10.3e  %10.3e  %10.1e   %10.3e   %10.3e   %10.3f\n",
            'H', info.status, info.time, info.nodecount, RET["BestObjVal"], RET["BestBound"],
            Gap, 0.0, 0.0, (time_ns()-tstart)/(1E9))

    # write a macro for this?
    TimeLeft = ops.TimeLimit-(time_ns()-tstart)/(1.0e9);
    if (TimeLeft<0)
        println("User time limit reached.")
        RET["TotalTime"] = (time_ns()-tstart)/(1.0e9);
        return RET;
    end
    (var_viol, obj_viol, cons_viol, info) = RelaxationRun(P, m_rlx, x, y, RET, env_rlx, TimeLeft, ops, false)

    Gap  = (RET["BestObjVal"]-RET["BestBound"])/(1E-6+abs(RET["BestObjVal"]));
    @printf("%4d  %12s  %7.3f  %7d    %10.3e  %10.3e  %10.1e   %10.3e   %10.3e   %10.3f\n",
            0, info.status, info.time, info.nodecount, RET["BestObjVal"], RET["BestBound"], Gap, obj_viol,
            cons_viol, (time_ns()-tstart)/(1E9) )

    info = HeuristicRun(m_heu, RET, xH, getvalue(x) );

    Gap  = (RET["BestObjVal"]-RET["BestBound"])/(1E-6+abs(RET["BestObjVal"]));
    @printf("%4c  %12s  %7.3f  %7d    %10.3e  %10.3e  %10.1e   %10.3e   %10.3e   %10.3f\n",
            'H', info.status, info.time,
            info.nodecount, RET["BestObjVal"], RET["BestBound"], Gap, 0.0, 0.0, (time_ns()-tstart)/(1E9))

    (TerminationCheck(RET["BestObjVal"], RET["BestBound"], ops))?(RET["TotalTime"] = (time_ns()-tstart)/(1E9);return RET):();
    for iter = 1:ops.maxIter
        dy = fill(NaN,P.n);
        for i in P.S
            dy[i] = getvalue(y[i])
        end
        HasModification = false;
        for i = 1:min(ops.NumRefine, length(var_viol))
            if (var_viol[i][2] > 1E-5) # variable violation is big
                j = var_viol[i][1] # variable index
                RefineApprox!(m_rlx, get(_dalist, j, nothing), x[j], y[j])
                HasModification = true
                #push!(_dalist, j=>SetUpApprox(m_rlx, x[j], y[j], inst.lb[j], inst.ub[j], nu));
            end
        end
        if (HasModification)
            TimeLeft = ops.TimeLimit-(time_ns()-tstart)/(1.0e9);
            if (TimeLeft<0)
                println("User time limit reached.")
                break;
            end

            # Relaxation
            (var_viol, obj_viol, cons_viol, info) = RelaxationRun(P, m_rlx, x, y, RET, env_rlx, TimeLeft, ops, true)
            Gap  = (RET["BestObjVal"]-RET["BestBound"])/(1E-6+abs(RET["BestObjVal"]));
            @printf("%4d  %12s  %7.3f  %7d    %10.3e  %10.3e  %10.1e   %10.3e   %10.3e   %10.3f\n",
                    iter, info.status, info.time,
                    info.nodecount, RET["BestObjVal"], RET["BestBound"],  Gap, obj_viol, cons_viol, (time_ns()-tstart)/(1E9))

            # Heuristics
            info = HeuristicRun(m_heu, RET, xH, getvalue(x));
            Gap  = (RET["BestObjVal"]-RET["BestBound"])/(1E-6+abs(RET["BestObjVal"]));
            @printf("%4c  %12s  %7.3f  %7d    %10.3e  %10.3e  %10.1e   %10.3e   %10.3e   %10.3f\n",
                    'H', info.status, info.time, info.nodecount, RET["BestObjVal"], RET["BestBound"], Gap,
                    0.0, 0.0, (time_ns()-tstart)/(1E9))

            (TerminationCheck(RET["BestObjVal"], RET["BestBound"], ops))?(RET["TotalTime"] = (time_ns()-tstart)/(1E9);return RET):();
        else
            println("no more refinement.")
            break;
        end
    end
    println("Max # of iterations reached.")
    RET["TotalTime"] = (time_ns()-tstart)/(1E9);
    return RET;
end

"Solve BoxQP instances with CPLEX, return an object of type CDA.SolverInfo.
"
function MIQCP_CplexDriver(P::CDA.MIQCP)
    m = Model(solver=CplexSolver(CPX_PARAM_OPTIMALITYTARGET=3, CPX_PARAM_TILIM=300))
    #m = Model(solver=GurobiSolver())
    lb = P.lb;
    ub = P.ub;
    n  = P.n;
    @variable(m, x[i=1:n], lowerbound = lb[i], upperbound = ub[i])
    for i=1:n
        if (inst.VarType[i]==CDA.BINARY)
            setcategory(x[i], :Bin)
        elseif (inst.VarType[i]==CDA.CONTINUOUS)
            setcategory(x[i], :Cont)
        elseif (inst.VarType[i]==CDA.INTEGER)
            setcategory(x[i], :Integer)
        else
            error("unsupported variable type.")
        end
    end
    @objective(m, Min, x'*P.Qobj*x + P.cobj'*x)

    for i=1:length(P.LHS)
        if (inst.LHS[i]>-Inf)
            @constraint( m, P.A[i,:]'*x >= P.LHS[i] );
        end
        if (inst.RHS[i]<Inf)
            @constraint( m, P.A[i,:]'*x <= P.RHS[i] );
        end
    end
    #varlist = SetUpApprox(m, x, y, lb, ub, 4);
    #println(m)
    println("Start solving full model: ")
    println(m)
    status = solve(m)
    println("Solution status: ", status)
    println("Objective value: ", getobjectivevalue(m))
    #ret = CDA.SolverInfo("CPLEX", status, getsolvetime(m), getnodecount(m),
    #                        getobjgap(m), getobjbound(m), getobjectivevalue(m));
    return 0;
end

# function GatherInfo(solver, m, status::Symbol)
#
#     if (solver=="Gurobi")
#         SolverInfo("Gurobi", status, Gurobi.getsolvetime(internalmodel(m)),
#         Gurobi.getnodecount(internalmodel(m)),
#         )
#     end
# end

#run(`clear`) # clear output console
#println("JuMP version: ", Pkg.installed("JuMP"))
#println("Gurobi package version: ", Pkg.installed("Gurobi"))



# ###
# ########## Code block: Run BoxQP instances with CPLEX ######################
# outputfile = open("./BoxQP_cplex_all.txt","w");
# folderlist = ( "/Users/hdong/workspace/Manuscripts/cda/data/boxqp/basic",
# "/Users/hdong/workspace/Manuscripts/cda/data/boxqp/extended",
# "/Users/hdong/workspace/Manuscripts/cda/data/boxqp/extended2"
# )
# for folder in folderlist
#     #folder = "/Users/hdong/workspace/Manuscripts/cda/data/boxqp/basic";
#     for fname in readdir(folder)
#         open(joinpath(folder,fname)) do f
#             (n,c,Q) = ParseBoxQP(f)
#             info = BoxQP_CplexDriver(n,c,Q)
#             write(outputfile, joinpath(folder,fname), "\n")
#             write(outputfile, @sprintf(" %7.3f, %10d, %10.5f, %7.3e, %7.3e\n", info.time,
#             info.nodecount, info.gap, info.bound, info.objval))
#             flush(outputfile)
#         end
#     end
# end
# close(outputfile)
# ########## End Code block: Run BoxQP instances with CPLEX ######################
# ###
