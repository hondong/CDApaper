module CDA

export DAlist, SolverInfo, MIQCP, CDAoption

mutable struct DAlist
    nu::Int
    l::Float64
    u::Float64
    xi  # list of variables of length nu
    eta # list of variables of length nu
    lam1 # list of variables of length nu
    lam2 # list of variables of length nu
    z # list of variables of length nu
    #DynCons::Array{Any,1}
end

mutable struct SolverInfo
    solvername::String # CPLEX, Gurobi, etc.
    status::String   # optimal, UserLimit, etc.
    time::Float64    # elapsed time
    nodecount::Int   # number of nodes explored
    gap::Float64     # gap in objective value
    bound::Float64   # lower bound
    objval::Float64  # objective value
    numiter::Int     # number of iterations
    SolverInfo() = (new("", "", NaN, -1, NaN, NaN, NaN, -1))
    SolverInfo(_sname, _status, _time, _nc, _gap, _bound, _objval, _niter) = new(_sname, _status, _time, _nc, _gap, _bound, _objval, _niter);
end

mutable struct CDAoption
    UseCQ::Int # do we model y>=x^2 explicitly,  = 1, yes. Default (=0)
    TimeLimit::Number  # time in seconds, Default (=1E10)
    maxIter::Int # Number of iterations in refine approx.
    NumRefine::Int # Max number of refined variables in an iteration.
    initialNu::Int
    display::Int
    GapTol::Float64
    CDAoption() = (new(0, 1E10, 100, 20, 3, 1, 1E-4);)
end

@enum VARTYPE CONTINUOUS BINARY INTEGER SEMICONT UNKNOWN

struct VectorPair
    left::AbstractVector{Float64}
    right::AbstractVector{Float64}
end

#
# min x'*Qobj*x + cobj'*x + constant
# s.t., LHS[i] <= x'*Ql[i]*x + A[i,:]'*x <= RHS[i]
#       lb <= x <= ub

mutable struct MIQCP
    n::Int  # number of total variables
    lb::Array{Float64,1}                 # lower bound
    ub::Array{Float64,1}                 # upper bound
    VarType::Array{VARTYPE,1}            #subset of binary variables
    Qobj::AbstractMatrix{Float64}        # quadratic objective, if exists
    cobj::AbstractVector{Float64}        # linear coefficients
    constant::Float64                    # constant term in objective
    Ql::Array{AbstractMatrix{Float64},1} # quadratic forms
    #cl::Array{AbstractVector{Float64},1} # linear vectors
    #dl::Array{Float64,1}                 # right hand side
    A::AbstractMatrix{Float64}           # LHS[i] <= x'*Q[i]*x + A[i,:]*x <= RHS[i]
    LHS::AbstractVector{Float64}         # Left hand side in QCs
    RHS::AbstractVector{Float64}         # Right hand side in QCs
    #AugList::Array{AbstractVector{Float64},1}
    AugList::Dict{Int, VectorPair}
    S_bin::AbstractVector{Int}           # set of binary variables
    S_cont::AbstractVector{Int}          # set of continuous variables
    S::AbstractVector{Int}               # set of variables appearing quads
    S_lin::AbstractVector{Int}           # set of linear constraints
    S_diag::AbstractVector{Int}          # set of diagonal quadratic constraints
    S_quad::AbstractVector{Int}          # set of general quadratic constraints
    QuadIdx::Dict{Int, AbstractVector{Int}} # quadratic indices for each constraint
    MIQCP(N::Int) = (new(N,Array{Float64,1}(0),Array{Float64,1}(0),Array{VARTYPE,1}(0),
                    Array{Float64,2}(0,0),Array{Float64,1}(0), 0.0,
                    [],Array{Float64,2}(0,0),Array{Float64}(0), Array{Float64}(0),
                    Dict{Int, VectorPair}(), Vector{Int}(),
                    Vector{Int}(), Vector{Int}(),
                    Vector{Int}(),Vector{Int}(),Vector{Int}(),
                    Dict{Int, AbstractVector{Int}}() ); )
end

function displayMIQCP(P::MIQCP)
    n = P.n;
    assert(length(P.lb)==P.n)
    assert(length(P.ub)==P.n)
    assert(length(P.VarType)==P.n)
    assert(size(P.Qobj)==(P.n,P.n))
    assert(length(P.cobj)==P.n)
    m = length(P.Ql)
    #assert(size(P.A)==(m,n))
    assert(length(P.A)==m*n)
    assert(length(P.LHS)==m)
    assert(length(P.RHS)==m)
    if (iszero(P.Qobj))
        if (iszero(P.cobj))
            println("Feasibility Problem")
        else
            println("Linear Objective Minimization")
        end
    else
        println("Quadratic Objective Minimization")
    end

    println("Variables: ", P.n);

    S_cont = find(P.VarType[i]==CDA.CONTINUOUS for i=1:n);
    @printf(" - Continuous: %4d, LINF: %4d, RINF: %4d, LRINF: %4d \n",
        length(S_cont),
        count(P.lb[i]== -Inf for i in S_cont),
        count(P.ub[i]== Inf for i in S_cont),
        count((P.lb[i]== -Inf)&&(P.ub[i]==Inf) for i in S_cont) )

    S_bin = find(P.VarType[i]==CDA.BINARY for i in 1:n);
    @printf(" - Binary:    %4d, Imp. bounds:%4d \n", length(S_bin),
        count((P.lb[i] < 0 || P.ub[i] > 1) for i in S_bin) )

    S_int = find(P.VarType[i]==CDA.INTEGER for i in 1:n);
    @printf(" - Integer:    %4d, LINF: %4d, RINF: %4d, LRINF: %4d \n",
        length(S_int), count(P.lb[i]== -Inf for i in S_int),
        count(P.ub[i]== Inf for i in S_int),
        count((P.lb[i]== -Inf)&&(P.ub[i]==Inf) for i in S_int) )

    @printf(" - SemiCont:   %4d \n",count(P.VarType[i]==CDA.SEMICONT for i in 1:n) )
    @printf(" - Unknown:    %4d \n",count(P.VarType[i]==CDA.UNKNOWN for i in 1:n) )

    println("Constraints: ", m)
    S_lin = find( iszero(P.Ql[i]) || isempty(P.Ql[i]) for i in 1:m)
    @printf(" - Linear:    %4d, LINF: %4d, RINF: %4d, LRINF: %4d \n", length(S_lin),
        count(P.LHS[i]== -Inf for i in S_lin), count(P.RHS[i]== Inf for i in S_lin),
        count(P.LHS[i]== -Inf && P.RHS[i]==Inf for i in S_lin) )
    S_quad = setdiff(1:m, S_lin);
    @printf(" - Quadratic: %4d, LINF: %4d, RINF: %4d, LRINF: %4d \n", length(S_quad),
        count(P.LHS[i]== -Inf for i in S_quad), count(P.RHS[i]== Inf for i in S_quad),
        count(P.LHS[i]== -Inf && P.RHS[i]==Inf for i in S_quad) )
    if (m>0)
        nnz_quad = sum(nnz(P.Ql[i]) for i in 1:length(P.Ql))
        nnz_A    = nnz(P.A)
    else
        nnz_quad = 0
        nnz_A    = 0
    end

    @printf("nnz(Quad) = %6d, nnz(A) = %6d\n", nnz_quad, nnz_A )
end

function setSetIndices!(P::MIQCP)
    P.S_bin  = find(P.VarType[i]==CDA.BINARY for i=1:P.n)
    P.S_cont = find(P.VarType[i]==CDA.CONTINUOUS for i=1:P.n)
    NumQC    = length(P.Ql)
    P.S_lin  = find( (isempty(P.Ql[i]) || iszero(P.Ql[i])) for i=1:NumQC );
    P.S_diag   = find( isdiag(P.Ql[i]) for i in setdiff(1:NumQC, P.S_lin) );
    P.S_quad   = setdiff(1:NumQC, union(P.S_lin, P.S_diag) );

    (Val, _) = findmax(abs.(P.Qobj), 2)
    P.S = find(Val)
    push!(P.QuadIdx, 0=>P.S);

    for i in union(P.S_diag, P.S_quad)
        (Val, _) = findmax(abs.(P.Ql[i]),2)
        Si = find(Val) # set of nonzero columns
        if (!isempty(Si))
            P.S = union(P.S, Si);
            push!(P.QuadIdx, i=>Si)
        end
    end
end


setObj!(obj::MIQCP, v::AbstractVector{Float64}, Q::AbstractMatrix{Float64}=[], constant::Float64=0.0) =
    (
        if (length(v) == obj.n)
            obj.cobj = v;
            obj.constant = constant;
            if (!isempty(Q))
                (size(Q)==(obj.n, obj.n))?(obj.Qobj=Q;):(error("Qobj incorrect length."))
            end
        else
            error("setObjVec! Incorrect length.")
        end
    )

setBounds!(obj::MIQCP, L::AbstractVector, R::AbstractVector) =
    (if (length(L)!=obj.n || length(R)!=obj.n)
        error("Bounds length does not match!")
    else
        obj.lb = L; obj.ub = R;
    end
    )

setVarType!(obj::MIQCP, Type::Array{VARTYPE,1}) =
    (
        if(length(Type)==obj.n)
            obj.VarType = Type;
        else
            error("setVarType! Incorrect dimension.")
        end
    )

pushQC!(obj::MIQCP, Q::AbstractMatrix, c::AbstractVector, d::Float64) =
        (if (size(Q)==(obj.n,obj.n) &&
             size(c)==(obj.n,)          )
             push!(obj.Ql, Q); push!(obj.cl, c); push!(obj.dl, d);
        else
            error("pushQC! Incorrect dimension.")
        end)

using JuMP
using AmplNLWriter

function GenNL(P::CDA.MIQCP, fname::String)
    m = Model(solver=AmplNLSolver("couenne"))

    @variable(m, x[i=1:P.n], lowerbound = P.lb[i], upperbound = P.ub[i]);
    @variable(m, obj);
    for i=1:P.n
        if (P.VarType[i]==CDA.CONTINUOUS)
            setcategory(x[i], :Cont)
        elseif (P.VarType[i]==CDA.BINARY)
            setcategory(x[i], :Bin)
        elseif (P.VarType[i]==CDA.INTEGER)
            setcategory(x[i], :Int)
        else
            error("Couenne_Driver: Unknown variable type?")
        end
    end

    for i in 1:length(P.Ql) # processing QCs
        isLeftBounded = (P.LHS[i]!=-Inf);
        isRightBounded = (P.RHS[i]!=Inf);

        #@NLexpression(m, QF, x'*full(P.Ql[i])*x + P.A[i,:]'*x)
        @NLexpression(m, QF, sum(P.Ql[i][j,k]*x[j]*x[k] for j=1:P.n, k=1:P.n) + sum(P.A[i,j]*x[j] for j=1:P.n) )
        if (isLeftBounded)
            #@NLconstraint( m, P.LHS[i] <= x'*full(P.Ql[i])*x + P.A[i,:]'*x  );
            @NLconstraint( m,  QF >= P.LHS[i] );
        end
        if (isRightBounded)
            #@NLconstraint( m,  x'*full(P.Ql[i])*x + P.A[i,:]'*x <= P.RHS[i] );
            @NLconstraint( m,  QF <= P.RHS[i] );
        end
    end

    # set objective
    #@objective(m, Min, x'*P.Qobj*x + P.cobj'*x + P.constant)
    @NLexpression(m, QExpr, sum(P.Qobj[i,j]*x[i]*x[j] for i=1:P.n, j=1:P.n) + sum(P.cobj[i]*x[i] for i=1:P.n) )
    @NLobjective(m,  Min,  obj)
    @NLconstraint(m, QExpr-obj <= 0 );

    f = open(fname,"w")
    JuMP.build(m)

    mi = internalmodel(m)
    mi.inner.status = :NotSolved
    mi.inner.solve_exitcode = -1
    mi.inner.solve_result_num = -1
    mi.inner.solve_result = "?"
    mi.inner.solve_message = ""

    AmplNLWriter.make_var_index!(mi.inner)
    AmplNLWriter.make_con_index!(mi.inner)

    AmplNLWriter.write_nl_file(f, mi.inner)
    close(f)
end

end
