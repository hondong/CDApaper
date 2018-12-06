using Ipopt

function Ipopt_Driver(P::CDA.MIQCP, xinit::AbstractVector{Float64})
    RET = Dict{String,AbstractVector{Float64}}();
    #S_cont = find(P.VarType[i]==CDA.CONTINUOUS for i=1:P.n)
    #assert(length(S_cont)==P.n)

    m = Model(solver=IpoptSolver())

    @variable(m, x[i=1:P.n], lowerbound = P.lb[i], upperbound = P.ub[i], start=xinit[i]);
    for i in 1:length(P.Ql) # processing QCs
        isLeftBounded = (P.LHS[i]!=-Inf);
        isRightBounded = (P.RHS[i]!=Inf);
        if (isLeftBounded)
            @constraint( m, x'*P.Ql[i]*x + P.A[i,:]'*x >= P.LHS[i] );
        end
        if (isRightBounded)
            @constraint( m, -x'*P.Ql[i]*x - P.A[i,:]'*x >= -P.RHS[i] );
        end
    end

    # set objective!
    @objective(m, Min, x'*P.Qobj*x + P.cobj'*x + P.constant)

    status = solve(m)
    println(m)
    return getvalue(x)
end
