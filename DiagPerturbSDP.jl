"Solve SDP problem
    min c'* x s.t. Q + diag(x) >=0.
The format of sdpa is desribed here: http://plato.asu.edu/ftp/sdpa_format.txt
whose primal problem is
(P)    min c1*x1+c2*x2+...+cn*xn
       st  F1*x1+F2*x2+...+Fn*xn-F0=X
                                 X >= 0
so F0 will be -Q, Fi will be a zero matrix with Fi[i,i]=1.
"
function writesdpa(fout::String, c::AbstractVector{Float64}, Q::AbstractMatrix{Float64})
    n = length(c)
    assert(size(Q)==(n,n))
    f = open(fout,"w")
    write(f, "$n == number of xi vars \n")
    write(f, "1 == number of blocks \n")
    write(f, "$n \n")
    for i=1:n
        write(f, "$(c[i]) ")
    end
    write(f, "\n");
    for i in 1:n
        for j in i:n
            if !iszero(Q[i,j])
                write(f, "0 1 $i $j $(-Q[i,j]) \n")
            end
        end
    end
    for i in 1:n
        write(f, "$i 1 $i $i 1 \n")
    end
    close(f)
end

function ComputeDiagPerturbByDSDP(c::AbstractVector{Float64}, Q::AbstractMatrix{Float64})
    fname = ".sdpa_file_4dsdp"
    writesdpa(fname, c, Q)
    solname = ".sdpa_sol_dsdp"
    #(si,pr) = Base.writesto(`dsdp $fname -save $solname`)
    #Base.wait(pr)
    #tic();
    tstart = time_ns();
    Base.run(pipeline(`dsdp $fname -save $solname -dloginfo 0 -dlogsummary 0`,stdout=DevNull))
    T = (time_ns() - tstart)/1E9;
    #toc();
    #@printf("DSDP solving SDP of size %4d ..., T=%5.3f (s)\n", length(c), T)
    f = open(solname)
    A = split(readline(f))
    dsol = parse.(Float64, A)
    close(f)
    Base.run(`rm $fname` & `rm $solname`)
    return (dsol,T)
end

#tic();
#n = 200;
#dsol = ComputeDiagPerturbByDSDP(rand(n), randn(n,n))
#toc();

# some warning:
# WARNING: The addition operator has been used on JuMP expressions a large number of times.
# This warning is safe to ignore but may indicate that model generation is slower than necessary.
# For performance reasons, you should not add expressions in a loop. Instead of x += y, use append!(x,y)
# to modify x in place. If y is a single variable, you may also use push!(x, coef, y) in place of x += coef*y.
