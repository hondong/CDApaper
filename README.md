# CDApaper
This repository holds Julia/JuMP code for the manuscript "Compact Disjunctive Approximations to Nonconvex Quadratically Constrained Programs" by H. Dong and Y. Luo 
submitted for publication in Nov. 2018.

This code has been tested with Julia 0.6.4, JuMP 0.18.2 on a Mac OS X system. To run these tests you need CPLEX, Gurobi, DSDP and 
IPOPT to be properly installed. The code currently does very little error checking. I plan to further clean the code and improve this 
README file hopefully in the near future.  

Assuming you have all the packages properly installed (e.g., JuMP.jl, CPLEX.jl, Gurobi.jl, Ipopt.jl) and can run DSDP solver by "dsdp" 
in command line, then the numerical experiments in our paper can be reproduced by running the following commands in Julia:

julia> include("BoxQP.jl");

julia> BoxQPTest(["./data/boxqp/small"],"boxqp\_output.txt")

julia> include("QCQP.jl"); 

julia> RandQCQPTest(["./data/qcqp/small"],"qcqp\_output.txt")
