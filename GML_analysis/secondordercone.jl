using LinearAlgebra
using JuMP
using Mosek
using MosekTools
m = Model(with_optimizer(Mosek.Optimizer, QUIET=true, MSK_DPAR_INTPNT_CO_TOL_REL_GAP=1e-3,
    MSK_DPAR_INTPNT_CO_TOL_MU_RED=1e-3))
@variable(m,x)
@variable(m,y)
# @variable(m,v)
@variable(m,u)
@constraint(m, 0.05*x>=y)
@constraint(m, 10>=u)
# @constraint(m, 10>=v)
@constraint(m, [0.5*u, u, x, y] in RotatedSecondOrderCone());
@objective(m, Min, -x)
status=optimize!(m);
cost=-JuMP.objective_value(m);
println(JuMP.value.(x))
println(JuMP.value.(y))
println(JuMP.value.(u))
# println(JuMP.value.(v))
