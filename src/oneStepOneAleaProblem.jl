#  Copyright 2017, V.Leclere, H.Gerard, F.Pacaud, T.Rigaut
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# Model and solve the One-Step One Alea problem in different settings
# - used to compute the optimal control (in forward phase / simulation)
# - used to compute the cuts in the Backward phase
#############################################################################

"""
Solve the Bellman equation at time t starting at state x under alea xi
with the current evaluation of Vt+1

# Description
The function solve
min_u current_cost(t,x,u,xi) + current_Bellman_Value_{t+1}(dynamic(t,x,u,xi))
and can return the optimal control and a subgradient of the value of the
problem with respect to the initial state x

# Arguments
* `model::SPmodel`:
the stochastic problem we want to optimize
* `param::SDDPparameters`:
the parameters of the SDDP algorithm
* `m::JuMP.Model`:
The linear problem to solve, in order to approximate the
current value functions
* `t::int`:
time step at which the problem is solved
* `xt::Array{Float}`:
current starting state
* `xi::Array{float}`:
current noise value
* `relaxation::Bool`: default is false
If problem is MILP, specify if it is needed to relax integer constraints.
* `init::Bool`:
If specified, approximate future cost as 0

# Returns
* `solved::Bool`:
True if the solution is feasible, false otherwise
* `NextStep`:
Store solution of the problem
* `ts::Float64`:
Solver's execution time
"""
function solve_one_step_one_alea(model,
                                 param,
                                 m::JuMP.Model,
                                 t::Int64,
                                 xt::Vector{Float64},
                                 xi::Vector{Float64};
                                 relaxation=false::Bool,
                                 init=false::Bool,
                                 verbosity::Int64=0)

    # Get var defined in JuMP.model:
    x = m[:x]
    u = m[:u]
    w = m[:w]
    alpha = m[:alpha]

    # Update value of w:
    JuMP.fix.(w, xi)

    #update objective
    if isa(model.costFunctions, Function)
        try
            @objective(m, Min, model.costFunctions(t, x, u, xi) + alpha)
        catch
            @objective(m, Min, model.costFunctions(m, t, x, u, xi) + alpha)
        end
    elseif isa(model.costFunctions, Vector{Function})
        cost = m[:cost]
        for i in 1:length(model.costFunctions)
            @constraint(m, cost >= model.costFunctions[i](t, x, u, xi))
        end
        @objective(m, Min, cost + alpha)
    end

    # Update constraint x == xt
    JuMP.fix.(m[:x_constant], xt)

    if false
        println("One step one alea problem at time t=", t)
        println("for x =", xt)
        println("and w=", xi)
        print(m)
    end

    JuMP.optimize!(m)

    # get time taken by the solver:
    solvetime = 0

    if true
        optimalControl = JuMP.value.(u)
        # Return object storing results:
        result = NLDSSolution(
                              true,
                              JuMP.objective_value(m),
                              model.dynamics(t, xt, optimalControl, xi),
                              optimalControl,
                              JuMP.dual.(m.ext[:cons]),
                              JuMP.value(alpha),
                              getcutsmultipliers(m))
    else
        println(m)
        println(status)
        error("Foo")
        # If no solution is found, then return nothing
        result = NLDSSolution()
    end

    return result, solvetime
end


"""Solve model in Decision-Hazard."""
function solve_dh(model, param, t, xt, m; verbosity::Int64=0)
    xf = m[:xf]
    u = m[:u]
    alpha = m[:alpha]
    JuMP.fix.(m[:x_constant], xt)

    (verbosity>5) && println("Decision Hazard model")
    (verbosity>5) && print(m)

    status = JuMP.optimize!(m)
    solved = JuMP.termination_status(m) == JuMP.MOI.OPTIMAL

    if ~solved
        println(m)
        println(getvalue(u))
        println(getvalue(alpha))
        println(getvalue(xf))
        warn("dh model not solved at time t=",t)
    end

    # TODO
    solvetime = 0

    if solved
        # Computation of subgradient:
        λ = Float64[JuMP.dual(m.ext[:cons][i]) for i in 1:model.dimStates]
        result = DHNLDSSolution(solved,
                                JuMP.objective_value(m),
                                JuMP.value.(xf),
                                JuMP.value.(u)[:, 1],
                                λ,
                                JuMP.value.(alpha),
                                getcutsmultipliers(m))
    else
        # If no solution is found, then return nothing
        result = NLDSSolution()
    end

    return result, solvetime
end


# Solve local problem with a quadratic penalization:
function regularize(model, param,
                    regularizer::AbstractRegularization,
                    m::JuMP.Model,
                    t::Int64,
                    xt::Vector{Float64}, xi::Vector{Float64}, xp::Vector{Float64};verbosity=0::Int64)
    # store current objective:
    pobj = m.obj
    xf = JuMP.index(m, :xf)
    qexp = getpenaltyexpr(regularizer, xf, xp)
    # and update model objective:
    @objective(m, :Min, m.obj + qexp)
    res = solve_one_step_one_alea(model, param, m, t, xt, xi,verbosity=verbosity)
    m.obj = pobj

    return res
end


"""Solve relaxed MILP problem."""
function solve_relaxed!(m, param,verbosity::Int64=0)
    setsolver(m, param.SOLVER)
    status = JuMP.optimize!(m)
    return status == JuMP.MOI.OPTIMAL
end


"""Solve original MILP problem."""
function solve_mip!(m, param,verbosity::Int64=0)
    setsolver(m, get(param.MIPSOLVER))
    status = solve(m, relaxation=false)
    return status == :Optimal
end


getcutsmultipliers(m::JuMP.Model)=zeros(m.ext[:ncuts])
#_getdual(m)[end-m.ext[:ncuts]+1:end]
#= _getdual(m::JuMP.Model) = JuMP.getdual(m) =#
