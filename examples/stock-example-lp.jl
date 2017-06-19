#  Copyright 2017, V.Leclere, H.Gerard, F.Pacaud, T.Rigaut
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# Compare different ways of solving a stock problem :
# Min   E [\sum_{t=1}^TF c_t u_t]
# s.t.    s_{t+1} = s_t + u_t - xi_t, s_0 given
#         0 <= s_t <= 1
#         u_min <= u_t <= u_max
#         u_t choosen knowing xi_1 .. xi_t
#############################################################################

using StochDynamicProgramming, CPLEX
println("library loaded")

run_exhaustive_sdp  = true # false if you don't want to run sdp
run_lp_sdp   = true # false if you don't want to run extensive formulation
test_simulation = true # false if you don't want to test your strategies

######## Optimization parameters  ########
# choose the LP solver used.
#SOLVER = ClpSolver() 			   # require "using Clp"
SOLVER = CplexSolver(CPX_PARAM_SIMDISPLAY=0) # require "using CPLEX"

# convergence test
MAX_ITER = 10 # number of iterations of SDDP
step = 0.01   # discretization step of SDP

######## Stochastic Model  Parameters  ########
N_STAGES = 6              # number of stages of the SP problem
COSTS = [sin(3*t)-1 for t in 1:N_STAGES-1]

CONTROL_MAX = 0.5         # bounds on the control
CONTROL_MIN = 0

XI_MAX = 0.3              # bounds on the noise
XI_MIN = 0
N_XI = 10                 # discretization of the noise

S0 = 0.5                  # initial stock

# create law of noises
proba = 1/N_XI*ones(N_XI) # uniform probabilities
xi_support = collect(linspace(XI_MIN,XI_MAX,N_XI))
xi_law = NoiseLaw(xi_support, proba)
xi_laws = NoiseLaw[xi_law for t in 1:N_STAGES-1]

# Define dynamic of the stock:
function dynamic(t, x, u, xi)
    return [x[1] + u[1] - xi[1]]
end

# Define cost corresponding to each timestep:
function cost_t(t, x, u, w)
    return COSTS[t] * u[1]
end

######## Setting up the SPmodel
s_bounds = [(0, 1)] 			# bounds on the state
u_bounds = [(CONTROL_MIN, CONTROL_MAX)] # bounds on controls
spmodel = StochDynModel(N_STAGES,u_bounds,[S0],cost_t,dynamic,xi_laws,
                        xbounds = s_bounds)
println("Model set up")

######### Solving the problem via Dynamic Programming
if run_exhaustive_sdp
    tic()
    println("Starting exhaustive resolution by SDP")
    stateSteps = [step] # discretization step of the state
    controlSteps = [step] # discretization step of the control
    paramSDP = ExhaustiveSdpParameters(stateSteps, controlSteps, infoStructure = :hd)
    Vs = solve_dp(spmodel,paramSDP, 1)
    value_sdp = StochDynamicProgramming.get_bellman_value(spmodel,paramSDP,Vs)
    println("Value obtained by SDP: "*string(round(value_sdp,4)))
    toc(); println();
end

######### Solving the problem via Extensive Formulation
if run_lp_sdp
    tic()
    println("Starting LP resolution by SDP")
    paramLPSDP = MathProgSdpParameters(stateSteps, SOLVER, infoStructure = :hd)
    Vslp = solve_dp(spmodel,paramLPSDP, 1)
    value_sdp = StochDynamicProgramming.get_bellman_value(spmodel,paramLPSDP,Vslp)
    println("Value obtained by SDP: "*string(round(value_sdp,4)))
    toc(); println();
end

######### Comparing the solutions on simulated scenarios.

#srand(1234) # to fix the random seed accross runs
# if run_sddp && run_sdp && test_simulation
#     scenarios = StochDynamicProgramming.simulate_scenarios(xi_laws,1000)
#     costsddp, stocks = forward_simulations(spmodel, paramSDDP, sddp.solverinterface, scenarios)
#     costsdp, states, controls = forward_simulations(spmodel,paramSDP, Vs, scenarios)
#     println("Simulated relative gain of sddp over sdp: "
#             *string(round(200*mean(costsdp-costsddp)/abs(mean(costsddp+costsdp)),3))*"%")
# end
