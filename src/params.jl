#  Copyright 2017, V.Leclere, H.Gerard, F.Pacaud, T.Rigaut
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
#  Definition of SDDP parameters
#############################################################################

type SDDPparameters
    # Solver used to solve LP
    SOLVER::MathProgBase.AbstractMathProgSolver
    # Solver used to solve MILP (default is nothing):
    MIPSOLVER::Nullable{MathProgBase.AbstractMathProgSolver}
    # number of scenarios in the forward pass
    forwardPassNumber::Int64
    # tolerance upon confidence interval:
    confidence_level::Float64
    # Estimate upper-bound every %% iterations:
    compute_ub::Int64
    # Number of MonteCarlo simulation to perform to estimate upper-bound:
    monteCarloSize::Int64
    # Number of MonteCarlo simulation to estimate the upper bound during one iteration
    in_iter_mc::Int64
    # Refresh JuMP Model:
    reload::Int

    function SDDPparameters(solver; passnumber=10, gap=0., confidence=.975,
                            max_iterations=20, prune_cuts=0,
                            pruning_algo="none",
                            compute_ub=-1, montecarlo_final=1000, montecarlo_in_iter=100,
                            mipsolver=nothing,
                            rho0=0., alpha=1., reload=-1)

        return new(solver, mipsolver, passnumber, confidence,
                   compute_ub, montecarlo_final, montecarlo_in_iter, reload)
    end
end


"""
Test compatibility of parameters.

# Arguments
* `model::SPModel`:
    Parametrization of the problem
* `param::SDDPparameters`:
    Parameters of SDDP
* `verbose:Int64`:

# Return
`Bool`
"""
function check_SDDPparameters(model::SPModel, param::SDDPparameters, verbose=0::Int64)
    if model.IS_SMIP && isnull(param.MIPSOLVER)
        error("MIP solver is not defined. Please set `param.MIPSOLVER`")
    end

    (verbose > 0) && (model.IS_SMIP) && println("SMIP SDDP")
end

abstract SdpParameters

type ExhaustiveSdpParameters <: SdpParameters

    stateSteps::Array
    controlSteps::Array
    infoStructure::Symbol
    expectation_computation::String
    monteCarloSize::Int
    buildSearchSpace::Nullable{Function}
    dynamicsType::Symbol

    function ExhaustiveSdpParameters(stateSteps, controlSteps; infoStructure = :dh,
                            expectation_computation="Exact" ,monteCarloSize=1000,
                            search_space_builder = Nullable{Function}(), dynamicsType = :classic)

        (expectation_computation != "Exact") && (expectation_computation != "MonteCarlo") && 
        warn("Expectation computation defaulted to Exact") && (expectation_computation="Exact")

        (dynamicsType == :classic) || (dynamicsType == :random) || warn("Dynamics type defaulted to classic")


        return new(stateSteps, controlSteps, infoStructure,
                    expectation_computation, monteCarloSize, 
                    search_space_builder, dynamicsType)
    end

end

type MathProgSdpParameters <: SdpParameters

    stateSteps::Array
    solver::MathProgBase.AbstractMathProgSolver
    programType::Symbol # LP or MILP (SOS2 constraints) or NLP
    infoStructure::Symbol
    expectation_computation::String
    monteCarloSize::Int
    buildSearchSpace::Nullable{Function}

    function MathProgSdpParameters(stateSteps::Array, solver; programType = :lp,
                                infoStructure = :dh,
                                expectation_computation="Exact",
                                monteCarloSize=1000, search_space_builder = Nullable{Function}())

        if (expectation_computation != "Exact") && (expectation_computation != "MonteCarlo")
            warn("Expectation computation defaulted to Exact")
            expectation_computation="Exact"
        end

        return new(stateSteps, solver, programType, infoStructure,
                    expectation_computation, monteCarloSize, search_space_builder)
    end

end