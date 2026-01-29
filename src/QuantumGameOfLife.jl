module QuantumGameOfLife

using ITensors
using ITensorMPS
using InteractiveUtils
using ProgressMeter
using CairoMakie
using SplitApplyCombine
using DefaultApplication
using TensorTimeSteps
using Random
using InteractiveUtils
using Combinatorics

import ArgParse: ArgParse, ArgParseSettings, parse_args, add_arg_group!, @add_arg_table!
include("parsing/types.jl")
include("parsing/initialStates.jl")
include("parsing/parsing.jl")

include("utils.jl")
include("hamiltonianMpoCreation.jl")
include("algorithms/exact.jl")
include("algorithms/tdvp.jl")
include("algorithms/sierpinski.jl")
include("algorithms/tebd.jl")
include("measuring.jl")
include("plotting.jl")
include("fragmentation_analysis/fragmentationAnalysis.jl")

export start

function start(args::String)
    empty!(ARGS)
    append!(ARGS, split(args, " "))
    start()
end

start() = start(get_args())

start(::Type{Nothing}) = return

function start(args::Args)
    ITensors.set_warn_order(args.num_cells * 2 + 1)
    H = build_MPO_hamiltonian(args.site_inds, args)

    if !isempty(args.initial_states) && !isempty(args.plots)
        psi_0_vec = if args.superposition
            [normalize(sum(args.initial_states))]
        else
            [normalize(state) for state in args.initial_states]
        end

        # If the only given plot is classical, no need to run quantum versions
        results = if length(args.plots) == 1 && in(Classical(), args.plots)
            [[psi_0 for _ in 1:args.num_steps] for psi_0 in psi_0_vec]
        else
            evolve(args.algorithm, psi_0_vec, H, args)
        end
        measurements = [measure(result, args) for result in results]
        plot(measurements, args)
    end # if

    if args.plot_eigval_vs_cbe
        eigval, cbe = eigval_vs_cbe(H)
        plot_eigval_vs_cbe(eigval, cbe, args)
    end

    if args.plot_fragment_sizes
        plot_fragment_sizes(
            fragment_sizes(
                H,
                args.periodic_boundaries
            ),
            args
        )
    end
end
end
