module QuantumGameOfLife

using ITensors
using ITensorMPS
using ProgressMeter
using CairoMakie
using SplitApplyCombine
using DefaultApplication
using TensorTimeSteps
using Random
using Logging

import ArgParse: ArgParse, ArgParseSettings, parse_args, add_arg_group!, @add_arg_table!

include("parsing/types.jl")
include("parsing/initial_states.jl")
include("parsing/parsing.jl")

include("utils.jl")
include("hamiltonian_mpo_creation.jl")
include("algorithms/exact.jl")
include("algorithms/tdvp.jl")
include("algorithms/sierpinski.jl")
include("algorithms/tebd.jl")
include("measuring.jl")
include("plotting.jl")
include("fragmentation_analysis/utils.jl")
include("fragmentation_analysis/fragmentation_analysis.jl")

export start

"""
    start()

Parse the command line arguments and run the simulation. Returns `nothing` without doing
any work if the arguments could not be parsed (for example when `--help` was requested).
"""
function start()
    args = parse_commandline()
    isnothing(args) && return nothing
    start(args)
end

"""
    start(args::AbstractString)

REPL convenience wrapper: split `args` on whitespace, treat the result as the command line
arguments and run the simulation. For example `start("--show --rule 150")`.
"""
function start(args::AbstractString)
    empty!(ARGS)
    append!(ARGS, split(args))
    start()
end

"""
    start(args::Args)

Run the simulation and produce the requested plots for the parsed configuration `args`.
"""
function start(args::Args)
    ITensors.set_warn_order(args.num_cells * 2 + 1)
    H = build_hamiltonian_mpo(args.site_inds, args)

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

end # module
