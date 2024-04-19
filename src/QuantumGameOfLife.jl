module QuantumGameOfLife

using ITensors
using InteractiveUtils
using ProgressMeter

include("parsing/parsing.jl")
include("utils.jl")
include("hamiltonianMpoCreation.jl")
include("initialStates.jl")
include("algorithms/exact.jl")
include("algorithms/tdvp.jl")
include("algorithms/sierpinski.jl")
include("measuring.jl")
include("plotting.jl")
include("fragmentation_analysis/fragmentationAnalysis.jl")

export start


start() = start(get_args())

function start(args::String)
    empty!(ARGS)
    append!(ARGS, split(args, " "))
    start()
end

function start(args::Args)
    ITensors.set_warn_order(args.num_cells * 2 + 1)
    site_inds = siteinds("Qubit", args.num_cells)
    H_MPO = build_MPO_hamiltonian(site_inds, args)

    if length(args.initial_state) > 0 && length(args.plots) > 0
        psi_0_mps = sum([getfield(QuantumGameOfLife, Symbol(state_name))(site_inds) for state_name in args.initial_state])
        normalize!(psi_0_mps)
        results = evolve(args.algorithm, H_MPO, psi_0_mps, args)
        measurements = measure(results, args)
        plot(measurements, args)
    end # if

    if args.plot_eigval_vs_cbe
        eigval, cbe = eigval_vs_cbe(H_MPO)
        plot_eigval_vs_cbe(eigval, cbe, args)
    end

    if args.plot_fragment_sizes
        plot_fragment_sizes(
            fragment_sizes(
                H_MPO,
                args.periodic_boundaries
            ),
            args
        )
    end
end
end
