module QuantumGameOfLife

using ITensors
using InteractiveUtils

include("types.jl")
include("utils.jl")
include("parsing.jl")
include("hamiltonianMpoCreation.jl")
include("initialStates.jl")
include("timeEvolution.jl")
include("measure.jl")
include("plotting.jl")
include("fragmentationAnalysis.jl")

export start

function start()
    start(get_args())
end

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
        results =
            if args.algorithm == "exact"
                evolve_exact(H_MPO, psi_0_mps, args)
            else
                evolve_serpinsky(psi_0_mps, args)
            end

        measurements = measure(results, args)
        plot(measurements, args)
    end # if

    if args.plot_eigval_vs_cbe
        eigval, cbe = eigval_vs_cbe(H_MPO)
        plot_eigval_vs_cbe(eigval, cbe, args)
    end

    if args.plot_fragment_sizes
        plot_fragment_sizes(
            fragment_sizes(H_MPO, args.periodic_boundaries),
            args
        )
    end
end
end
