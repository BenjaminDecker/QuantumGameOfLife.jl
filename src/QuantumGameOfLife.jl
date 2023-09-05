module QuantumGameOfLife

using ITensors

include("constants.jl")
include("Parser.jl")
include("InitialStates.jl")
include("HamiltonianMpoCreator.jl")
include("timeEvolution.jl")
include("measure.jl")
include("plotting.jl")

using .Parser
using InteractiveUtils
using ITensors
using .HamiltonianMpoCreator

function start()
    my_args = Parser.parse_args()
    start(my_args)
end

function start(args::String)
    empty!(ARGS)
    append!(ARGS, split(args, " "))
    start()
end

function start(args::Dict{Symbol,Any})
    site_inds = siteinds("Qubit", args[:num_cells])
    H_mpo = build_hamiltonian_MPO(site_inds, args[:distance], args[:activation_interval][1], args[:activation_interval][2], args[:periodic_boundaries])

    for state in args[:initial_states]
        psi_0_mps = getfield(InitialStates, Symbol(state))(site_inds)

        results = evolve(H_mpo, psi_0_mps, args[:num_steps], args[:step_size])

        planned_measurements = Set{MeasurementType}([ExpectationValue()])
        if args[:plot_sse]
            push!(planned_measurements, SingleSiteEntropy())
        end
        if args[:plot_rounded]
            push!(planned_measurements, Rounded())
        end
        if args[:plot_bond_dims]
            push!(planned_measurements, BondDimensions())
        end
        measurements = measure(results, planned_measurements)

        file_name = "$(state)$(args[:num_cells])_$(args[:distance])$(args[:activation_interval][1])$(args[:activation_interval][2])"
        continuous_measurements = filter(x -> haskey(measurements, x()), subtypes(Continuous))
        discrete_measurements = filter(x -> haskey(measurements, x()), subtypes(Discrete))
        continuous_plots = map(x -> LabeledPlot(PlotLabels[x()], measurements[x()]), continuous_measurements)
        discrete_plots = map(x -> LabeledPlot(PlotLabels[x()], measurements[x()]), discrete_measurements)
        plot_results(continuous_plots=continuous_plots, discrete_plots=discrete_plots, path="plots/$(file_name)", file_formats=args[:file_formats])
    end
end
end