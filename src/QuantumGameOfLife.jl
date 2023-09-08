module QuantumGameOfLife

using ITensors

include("utils.jl")
include("constants.jl")
include("Parser.jl")
include("InitialStates.jl")
include("HamiltonianMpoCreator.jl")
include("timeEvolution.jl")
include("measure.jl")
include("fragmentationAnalysis.jl")
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
    rule_filename = "$(args[:distance])$(args[:activation_interval][1])$(args[:activation_interval][2])"
    for state_name in args[:initial_states]
        psi_0_mps = getfield(InitialStates, Symbol(state_name))(site_inds)

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
        if args[:plot_cbe]
            push!(planned_measurements, CenterBipartiteEntropy())
        end
        measurements = measure(results, planned_measurements)

        state_filename = "$(state_name)$(args[:num_cells])"
        heatmap_continuous_measurements = filter(x -> haskey(measurements, x()), subtypes(HeatmapContinuous))
        heatmap_discrete_measurements = filter(x -> haskey(measurements, x()), subtypes(HeatmapDiscrete))
        line_plot_measurements = filter(x -> haskey(measurements, x()), subtypes(LinePlot))
        heatmaps_continuous = map(x -> LabeledPlot(PlotLabels[x()], measurements[x()]), heatmap_continuous_measurements)
        heatmaps_discrete = map(x -> LabeledPlot(PlotLabels[x()], measurements[x()]), heatmap_discrete_measurements)
        line_plots = map(x -> LabeledPlot(PlotLabels[x()], measurements[x()]), line_plot_measurements)
        plot_results(
            heatmaps_continuous=heatmaps_continuous,
            heatmaps_discrete=heatmaps_discrete,
            line_plots=line_plots,
            path="$(args[:plot_file_path])/$(state_filename)_$(rule_filename)",
            file_formats=args[:file_formats],
            show=args[:show]
        )

    end
    if args[:plot_eigval_vs_cbe]
        eigval, cbe = eigval_vs_cbe(H_mpo)
        plot_eigval_vs_cbe(
            eigval=eigval,
            cbe=cbe,
            path="$(args[:plot_file_path])/$(rule_filename)",
            file_formats=args[:file_formats],
            show=args[:show]
        )
    end
end
end