module QuantumGameOfLife

using ITensors
using InteractiveUtils
using FromFile

@from "Parser.jl" using Parser
@from "Types.jl" using Types
@from "HamiltonianMpoCreator.jl" using HamiltonianMpoCreator
@from "InitialStates.jl" using InitialStates
@from "TimeEvolution.jl" using TimeEvolution
@from "Measure.jl" using Measure
@from "Plotting.jl" using Plotting
@from "FragmentationAnalysis.jl" using FragmentationAnalysis

export start

function start()
    my_args = get_args()
    start(my_args)
end

function start(args::String)
    empty!(ARGS)
    append!(ARGS, split(args, " "))
    start()
end

function start(args::Dict{Symbol,Any})
    for num_cells in args[:num_cells] # args[:num_cells] is an array of num_cells values
        ITensors.set_warn_order(num_cells * 2 + 1)
        site_inds = siteinds("Qubit", num_cells)
        H_mpo = build_MPO_hamiltonian(site_inds, args[:distance], args[:activation_interval][1], args[:activation_interval][2], args[:periodic_boundaries])
        rule_filename = "$(num_cells)-$(args[:distance])-$(args[:activation_interval][1])$(args[:activation_interval][2])"
        for state_name in args[:initial_states]
            psi_0_mps = getfield(InitialStates, Symbol(state_name))(site_inds)
            for algorithm in args[:algorithm]
                results = evolve(H_mpo, psi_0_mps, args[:num_steps], args[:step_size], algorithm)
                measure_and_plot(results, "$(args[:plot_file_path])/$(rule_filename)-$(state_name)", args)
            end
        end # for

        if args[:plot_eigval_vs_cbe]
            eigval, cbe = eigval_vs_cbe(H_mpo)
            plot_eigval_vs_cbe(eigval=eigval, cbe=cbe, path="$(args[:plot_file_path])/$(rule_filename)-eigval_vs_cbe", file_formats=args[:file_formats], show=args[:show])
        end

        if args[:plot_fragment_sizes]
            plot_fragment_sizes(fragment_sizes=fragment_sizes(H_mpo, args[:periodic_boundaries]), path="$(args[:plot_file_path])/$(rule_filename)-fragment_sizes", file_formats=args[:file_formats], show=args[:show])
        end
    end # for
end

function measure_and_plot(results::Vector{MPS}, path::String, args::Dict{Symbol,Any})
    planned_measurements = Set{Types.MeasurementType}([Types.ExpectationValue()])
    if args[:plot_sse]
        push!(planned_measurements, Types.SingleSiteEntropy())
    end
    if args[:plot_rounded]
        push!(planned_measurements, Types.Rounded())
    end
    if args[:plot_bond_dims]
        push!(planned_measurements, Types.BondDimensions())
    end
    if args[:plot_cbe]
        push!(planned_measurements, Types.CenterBipartiteEntropy())
    end
    measurements = measure(results, planned_measurements)

    heatmap_continuous_measurements = filter(x -> haskey(measurements, x()), subtypes(Types.HeatmapContinuous))
    heatmap_discrete_measurements = filter(x -> haskey(measurements, x()), subtypes(Types.HeatmapDiscrete))
    line_plot_measurements = filter(x -> haskey(measurements, x()), subtypes(Types.LinePlot))
    heatmaps_continuous = map(x -> LabeledPlot(Types.label(x()), measurements[x()]), heatmap_continuous_measurements)
    heatmaps_discrete = map(x -> LabeledPlot(Types.label(x()), measurements[x()]), heatmap_discrete_measurements)
    line_plots = map(x -> LabeledPlot(Types.label(x()), measurements[x()]), line_plot_measurements)
    plot_results(heatmaps_continuous=heatmaps_continuous, heatmaps_discrete=heatmaps_discrete, line_plots=line_plots, path=path, file_formats=args[:file_formats], show=args[:show])
end
end