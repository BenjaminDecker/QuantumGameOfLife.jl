module QuantumGameOfLife

using ITensors
using InteractiveUtils

include("utils.jl")
include("types.jl")
include("parsing.jl")
include("hamiltonianMpoCreation.jl")
include("initialStates.jl")
include("timeEvolution.jl")
include("measure.jl")
include("plotting.jl")
include("fragmentationAnalysis.jl")

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
    num_steps::Int = args[:num_steps]
    distance::Int = args[:distance]
    activation_interval::UnitRange = (args[:activation_interval][1]):(args[:activation_interval][2])
    periodic::Bool = args[:periodic_boundaries]
    step_size::Float64 = args[:step_size]
    initial_state::Vector{String} = args[:initial_state]
    algorithm::String = args[:algorithm]
    num_cells::Int = args[:num_cells]

    ITensors.set_warn_order(num_cells * 2 + 1)
    site_inds = siteinds("Qubit", num_cells)
    H_MPO = build_MPO_hamiltonian(site_inds, distance, activation_interval, periodic)
    rule_filename = "$(num_cells)-$(distance)-$(activation_interval.start)$(activation_interval.stop)-$(num_steps)-$(replace(string(step_size), "." => ""))-$(periodic ? "periodic" : "open")"

    if length(initial_state) > 0
        psi_0_mps = sum([getfield(QuantumGameOfLife, Symbol(state_name))(site_inds) for state_name in initial_state])
        normalize!(psi_0_mps)
        results =
            if algorithm == "exact"
                evolve_exact(H_MPO, psi_0_mps, num_steps, args[:step_size])
            else
                evolve_serpinsky(psi_0_mps, args[:num_steps], args[:step_size], args[:sweeps_per_time_step], args[:max_bond_dim])
            end
        measure_and_plot(results, "$(args[:plot_file_path])/$(rule_filename)-$(replace(string(args[:initial_state]), "\"" => ""))$(args[:max_bond_dim])", args)
    end # if

    if args[:plot_eigval_vs_cbe]
        eigval, cbe = eigval_vs_cbe(H_MPO)
        plot_eigval_vs_cbe(eigval=eigval, cbe=cbe, path="$(args[:plot_file_path])/$(rule_filename)-eigval_vs_cbe", file_formats=args[:file_formats], show=args[:show], num_cells=num_cells)
    end

    if args[:plot_fragment_sizes]
        plot_fragment_sizes(fragment_sizes=fragment_sizes(H_MPO, args[:periodic_boundaries]), path="$(args[:plot_file_path])/$(rule_filename)-fragment_sizes", file_formats=args[:file_formats], show=args[:show])
    end
end

function measure_and_plot(results::Vector{MPS}, path::String, args::Dict{Symbol,Any})
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

    heatmap_continuous_measurements = filter(x -> haskey(measurements, x()), subtypes(HeatmapContinuous))
    heatmap_discrete_measurements = filter(x -> haskey(measurements, x()), subtypes(HeatmapDiscrete))
    line_plot_measurements = filter(x -> haskey(measurements, x()), subtypes(LinePlot))
    heatmaps_continuous = map(x -> LabeledPlot(label(x()), measurements[x()]), heatmap_continuous_measurements)
    heatmaps_discrete = map(x -> LabeledPlot(label(x()), measurements[x()]), heatmap_discrete_measurements)
    line_plots = map(x -> LabeledPlot(label(x()), measurements[x()]), line_plot_measurements)
    plot_results(heatmaps_continuous=heatmaps_continuous, heatmaps_discrete=heatmaps_discrete, line_plots=line_plots, path=path, file_formats=args[:file_formats], show=args[:show])
end
end
