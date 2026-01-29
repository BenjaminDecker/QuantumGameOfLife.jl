abstract type PlotType end

abstract type HeatmapContinuous <: PlotType end
struct ExpectationValue <: HeatmapContinuous end
filename_identifier(::ExpectationValue) = "expect"
name(::ExpectationValue) = "Expectation Value"
label(::ExpectationValue) = "Expectation\nValue"
ordering_index(::ExpectationValue) = 2
struct SingleSiteEntropy <: HeatmapContinuous end
filename_identifier(::SingleSiteEntropy) = "sse"
name(::SingleSiteEntropy) = "Single Site Entropy"
label(::SingleSiteEntropy) = "Single Site\nEntropy"
ordering_index(::SingleSiteEntropy) = 3


abstract type HeatmapDiscrete <: PlotType end
struct Rounded <: HeatmapDiscrete end
filename_identifier(::Rounded) = "rounded"
name(::Rounded) = "Rounded"
label(::Rounded) = "Rounded"
ordering_index(::Rounded) = 4
struct BondDimensions <: HeatmapDiscrete end
filename_identifier(::BondDimensions) = "bond_dims"
name(::BondDimensions) = "Bond Dimensions"
label(::BondDimensions) = "Bond\nDimension"
ordering_index(::BondDimensions) = 5
struct Classical <: HeatmapDiscrete end
filename_identifier(::Classical) = "classical"
name(::Classical) = "Classical"
label(::Classical) = "Classical"
ordering_index(::Classical) = 1

abstract type LinePlot <: PlotType end
struct CenterBipartiteEntropy <: LinePlot end
filename_identifier(::CenterBipartiteEntropy) = "cbe"
name(::CenterBipartiteEntropy) = "Center Bipartite Entropy"
label(::CenterBipartiteEntropy) = "Center\nBipartite\nEntropy"
ordering_index(::CenterBipartiteEntropy) = 6
struct Autocorrelation <: LinePlot end
filename_identifier(::Autocorrelation) = "autocorrelation"
name(::Autocorrelation) = "Autocorrelation"
label(::Autocorrelation) = "Autocorrelation"
ordering_index(::Autocorrelation) = 7


Base.isless(a::PlotType, b::PlotType) = ordering_index(a) < ordering_index(b)


abstract type Algorithm end

struct Exact <: Algorithm end
name(::Exact) = "Exact"
struct TDVP1 <: Algorithm end
name(::TDVP1) = "TDVP1"
struct TDVP2 <: Algorithm end
name(::TDVP2) = "TDVP2"
struct Sierpinski <: Algorithm end
name(::Sierpinski) = "Sierpinski"
struct TEBD <: Algorithm end
name(::TEBD) = "TEBD"

struct Args
    num_steps::Int
    distance::Int
    # activation_interval::UnitRange
    rule::Int
    periodic::Bool
    step_size::Float64
    site_inds::Vector{Index{Int64}}
    initial_states_names::Vector{String}
    initial_states::Vector{MPS}
    superposition::Bool
    algorithm::Algorithm
    num_cells::Int
    sweeps_per_time_step::Int
    max_bond_dim::Int
    svd_epsilon::Float64
    operator_set::Int
    periodic_boundaries::Bool
    plots::Set{PlotType}
    plotting_file_path::String
    file_formats::Set{String}
    width::Union{Nothing,Int}
    page_entropy::Bool
    px_per_unit::Float64
    plot_eigval_vs_cbe::Bool
    plot_fragment_sizes::Bool
    show::Bool
end

function get_rule(args::Dict{Symbol,Any})::Int
    if isempty(args[:activation_interval])
        return args[:rule]
    end
    # TODO
    return 150
end

function Args(args::Dict{Symbol,Any})::Args
    site_inds = siteinds("Qubit", args[:num_cells])
    return Args(
        args[:num_steps],
        args[:distance],
        get_rule(args),
        args[:periodic_boundaries],
        args[:step_size],
        site_inds,
        args[:initial_states],
        [getfield(QuantumGameOfLife, Symbol(state_name))(site_inds) for state_name in args[:initial_states]],
        args[:superposition],
        args[:algorithm],
        args[:num_cells],
        args[:sweeps_per_time_step],
        args[:max_bond_dim],
        args[:svd_epsilon],
        args[:operator_set],
        args[:periodic_boundaries],
        Set{PlotType}(args[:plot]),
        args[:plotting_file_path],
        Set{String}(args[:file_formats]),
        args[:width],
        args[:page_entropy],
        args[:px_per_unit],
        args[:plot_eigval_vs_cbe],
        args[:plot_fragment_sizes],
        args[:show]
    )
end
