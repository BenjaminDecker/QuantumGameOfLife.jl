abstract type PlotType end

abstract type HeatmapContinuous <: PlotType end
struct ExpectationValue <: HeatmapContinuous end
struct SingleSiteEntropy <: HeatmapContinuous end

abstract type HeatmapDiscrete <: PlotType end
struct Rounded <: HeatmapDiscrete end
struct BondDimensions <: HeatmapDiscrete end
struct Classical <: HeatmapDiscrete end

abstract type LinePlot <: PlotType end
struct CenterBipartiteEntropy <: LinePlot end
struct Autocorrelation <: LinePlot end

"""
    PlotSpec

Everything the rest of the package needs to know about a `PlotType` that is not encoded in
its dispatch behavior: the file-name identifier (`id`), the human-readable title (`name`) and
axis label (`label`), and the alternative spellings accepted on the command line and in a
config file (`aliases`).
"""
struct PlotSpec
    id::String
    name::String
    label::String
    aliases::Vector{String}
end

# Single source of truth for the available plot types. The order of this vector is the order
# plots are drawn in (`Base.isless`); the metadata accessors, `PLOTS_CHOICES`, the help
# strings and command line / config parsing are all derived from it. Add a plot type by
# adding one row here and a matching `measure` method (and a `render_row!` method in
# plotting.jl if it needs a new kind of axis).
const PLOT_TYPES = [
    Classical => PlotSpec(
        "classical", "Classical", "Classical",
        ["classic"]),
    ExpectationValue => PlotSpec(
        "expect", "Expectation Value", "Expectation\nValue",
        ["expectation", "expectation_value", "expectation-value"]),
    SingleSiteEntropy => PlotSpec(
        "sse", "Single Site Entropy", "Single Site\nEntropy",
        ["single_site_entropy", "single-site-entropy"]),
    Rounded => PlotSpec(
        "rounded", "Rounded", "Rounded",
        ["round"]),
    BondDimensions => PlotSpec(
        "bond_dims", "Bond Dimensions", "Bond\nDimension",
        ["bond_dim", "bond_dimension", "bond_dimensions",
         "bond-dim", "bond-dims", "bond-dimension", "bond-dimensions"]),
    CenterBipartiteEntropy => PlotSpec(
        "cbe", "Center Bipartite Entropy", "Center\nBipartite\nEntropy",
        ["center_bipartite_entropy", "center-bipartite-entropy"]),
    Autocorrelation => PlotSpec(
        "autocorrelation", "Autocorrelation", "Autocorrelation",
        String[]),
]

const _PLOT_INDEX = Dict{DataType,Int}(T => i for (i, (T, _)) in enumerate(PLOT_TYPES))
_plot_spec(p::PlotType) = PLOT_TYPES[_PLOT_INDEX[typeof(p)]].second

filename_identifier(p::PlotType) = _plot_spec(p).id
name(p::PlotType) = _plot_spec(p).name
label(p::PlotType) = _plot_spec(p).label
ordering_index(p::PlotType) = _PLOT_INDEX[typeof(p)]

# Maps every accepted spelling (canonical id or alias) to the corresponding PlotType.
const PLOT_LOOKUP = Dict{String,PlotType}(
    alias => T() for (T, spec) in PLOT_TYPES for alias in [spec.id; spec.aliases])

const PLOTS_CHOICES = [spec.id for (_, spec) in PLOT_TYPES]

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

"""
    AlgorithmSpec

The command line / config file spelling(s) of an `Algorithm` (`id` and `aliases`) and whether
it is `advertised`, i.e. listed in the help and choices.
"""
struct AlgorithmSpec
    id::String
    aliases::Vector{String}
    advertised::Bool
end

# Single source of truth for the algorithms that can be named; `parse_algorithm` and
# `ALGORITHM_CHOICES` are derived from it. TEBD parses but is not advertised because it is
# not implemented yet (see `algorithms/tebd.jl`).
const ALGORITHMS = [
    Exact()      => AlgorithmSpec("exact", String[], true),
    TDVP1()      => AlgorithmSpec("tdvp1", String[], true),
    TDVP2()      => AlgorithmSpec("tdvp2", String[], true),
    Sierpinski() => AlgorithmSpec("sierpinski", ["sierpiński"], true),
    TEBD()       => AlgorithmSpec("tebd", String[], false),
]

const ALGORITHM_LOOKUP = Dict{String,Algorithm}(
    alias => a for (a, spec) in ALGORITHMS for alias in [spec.id; spec.aliases])

const ALGORITHM_CHOICES = [spec.id for (_, spec) in ALGORITHMS if spec.advertised]

"""
    Args

The parsed and derived simulation configuration. Besides the values coming from the command
line or a TOML config file, it also holds derived data such as the site indices
(`site_inds`) and the constructed initial state `MPS` objects (`initial_states`). This is
the single object that is threaded through Hamiltonian construction, time evolution,
measuring and plotting.
"""
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
        [INITIAL_STATE_REGISTRY[state_name](site_inds) for state_name in args[:initial_states]],
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
