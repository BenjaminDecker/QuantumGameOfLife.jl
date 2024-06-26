abstract type PlotType end

abstract type HeatmapContinuous <: PlotType end
struct ExpectationValue <: HeatmapContinuous end
name(::ExpectationValue) = "Expectation Value"
label(::ExpectationValue) = "Expectation\nValue"
ordering_index(::ExpectationValue) = 2
struct SingleSiteEntropy <: HeatmapContinuous end
name(::SingleSiteEntropy) = "Single Site Entropy"
label(::SingleSiteEntropy) = "Single Site\nEntropy"
ordering_index(::SingleSiteEntropy) = 3

abstract type HeatmapDiscrete <: PlotType end
struct Rounded <: HeatmapDiscrete end
name(::Rounded) = "Rounded"
label(::Rounded) = "Rounded"
ordering_index(::Rounded) = 4
struct BondDimensions <: HeatmapDiscrete end
name(::BondDimensions) = "Bond Dimensions"
label(::BondDimensions) = "Bond\nDimensions"
ordering_index(::BondDimensions) = 5
struct Classical <: HeatmapDiscrete end
name(::Classical) = "Classical"
label(::Classical) = "Classical"
ordering_index(::Classical) = 1

abstract type LinePlot <: PlotType end
struct CenterBipartiteEntropy <: LinePlot end
name(::CenterBipartiteEntropy) = "Center Bipartite Entropy"
label(::CenterBipartiteEntropy) = "Center\nBipartite\nEntropy"
ordering_index(::CenterBipartiteEntropy) = 6

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

struct Args
    num_steps::Int
    distance::Int
    activation_interval::UnitRange
    periodic::Bool
    step_size::Float64
    initial_state::Vector{String}
    algorithm::Algorithm
    num_cells::Int
    sweeps_per_time_step::Int
    max_bond_dim::Int
    svd_epsilon::Float64
    periodic_boundaries::Bool
    plots::Set{PlotType}
    plotting_file_path::String
    file_formats::Set{String}
    dpi::Float64
    plot_eigval_vs_cbe::Bool
    plot_fragment_sizes::Bool
    show::Bool
end

Args(args::Dict{Symbol,Any})::Args = Args(
    args[:num_steps],
    args[:distance],
    (args[:activation_interval][1]):(args[:activation_interval][2]),
    args[:periodic_boundaries],
    args[:step_size],
    args[:initial_state],
    args[:algorithm],
    args[:num_cells],
    args[:sweeps_per_time_step],
    args[:max_bond_dim],
    args[:svd_epsilon],
    args[:periodic_boundaries],
    Set{PlotType}(args[:plot]),
    args[:plotting_file_path],
    Set{String}(args[:file_formats]),
    args[:dpi],
    args[:plot_eigval_vs_cbe],
    args[:plot_fragment_sizes],
    args[:show]
)
