abstract type PlotType end

abstract type HeatmapContinuous <: PlotType end
struct ExpectationValue <: HeatmapContinuous end
name(_::ExpectationValue) = "Expectation Value"
label(_::ExpectationValue) = "Expectation\nValue"
ordering_index(_::ExpectationValue) = 2
struct SingleSiteEntropy <: HeatmapContinuous end
name(_::SingleSiteEntropy) = "Single Site Entropy"
label(_::SingleSiteEntropy) = "Single Site\nEntropy"
ordering_index(_::SingleSiteEntropy) = 3

abstract type HeatmapDiscrete <: PlotType end
struct Rounded <: HeatmapDiscrete end
name(_::Rounded) = "Rounded"
label(_::Rounded) = "Rounded"
ordering_index(_::Rounded) = 4
struct BondDimensions <: HeatmapDiscrete end
name(_::BondDimensions) = "Bond Dimensions"
label(_::BondDimensions) = "Bond\nDimensions"
ordering_index(_::BondDimensions) = 5
struct Classical <: HeatmapDiscrete end
name(_::Classical) = "Classical"
label(_::Classical) = "Classical"
ordering_index(_::Classical) = 1

abstract type LinePlot <: PlotType end
struct CenterBipartiteEntropy <: LinePlot end
name(_::CenterBipartiteEntropy) = "Center Bipartite Entropy"
label(_::CenterBipartiteEntropy) = "Center\nBipartite\nEntropy"
ordering_index(_::CenterBipartiteEntropy) = 6

Base.isless(a::PlotType, b::PlotType) = ordering_index(a) < ordering_index(b)

struct Args
    num_steps::Int
    distance::Int
    activation_interval::UnitRange
    periodic::Bool
    step_size::Float64
    initial_state::Vector{String}
    algorithm::String
    num_cells::Int
    sweeps_per_time_step::Int
    max_bond_dim::Int
    periodic_boundaries::Bool
    plots::Set{PlotType}
    plotting_file_path::String
    file_formats::Set{String}
    plot_eigval_vs_cbe::Bool
    plot_fragment_sizes::Bool
    show::Bool
end

function Args(args::Dict{Symbol,Any})
    return Args(
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
        args[:periodic_boundaries],
        Set{PlotType}(args[:plot]),
        args[:plotting_file_path],
        Set{String}(args[:file_formats]),
        args[:plot_eigval_vs_cbe],
        args[:plot_fragment_sizes],
        args[:show]
    )
end
