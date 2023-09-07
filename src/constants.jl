abstract type MeasurementType end

abstract type HeatmapContinuous <: MeasurementType end
struct ExpectationValue <: HeatmapContinuous end
struct SingleSiteEntropy <: HeatmapContinuous end

abstract type HeatmapDiscrete <: MeasurementType end
struct Rounded <: HeatmapDiscrete end
struct BondDimensions <: HeatmapDiscrete end

abstract type LinePlot <: MeasurementType end
struct CenterBipartiteEntropy <: LinePlot end

const PlotLabels = Dict{MeasurementType,String}(
    ExpectationValue() => "Expectation\nValue",
    SingleSiteEntropy() => "Single Site\nEntropy",
    Rounded() => "Rounded",
    BondDimensions() => "Bond\nDimensions",
    CenterBipartiteEntropy() => "Center\nBipartite\nEntropy"
)