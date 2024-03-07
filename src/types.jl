abstract type MeasurementType end

abstract type HeatmapContinuous <: MeasurementType end
struct ExpectationValue <: HeatmapContinuous end
name(_::ExpectationValue) = "ExpectationValue"
label(_::ExpectationValue) = "Expectation\nValue"
struct SingleSiteEntropy <: HeatmapContinuous end
name(_::SingleSiteEntropy) = "SingleSiteEntropy"
label(_::SingleSiteEntropy) = "Single Site\nEntropy"

abstract type HeatmapDiscrete <: MeasurementType end
struct Rounded <: HeatmapDiscrete end
name(_::Rounded) = "Rounded"
label(_::Rounded) = "Rounded"
struct BondDimensions <: HeatmapDiscrete end
name(_::BondDimensions) = "BondDimensions"
label(_::BondDimensions) = "Bond\nDimensions"

abstract type LinePlot <: MeasurementType end
struct CenterBipartiteEntropy <: LinePlot end
name(_::CenterBipartiteEntropy) = "CenterBipartiteEntropy"
label(_::CenterBipartiteEntropy) = "Center\nBipartite\nEntropy"
