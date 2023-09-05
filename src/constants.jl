abstract type MeasurementType end

abstract type Continuous <: MeasurementType end
struct ExpectationValue <: Continuous end
struct SingleSiteEntropy <: Continuous end

abstract type Discrete <: MeasurementType end
struct Rounded <: Discrete end
struct BondDimensions <: Discrete end

const PlotLabels = Dict{MeasurementType,String}(
    ExpectationValue() => "Expectation\nValue",
    SingleSiteEntropy() => "Single Site\nEntropy",
    Rounded() => "Rounded",
    BondDimensions() => "Bond\nDimensions"
)