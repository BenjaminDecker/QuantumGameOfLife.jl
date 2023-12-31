module Measure

using ITensors
using ProgressMeter
using FromFile

@from "Types.jl" using Types
@from "Utils.jl" using Utils

export measure

measure(state::MPS, measurement_type::Types.ExpectationValue)::Vector{Float64} = expect(state, "Proj1")
measure(state::MPS, measurement_type::Types.Rounded)::Vector{Float64} = round.(measure(state, Types.ExpectationValue()))
measure(state::MPS, measurement_type::Types.BondDimensions)::Vector{Float64} = [dim(linkind(state, i)) for i in 1:(length(state)-1)]

function measure(state::MPS, measurement_type::Types.SingleSiteEntropy)::Vector{Float64}
    site_inds = siteinds(state)
    num_sites = length(site_inds)
    sse = Vector{Float64}(undef, num_sites)
    for i in 1:num_sites
        orthogonalize!(state, i)
        ket = state[i]
        bra = conj(prime(ket, tags="Site"))
        ket_bra = ket * bra
        density_matrix = Matrix(ket_bra, inds(ket, tags="Site"), inds(bra, tags="Site"))
        result = real(-tr(density_matrix * log(density_matrix) / log(2)))
        sse[i] = isnan(result) ? zero(result) : result
    end
    return sse
end

function measure(state::MPS, measurement_type::Types.CenterBipartiteEntropy)::Vector{Float64}
    Vector{Float64}([Utils.center_bipartite_entropy(state)])
end

function measure(states::Vector{MPS}, measurements_types::Set{Types.MeasurementType})::Dict{Types.MeasurementType,Vector{Vector{Float64}}}
    measurements = Dict{Types.MeasurementType,Vector{Vector{Float64}}}()
    p = Progress(length(states) * length(measurements_types); desc="Measuring")
    for measurement_type in measurements_types
        result_vector = get!(measurements, measurement_type, [])
        sizehint!(result_vector, length(states))
        for state in states
            measurement = measure(state, measurement_type)
            push!(result_vector, measurement)
            next!(p)
        end
    end
    finish!(p)
    return measurements
end
end