using ITensors
using ProgressMeter

measure(state::MPS, measurement_type::ExpectationValue)::Vector{Float64} = expect(state, "Proj1")
measure(state::MPS, measurement_type::Rounded)::Vector{Float64} = round.(measure(state, ExpectationValue()))
measure(state::MPS, measurement_type::BondDimensions)::Vector{Float64} = [dim(linkind(state, i)) for i in 1:(length(state)-1)]

function measure(state::MPS, measurement_type::SingleSiteEntropy)::Vector{Float64}
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

function measure(state::MPS, measurement_type::CenterBipartiteEntropy)::Vector{Float64}
    Vector{Float64}([center_bipartite_entropy(state)])
end

function measure(states::Vector{MPS}, measurements_types::Set{MeasurementType})::Dict{MeasurementType,Vector{Vector{Float64}}}
    measurements = Dict{MeasurementType,Vector{Vector{Float64}}}()
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
