using ITensors
using ProgressMeter

measure(state::MPS, measurement_type::ExpectationValue) = expect(state, "Proj1")
measure(state::MPS, measurement_type::Rounded) = round.(measure(state, ExpectationValue()))
measure(state::MPS, measurement_type::BondDimensions) = [dim(linkind(state, i)) for i in 1:(length(state)-1)]

function measure(state::MPS, measurement_type::SingleSiteEntropy)
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

function measure(state::MPS, measurement_type::CenterBipartiteEntropy)
    site_inds = siteinds(state)
    Vector{Float64}([center_bipartite_entropy(state)])
end

function measure(states::Vector{MPS}, measurements_types::Set{MeasurementType})
    measurements = Dict{MeasurementType,Vector{Vector{Float64}}}()
    @showprogress "Measuring..." for state in states
        for measurement_type in measurements_types
            measurement = measure(state, measurement_type)
            result_vector = get!(measurements, measurement_type, [])
            push!(result_vector, measurement)
        end
    end
    return measurements
end