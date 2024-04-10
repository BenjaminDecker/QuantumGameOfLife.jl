using ITensors
using ProgressMeter

measure(state::MPS, _::ExpectationValue)::Vector{Float64} =
    expect(state, "Proj1")

measure(state::MPS, _::Rounded)::Vector{Float64} =
    round.(measure(state, ExpectationValue()))

measure(state::MPS, _::BondDimensions)::Vector{Float64} =
    [dim(linkind(state, i)) for i in 1:(length(state)-1)]

measure(state::MPS, _::CenterBipartiteEntropy)::Vector{Float64} =
    [center_bipartite_entropy(state)]

function measure(state::MPS, _::SingleSiteEntropy)::Vector{Float64}
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

function measure(
    states::Vector{MPS},
    _::Args,
    plot::PlotType
)::Vector{Vector{Float64}}
    return @showprogress desc = "Measuring $(name(plot))" map(state -> measure(state, plot), states)
end

function measure(
    states::Vector{MPS},
    args::Args,
    _::Classical
)::Vector{Vector{Float64}}
    states = Vector{Vector{Bool}}([measure(states[1], Rounded())])
    for step in 2:args.num_steps
        last_state = states[step-1]
        next_state = fill(false, length(states[1]))
        for index in eachindex(next_state)
            alive_neighbors = 0
            for offset in 1:args.distance
                alive_neighbors += get_value(
                    last_state,
                    index - offset,
                    args.periodic
                )
                alive_neighbors += get_value(
                    last_state,
                    index + offset,
                    args.periodic
                )
            end
            if alive_neighbors in args.activation_interval
                next_state[index] = !last_state[index]
            end
        end
        push!(states, next_state)
    end
    return Vector{Vector{Float64}}(states)
end

function measure(
    states::Vector{MPS},
    args::Args
)::Dict{PlotType,Vector{Vector{Float64}}}
    measurements = Dict{PlotType,Vector{Vector{Float64}}}()
    for plot in args.plots
        measurements[plot] = measure(states, args, plot)
    end
    return measurements
end

function get_value(state::Vector{Bool}, index::Int, periodic::Bool)::Int
    if !(index in eachindex(state)) && !periodic
        return 0
    end
    num_cells = length(state)
    if index < 1
        index += num_cells
    elseif index > num_cells
        index -= num_cells
    end
    return Int(state[index])
end
