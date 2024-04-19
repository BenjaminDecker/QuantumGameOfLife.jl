measure(::ExpectationValue, state::MPS)::Vector{Float64} =
    expect(state, "Proj1")

measure(::Rounded, state::MPS)::Vector{Float64} =
    round.(measure(ExpectationValue(), state))

measure(::BondDimensions, state::MPS)::Vector{Float64} =
    [dim(linkind(state, i)) for i in 1:(length(state)-1)]

measure(::CenterBipartiteEntropy, state::MPS)::Vector{Float64} =
    [center_bipartite_entropy(state)]

function measure(::SingleSiteEntropy, state::MPS)::Vector{Float64}
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
    ::Classical,
    states::Vector{MPS},
    args::Args
)::Vector{Vector{Float64}}
    states = Vector{Vector{Bool}}([measure(Rounded(), states[1])])
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
    plot_type::PlotType,
    states::Vector{MPS},
    ::Args
)::Vector{Vector{Float64}}
    return @showprogress desc = "Measuring $(name(plot_type))" map(state -> measure(plot_type, state), states)
end

function measure(
    states::Vector{MPS},
    args::Args
)::Dict{PlotType,Vector{Vector{Float64}}}
    measurements = Dict{PlotType,Vector{Vector{Float64}}}()
    for plot_type in args.plots
        measurements[plot_type] = measure(plot_type, states, args)
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
