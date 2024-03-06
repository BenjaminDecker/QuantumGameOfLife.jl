module Utils

using ITensors
using Combinatorics

function bipartite_entropy(psi::MPS, seperator_index::Int)::Float64
    orthogonalize!(psi, seperator_index)
    _, S, _ = svd(
        psi[seperator_index],
        (linkind(psi, seperator_index - 1), siteind(psi, seperator_index))
    )
    SvN = 0.0
    for n = 1:dim(S, 1)
        p = S[n, n]^2
        if p > 0.0
            SvN -= p * log(p)
        end
    end
    real(SvN)
end

function center_bipartite_entropy(psi::MPS)::Float64
    center = trunc(Int, length(psi) / 2)
    bipartite_entropy(psi, center)
end


function get_combinations(num_cells::Int, index::Int, distance::Int, activation_interval::UnitRange{Int64}, periodic::Bool)::Vector{Tuple{Vector{Int},Vector{Int}}}
    valid_indices::Vector{Int} = []
    append!(valid_indices, (max(1, index - distance):(index-1)))
    append!(valid_indices, (index+1):(min(num_cells, index + distance)))
    if periodic
        append!(valid_indices, ((num_cells+index-distance):num_cells))
        append!(valid_indices, (1:(index-num_cells+distance)))
    end
    ket_1_indices::Vector{Vector{Int}} = []
    for num_ket_1 in activation_interval
        append!(ket_1_indices, combinations(valid_indices, num_ket_1))
    end
    return map(combination -> (combination, setdiff(valid_indices, combination)), ket_1_indices)
end
end
