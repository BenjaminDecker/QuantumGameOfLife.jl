module HamiltonianMpoCreator

using ITensors
using Combinatorics

export build_hamiltonian_MPO

function get_combinations(num_cells::Int, index::Int, distance::Int, lower_bound::Int, upper_bound::Int, periodic::Bool)
    valid_indices::Array{Int,1} = []
    append!(valid_indices, (max(1, index - distance):(index-1)))
    append!(valid_indices, (index+1):(min(num_cells, index + distance)))
    if periodic
        append!(valid_indices, ((num_cells+index-distance):num_cells))
        append!(valid_indices, (1:(index-num_cells+distance)))
    end
    ket_1_indices::Array{Array{Int,1},1} = []
    for num_ket_1 in lower_bound:(upper_bound-1)
        append!(ket_1_indices, combinations(valid_indices, num_ket_1))
    end
    return map(combination -> (combination, setdiff(valid_indices, combination)), ket_1_indices)
end

function build_hamiltonian_MPO(site_inds::Vector{ITensors.Index{Int64}}, distance::Int, lower_activation_bound::Int, upper_activation_bound::Int, periodic::Bool)::MPO
    num_cells = length(site_inds)
    os = OpSum()
    for site_index in 1:num_cells
        for (ket_1_indices, ket_0_indices) in get_combinations(num_cells, site_index, distance, lower_activation_bound, upper_activation_bound, periodic)
            op::Array{Union{Int,String},1} = []
            push!(op, "X")
            push!(op, site_index)
            for ket_1_index in ket_1_indices
                push!(op, "Proj1")
                push!(op, ket_1_index)
            end
            for ket_0_index in ket_0_indices
                push!(op, "Proj0")
                push!(op, ket_0_index)
            end
            os += Tuple(op)
        end
    end
    return MPO(os, site_inds)
end
end