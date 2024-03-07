using ITensors

"""
    build_MPO_hamiltonian(
    site_inds::Vector{ITensors.Index{Int64}},
    distance::Int,
    activation_interval::UnitRange{Int64},
    periodic::Bool
)::MPO

TBW
"""
function build_MPO_hamiltonian(
    site_inds::Vector{ITensors.Index{Int64}},
    distance::Int,
    activation_interval::UnitRange{Int64},
    periodic::Bool
)::MPO
    num_cells = length(site_inds)
    os = OpSum()
    for site_index in 1:num_cells
        for (ket_1_indices, ket_0_indices) in get_combinations(num_cells, site_index, distance, activation_interval, periodic)
            op::Vector{Union{Int,String}} = []
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
