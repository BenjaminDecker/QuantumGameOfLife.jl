"""
    evolve_tebd(distance::Int, lower_activation_bound::Int, upper_activation_bound::Int, periodic::Bool, psi_0_mps::MPS, num_steps::Int, step_size::Float64)::Vector{MPS}

Calculate a time evolution of `psi_0_mps` with the TEBD algorithm.

Inspired by code from [ITensors.jl](https://itensor.github.io/ITensors.jl/dev/tutorials/MPSTimeEvolution.html)
"""
function evolve_tebd(distance::Int, lower_activation_bound::Int, upper_activation_bound::Int, periodic::Bool, psi_0_mps::MPS, num_steps::Int, step_size::Float64, sweeps::Int, max_bond_dim::Int)::Vector{MPS}
    site_inds = siteinds(psi_0_mps)
    trotter_gates = ITensor[]

    for i in eachindex(site_inds)
        summands = ITensor[]
        for (ket_1_indices, ket_0_indices) in get_combinations(length(site_inds), i, distance, lower_activation_bound, upper_activation_bound, periodic)
            opj = op("X", site_inds[i])
            for ket_1_index in ket_1_indices
                opj *= op("Proj1", site_inds[ket_1_index])
            end
            for ket_0_index in ket_0_indices
                opj *= op("Proj0", site_inds[ket_0_index])
            end
            # op1 = mapreduce(ket_1_index -> op("Proj1", site_inds[ket_1_index]), *, ket_1_indices; init=1.0)
            # op0 = mapreduce(ket_0_index -> op("Proj0", site_inds[ket_0_index]), *, ket_0_indices; init=1.0)
            push!(summands, opj)
        end
        hj = reduce(+, summands)
        uj = let
            t = step_size * pi / (4 * sweeps)
            exp(-im * t * hj)
        end
        push!(trotter_gates, uj)
    end

    append!(trotter_gates, reverse(trotter_gates))

    results = [psi_0_mps]
    sizehint!(results, num_steps)
    psi_mps = psi_0_mps

    cutoff = 1E-8

    @showprogress "Calculating Time Evolution" for _ in 2:num_steps
        for _ in 1:sweeps
            psi_mps = apply(trotter_gates, psi_mps; maxdim=max_bond_dim, cutoff=0.0)
            normalize!(psi_mps)
        end
        push!(results, psi_mps)
    end
    return results
end
