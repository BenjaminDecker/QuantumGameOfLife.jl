using ITensors
using ProgressMeter

function evolve_exact(H_mpo::MPO, psi_0_mps::MPS, args::Args)::Vector{MPS}
    U_tensor = let
        print("Calculating Time Evolution Operator...")
        t = args.step_size * pi / 2
        exp(-im * contract(H_mpo) * t)
    end
    println("done")
    results = [psi_0_mps]
    sizehint!(results, args.num_steps)
    psi_tensor = contract(psi_0_mps)
    @showprogress "Calculating Time Evolution" for _ in 2:args.num_steps
        psi_tensor = noprime(U_tensor * psi_tensor)
        push!(results, MPS(psi_tensor, siteinds(psi_0_mps)))
    end
    return results
end

function evolve_serpinsky(psi_0_mps::MPS, args::Args)::Vector{MPS}
    site_inds = siteinds(psi_0_mps)
    trotter_gates = ITensor[]
    for i in 1:(length(site_inds)-1)
        h_i = op("Proj1", site_inds[i]) * op("X", site_inds[i+1])
        t = args.step_size * pi / (2 * args.sweeps_per_time_step)
        push!(trotter_gates, exp(-im * t * h_i))
    end

    results = [psi_0_mps]
    sizehint!(results, args.num_steps)
    psi_mps = psi_0_mps

    cutoff = 1E-8
    trotter_gates = reverse(trotter_gates)
    # append!(trotter_gates, reverse(trotter_gates))

    @showprogress "Calculating Time Evolution" for step in 2:args.num_steps

        for _ in 1:args.sweeps_per_time_step
            psi_mps = apply(trotter_gates, psi_mps; maxdim=args.max_bond_dim, cutoff=0.0)
            normalize!(psi_mps)
        end
        # if (step + 7) % 16 <= 7
        # measurements = expect(psi_mps, "Proj1")
        # collapse = rand(Float64) < measurements[9] ? "Proj1" : "Proj0"
        # newA = psi_mps[9] * op(collapse, site_inds[9])
        # noprime!(newA)
        # psi_mps[9] = newA
        # normalize!(psi_mps)
        # end
        push!(results, psi_mps)
        # for trotter_gate in trotter_gates
        #     psi_mps = apply(trotter_gate, psi_mps)
        #     normalize!(psi_mps)
        #     push!(results, psi_mps)
        # end
    end
    return results
end

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

# function evolve(H_mpo::MPO, psi_0_mps::MPS, num_steps::Int, step_size::Float64, algorithm::Types.Algorithm)
# return evolve_exact(H_mpo, psi_0_mps, num_steps, step_size)
# end
