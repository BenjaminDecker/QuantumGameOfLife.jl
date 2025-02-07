function evolve(::Sierpinski, psi_0_vec::Vector{MPS}, ::MPO, args::Args)::Vector{Vector{MPS}}
    site_inds = siteinds(psi_0_vec[1])
    trotter_gates = ITensor[]
    for i in 1:(length(site_inds)-1)
        h_i = op("Proj1", site_inds[i]) * op("X", site_inds[i+1])
        t = args.step_size * pi / (2 * args.sweeps_per_time_step)
        push!(trotter_gates, exp(-im * t * h_i))
    end
    trotter_gates = reverse(trotter_gates)
    # append!(trotter_gates, reverse(trotter_gates))

    results_vec = []

    for psi_0 in psi_0_vec
        results = [psi_0]
        sizehint!(results, args.num_steps)
        psi_mps = psi_0
        @showprogress "Calculating Time Evolution" for _ in 2:args.num_steps

            for _ in 1:args.sweeps_per_time_step
                psi_mps = apply(trotter_gates, psi_mps; maxdim=args.max_bond_dim, cutoff=args.svd_epsilon)
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
        push!(results_vec, results)
    end
    return results_vec
end
