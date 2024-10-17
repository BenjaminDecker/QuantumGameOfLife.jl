function evolve(::Exact, psi_0_vec::Vector{MPS}, H::MPO, args::Args)::Vector{Vector{MPS}}
    U_tensor = let
        print("Calculating Time Evolution Operator...")
        t = args.step_size * pi / 2
        exp(-im * contract(H) * t)
    end
    println("done")

    results_vec = []
    for psi_0 in psi_0_vec
        results = [psi_0]
        sizehint!(results, args.num_steps)
        psi_tensor = contract(psi_0)

        @showprogress "Calculating Time Evolution" for _ in 2:args.num_steps
            psi_tensor = noprime(U_tensor * psi_tensor)
            normalize!(psi_tensor)
            push!(results, MPS(psi_tensor, siteinds(psi_0)))
        end
        push!(results_vec, results)
    end
    return results_vec
end
