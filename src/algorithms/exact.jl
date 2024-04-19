function evolve(::Exact, H_mpo::MPO, psi_0_mps::MPS, args::Args)::Vector{MPS}
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
