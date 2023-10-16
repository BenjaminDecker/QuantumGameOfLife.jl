module TimeEvolution

using ITensors
using ProgressMeter
using FromFile

@from "Types.jl" using Types

export evolve

function evolve_exact(H_mpo::MPO, psi_0_mps::MPS, num_steps::Int, step_size::Float64)
    U_tensor = let
        print("Calculating Time Evolution Operator...")
        t = step_size * pi / 2
        exp(-im * contract(H_mpo) * t)
    end
    println("done")
    results = [psi_0_mps]
    sizehint!(results, num_steps)
    psi_tensor = contract(psi_0_mps)
    @showprogress "Calculating Time Evolution" for _ in 2:num_steps
        psi_tensor = noprime(U_tensor * psi_tensor)
        push!(results, MPS(psi_tensor, siteinds(psi_0_mps)))
    end
    return results
end

function evolve(H_mpo::MPO, psi_0_mps::MPS, num_steps::Int, step_size::Float64, algorithm::Types.Algorithm)
    return evolve_exact(H_mpo, psi_0_mps, num_steps, step_size)
end
end