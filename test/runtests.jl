using Test
using QuantumGameOfLife
using ITensors
using ITensorMPS

const QGL = QuantumGameOfLife

"""
Build an `Args` object by running the real command line parser on `flags`.
This exercises argument parsing, the initial-state registry and site-index creation.
"""
function parse_args_for_test(flags::Vector{String})::QGL.Args
    saved = copy(ARGS)
    try
        empty!(ARGS)
        append!(ARGS, flags)
        return QGL.parse_commandline()
    finally
        empty!(ARGS)
        append!(ARGS, saved)
    end
end

@testset "QuantumGameOfLife.jl" begin

    @testset "necklace" begin
        @test QGL.necklace(BitVector([1, 0, 0])) == BitVector([0, 0, 1])
        @test QGL.necklace(BitVector([1, 0, 1])) == BitVector([0, 1, 1])
        # rotating an already minimal representative is a no-op
        @test QGL.necklace(BitVector([0, 0, 0])) == BitVector([0, 0, 0])
        # all elements of a rotation class map to the same necklace
        @test QGL.necklace(BitVector([0, 1, 1, 0])) == QGL.necklace(BitVector([1, 1, 0, 0]))
    end

    @testset "bitvector_to_mps" begin
        bits = BitVector([1, 0, 1, 1])
        psi = QGL.bitvector_to_mps(bits, siteinds("Qubit", 4))
        @test expect(psi, "Proj1") ≈ Float64.(bits)
    end

    @testset "bipartite entropy" begin
        # A product state has zero entanglement across any cut.
        product_state = QGL.bitvector_to_mps(BitVector([1, 0, 1, 0]), siteinds("Qubit", 4))
        @test QGL.bipartite_entropy(product_state, 2) ≈ 0.0 atol = 1e-10
        @test QGL.center_bipartite_entropy(product_state) ≈ 0.0 atol = 1e-10
    end

    @testset "configuration_id" begin
        args = parse_args_for_test(["--num-cells", "3"])  # distance=1, open boundaries
        # Neighbors (1, 2, 3) of cell 2 are all "alive" -> binary 111 == 7
        @test QGL.configuration_id(fill(true, 3), 2, args) == 0b111
        # Only the cell itself is "alive" -> binary 010 == 2
        @test QGL.configuration_id([false, true, false], 2, args) == 0b010
    end

    @testset "initial state registry" begin
        # The command line choices are derived from the registry keys.
        @test Set(QGL.INITIAL_STATE_CHOICES) == Set(keys(QGL.INITIAL_STATE_REGISTRY))
        site_inds = siteinds("Qubit", 6)
        for (name, make_state) in QGL.INITIAL_STATE_REGISTRY
            psi = make_state(site_inds)
            @test psi isa MPS
            @test length(psi) == 6
        end
    end

    @testset "hamiltonian construction" begin
        args = parse_args_for_test(["--num-cells", "4"])
        H = QGL.build_hamiltonian_mpo(args.site_inds, args)
        @test H isa MPO
        @test length(firstsiteinds(H)) == 4
        # For the default rule the Hamiltonian is hermitian, so the evolution is unitary.
        @test norm(H - dag(swapprime(H, 0, 1))) < 1e-5
    end

    @testset "exact time evolution" begin
        args = parse_args_for_test([
            "--num-cells", "3", "--num-steps", "2", "--algorithm", "exact",
            "--initial-states", "blinker",
        ])
        H = QGL.build_hamiltonian_mpo(args.site_inds, args)
        results = QGL.evolve(args.algorithm, args.initial_states, H, args)
        @test results isa Vector{Vector{MPS}}
        @test length(results) == length(args.initial_states)
        trajectory = only(results)
        @test length(trajectory) == args.num_steps + 1  # includes psi_0
        for psi in trajectory
            @test norm(psi) ≈ 1.0 atol = 1e-8
        end
    end

    @testset "measuring" begin
        args = parse_args_for_test([
            "--num-cells", "4", "--num-steps", "2", "--initial-states", "blinker",
        ])
        psi = QGL.bitvector_to_mps(BitVector([0, 1, 0, 0]), args.site_inds)

        expectation = QGL.measure(QGL.ExpectationValue(), psi)
        @test expectation ≈ [0.0, 1.0, 0.0, 0.0]

        # A single state measured over the classical rule stays deterministic.
        classical = QGL.measure(QGL.Classical(), [psi], args)
        @test length(classical) == args.num_steps + 1
        @test all(step -> all(x -> x in (0.0, 1.0), step), classical)
    end

    @testset "TEBD is not implemented" begin
        args = parse_args_for_test(["--num-cells", "3", "--algorithm", "tebd"])
        H = QGL.build_hamiltonian_mpo(args.site_inds, args)
        @test_throws ErrorException QGL.evolve(args.algorithm, args.initial_states, H, args)
    end

end
