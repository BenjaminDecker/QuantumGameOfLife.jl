using ITensors
using ProgressMeter

"""
    eigval_vs_cbe(H::MPO)::Tuple{Vector{Float64},Vector{Float64}}

For hamiltonian `H`, calculate its eigenvalues divided by the system size and their
respective center bipartite entropies.

# Returns
- `Tuple{Vector{Float64},Vector{Float64}}`: Tuple of eigenvalues and center bipartite
 entropies (eigvals, cbes)

"""
function eigval_vs_cbe(H::MPO)::Tuple{Vector{Float64},Vector{Float64}}
    site_inds = firstsiteinds(H)
    print("Calculating eigen decomposition...")
    D, U = eigen(contract(H); ishermitian=true)
    println("done")

    eigen_ind = inds(U, tags="eigen")[1]
    eigvals = [D[eigen_ind=>i, eigen_ind'=>i] for i in eachval(eigen_ind)]
    cbe = @showprogress "Calculating center bipartite entropy of eigenvectors" [
        center_bipartite_entropy(MPS(U * onehot(eigen_ind => i), site_inds))
        for i in eachval(eigen_ind)
    ]

    (eigvals ./ length(site_inds), cbe)
end


"""
    necklace(bits::BitVector)::BitVector

Calculate the numerically smallest bitstring from the binary necklace equivalence class of
`bits`.

A binary necklace is an equivalence class of bitstrings of fixed length that are equivalent
up to some rotation.
https://en.wikipedia.org/wiki/Necklace_(combinatorics)

# Returns
- `BitVector`: The numerically smallest bitstring of length `length(bits)` in the necklace
 of `bits`

# Examples
```julia
julia> necklace(BitVector([1, 0, 0]))
3-element BitVector:
 0
 0
 1

julia> necklace(BitVector([1, 0, 1]))
3-element BitVector:
 0
 1
 1
```
"""
function necklace(bits::BitVector)::BitVector
    best_candidate = bits
    for i in eachindex(bits)
        candidate = circshift(bits, i)
        if candidate < best_candidate
            best_candidate = candidate
        end
    end
    return best_candidate
end



"""
    bitvector_to_mps(bitvector::BitVector, site_inds::Vector{Index{Int64}})::MPS

Create a computational basis state in MPS form from `bitvector`.

The state is equivalent to a product state of Qubits |0> and |1>, corresponding to the
entries of `bitvector`.

# Examples
```julia
julia> mps1 = bitvector_to_mps(BitVector([1, 0, 0]), siteinds("Qubit", 3));
julia> expect(mps1, "Proj1")
3-element Vector{Float64}:
 1.0
 0.0
 0.0

julia> mps2 = bitvector_to_mps(BitVector([0, 1, 1, 0]), siteinds("Qubit", 4));
julia> expect(mps2, "Proj1")
4-element Vector{Float64}:
 0.0
 1.0
 1.0
 0.0
```
"""
function bitvector_to_mps(bitvector::BitVector, site_inds::Vector{Index{Int64}})::MPS
    return MPS(
        site_inds,
        map(string, filter(x -> x == '1' || x == '0', collect(bitstring(bitvector))))
    )
end


"""
    fragment_sizes(H::MPO, periodic::Bool)::Vector{Int}

Calculate the sizes of the Hilbert space fragments of `H`.

A fragment is a Krylov subspace of dynamically connected basis states of `H`. Two basis
states ψ_1 and ψ_2 are in the same fragment if their Krylov subspaces generated with `H` are
equal. The size of a fragment is equal to the dimension of the corresponding Krylov
subspace, which coincides with the number of basis states in it. If `periodic` is `true`,
only necklaces of basis states are considered.
https://en.wikipedia.org/wiki/Krylov_subspace
https://en.wikipedia.org/wiki/Necklace_(combinatorics)

# Returns
- `Vector{Int}`: Vector of fragment sizes, sorted in ascending order
"""
function fragment_sizes(H::MPO, periodic::Bool)::Vector{Int}
    site_inds = firstsiteinds(H)
    num_sites = length(site_inds)
    C = combiner(site_inds...)
    ci = combinedind(C)
    basis_states_with_fragment_id = Dict{BitVector,Int}()
    fragments = Dict{Int,Set{BitVector}}()
    fragment_id_counter = 0
    prog = ProgressUnknown(desc="Calculating fragment sizes:", spinner=true)
    for root_state in 0:(2^num_sites-1)
        root_state = BitVector(digits(root_state; base=2, pad=num_sites))
        if periodic
            root_state = necklace(root_state)
        end
        if !haskey(basis_states_with_fragment_id, root_state)
            working_set = Set{BitVector}([root_state])
            finished_set = Set{BitVector}()
            while length(working_set) > 0
                next!(prog)
                next_state = pop!(working_set)
                push!(finished_set, next_state)
                next_state_mps = bitvector_to_mps(next_state, site_inds)
                superposition = C * contract(apply(H, next_state_mps))
                for index in eachval(ci)
                    if !isapprox(superposition[index], 0.0; atol=1e-5)
                        state = BitVector(digits(index - 1; base=2, pad=num_sites))
                        if periodic
                            state = necklace(state)
                        end
                        if state ∉ finished_set
                            push!(working_set, state)
                        end
                    end
                end
            end
            for state in finished_set
                basis_states_with_fragment_id[state] = fragment_id_counter
                fragments[fragment_id_counter] = finished_set
            end
            fragment_id_counter += 1
        end
    end
    finish!(prog)

    return sort([length(fragments[key]) for key in eachindex(fragments)])
end
