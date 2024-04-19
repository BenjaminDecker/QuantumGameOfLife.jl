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
