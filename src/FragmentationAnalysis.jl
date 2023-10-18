module FragmentationAnalysis

using ITensors
using ProgressMeter
using FromFile

@from "Utils.jl" using Utils

export eigval_vs_cbe, fragment_sizes

function eigval_vs_cbe(H::MPO)::Tuple{Vector{Float64},Vector{Float64}}
    site_inds = firstsiteinds(H)
    print("Calculating eigen decomposition...")
    D, U = eigen(contract(H); ishermitian=true)
    println("done")
    eigen_ind = inds(U, tags="eigen")[1]
    eigvals = [D[eigen_ind=>i, eigen_ind'=>i] for i in eachval(eigen_ind)]
    cbe = @showprogress "Calculating center bipartite entropy of eigenvectors" [
        Utils.center_bipartite_entropy(MPS(U * onehot(eigen_ind => i), site_inds))
        for i in eachval(eigen_ind)
    ]
    (eigvals ./ length(site_inds), cbe)
end

function eigval_vs_cbe(H::MPO, eigenvectors::Set{BitVector})::Tuple{Vector{Float64},Vector{Float64}}
    site_inds = firstsiteinds(H)
    eigvals = Vector{Float64}()
    cbe = Vector{Float64}()
    for eigvec in eigenvectors
        eigvec_mps = bitvector_to_mps(eigvec, site_inds)
        push!(cbe, center_bipartite_entropy(eigvec_mps))
        push!(eigvals, sqrt(inner(eigvec_mps, eigvec_mps)))
    end
    (eigvals ./ length(site_inds), cbe)
end

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

# function is_in_subspace(bits::BitVector)::Bool
#     return true
#     last_value = bits[end]
#     for value in bits
#         if last_value && value
#             return false
#         end
#         last_value = value
#     end
#     return true
# end

bitvector_to_mps(bitvector::BitVector, site_inds::Vector{Index{Int64}})::MPS = MPS(site_inds, map(string, filter(x -> x == '1' || x == '0', collect(bitstring(bitvector)))))

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
        # if !is_in_subspace(root_state)
        # continue
        # end
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
                superposition = C * contract(apply(H, next_state_mps; alg="naive", truncate=false))
                for index in eachval(ci)
                    if superposition[index] != 0.0
                        state = BitVector(digits(index - 1; base=2, pad=num_sites))
                        # if !is_in_subspace(state)
                        # continue
                        # end
                        if periodic
                            state = necklace(state)
                        end
                        if !(state in finished_set)
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
end