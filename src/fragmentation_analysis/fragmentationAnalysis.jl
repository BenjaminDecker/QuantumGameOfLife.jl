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
    D, U = eigen(contract(H), ishermitian=true) # TODO what about non-hermitian matrices?
    println("done")

    eigen_ind = inds(U, tags="eigen")[1]
    eigvals = [D[eigen_ind=>i, eigen_ind'=>i] for i in eachval(eigen_ind)]


    # rat_eigvals = [rationalize(BigInt, eigval) for eigval in eigvals]
    # filter!(x -> x != 0 // 1, rat_eigvals)
    # # println([rationalize(BigInt, eigval) for eigval in eigvals])
    # for i in eachindex(rat_eigvals)
    #     lcm_result = lcm(rat_eigvals[1:i])
    #     if (lcm_result == 0 // 1)
    #         println(rat_eigvals[i])
    #         asd()
    #     end
    #     println(lcm(rat_eigvals[1:i]))
    # end


    cbe = @showprogress "Calculating center bipartite entropy of eigenvectors" [
        center_bipartite_entropy(MPS(U * onehot(eigen_ind => i), site_inds))
        for i in eachval(eigen_ind)
    ]

    (eigvals ./ length(site_inds), cbe)
end


function bitvec_to_int(arr)
    v = 1
    acc = 0
    for bit in view(arr, 1:length(arr))
        acc += v * bit
        v <<= 1
    end
    acc
end


function is_frozen(arr)
    last_bit = arr[end]
    for bit in view(arr, 1:length(arr))
        if last_bit && bit
            return true
        end
        last_bit = bit
    end
    return false
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
    basis_states::Vector{Union{Nothing,Int}} = [nothing for _ in 1:2^num_sites]
    fragments = Dict{Int,Set{BitVector}}()
    fragment_id_counter = 0
    prog = ProgressUnknown(desc="Calculating fragment sizes:", spinner=true)
    for root_state_idx in 0:(2^num_sites-1)
        root_state = BitVector(digits(root_state_idx; base=2, pad=num_sites))
        if periodic
            # root_state = necklace(root_state)
        end
        if isnothing(basis_states[root_state_idx+1])
            fragment_id_counter += 1
            # println(root_state)
            working_set = Set{BitVector}([root_state])
            finished_set = Set{BitVector}()
            while length(working_set) > 0
                next!(prog)
                next_state = pop!(working_set)
                push!(finished_set, next_state)
                next_state_mps = bitvector_to_mps(next_state, site_inds)
                superposition = C * contract(apply(H, next_state_mps))
                for index in eachval(ci)
                    if !isapprox(superposition[index], 0.0; atol=1e-10)
                        if !isnothing(basis_states[index])
                            # throw(ErrorException("Uff"))
                        end
                        state = BitVector(digits(index - 1; base=2, pad=num_sites))
                        if periodic
                            # state = necklace(state)
                        end
                        if state ∉ finished_set
                            push!(working_set, state)
                        end
                    end
                end
            end
            for state in finished_set
                basis_states[bitvec_to_int(state)+1] = fragment_id_counter
                fragments[fragment_id_counter] = finished_set
            end
        end
    end
    finish!(prog)

    println(sum(length(fragments[key]) for key in eachindex(fragments)))
    # sorted_fragments = sort([(fragments[key]) for key in eachindex(fragments)], by=x -> length(x))
    # for f in sorted_fragments
    #     println(f)
    # end
    # println(sort(collect(sorted_fragments[1])))

    println(sort([length(fragments[key]) for key in eachindex(fragments)]))
    return sort([length(fragments[key]) for key in eachindex(fragments)])
end
