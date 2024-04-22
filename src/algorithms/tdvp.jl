using ITensors

const TENSOR_1::ITensor = ITensor(1)

evolve(A::ITensor, H_eff::ITensor, T::ComplexF64)::ITensor =
    noprime(A * exp(H_eff * T))

push_layer!(layers::Vector{ITensor}, site::ITensor, H_site::ITensor) =
    push!(layers, layers[end] * site * H_site * dag(prime(site)))

pop_layer!(layers::Vector{ITensor}) = pop!(layers)


function evolve(::TDVP1, psi_0::MPS, H::MPO, args::Args)::Vector{MPS}
    T = (-im) * (args.step_size / args.sweeps_per_time_step) * (pi / 2) / 2
    layers_left::Vector{ITensor} = [TENSOR_1]
    layers_right::Vector{ITensor} = [TENSOR_1]
    results = [psi_0]
    sizehint!(results, args.num_steps)
    psi = copy(psi_0)

    # Fix the bond dimensions
    # https://itensor.discourse.group/t/how-do-i-set-an-mps-bond-dimension-that-is-higher-than-needed/1637
    psi = +(
        psi,
        0 * randomMPS(siteinds(psi); linkdims=args.max_bond_dim - maximum(linkdims(psi)));
        alg="directsum"
    )
    @assert maximum(linkdims(psi)) == args.max_bond_dim
    @assert isapprox(inner(psi, psi_0), 1)


    orthogonalize!(psi, 1)
    for site_idx in args.num_cells:-1:2
        push_layer!(
            layers_right,
            psi[site_idx],
            H[site_idx]
        )
    end

    @showprogress "Calculating Time Evolution" for _ in 2:args.num_steps
        for _ in 1:args.sweeps_per_time_step
            for site_idx in 1:args.num_cells
                psi[site_idx] = evolve(
                    psi[site_idx],
                    layers_left[end] * H[site_idx] * layers_right[end],
                    T
                )
                if site_idx != args.num_cells
                    inds_left = uniqueinds(psi[site_idx], psi[site_idx+1])
                    site_orthogonal, bond = qr(psi[site_idx], inds_left)
                    psi[site_idx] = site_orthogonal
                    push_layer!(
                        layers_left,
                        psi[site_idx],
                        H[site_idx]
                    )
                    new_bond = evolve(
                        bond,
                        layers_left[end] * layers_right[end],
                        -T
                    )
                    pop_layer!(layers_right)
                    psi[site_idx+1] *= new_bond
                end
            end
            for site_idx in args.num_cells:-1:1
                psi[site_idx] = evolve(
                    psi[site_idx],
                    layers_left[end] * H[site_idx] * layers_right[end],
                    T
                )
                if site_idx != 1
                    inds_right = uniqueinds(psi[site_idx], psi[site_idx-1])
                    site_orthogonal, bond = qr(psi[site_idx], inds_right)
                    psi[site_idx] = site_orthogonal
                    push_layer!(
                        layers_right,
                        psi[site_idx],
                        H[site_idx]
                    )
                    new_bond = evolve(
                        bond,
                        layers_left[end] * layers_right[end],
                        -T
                    )
                    pop_layer!(layers_left)
                    psi[site_idx-1] *= new_bond
                end
            end
        end
        push!(results, copy(psi))
    end
    return results
end


function evolve(::TDVP2, psi_0::MPS, H::MPO, args::Args)::Vector{MPS}
    T = (-im) * (args.step_size / args.sweeps_per_time_step) * (pi / 2) / 2
    layers_left::Vector{ITensor} = [TENSOR_1]
    layers_right::Vector{ITensor} = [TENSOR_1]
    max_bond_dims::Vector{Int} = [
        min(2^i, 2^(args.num_cells - i), args.max_bond_dim)
        for i in 1:(args.num_cells-1)
    ]
    results = [psi_0]
    sizehint!(results, args.num_steps)
    psi = copy(psi_0)

    orthogonalize!(psi, 1)
    for site_idx in args.num_cells:-1:3
        push_layer!(
            layers_right,
            psi[site_idx],
            H[site_idx]
        )
    end

    @showprogress "Calculating Time Evolution" for _ in 2:args.num_steps
        for _ in 1:args.sweeps_per_time_step
            for site_idx in 1:(args.num_cells-1)
                two_site_tensor = evolve(
                    psi[site_idx] * psi[site_idx+1],
                    layers_left[end] * H[site_idx] * H[site_idx+1] * layers_right[end],
                    T
                )
                inds_left = uniqueinds(psi[site_idx], psi[site_idx+1])
                left, S, right = svd(
                    two_site_tensor,
                    inds_left;
                    maxdim=max_bond_dims[site_idx],
                    cutoff=args.svd_epsilon
                )
                psi[site_idx] = left
                psi[site_idx+1] = S * right
                if site_idx != (args.num_cells - 1)
                    push_layer!(
                        layers_left,
                        psi[site_idx],
                        H[site_idx]
                    )
                    psi[site_idx+1] = evolve(
                        psi[site_idx+1],
                        layers_left[end] * H[site_idx+1] * layers_right[end],
                        -T
                    )
                    pop_layer!(layers_right)
                end
            end
            for site_idx in args.num_cells:-1:2
                two_site_tensor = evolve(
                    psi[site_idx-1] * psi[site_idx],
                    layers_left[end] * H[site_idx-1] * H[site_idx] * layers_right[end],
                    T
                )
                inds_right = uniqueinds(psi[site_idx], psi[site_idx-1])
                right, S, left = svd(
                    two_site_tensor,
                    inds_right;
                    maxdim=max_bond_dims[site_idx-1],
                    cutoff=args.svd_epsilon
                )
                psi[site_idx] = right
                psi[site_idx-1] = left * S
                if site_idx != 2
                    push_layer!(
                        layers_right,
                        psi[site_idx],
                        H[site_idx]
                    )
                    psi[site_idx-1] = evolve(
                        psi[site_idx-1],
                        layers_left[end] * H[site_idx-1] * layers_right[end],
                        -T
                    )
                    pop_layer!(layers_left)
                end
            end
        end
        push!(results, copy(psi))
    end
    return results
end
