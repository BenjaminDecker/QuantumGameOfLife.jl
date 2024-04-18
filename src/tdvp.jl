using ITensors

const TENSOR_1::ITensor = ITensor(1)

struct Tdvp_data
    H::MPO
    psi::MPS

    H_eff_layers_left::Vector{ITensor}
    H_eff_layers_right::Vector{ITensor}
    max_bond_dims::Vector{Int}
end

function evolve_tdvp(H::MPO, psi_0::MPS, args::Args)::Vector{MPS}

    # ---------- Data Preparation ----------
    T = (-im) * (args.step_size / args.sweeps_per_time_step) * (pi / 2) / 2
    H_eff_layers_left::Vector{ITensor} = [TENSOR_1]
    H_eff_layers_right::Vector{ITensor} = [TENSOR_1]
    # max_bond_dims::Vector{Int} = [
    # min(2^i, 2^(args.num_cells - i), args.max_bond_dim)
    # for i in 1:(args.num_cells-1)
    # ]
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

    # ---------- Function Definitions ------
    function evolve(A::ITensor, H_eff::ITensor, backwards::Bool)::ITensor
        noprime(A * exp(H_eff * T * (backwards ? -1 : 1)))
    end


    evolve_site(site::ITensor, H_site::ITensor)::ITensor = evolve(
        site,
        H_eff_layers_left[end] * H_site * H_eff_layers_right[end],
        false
    )

    evolve_bond(bond::ITensor)::ITensor = evolve(
        bond,
        H_eff_layers_left[end] * H_eff_layers_right[end],
        true
    )

    create_layer(site_idx::Int)::ITensor = contract(
        psi[site_idx],
        H[site_idx],
        dag(prime(psi[site_idx]))
    )

    push_layer_left!(site_idx::Int) = push!(
        H_eff_layers_left,
        H_eff_layers_left[end] * create_layer(site_idx)
    )

    push_layer_right!(site_idx::Int) = push!(
        H_eff_layers_right,
        H_eff_layers_right[end] * create_layer(site_idx)
    )


    # ---------- Main Code -----------------
    orthogonalize!(psi, 1)
    for site_idx in args.num_cells:-1:2
        push_layer_right!(site_idx)
    end

    @showprogress "Calculating Time Evolution" for _ in 2:args.num_steps
        for _ in 1:args.sweeps_per_time_step
            for site_idx in 1:args.num_cells
                psi[site_idx] = evolve_site(psi[site_idx], H[site_idx])
                if site_idx != args.num_cells
                    inds_left = uniqueinds(psi[site_idx], psi[site_idx+1])
                    site_orthogonal, bond = qr(psi[site_idx], inds_left)
                    psi[site_idx] = site_orthogonal
                    push_layer_left!(site_idx)
                    new_bond = evolve_bond(bond)
                    pop!(H_eff_layers_right)
                    psi[site_idx+1] = psi[site_idx+1] * new_bond
                end
            end
            for site_idx in args.num_cells:-1:1
                psi[site_idx] = evolve_site(psi[site_idx], H[site_idx])
                if site_idx != 1
                    inds_right = uniqueinds(psi[site_idx], psi[site_idx-1])
                    site_orthogonal, bond = qr(psi[site_idx], inds_right)
                    psi[site_idx] = site_orthogonal
                    push_layer_right!(site_idx)
                    new_bond = evolve_bond(bond)
                    pop!(H_eff_layers_left)
                    psi[site_idx-1] = psi[site_idx-1] * new_bond
                end
            end
        end
        push!(results, copy(psi))
    end
    return results
end
