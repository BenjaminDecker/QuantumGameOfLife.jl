using TensorTimeSteps

function evolve(::TDVP1, psi_0_vec::Vector{MPS}, H::MPO, args)::Vector{Vector{MPS}}
    return [
        tdvp1(
            H,
            psi_0;
            step_size=-im * args.step_size * (pi / 2),
            num_steps=args.num_steps,
            sweeps_per_time_step=args.sweeps_per_time_step,
            max_bond_dim=args.max_bond_dim
        )
        for psi_0 in psi_0_vec
    ]
end


function evolve(::TDVP2, psi_0_vec::Vector{MPS}, H::MPO, args)::Vector{Vector{MPS}}
    return [
        tdvp2(
            H,
            psi_0;
            step_size=-im * args.step_size * (pi / 2),
            num_steps=args.num_steps,
            sweeps_per_time_step=args.sweeps_per_time_step,
            max_bond_dim=args.max_bond_dim,
            svd_epsilon=args.svd_epsilon
        )
        for psi_0 in psi_0_vec
    ]
end
