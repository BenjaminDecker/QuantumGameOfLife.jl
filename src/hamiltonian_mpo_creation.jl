ITensors.op(::OpName"imS-im1", ::SiteType"Qubit") = [0 im; 0 -im]

ITensors.op(::OpName"imS-im0", ::SiteType"Qubit") = [-im 0; im 0]

ITensors.op(::OpName"0", ::SiteType"Qubit") = [0 0; 0 0]

ITensors.op(::OpName"1", ::SiteType"Qubit") = [1 0; 0 1]

ITensors.op(::OpName"P0S+", ::SiteType"Qubit") = [1 1; 0 0]

ITensors.op(::OpName"P1S-", ::SiteType"Qubit") = [0 0; 1 1]

ITensors.op(::OpName"IS+", ::SiteType"Qubit") = [1 1; 0 1]

ITensors.op(::OpName"IS-", ::SiteType"Qubit") = [1 0; 1 1]

const operators::Vector{Dict{Tuple{Int,Int},String}} = [
    Dict{Tuple{Int,Int},String}(
        (0, 0) => "imS-im1",
        (0, 1) => "0",
        (1, 0) => "X",
        (1, 1) => "imS-im0"
    ),
    Dict{Tuple{Int,Int},String}(
        (0, 0) => "S+",
        (0, 1) => "0",
        (1, 0) => "X",
        (1, 1) => "S-"
    ),
    Dict{Tuple{Int,Int},String}(
        (0, 0) => "P0S+",
        (0, 1) => "0",
        (1, 0) => "X",
        (1, 1) => "P1S-"
    ),
    Dict{Tuple{Int,Int},String}(
        (0, 0) => "IS+",
        (0, 1) => "0",
        (1, 0) => "X",
        (1, 1) => "IS-"
    ),
]

function get_correct_index(index::Int, args::Args)::Union{Nothing,Int}
    if !checkbounds(Bool, 1:args.num_cells, index) && !args.periodic
        return nothing
    end
    return (index - 1 + args.num_cells) % args.num_cells + 1
end

function local_opsum_tuple_list(index::Int, args::Args)::Vector{Tuple}
    op = []
    op_dict = operators[args.operator_set]
    for left_side in 0:(2^(args.distance)-1)
        for right_side in 0:(2^(args.distance)-1)
            ops = []
            conf_id = (left_side << (args.distance + 1)) + (right_side)
            rule_result_0 = args.rule & 1 << conf_id != 0
            rule_result_1 = args.rule & 1 << (conf_id + 1 << args.distance) != 0
            operator = op_dict[(rule_result_0, rule_result_1)]
            push!(ops, operator)
            push!(ops, index)

            bits_left = digits(left_side, base=2, pad=args.distance)
            bits_right = reverse(digits(right_side, base=2, pad=args.distance))
            for i in 1:args.distance
                correct_index_left = get_correct_index(index - i, args)
                if !isnothing(correct_index_left)
                    push!(ops, Bool(bits_left[i]) ? "Proj1" : "Proj0")
                    push!(ops, correct_index_left)
                elseif Bool(bits_left[i])
                    @goto escape_label
                end
                correct_index_right = get_correct_index(index + i, args)
                if !isnothing(correct_index_right)
                    push!(ops, Bool(bits_right[i]) ? "Proj1" : "Proj0")
                    push!(ops, correct_index_right)
                elseif Bool(bits_right[i])
                    @goto escape_label
                end
            end

            push!(op, Tuple(ops))
            @label escape_label
        end
    end
    return op
end


"""
    build_hamiltonian_mpo(site_inds::Vector{ITensors.Index{Int64}}, args::Args)::MPO

Build the Hamiltonian of the quantum game of life as an `MPO` acting on `site_inds`.

For every cell, `local_opsum_tuple_list` generates the local operator summand determined by
the Wolfram `rule`, the interaction `distance` and the selected `operator_set`. The
resulting `OpSum` is converted to an `MPO`. If the resulting Hamiltonian is (numerically)
hermitian, the time evolution is unitary and this is reported via the `Logging` macro.
"""
function build_hamiltonian_mpo(
    site_inds::Vector{ITensors.Index{Int64}},
    args::Args
)::MPO
    @assert length(site_inds) == args.num_cells

    os = OpSum()
    for idx in 1:args.num_cells
        for local_opsum_tuple in local_opsum_tuple_list(idx, args)
            os .+= local_opsum_tuple
        end
    end
    mpo = MPO(os, site_inds)
    hermiticity_residual = norm(mpo - dag(swapprime(mpo, 0, 1)))
    if hermiticity_residual < 1e-5
        @info "The Hamiltonian is hermitian, the time evolution is unitary."
    end
    return mpo
end
