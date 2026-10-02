"""
    evolve(::TEBD, psi_0_vec::Vector{MPS}, H::MPO, args::Args)

Time evolution via the TEBD algorithm.

!!! warning "Not implemented"
    TEBD is not implemented yet. Calling this method throws an
    `ErrorException`. It is kept as a placeholder so the `TEBD` algorithm type stays
    selectable while the implementation is being written.
"""
function evolve(::TEBD, ::Vector{MPS}, ::MPO, ::Args)::Vector{Vector{MPS}}
    error("The TEBD algorithm is not implemented yet.")
end
