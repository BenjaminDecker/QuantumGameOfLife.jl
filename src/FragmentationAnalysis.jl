module FragmentationAnalysis

using ITensors
using Plots

import .Utils

function eigval_vs_entropy(H::MPO)
    site_inds = firstsiteinds(H)
    D, U = @disable_warn_order eigen(contract(H))
    eigen_ind = inds(U, tags="eigen")[1]
    eigvals = [D[eigen_ind=>i, eigen_ind'=>i] for i in eachval(eigen_ind)]
    cbe = [Utils.bipartite_entropy(MPS(U * onehot(eigen_ind => i), site_inds)) for i in eachval(eigen_ind)]
    (real(eigvals), real(cbe))
end
end