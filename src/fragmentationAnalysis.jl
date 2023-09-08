using ITensors

function eigval_vs_cbe(H::MPO)
    site_inds = firstsiteinds(H)
    print("Performing eigen decomposition...")
    D, U = @disable_warn_order eigen(contract(H))
    println("done")
    eigen_ind = inds(U, tags="eigen")[1]
    eigvals = [D[eigen_ind=>i, eigen_ind'=>i] for i in eachval(eigen_ind)]
    print("Calculating center bipartite entropies of eigenvectors...")
    cbe = [center_bipartite_entropy(MPS(U * onehot(eigen_ind => i), site_inds)) for i in eachval(eigen_ind)]
    println("done")
    (real(eigvals) ./ length(site_inds), real(cbe))
end