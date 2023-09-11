using ITensors
using ProgressMeter

function eigval_vs_cbe(H::MPO)
    site_inds = firstsiteinds(H)
    print("Calculating eigen decomposition...")
    D, U = @disable_warn_order eigen(contract(H); ishermitian=true)
    println("done")
    eigen_ind = inds(U, tags="eigen")[1]
    eigvals = [D[eigen_ind=>i, eigen_ind'=>i] for i in eachval(eigen_ind)]
    cbe = @showprogress "Calculating center bipartite entropy of eigenvectors" [
        center_bipartite_entropy(MPS(U * onehot(eigen_ind => i), site_inds))
        for i in eachval(eigen_ind)
    ]
    (eigvals ./ length(site_inds), cbe)
end