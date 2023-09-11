using ITensors

function bipartite_entropy(psi::MPS, seperator_index::Int)
    orthogonalize!(psi, seperator_index)
    _, S, _ = svd(
        psi[seperator_index],
        (linkind(psi, seperator_index - 1), siteind(psi, seperator_index))
    )
    SvN = 0.0
    for n = 1:dim(S, 1)
        p = S[n, n]^2
        if p > 0.0
            SvN -= p * log(p)
        end
    end
    real(SvN)
end

function center_bipartite_entropy(psi::MPS)
    center = trunc(Int, length(psi) / 2)
    bipartite_entropy(psi, center)
end