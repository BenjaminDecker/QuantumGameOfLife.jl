using Combinatorics

function bipartite_entropy(psi::MPS, separator_index::Int)::Float64
    orthogonalize!(psi, separator_index)
    try
        _, S, _ = svd(
            psi[separator_index],
            (linkind(psi, separator_index - 1), siteind(psi, separator_index))
        )
        SvN = 0.0
        for n = 1:dim(S, 1)
            p = S[n, n]^2
            if p > 0.0
                SvN -= p * log(p)
            end
        end
        return real(SvN)
    catch _
        return NaN
    end
end

function center_bipartite_entropy(psi::MPS)::Float64
    center = trunc(Int, length(psi) / 2)
    bipartite_entropy(psi, center)
end
