module InitialStates

using ITensors
using Random

function blinker(site_inds::Vector{ITensors.Index{Int64}}, width::Int=1)
    mid = floor(length(site_inds) / 2) + 1
    plist = [i == (mid - width) || i == (mid + width) ? "1" : "0" for i in 1:length(site_inds)]
    MPS(site_inds, plist)
end

function triple_blinker(site_inds::Vector{ITensors.Index{Int64}})
    mid = floor(length(site_inds) / 2) + 1
    plist = [i == (mid - 2) || i == mid || i == (mid + 2) ? "1" : "0" for i in 1:length(site_inds)]
    MPS(site_inds, plist)
end

function full_blinker(site_inds::Vector{ITensors.Index{Int64}})
    plist = [i % 2 == 0 ? "1" : "0" for i in 1:length(site_inds)]
    MPS(site_inds, plist)
end

function single(site_inds::Vector{ITensors.Index{Int64}}, position::Int=-1)
    if !(position in 1:length(site_inds))
        position = floor(length(site_inds) / 2) + 1
    end
    plist = [i == position ? "1" : "0" for i in 1:length(site_inds)]
    MPS(site_inds, plist)
end

function single_bottom(site_inds::Vector{ITensors.Index{Int64}})
    single(site_inds, 1)
end

function all_ket_0(site_inds::Vector{ITensors.Index{Int64}})
    plist = fill("0", length(site_inds))
    MPS(site_inds, plist)
end

function all_ket_1(site_inds::Vector{ITensors.Index{Int64}})
    plist = fill("1", length(site_inds))
    MPS(site_inds, plist)
end

function only_outer(site_inds::Vector{ITensors.Index{Int64}})
    plist = [i == 1 || i == length(site_inds) ? "1" : "0" for i in 1:length(site_inds)]
    MPS(site_inds, plist)
end

function all_ket_0_but_outer(site_inds::Vector{ITensors.Index{Int64}})
    plist = [i == 1 || i == length(site_inds) ? "1" : "0" for i in 1:length(site_inds)]
    MPS(site_inds, plist)
end

function all_ket_1_but_outer(site_inds::Vector{ITensors.Index{Int64}})
    plist = [i == 1 || i == length(site_inds) ? "0" : "1" for i in 1:length(site_inds)]
    MPS(site_inds, plist)
end

function equal_superposition(site_inds::Vector{ITensors.Index{Int64}})
    plist = fill("+", length(site_inds))
    MPS(site_inds, plist)
end

function equal_superposition_but_outer_ket_0(site_inds::Vector{ITensors.Index{Int64}})
    plist = fill("+", length(site_inds))
    plist[1] = plist[length(site_inds)] = "0"
    MPS(site_inds, plist)
end

function equal_superposition_but_outer_ket_1(site_inds::Vector{ITensors.Index{Int64}})
    plist = fill("+", length(site_inds))
    plist[1] = plist[length(site_inds)] = "1"
    MPS(site_inds, plist)
end

function random(site_inds::Vector{ITensors.Index{Int64}})
    plist = map(x -> x ? "1" : "0", bitrand(length(site_inds)))
    MPS(site_inds, plist)
end
end