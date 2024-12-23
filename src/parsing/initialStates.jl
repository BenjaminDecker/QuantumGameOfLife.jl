using Random

function blinker(site_inds::Vector{ITensors.Index{Int64}}, width::Int=1)::MPS
    mid = floor(length(site_inds) / 2) + 1
    plist = [i == (mid - width) || i == (mid + width) ? "1" : "0" for i in 1:length(site_inds)]
    MPS(site_inds, plist)
end

blinker_wide(site_inds::Vector{ITensors.Index{Int64}})::MPS = blinker(site_inds, 2)

function triple_blinker(site_inds::Vector{ITensors.Index{Int64}})::MPS
    mid = floor(length(site_inds) / 2) + 1
    plist = [i == (mid - 2) || i == mid || i == (mid + 2) ? "1" : "0" for i in 1:length(site_inds)]
    MPS(site_inds, plist)
end

function alternating(site_inds::Vector{ITensors.Index{Int64}}, reversed::Bool=false)::MPS
    plist = [xor(i % 2 == 0, reversed) ? "1" : "0" for i in 1:length(site_inds)]
    MPS(site_inds, plist)
end

alternating_reversed(site_inds::Vector{ITensors.Index{Int64}})::MPS = alternating(site_inds, true)

function single(site_inds::Vector{ITensors.Index{Int64}}, position::Int=-1)::MPS
    if !(position in 1:length(site_inds))
        position = floor(length(site_inds) / 2) + 1
    end
    plist = [i == position ? "1" : "0" for i in 1:length(site_inds)]
    MPS(site_inds, plist)
end

function single_wide(site_inds::Vector{ITensors.Index{Int64}}, position::Int=-1)::MPS
    if !(position in 1:length(site_inds))
        position = floor(length(site_inds) / 2) + 1
    end
    plist = [i == position || i == (position - 1) ? "1" : "0" for i in 1:length(site_inds)]
    MPS(site_inds, plist)
end

single_bottom(site_inds::Vector{ITensors.Index{Int64}})::MPS = single(site_inds, firstindex(site_inds))

single_top(site_inds::Vector{ITensors.Index{Int64}})::MPS = single(site_inds, lastindex(site_inds))

function single_bottom_half(site_inds::Vector{ITensors.Index{Int64}})::MPS
    single(site_inds, Int(length(site_inds) / 4))
end

function single_top_half(site_inds::Vector{ITensors.Index{Int64}})::MPS
    single(site_inds, Int(length(site_inds) - length(site_inds) / 4))
end

function all_ket_0(site_inds::Vector{ITensors.Index{Int64}})::MPS
    plist = fill("0", length(site_inds))
    MPS(site_inds, plist)
end

function all_ket_1(site_inds::Vector{ITensors.Index{Int64}})::MPS
    plist = fill("1", length(site_inds))
    MPS(site_inds, plist)
end

function all_ket_0_but_outer(site_inds::Vector{ITensors.Index{Int64}})::MPS
    plist = [i == 1 || i == length(site_inds) ? "1" : "0" for i in 1:length(site_inds)]
    MPS(site_inds, plist)
end

function all_ket_1_but_outer(site_inds::Vector{ITensors.Index{Int64}})::MPS
    plist = [i == 1 || i == length(site_inds) ? "0" : "1" for i in 1:length(site_inds)]
    MPS(site_inds, plist)
end

function equal_superposition(site_inds::Vector{ITensors.Index{Int64}})::MPS
    plist = fill("+", length(site_inds))
    MPS(site_inds, plist)
end

function equal_superposition_but_outer_ket_0(site_inds::Vector{ITensors.Index{Int64}})::MPS
    plist = fill("+", length(site_inds))
    plist[1] = plist[length(site_inds)] = "0"
    MPS(site_inds, plist)
end

function equal_superposition_but_outer_ket_1(site_inds::Vector{ITensors.Index{Int64}})::MPS
    plist = fill("+", length(site_inds))
    plist[1] = plist[length(site_inds)] = "1"
    MPS(site_inds, plist)
end

function single_bottom_blinker_top(site_inds::Vector{ITensors.Index{Int64}})::MPS
    plist = ["1", "0", "1"]
    append!(plist, fill("0", length(site_inds) - 4))
    push!(plist, "1")
    MPS(site_inds, plist)
end

random(site_inds::Vector{ITensors.Index{Int64}})::MPS = randomMPS(site_inds)

function random_product(site_inds::Vector{ITensors.Index{Int64}})::MPS
    plist = map(x -> x ? "1" : "0", bitrand(length(site_inds)))
    MPS(site_inds, plist)
end
