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
    single(site_inds, length(site_inds) ÷ 4)
end

function single_top_half(site_inds::Vector{ITensors.Index{Int64}})::MPS
    single(site_inds, length(site_inds) - length(site_inds) ÷ 4)
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

"""
    INITIAL_STATE_REGISTRY::Dict{String,Function}

Maps the name of each supported initial state (as accepted by the `--initial-states`
command line option) to the function that builds the corresponding `MPS`. Each function
takes a `Vector{ITensors.Index{Int64}}` of site indices as its only argument.
"""
const INITIAL_STATE_REGISTRY = Dict{String, Function}(
    "blinker" => blinker,
    "blinker_wide" => blinker_wide,
    "triple_blinker" => triple_blinker,
    "alternating" => alternating,
    "alternating_reversed" => alternating_reversed,
    "single" => single,
    "single_wide" => single_wide,
    "single_bottom" => single_bottom,
    "single_top" => single_top,
    "single_bottom_half" => single_bottom_half,
    "single_top_half" => single_top_half,
    "all_ket_0" => all_ket_0,
    "all_ket_1" => all_ket_1,
    "all_ket_0_but_outer" => all_ket_0_but_outer,
    "all_ket_1_but_outer" => all_ket_1_but_outer,
    "equal_superposition" => equal_superposition,
    "equal_superposition_but_outer_ket_0" => equal_superposition_but_outer_ket_0,
    "equal_superposition_but_outer_ket_1" => equal_superposition_but_outer_ket_1,
    "single_bottom_blinker_top" => single_bottom_blinker_top,
    "random" => random,
    "random_product" => random_product,
)
