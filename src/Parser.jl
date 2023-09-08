module Parser

import ArgParse: ArgParseSettings, parse_args, add_arg_group!, @add_arg_table!

const INITIAL_STATE_CHOICES = ["blinker", "triple_blinker", "full_blinker", "single", "single_bottom", "all_ket_0", "all_ket_1", "only_outer", "all_ket_0_but_outer", "all_ket_1_but_outer", "equal_superposition", "equal_superposition_but_outer_ket_0", "equal_superposition_but_outer_ket_1", "random"]
const ALGORITHM_CHOICES = ["exact", "1tdvp", "2tdvp", "a1tdvp"]
const FILE_FORMAT_CHOICES = ["html", "eps", "jpeg", "jpg", "pdf", "pgf",
    "png", "ps", "raw", "rgba", "svg", "svgz", "tif", "tiff"]
const PLOTS_CHOICES = ["expect", "sse"]

s = ArgParseSettings(
    prog="main.jl",
    description="A classical simulation of the quantum game of life",
    autofix_names=true,
    error_on_conflict=false
)
add_arg_group!(s, "Rules")
@add_arg_table! s begin
    "--num-cells"
    arg_type = Int
    default = 9
    help = "The number of cells to use in the simulation. Computation running time scales exponentially with NUM_CELLS. Anything higher than 11 takes a lot of time and memory to compute."

    "--distance"
    arg_type = Int
    default = 1
    help = "The distance each cell looks for alive or dead neighbours"

    "--activation-interval"
    arg_type = Int
    nargs = 2
    default = Int[1, 2]
    metavar = ["LowerBound", "UpperBound"]
    help = "Range of alive neighbours required for a flip, upper bound is excluded"

    "--num-steps"
    arg_type = Int
    default = 100
    help = "Number of time steps to simulate"

    "--periodic-boundaries"
    action = :store_true
    help = "Use periodic instead of open boundary conditions"
end

add_arg_group!(s, "Initial States")
@add_arg_table! s begin
    "--initial-states"
    arg_type = String
    nargs = '*'
    default = ["blinker"]
    range_tester = x -> x in INITIAL_STATE_CHOICES
    help = "Initial States. Choices are: " * string(INITIAL_STATE_CHOICES)
end

add_arg_group!(s, "Algorithm")
@add_arg_table! s begin
    "--algorithm"
    arg_type = String
    default = "exact"
    range_tester = x -> x in ALGORITHM_CHOICES
    help = "The algorithm used for the time evolution. Use 'exact' for a small number of cells, otherwise some version of 'tdvp'. Choices are: " * string(ALGORITHM_CHOICES)

    "--step-size"
    arg_type = Float64
    default = 1.0
    help = "Size of one time step. The time step size is calculated as (STEP_SIZE * pi/2)"

    "--max-bond-dim"
    arg_type = Int
    default = 32
    help = "The maximum that a bond of the MPS is allowed to grow to during simulation. Is ignored if the chosen algorithm is not 'tdvp'."

    "--svd-epsilon"
    arg_type = Float64
    default = 0.00005
    help = "A measure of accuracy for the truncation step after splitting a mps tensor. This parameter controls how quickly the bond dimension of the mps grows during the simulation. Lower means more accurate, but slower."
end

add_arg_group!(s, "Plot")
@add_arg_table! s begin
    "--plotting-frequency"
    arg_type = Float64
    default = 1.0
    help = "Frequency at which time steps are plotted. Time between plot steps is calculated as (pi/2 * 1/PLOT_FREQUENCY * 1/STEP_SIZE)"

    "--plots"
    arg_type = String
    nargs = '*'
    default = String["expect"]
    range_tester = x -> x in PLOTS_CHOICES
    help = "Plots to create. Choices are: " * string(FILE_FORMAT_CHOICES)

    "--plot-sse"
    action = :store_true
    help = "Plot the single site entropy"

    "--plot-bond-dims"
    action = :store_true
    help = "Plot the bond dimensions of the mps"

    "--plot-rounded"
    action = :store_true
    help = "Plot a rounded version of the probability"

    "--plot-cbe"
    action = :store_true
    help = "Plot the center bipartite entropy"

    "--show"
    action = :store_true
    help = "Open plots in their respective default applications"

    "--plot-file-path"
    arg_type = String
    default = "plots"
    help = "Write files to a directory at the specified relative location"

    "--file-formats"
    arg_type = String
    nargs = '*'
    default = String["pdf"]
    range_tester = x -> x in FILE_FORMAT_CHOICES
    help = "File formats for plots. Choices are: " * string(FILE_FORMAT_CHOICES)

    "--plot-fragmentation"
    action = :store_true
    help = "Plot the the fragment sizes of the state fragmentation of the Hamiltonian"
end

add_arg_group!(s, "Fragmentation Analysis")
@add_arg_table! s begin
    "--plot-eigval-vs-cbe"
    action = :store_true
    help = "Plot the eigenvalues vs. the center bipartite entropy of the hamiltonian's eigenvectors"

    "--plot-fragment-sizes"
    action = :store_true
    help = "Plot the fragment sizes of the Hamiltonian"
end

parse_args() = parse_args(s; as_symbols=true)

end