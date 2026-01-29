const INITIAL_STATE_CHOICES = ["blinker", "blinker_wide", "triple_blinker", "alternating", "alternating_reversed", "single", "single_wide", "single_bottom", "single_top", "single_bottom_half", "single_top_half", "all_ket_0", "all_ket_1", "all_ket_0_but_outer", "all_ket_1_but_outer", "equal_superposition", "equal_superposition_but_outer_ket_0", "equal_superposition_but_outer_ket_1", "single_bottom_blinker_top", "random", "random_product"]
const FILE_FORMAT_CHOICES = ["pdf", "png", "svg", "eps"]
const PLOTS_CHOICES = ["classical", "expect", "sse", "rounded", "bond_dims", "cbe", "autocorrelation"]
const ALGORITHM_CHOICES = ["exact", "tdvp1", "tdvp2", "sierpinski"] #TODO tebd

s = ArgParseSettings(
    prog="main.jl",
    description="A classical simulation of the quantum game of life",
    autofix_names=true,
    error_on_conflict=false,
    exit_after_help=false
)
add_arg_group!(s, "Setup")
@add_arg_table! s begin
    "--num-cells"
    arg_type = Int
    default = 9
    help = "The number of cells to use in the simulation. Depending on the algorithm used, the running time can scale exponentially(exact) or linearly(tdvp) with the number of cells."

    "--initial-states"
    arg_type = String
    nargs = '*'
    default = String["blinker"]
    range_tester = x -> x in INITIAL_STATE_CHOICES
    help = "Initial State. If more than one is given, an equal superposition of the states is used. Choices are: " * string(INITIAL_STATE_CHOICES)

    "--superposition"
    action = :store_true
    help = "Create a superposition of all states given in --initial-states instead of calculating a time evolution for each one separately"
end

add_arg_group!(s, "Rule")
@add_arg_table! s begin
    "--distance"
    arg_type = Int
    default = 1
    help = "The interaction distance to the left and right of a cell. This controls the range of local operators in the Hamiltonian. For only nearest neighbor interactions, use 1."

    "--rule"
    arg_type = Int
    default = 150

    "--activation-interval"
    arg_type = Int
    nargs = 2
    # default = Int[1, 1]
    metavar = ["LowerBound", "UpperBound"]
    help = "Range of alive neighbors required for a flip, upper bound is included"
end

add_arg_group!(s, "Algorithm")
@add_arg_table! s begin
    "--algorithm"
    arg_type = Algorithm
    default = Exact()
    help = "The algorithm used for the time evolution. 'exact' is fast and most accurate for a small numbers of cells. Choices are: " * string(ALGORITHM_CHOICES)

    "--num-steps"
    arg_type = Int
    default = 100
    help = "Number of time steps to simulate"

    "--periodic-boundaries"
    action = :store_true
    help = "Use periodic instead of open boundary conditions"

    "--step-size"
    arg_type = Float64
    default = 1.0
    help = "Size of one time step. The time step size is calculated as (STEP_SIZE * pi/2)"

    "--sweeps-per-time-step"
    arg_type = Int
    default = 100
    help = "The number of sweeps to perform per time step. Is ignored if the chosen algorithm is 'exact'."

    "--max-bond-dim"
    arg_type = Int
    default = 32
    help = "The maximum that a bond of the MPS is allowed to grow to during simulation. Is ignored if the chosen algorithm is not 'tdvp'."

    "--svd-epsilon"
    arg_type = Float64
    default = 1e-10
    help = "A measure of accuracy for the truncation step after splitting a mps tensor. This parameter controls how quickly the bond dimension of the mps grows during the simulation. Lower means more accurate, but slower."

    "--operator-set"
    arg_type = Int
    default = 1
    range_tester = x -> x in 1:4
    help = "Set of operators used to build up the hamiltonian in the non-hermitian case."
end

add_arg_group!(s, "Plot")
@add_arg_table! s begin
    "--show"
    action = :store_true
    help = "Open plots in their respective default applications"

    "--plot"
    arg_type = PlotType
    nargs = '*'
    default = [ExpectationValue()]
    help = "Plots to create. Choices are: " * string(PLOTS_CHOICES)

    "--plotting-file-path"
    arg_type = String
    default = "plots"
    help = "Write files to a directory at the specified relative location"

    "--file-formats"
    arg_type = String
    nargs = '*'
    default = String["pdf"]
    range_tester = x -> x in FILE_FORMAT_CHOICES
    help = "File formats for plots. Choices are: " * string(FILE_FORMAT_CHOICES)

    "--width"
    arg_type = Int
    help = "Plot width. If omitted, plot width will grow with num-steps such that heatmap datapoints look like squares. If you are unsure, a good default value is 600."

    "--page-entropy"
    action = :store_true
    help = "Show the page entropy value in cbe plots"

    "--px-per-unit"
    arg_type = Float64
    default = 2.0
    help = "The size of one unit length of the plot in px"
end

add_arg_group!(s, "Fragmentation Analysis")
@add_arg_table! s begin
    "--plot-eigval-vs-cbe"
    action = :store_true
    help = "Plot the eigenvalues vs. the center bipartite entropy of the hamiltonian's eigenvectors"

    "--plot-fragment-sizes"
    action = :store_true
    help = "Plot the fragment sizes of the Hamiltonian. Sizes are plotted up to periodic symmetry if periodic boundary conditions are used."

    # "--include-frozen-states"
    # action = :store_true
    # help = "Include frozen states. Frozen states are eigenstates of the Hamiltonian that are also product states in the z-basis. This option is ignored if none of the othr Fragmentation Analysis options is set."
end

function ArgParse.parse_item(::Type{PlotType}, x::AbstractString)
    x = lowercase(x)
    if x in ["classic", "classical"]
        return Classical()
    end
    if x in ["expect", "expectation", "expectation_value", "expectation-value"]
        return ExpectationValue()
    end
    if x in ["sse", "single_site_entropy", "single-site-entropy"]
        return SingleSiteEntropy()
    end
    if x in ["round", "rounded"]
        return Rounded()
    end
    if x in ["bond_dim", "bond_dims", "bond_dimension", "bond_dimensions", "bond-dim", "bond-dims", "bond-dimension", "bond-dimensions"]
        return BondDimensions()
    end
    if x in ["cbe", "center_bipartite_entropy", "center-bipartite-entropy"]
        return CenterBipartiteEntropy()
    end
    if x in ["autocorrelation"]
        return Autocorrelation()
    end
    throw(ArgumentError("Not a valid plot type"))
end

function ArgParse.parse_item(::Type{Algorithm}, x::AbstractString)
    x = lowercase(x)
    if x == "exact"
        return Exact()
    end
    if x == "tdvp1"
        return TDVP1()
    end
    if x == "tdvp2"
        return TDVP2()
    end
    if x == "sierpinski" || x == "sierpiński"
        return Sierpinski()
    end
    if x == "tebd"
        return TEBD()
    end
    throw(ArgumentError("Not a valid Algorithm"))
end

function get_args()
    args = parse_args(s; as_symbols=true)
    return isnothing(args) ? Nothing : Args(args)
end
