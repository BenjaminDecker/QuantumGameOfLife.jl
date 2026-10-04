const INITIAL_STATE_CHOICES = String.(keys(INITIAL_STATE_REGISTRY))
const FILE_FORMAT_CHOICES = ["pdf", "png", "svg", "eps"]
# PLOTS_CHOICES and ALGORITHM_CHOICES are derived from the registries in types.jl.

# Single source of truth for all option defaults, shared with the TOML config loader.
const DEFAULTS = default_config()

const settings = ArgParseSettings(
    prog="cli.jl",
    description="A classical simulation of the quantum game of life. Start a run by giving a TOML config file (julia cli.jl run.toml). If command line options are used instead, the effective configuration is written to a timestamped config file so the run can be reproduced exactly.",
    autofix_names=true,
    error_on_conflict=false,
    exit_after_help=false
)
add_arg_group!(settings, "Config")
@add_arg_table! settings begin
    "config"
    help = "Path to a TOML config file. If given, it fully determines the run and all other command line options are ignored."

    "--write-config"
    arg_type = String
    default = ""
    help = "Write the effective configuration to this path (only relevant when starting the run with command line options). Defaults to a timestamped file in the plot directory."
end

add_arg_group!(settings, "Setup")
@add_arg_table! settings begin
    "--num-cells"
    arg_type = Int
    default = DEFAULTS[:num_cells]
    help = "The number of cells to use in the simulation. Depending on the algorithm used, the running time can scale exponentially(exact) or linearly(tdvp) with the number of cells."

    "--initial-states"
    arg_type = String
    nargs = '*'
    default = DEFAULTS[:initial_states]
    range_tester = x -> x in INITIAL_STATE_CHOICES
    help = "Initial State. If more than one is given, an equal superposition of the states is used. Choices are: " * string(INITIAL_STATE_CHOICES)

    "--superposition"
    action = :store_true
    help = "Create a superposition of all states given in --initial-states instead of calculating a time evolution for each one separately"
end

add_arg_group!(settings, "Rule")
@add_arg_table! settings begin
    "--distance"
    arg_type = Int
    default = DEFAULTS[:distance]
    help = "The interaction distance to the left and right of a cell. This controls the range of local operators in the Hamiltonian. For only nearest neighbor interactions, use 1."

    "--rule"
    arg_type = Int
    default = DEFAULTS[:rule]

    "--activation-interval"
    arg_type = Int
    nargs = 2
    # default = Int[1, 1]
    metavar = ["LowerBound", "UpperBound"]
    help = "Range of alive neighbors required for a flip, upper bound is included"
end

add_arg_group!(settings, "Algorithm")
@add_arg_table! settings begin
    "--algorithm"
    arg_type = Algorithm
    default = DEFAULTS[:algorithm]
    help = "The algorithm used for the time evolution. 'exact' is fast and most accurate for a small numbers of cells. Choices are: " * string(ALGORITHM_CHOICES)

    "--num-steps"
    arg_type = Int
    default = DEFAULTS[:num_steps]
    help = "Number of time steps to simulate"

    "--periodic-boundaries"
    action = :store_true
    help = "Use periodic instead of open boundary conditions"

    "--step-size"
    arg_type = Float64
    default = DEFAULTS[:step_size]
    help = "Size of one time step. The time step size is calculated as (STEP_SIZE * pi/2)"

    "--sweeps-per-time-step"
    arg_type = Int
    default = DEFAULTS[:sweeps_per_time_step]
    help = "The number of sweeps to perform per time step. Is ignored if the chosen algorithm is 'exact'."

    "--max-bond-dim"
    arg_type = Int
    default = DEFAULTS[:max_bond_dim]
    help = "The maximum that a bond of the MPS is allowed to grow to during simulation. Is ignored if the chosen algorithm is not 'tdvp'."

    "--svd-epsilon"
    arg_type = Float64
    default = DEFAULTS[:svd_epsilon]
    help = "A measure of accuracy for the truncation step after splitting a mps tensor. This parameter controls how quickly the bond dimension of the mps grows during the simulation. Lower means more accurate, but slower."

    "--operator-set"
    arg_type = Int
    default = DEFAULTS[:operator_set]
    range_tester = x -> x in 1:4
    help = "Set of operators used to build up the hamiltonian in the non-hermitian case."
end

add_arg_group!(settings, "Plot")
@add_arg_table! settings begin
    "--show"
    action = :store_true
    help = "Open plots in their respective default applications"

    "--plot"
    arg_type = PlotType
    nargs = '*'
    default = DEFAULTS[:plot]
    help = "Plots to create. Choices are: " * string(PLOTS_CHOICES)

    "--plotting-file-path"
    arg_type = String
    default = DEFAULTS[:plotting_file_path]
    help = "Write files to a directory at the specified relative location"

    "--file-formats"
    arg_type = String
    nargs = '*'
    default = DEFAULTS[:file_formats]
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
    default = DEFAULTS[:px_per_unit]
    help = "The size of one unit length of the plot in px"
end

add_arg_group!(settings, "Fragmentation Analysis")
@add_arg_table! settings begin
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
    return parse_plot_type(x)
end

function ArgParse.parse_item(::Type{Algorithm}, x::AbstractString)
    return parse_algorithm(x)
end

"""
    parse_commandline()::Union{Dict{Symbol,Any},Nothing}

Parse `ARGS` using the configured ArgParse settings and return the effective configuration
as a normalized dictionary that can be passed to the `Args` constructor or written to a
config file with `write_run_config`. Returns `nothing` when parsing was interrupted (for
example after printing the `--help` message).

If a config file is given on the command line, it fully determines the run and any other
command line options are ignored (with a warning). Otherwise, the effective configuration is
written to a config file so that the exact same run can be started later.
"""
function parse_commandline()::Union{Dict{Symbol,Any},Nothing}
    parsed = parse_args(settings; as_symbols=true)
    isnothing(parsed) && return nothing

    config_file = pop!(parsed, :config)
    config_file = isnothing(config_file) ? "" : String(config_file)
    write_config_path = String(pop!(parsed, :write_config))

    cfg = if !isempty(config_file)
        ignored = _explicit_cli_options(parsed)
        if !isempty(ignored)
            @warn "Ignoring command line options because a config file was given: $(join(ignored, ", "))"
        end
        load_config(config_file)
    else
        path = write_run_config(parsed; path=isempty(write_config_path) ? nothing : write_config_path)
        @info "Wrote the effective configuration to '$path'. Start the exact same run with: julia cli.jl $path"
        normalize_config(parsed)
    end
    return cfg
end

"""
    _explicit_cli_options(parsed::Dict{Symbol,Any})::Vector{Symbol}

Return the names of all options in `parsed` (as returned by ArgParse) whose value differs
from the default value, i.e. the ones the user explicitly set.
"""
function _explicit_cli_options(parsed::Dict{Symbol,Any})::Vector{Symbol}
    return [key for (key, value) in parsed if get(DEFAULTS, key, nothing) != value]
end
