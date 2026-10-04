"""
TOML based run configuration.

A run can be fully described by a TOML config file which is passed to `cli.jl` as the first
argument. For reproducibility, runs that are started with command line options instead write
their effective configuration (given options plus defaults) to such a file.
"""

# The sections and options of a config file. Option keys are unique across all sections and
# map directly to the symbol keys that the `Args` constructor expects.
const CONFIG_SECTIONS = [
    "setup" => [:num_cells, :initial_states, :superposition],
    "rule" => [:distance, :rule, :activation_interval],
    "algorithm" =>
        [:algorithm, :num_steps, :periodic_boundaries, :step_size, :sweeps_per_time_step, :max_bond_dim, :svd_epsilon, :operator_set],
    "plot" => [:show, :plot, :plotting_file_path, :file_formats, :width, :page_entropy, :px_per_unit],
    "fragmentation_analysis" => [:plot_eigval_vs_cbe, :plot_fragment_sizes],
]

"""
    default_config()::Dict{Symbol,Any}

The effective default value for every option. This is the single source of truth for the
defaults used by the command line parser and the TOML config loader.
"""
function default_config()::Dict{Symbol,Any}
    return Dict{Symbol,Any}(
        :num_cells => 9,
        :initial_states => String["blinker"],
        :superposition => false,
        :distance => 1,
        :rule => 150,
        :activation_interval => Int[],
        :algorithm => Exact(),
        :num_steps => 100,
        :periodic_boundaries => false,
        :step_size => 1.0,
        :sweeps_per_time_step => 100,
        :max_bond_dim => 32,
        :svd_epsilon => 1e-10,
        :operator_set => 1,
        :show => false,
        :plot => PlotType[ExpectationValue()],
        :plotting_file_path => "plots",
        :file_formats => String["pdf"],
        :width => nothing,
        :page_entropy => false,
        :px_per_unit => 2.0,
        :plot_eigval_vs_cbe => false,
        :plot_fragment_sizes => false
    )
end

function parse_plot_type(x::AbstractString)::PlotType
    key = lowercase(x)
    return get(PLOT_LOOKUP, key) do
        throw(ArgumentError("Not a valid plot type: '$key'. Choices are: " * string(PLOTS_CHOICES)))
    end
end

function parse_algorithm(x::AbstractString)::Algorithm
    key = lowercase(x)
    return get(ALGORITHM_LOOKUP, key) do
        throw(ArgumentError("Not a valid algorithm: '$key'. Choices are: " * string(ALGORITHM_CHOICES)))
    end
end

# Value checking and conversion for the values coming from a TOML file.
_config_type_error(key::Symbol, expected::String, value) =
    error("Invalid value for '$key': expected $expected, got $(repr(value))")

function _as_int(key::Symbol, value)::Int
    value isa Integer || _config_type_error(key, "an integer", value)
    return Int(value)
end

function _as_float(key::Symbol, value)::Float64
    value isa Real || _config_type_error(key, "a number", value)
    return Float64(value)
end

function _as_bool(key::Symbol, value)::Bool
    value isa Bool || _config_type_error(key, "a boolean", value)
    return value
end

function _as_string(key::Symbol, value)::String
    value isa AbstractString || _config_type_error(key, "a string", value)
    return String(value)
end

function _as_int_vector(key::Symbol, value)::Vector{Int}
    value isa AbstractVector || _config_type_error(key, "an array of integers", value)
    return [_as_int(key, v) for v in value]
end

function _as_string_vector(key::Symbol, value)::Vector{String}
    value isa AbstractVector || _config_type_error(key, "an array of strings", value)
    return [_as_string(key, v) for v in value]
end

"""
    normalize_config(cfg::Dict{Symbol,Any})::Dict{Symbol,Any}

Convert and validate the raw option values in `cfg` (TOML types) into the Julia types and
values that the `Args` constructor expects. Throws an informative error on invalid input.
"""
function normalize_config(cfg::Dict{Symbol,Any})::Dict{Symbol,Any}
    cfg[:num_cells] = _as_int(:num_cells, cfg[:num_cells])
    cfg[:num_cells] >= 1 || error("Invalid value for 'num_cells': must be at least 1")

    cfg[:initial_states] = _as_string_vector(:initial_states, cfg[:initial_states])
    isempty(cfg[:initial_states]) && error("Invalid value for 'initial_states': must not be empty")
    for state in cfg[:initial_states]
        state in INITIAL_STATE_CHOICES ||
            error("Invalid value for 'initial_states': '$state' is not one of " * string(INITIAL_STATE_CHOICES))
    end

    cfg[:superposition] = _as_bool(:superposition, cfg[:superposition])

    cfg[:distance] = _as_int(:distance, cfg[:distance])
    cfg[:distance] >= 1 || error("Invalid value for 'distance': must be at least 1")
    cfg[:rule] = _as_int(:rule, cfg[:rule])

    activation_interval = _as_int_vector(:activation_interval, cfg[:activation_interval])
    isempty(activation_interval) || length(activation_interval) == 2 ||
        error("Invalid value for 'activation_interval': expected [lower_bound, upper_bound] or []")
    cfg[:activation_interval] = activation_interval

    cfg[:algorithm] = cfg[:algorithm] isa Algorithm ? cfg[:algorithm] : parse_algorithm(_as_string(:algorithm, cfg[:algorithm]))
    cfg[:num_steps] = _as_int(:num_steps, cfg[:num_steps])
    cfg[:periodic_boundaries] = _as_bool(:periodic_boundaries, cfg[:periodic_boundaries])
    cfg[:step_size] = _as_float(:step_size, cfg[:step_size])
    cfg[:sweeps_per_time_step] = _as_int(:sweeps_per_time_step, cfg[:sweeps_per_time_step])
    cfg[:max_bond_dim] = _as_int(:max_bond_dim, cfg[:max_bond_dim])
    cfg[:svd_epsilon] = _as_float(:svd_epsilon, cfg[:svd_epsilon])

    cfg[:operator_set] = _as_int(:operator_set, cfg[:operator_set])
    cfg[:operator_set] in 1:4 || error("Invalid value for 'operator_set': must be one of 1, 2, 3, 4")

    cfg[:show] = _as_bool(:show, cfg[:show])
    plots = cfg[:plot] isa AbstractVector ? cfg[:plot] : [cfg[:plot]]
    cfg[:plot] = PlotType[p isa PlotType ? p : parse_plot_type(_as_string(:plot, p)) for p in plots]
    cfg[:plotting_file_path] = _as_string(:plotting_file_path, cfg[:plotting_file_path])

    cfg[:file_formats] = _as_string_vector(:file_formats, cfg[:file_formats])
    for format in cfg[:file_formats]
        format in FILE_FORMAT_CHOICES ||
            error("Invalid value for 'file_formats': '$format' is not one of " * string(FILE_FORMAT_CHOICES))
    end

    cfg[:width] = isnothing(cfg[:width]) ? nothing : _as_int(:width, cfg[:width])
    cfg[:page_entropy] = _as_bool(:page_entropy, cfg[:page_entropy])
    cfg[:px_per_unit] = _as_float(:px_per_unit, cfg[:px_per_unit])

    cfg[:plot_eigval_vs_cbe] = _as_bool(:plot_eigval_vs_cbe, cfg[:plot_eigval_vs_cbe])
    cfg[:plot_fragment_sizes] = _as_bool(:plot_fragment_sizes, cfg[:plot_fragment_sizes])

    return cfg
end

"""
    load_config(path::AbstractString)::Dict{Symbol,Any}

Read the TOML config file at `path`, fill in the default values for all options that are not
set in the file and normalize the result such that it can be passed to the `Args`
constructor. Unknown sections or options are reported as errors to catch typos.
"""
function load_config(path::AbstractString)::Dict{Symbol,Any}
    isfile(path) || error("Config file not found: $path")
    raw = try
        TOML.parse(read(path, String))
    catch err
        error("Could not parse '$path' as TOML: $(sprint(showerror, err))")
    end

    section_keys = Dict(section => keys for (section, keys) in CONFIG_SECTIONS)
    cfg = default_config()
    for (section, values) in raw
        haskey(section_keys, section) ||
            error("Unknown section '[$section]' in $path. Sections are: " * string(collect(keys(section_keys))))
        values isa Dict || error("Section '[$section]' in $path must be a table")
        for (key, value) in values
            option = Symbol(replace(key, '-' => '_'))
            option in section_keys[section] ||
                error("Unknown option '$key' in section '[$section]' of $path. Options are: " * string(section_keys[section]))
            cfg[option] = value
        end
    end
    return normalize_config(cfg)
end

function _toml_value(key::Symbol, value)
    if key === :algorithm
        return lowercase(name(value))
    elseif key === :plot
        return sort([filename_identifier(p) for p in value])
    elseif value isa AbstractVector
        return collect(value)
    else
        return value
    end
end

"""
    to_toml(cfg::Dict{Symbol,Any})::Dict{String,Any}

Convert a normalized configuration into a nested `String`-keyed dictionary that can be
written with `TOML.print`. Julia-only values (algorithm and plot types) are converted to
their string identifiers; options with no value (`nothing`) are omitted.
"""
function to_toml(cfg::Dict{Symbol,Any})::Dict{String,Any}
    document = Dict{String,Any}()
    for (section, keys) in CONFIG_SECTIONS
        table = Dict{String,Any}()
        for key in keys
            value = get(cfg, key, nothing)
            isnothing(value) && continue
            table[String(key)] = _toml_value(key, value)
        end
        document[section] = table
    end
    return document
end

"""
    write_run_config(cfg::Dict{Symbol,Any}; path::Union{AbstractString,Nothing}=nothing)::String

Write the effective configuration `cfg` to a TOML file and return the path. If no `path` is
given, a timestamped file is created in the configured plot directory. Running `cli.jl` with
the resulting file reproduces the same run exactly.
"""
function write_run_config(cfg::Dict{Symbol,Any}; path::Union{AbstractString,Nothing}=nothing)::String
    if isnothing(path)
        dir = String(cfg[:plotting_file_path])
        mkpath(dir)
        stamp = Dates.format(Dates.now(), "yyyymmdd-HHMMSS")
        path = joinpath(dir, "run_config_$(stamp).toml")
        n = 1
        while isfile(path)
            path = joinpath(dir, "run_config_$(stamp)_$(n).toml")
            n += 1
        end
    else
        dir = dirname(String(path))
        isempty(dir) || mkpath(dir)
    end

    open(path, "w") do io
        println(io, "# Run configuration written by QuantumGameOfLife.jl on $(Dates.now()).")
        println(io, "# Start the exact same run with: julia cli.jl $(path)")
        println(io)
        TOML.print(io, to_toml(cfg); sorted=true)
    end
    return String(path)
end
