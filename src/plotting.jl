struct LabeledPlot
    label::String
    data::Vector{Vector{Float64}}
end

function plot(
    measurements_vector::Vector{Dict{PlotType,Vector{Vector{Float64}}}},
    args::Args
)
    S_max = (args.num_cells * log(2) - 1) / 2

    for (j, measurements) in enumerate(measurements_vector)

        measurements_sorted = sort(collect(measurements), by=x -> x[1])

        width = @something args.width 600
        f = Figure(size=(width, width))


        axis_ids_of_HeatmapContinuous = []
        axes::Vector{Makie.Axis} = []

        for (i, (type, data)) in enumerate(measurements_sorted)

            if isa(type, HeatmapDiscrete)
                ax = Axis(f[i, 1], ylabel=label(type), yticks=[1; length(data[1])])
                num_categories = max(2, Int(maximum(maximum.(data))) + 1)
                cmap = cgrad(:inferno, num_categories; categorical=true)
                push!(axes, ax)
                heatmap!(
                    ax,
                    0:(length(data)-1),
                    1:(length(data[1])),
                    transpose(reduce(hcat, data)),
                    colormap=cmap,
                )
                pos = 1 / (2 * num_categories)
                Colorbar(
                    f[i, 2],
                    colormap=cmap,
                    ticks=([pos, 1 - pos], ["0", string(num_categories - 1)]),
                )
            elseif isa(type, LinePlot)
                ax = Axis(f[i, 1], limits=((-0.5, length(data) - 0.5), (0, nothing)), ylabel=label(type), yticks=0:ceil(S_max))
                if args.page_entropy
                    hlines!(S_max; color=:red, label="page entropy")
                    axislegend(ax, position=:rb)
                end
                push!(axes, ax)
                data = flatten(data)
                _ = lines!(
                    ax,
                    0:args.num_steps,
                    data,
                    linewidth=2
                )
            else
                ax = Axis(f[i, 1], ylabel=label(type), yticks=[1; length(data[1])])
                push!(axes, ax)
                _ = heatmap!(
                    ax,
                    0:(length(data)-1),
                    1:(length(data[1])),
                    transpose(reduce(hcat, data)),
                    colormap=:inferno,
                )
                push!(axis_ids_of_HeatmapContinuous, i)
            end
        end
        if !isempty(axis_ids_of_HeatmapContinuous)
            Colorbar(
                f[minimum(axis_ids_of_HeatmapContinuous):maximum(axis_ids_of_HeatmapContinuous), 2],
                colormap=:inferno
            )
        end

        axes[end].xlabel = "Time Steps"

        linkxaxes!(axes)
        for ax in axes[1:(end-1)]
            hidexdecorations!(
                ax,
                ticks=false,
                grid=false,
                minorgrid=false,
                minorticks=false,
            )
        end

        axis_height = 80
        for r in 1:nrows(f.layout)
            rowsize!(f.layout, r, Fixed(axis_height))
        end

        aspect_ratio = length(measurements_sorted[1][2]) / length(measurements_sorted[1][2][1])
        if isnothing(args.width)
            colsize!(f.layout, 1, Fixed(aspect_ratio * axis_height))
        else
            f.layout.width = args.width
        end

        resize_to_layout!(f)
        save_and_show(f, args, "$j")
    end

end


function plot_eigval_vs_cbe(
    eigval::Vector{Float64},
    cbe::Vector{Float64},
    args::Args
)

    f = Figure()
    S_max = (args.num_cells * log(2) - 1) / 2
    ax = Axis(f[1, 1], xlabel="Energy Density E/L", ylabel="Center Bipartite Entropy")
    scatter!(ax, eigval, cbe)
    hlines!(S_max, linestyle=:dash, color=:red, label=L"Page Entropy $\frac{L\cdot log(2)-1}{2}$")
    axislegend(ax, position=:lt)
    save_and_show(f, args, "eigval_vs_cbe")
end

function plot_fragment_sizes(fragment_sizes::Vector{Int}, args::Args)
    f = Figure()
    ax = Axis(f[1, 1], xlabel="Fragment", ylabel="Fragment Size", xticks=eachindex(fragment_sizes))
    barplot!(ax, eachindex(fragment_sizes), fragment_sizes)
    # hidedecorations!(
    #     ax,
    #     label=false,
    #     ticklabels=false,
    #     ticks=false,
    #     # grid=false,
    #     minorgrid=false,
    #     minorticks=false,
    # )
    save_and_show(f, args, "fragment_sizes")
end

function save_and_show(f, args::Args, id::String="")
    rule_filename =
        "$(args.num_cells)" *
        "-$(args.distance)" *
        "-$(args.rule)" *
        "-$(args.num_steps)" *
        "-$(replace(string(args.step_size), "." => ""))" *
        "-$(args.periodic ? "periodic" : "open")" *
        "-$(lowercase(name(args.algorithm)))"

    path =
        args.plotting_file_path *
        "/$(rule_filename)" *
        "-$(join(args.initial_states_names, "-"))" *
        "-$(join(filename_identifier.(args.plots), "-"))" *
        "-$(args.max_bond_dim)-" *
        id

    for suffix in args.file_formats
        file_path = "$(path).$(suffix)"
        save(file_path, f; px_per_unit=args.px_per_unit)
        if args.show
            DefaultApplication.open(file_path)
        end
    end
end
