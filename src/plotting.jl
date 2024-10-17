using PyPlot
using SplitApplyCombine
using DefaultApplication

struct LabeledPlot
    label::String
    data::Vector{Vector{Float64}}
end

function plot(
    measurements_vector::Vector{Dict{PlotType,Vector{Vector{Float64}}}},
    args::Args
)
    cmap = "inferno"
    magic_number = 1.2
    colorbar_width = 0.15 # As a fraction of the subplot height
    S_max = (args.num_cells * log(2) - 1) / 2



    for (j, measurements) in enumerate(measurements_vector)

        measurements_sorted = sort(collect(measurements), by=x -> x[1])
        num_plots = length(measurements_sorted)

        figsize = (magic_number / args.num_cells) .*
                  [args.num_steps, args.num_cells * num_plots]
        gs_kw = Dict(
            "width_ratios" => [1 * args.num_steps / args.num_cells, colorbar_width],
            "wspace" => 0.5 * magic_number / figsize[1],
            "hspace" => 0.3
        )

        mosaic = [
            [
                "plot_$(i)",
                "colorbar_$(isa(measurements_sorted[i][1], HeatmapContinuous) ? "continuous" : i)"
            ]
            for i in eachindex(measurements_sorted)
        ]
        fig, axs = plt.subplot_mosaic(mosaic, gridspec_kw=gs_kw)
        for i in 2:num_plots
            axs["plot_$(i)"].sharex(axs["plot_1"])
        end
        for i in 1:(num_plots-1)
            axs["plot_$(i)"].tick_params(labelbottom=false)
        end

        colorbar_created::Bool = false

        for (i, (type, data)) in enumerate(measurements_sorted)
            ax = axs["plot_$(i)"]
            ax.figure.set_size_inches(figsize...)
            ax.set_ylabel(label(type))
            if isa(type, HeatmapContinuous)
                img = ax.pcolormesh(
                    real(combinedims(data)),
                    cmap=cmap,
                    vmin=0.0,
                    vmax=1.0
                )
                if !colorbar_created
                    fig.colorbar(
                        img,
                        cax=axs["colorbar_continuous"]
                    )
                    colorbar_created = true
                end
            elseif isa(type, HeatmapDiscrete)
                vmin = trunc(Int, minimum(minimum.(data)))
                vmax = trunc(Int, maximum(maximum.(data)))
                img = ax.pcolormesh(
                    real(combinedims(data)),
                    cmap=get_cmap(cmap, vmax - vmin + 1),
                    vmin=vmin,
                    vmax=vmax
                )
                cbar = fig.colorbar(img, cax=axs["colorbar_$(i)"])
                cbar.ax.locator_params(nbins=min(vmax + 1, 6))
                img.set_clim(-0.5 + vmin, (vmax + 0.5))
            elseif isa(type, LinePlot)
                axs["colorbar_$(i)"].axis("off")
                ax.plot(data)
                if (isa(type, CenterBipartiteEntropy))
                    ax.hlines(
                        S_max,
                        0,
                        length(data) - 1,
                        colors="red",
                        linestyles="--",
                        label="page entropy"
                    )
                    ax.legend(loc="lower right")
                end
            end
        end

        axs["plot_$(num_plots)"].set_xlabel("Time Steps")
        write_and_show(args, "$j")
        PyPlot.close()
    end
end

function plot_eigval_vs_cbe(
    eigval::Vector{Float64},
    cbe::Vector{Float64},
    args::Args
)
    scatter(eigval, cbe)
    xlabel("Energy Density E/L")
    ylabel("Center Bipartite Entropy")
    S_max = (args.num_cells * log(2) - 1) / 2
    PyPlot.axhline(
        S_max,
        color="red",
        linestyle="--",
        label="Page Entropy \$\\frac{L\\cdot log(2)-1}{2}\$"
    )
    PyPlot.legend(loc="upper left", framealpha=1.0)
    write_and_show(args, "eigval_vs_cbe")
    PyPlot.close()
end

function plot_fragment_sizes(fragment_sizes::Vector{Int}, args::Args)
    bar(eachindex(fragment_sizes), fragment_sizes)
    xlabel("Fragment")
    ylabel("Fragment Size")
    PyPlot.yticks(
        Vector(1:ceil(Int, fragment_sizes[end] / 20):fragment_sizes[end])
    )
    PyPlot.xticks(
        Vector(1:ceil(Int, length(fragment_sizes) / 10):length(fragment_sizes))
    )
    write_and_show(args, "fragment_sizes")
    PyPlot.close()
end

function write_and_show(args::Args, id::String="")
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
        savefig(file_path, bbox_inches="tight"; dpi=args.dpi)
        if args.show
            DefaultApplication.open(file_path)
        end
    end
end
