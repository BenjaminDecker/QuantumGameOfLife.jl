using ITensors
using PyPlot
using PyCall
using SplitApplyCombine
using DefaultApplication

struct LabeledPlot
    label::String
    data::Vector{Vector{Float64}}
end

function plot_results(
    ;
    heatmaps_continuous::Vector{LabeledPlot},
    heatmaps_discrete::Vector{LabeledPlot}=Vector{LabeledPlot}(),
    line_plots::Vector{LabeledPlot}=Vector{LabeledPlot}(),
    path::String,
    file_formats::Vector{String}
)

    plt = pyimport("matplotlib.pyplot")

    cmap = "inferno"
    magic_number = 1.4

    num_plots = length(heatmaps_continuous) + length(heatmaps_discrete) + length(line_plots)
    num_cells = length(heatmaps_continuous[1].data[1])
    num_steps = length(heatmaps_continuous[1].data)

    colorbar_width = 0.15 # As a fraction of the subplot height

    figsize = (magic_number / num_cells) .* [num_steps, num_cells * num_plots]
    gs_kw = Dict("width_ratios" => [1 * num_steps / num_cells, colorbar_width], "wspace" => 0.5 * magic_number / figsize[1])

    mosaic = [
        ["plot_$(i)", "colorbar_$(i <= length(heatmaps_continuous) ? "continuous" : i)"]
        for i in 1:num_plots
    ]
    fig, axs = plt.subplot_mosaic(mosaic, gridspec_kw=gs_kw)
    for i in 2:num_plots
        axs["plot_$(i)"].sharex(axs["plot_1"])
    end
    for i in 1:(num_plots-1)
        axs["plot_$(i)"].tick_params(labelbottom=false)
    end

    heatmap_continuous_imgs = [
        let plot = heatmaps_continuous[i]
            ax = axs["plot_$(i)"]
            ax.figure.set_size_inches(figsize...)
            ax.set_ylabel(plot.label)
            ax.pcolormesh(real(combinedims(plot.data)), cmap=cmap, vmin=0.0, vmax=1.0)
        end
        for i in eachindex(heatmaps_continuous)
    ]
    cbar = fig.colorbar(heatmap_continuous_imgs[end], cax=axs["colorbar_continuous"])


    heatmap_discrete_imgs = [
        let plot = heatmaps_discrete[i]
            axs_index = i + length(heatmaps_continuous)
            ax = axs["plot_$(axs_index)"]
            ax.figure.set_size_inches(figsize...)
            ax.set_ylabel(plot.label)
            vmin = trunc(Int, minimum(minimum.(plot.data)))
            vmax = trunc(Int, maximum(maximum.(plot.data)))
            img = ax.pcolormesh(real(combinedims(plot.data)), cmap=get_cmap(cmap, vmax - vmin + 1), vmin=vmin, vmax=vmax)
            cbar = fig.colorbar(img, cax=axs["colorbar_$(axs_index)"])
            cbar.ax.locator_params(nbins=min(vmax + 1, 6))
            img.set_clim(-0.5 + vmin, (vmax + 0.5))
            img
        end
        for i in eachindex(heatmaps_discrete)
    ]

    line_plot_imgs = [
        let plot = line_plots[i]
            axs_index = i + length(heatmaps_continuous) + length(heatmaps_discrete)
            ax = axs["plot_$(axs_index)"]
            ax.figure.set_size_inches(figsize...)
            ax.set_ylabel(plot.label)
            axs["colorbar_$(axs_index)"].axis("off")
            ax.plot(plot.data)
        end
        for i in eachindex(line_plots)
    ]

    axs["plot_$(num_plots)"].set_xlabel("time steps")
    for suffix in file_formats
        file_path = "$(path).$(suffix)"
        savefig(file_path, bbox_inches="tight")
        DefaultApplication.open(file_path)
    end
    nothing
end