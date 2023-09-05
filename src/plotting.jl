using ITensors
using PyPlot
using SplitApplyCombine

struct LabeledPlot
    label::String
    data::Vector{Vector{Float64}}
end

function plot_results(; continuous_plots::Vector{LabeledPlot}, discrete_plots::Vector{LabeledPlot}=Vector{LabeledPlot}(), path::String, file_formats::Vector{String})

    cmap = "inferno"

    num_plots = length(continuous_plots) + length(discrete_plots)

    fig, axs = subplots(num_plots, 1, sharex=true)
    if num_plots == 1
        axs = [axs]
    end

    reference_aspect_ratio = (length(continuous_plots[1].data) / length(continuous_plots[1].data[1]))

    continuous_imgs = [
        let plot = continuous_plots[i], ax = axs[i]
            ax.set_aspect(length(plot.data) / length(plot.data[1]) / reference_aspect_ratio)
            ax.set_ylabel(plot.label)
            ax.pcolormesh(real(combinedims(plot.data)), cmap=cmap, vmin=0.0, vmax=1.0)
        end
        for i in eachindex(continuous_plots)
    ]
    fig.colorbar(continuous_imgs[end], ax=axs[eachindex(continuous_plots)], aspect=9 * length(continuous_plots))

    discrete_imgs = [
        let plot = discrete_plots[i], ax = axs[i+length(continuous_plots)]
            ax.set_aspect(length(plot.data) / length(plot.data[1]) / reference_aspect_ratio)
            ax.set_ylabel(plot.label)
            vmin = trunc(Int, minimum(minimum.(plot.data)))
            vmax = trunc(Int, maximum(maximum.(plot.data)))
            img = ax.pcolormesh(real(combinedims(plot.data)), cmap=get_cmap(cmap, vmax - vmin + 1), vmin=vmin, vmax=vmax)
            cbar = fig.colorbar(img, ax=ax, aspect=9)
            cbar.ax.locator_params(nbins=min(vmax + 1, 6))
            img.set_clim(-0.5 + vmin, (vmax + 0.5))
        end
        for i in eachindex(discrete_plots)
    ]

    axs[end].set_xlabel("time steps")
    for suffix in file_formats
        savefig("$(path).$(suffix)", bbox_inches="tight")
    end
    PyPlot.close()
    nothing
end