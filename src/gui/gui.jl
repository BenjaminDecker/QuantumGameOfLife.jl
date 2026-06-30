using Gtk4
using QuantumGameOfLife

gui_done = Condition()

win = GtkWindow("Parameters")
vbox = GtkBox(:v)
nb = GtkNotebook()
include("general.jl")
include("initialState.jl")
include("plot.jl")
push!(vbox, nb)
run = GtkButton("Run")
push!(vbox, run)
push!(win, vbox)

signal_connect(run, "clicked") do widget
    global str
    run.sensitive = false
    run.label = "Running..."

    sleep(0.1)
    @async begin
        sleep(0.1)
        include("src/gui/buildString.jl") 
        println("Starting simulation: ", str)

        try
            QuantumGameOfLife.start(str)
        catch err
            bt = catch_backtrace()
            println()
            showerror(stderr, err, bt)
        end
        run.sensitive = true
        run.label = "Run"
    end
end

signal_connect(win, "close-request") do widget
    notify(gui_done)
    return false
end

show(win)
wait(gui_done)
destroy(win)
