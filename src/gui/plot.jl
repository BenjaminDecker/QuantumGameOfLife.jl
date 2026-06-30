g3 = GtkGrid()

label31 = GtkLabel("Show")
s2 = GtkSwitch()
s2.halign = Gtk4.Align_START

label32 = GtkLabel("File Path")
entry2 = GtkEntry()
entry2.text = "plots"

label33 = GtkLabel("File Formats")
b30 = GtkCheckButton("pdf")
b30.active = true
b31 = GtkCheckButton("png")
b32 = GtkCheckButton("svg")
b33 = GtkCheckButton("eps")

label34 = GtkLabel("Px per Unit")
spin34 = GtkSpinButton(1.0, prevfloat(Inf), 0.1)
Gtk4.value(spin34, 2.0)
spin.hexpand = true

label35 = GtkLabel("Width")
spin35 = GtkSpinButton(0, prevfloat(Inf), 100)
Gtk4.value(spin35, 0)
spin.hexpand = true

label36 = GtkLabel("Page Entropy")
s3 = GtkSwitch()
s3.halign = Gtk4.Align_START

label37 = GtkLabel("Eigenvalue vs SBE")
s4 = GtkSwitch()
s4.halign = Gtk4.Align_START

label38 = GtkLabel("Fragment Sizes")
s5 = GtkSwitch()
s5.halign = Gtk4.Align_START

g3[1,1] = label31
g3[2,1] = s2
g3[1,2] = label32
g3[2,2] = entry2

g3[1,3] = label33
hbox = GtkBox(:h)
push!(hbox, b30)
push!(hbox, b31)
push!(hbox, b32)
push!(hbox, b33)
hbox.halign = Gtk4.Align_START
g3[2,3] = hbox

g3[1,4] = label34
g3[2,4] = spin34
g3[1,5] = label35
g3[2,5] = spin35
g3[1,6] = label36
g3[2,6] = s3
g3[1,7] = label37
g3[2,7] = s4
g3[1,8] = label38
g3[2,8] = s5

plots = ["classical", "expect", "sse", "rounded", "bond_dims", "cbe", "autocorrelation", "clustering", "avg_concurrence"]
plot_buttons = []
for (i, plot) in enumerate(plots)
    button = GtkCheckButton(plot)
    push!(plot_buttons, button)
    g3[((i + 17) % 2) + 1, div(i + 17, 2) + 1] = button
end
plot_buttons[1].active = true

g3.column_homogeneous = true
g3.column_spacing = 15
push!(nb, g3, "Plot")
