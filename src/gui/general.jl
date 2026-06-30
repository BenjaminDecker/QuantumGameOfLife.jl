g = GtkGrid()

label = GtkLabel("Cells")
spin = GtkSpinButton(1.0, prevfloat(Inf), 1.0)
Gtk4.value(spin, 9.0)
spin.hexpand = true

label2 = GtkLabel("Distance")
spin2 = GtkSpinButton(1.0, prevfloat(Inf), 1.0)
spin2.hexpand = true

label3 = GtkLabel("Rule")
spin3 = GtkSpinButton(0, prevfloat(Inf), 1.0)
spin3.hexpand = true
Gtk4.value(spin3, 150.0)

label4 = GtkLabel("Activation Interval")
spin4 = GtkSpinButton(0, prevfloat(Inf), 1.0)
spin4.hexpand = true
spin5 = GtkSpinButton(0, prevfloat(Inf), 1.0)
spin5.hexpand = true
hbox2 = GtkBox(:h)
push!(hbox2, spin4)
push!(hbox2, spin5)

label5 = GtkLabel("Algorithm")
algs = ["exact", "tdvp1", "tdvp2", "sierpinski"]
dd = GtkDropDown(algs)

label6 = GtkLabel("Steps")
spin6 = GtkSpinButton(0, prevfloat(Inf), 10.0)
spin6.hexpand = true
Gtk4.value(spin6, 100.0)

label7 = GtkLabel("Step Size")
spin7 = GtkSpinButton(0, prevfloat(Inf), 0.25)
spin7.hexpand = true
spin7.digits = 2
Gtk4.value(spin7, 1.0)

label8 = GtkLabel("Sweeps per Time Step")
spin8 = GtkSpinButton(1, prevfloat(Inf), 10.0)
spin8.hexpand = true
Gtk4.value(spin8, 100.0)

label9 = GtkLabel("Max Bond Dimension")
spin9 = GtkSpinButton(1, prevfloat(Inf), 1.0)
spin9.hexpand = true
Gtk4.value(spin9, 32.0)

label10 = GtkLabel("Periodic Boundaries")
s = GtkSwitch()
s.halign = Gtk4.Align_START

label11 = GtkLabel("Addional Parameters:")
label11.halign = Gtk4.Align_START
entry = GtkEntry()

g[1,1] = label
g[2,1] = spin
g[1,2] = label2
g[2,2] = spin2
g[1,3] = label3
g[2,3] = spin3
g[1,4] = label4
g[2,4] = hbox2
g[1,5] = label5
g[2,5] = dd
g[1,6] = label6
g[2,6] = spin6
g[1,7] = label7
g[2,7] = spin7
g[1,8] = label8
g[2,8] = spin8
g[1,9] = label9
g[2,9] = spin9
g[1,10] = label10
g[2,10] = s
g[1:2,11] = label11
g[1:2,12] = entry


g.column_homogeneous = true
g.column_spacing = 15
push!(nb, g, "General")
