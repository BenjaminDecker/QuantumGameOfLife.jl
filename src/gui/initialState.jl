g2 = GtkGrid()

label21 = GtkLabel("Superposition")
s21 = GtkSwitch()
s21.halign = Gtk4.Align_CENTER
g2[1:2,1] = label21
g2[2:3,1] = s21

states = ["blinker", "blinker_wide",
            "triple_blinker", "alternating",
            "alternating_reversed", "single",
            "single_wide", "single_bottom", "single_top",
            "single_bottom_half", "single_top_half",
            "all_ket_0", "all_ket_1",
            "all_ket_0_but_outer", "all_ket_1_but_outer",
            "equal_superposition",
            "equal_superposition_but_outer_ket_0",
            "equal_superposition_but_outer_ket_1",
            "single_bottom_blinker_top", "random",
            "random_product"
]
buttons = []
for (i, state) in enumerate(states)
    button = GtkCheckButton(state)
    push!(buttons, button)
    g2[((i + 2) % 3) + 1, div(i + 2, 3) + 1] = button
end

buttons[1].active = true
g2.column_homogeneous = true
g2.column_spacing = 15
push!(nb, g2, "Initial State")
