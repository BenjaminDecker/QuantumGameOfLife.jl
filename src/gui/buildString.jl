str = "--num-cells $(Int(spin.value)) "
str *= "--distance $(Int(spin2.value)) "
str *= "--rule $(Int(spin3.value)) "
str *= "--algorithm $(Gtk4.selected_string(dd)) "
str *= "--num-steps $(Int(spin6.value)) "
str *= "--step-size $(spin7.value) "
str *= "--sweeps-per-time-step $(Int(spin8.value)) "
str *= "--max-bond-dim $(Int(spin9.value)) "
str *= "--plotting-file-path $(entry2.text) "
str *= "--px-per-unit $(spin34.value) "

if spin5.value > 0
    str *= "--activation-interval $(Int(spin4.value)) $(Int(spin5.value)) "
end
if spin35.value > 0
    str *= "--width $(Int(spin35.value)) "
end
if s.active
    str *= "--periodic-boundaries "
end
if s2.active
    str *= "--show "
end
if s3.active
    str *= "--page-entropy "
end
if s4.active
    str *= "--plot-eigval-vs-cbe "
end
if s5.active
    str *= "--plot-fragment-sizes "
end
if s21.active
    str *= "--superposition "
end

str *= "--initial-states "
for button in buttons
    if button.active
        global str
        str *= button.label
        str *= " "
    end
end

str *= "--plot "
for button in plot_buttons
    if button.active
        global str
        str *= button.label
        str *= " "
    end
end

str *= "--file-formats"
for file_format_button in hbox
    global str
    if file_format_button.active
        str *= " $(file_format_button.label)"
    end
end

if entry.text != ""
    str *= " $(entry.text)"
end
