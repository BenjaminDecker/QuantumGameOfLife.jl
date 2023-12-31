# QuantumGameOfLife.jl
![](plots/plot.svg)

## Usage
List all available parameters with
```bash
julia --project cli.jl --help
```
<br/>

Create and show a plot
```bash
julia --project cli.jl --show --initial-states blinker
```
<br/>

Use different initial state vectors
```bash
julia --project cli.jl --show --initial-states triple_blinker random
```
<br/>

Plot different measurement information
```bash
julia --project cli.jl --show --initial-states blinker --plot-sse --plot-rounded --plot-cbe
```
<br/>

Use different rules
```bash
--julia --project cli.jl --show --initial-states blinker --distance 2 --activation-interval 2 4
```
<br/>

Write to different file formats
```bash
julia --project cli.jl --show --initial-states blinker --file-formats svg png pdf
```
<br/>

<!-- Plot the classical evolution and mps bond dimension
```bash
julia --project cli.jl --show --initial-states blinker --plot-classical --plot-bond-dims
```
<br/>

Try the TDVP algorithm (This can take a while)
```bash
julia --project cli.jl --show --initial-states blinker --algorithm 2tdvp --num-steps 1000 --plotting-frequency 10 --plot-bond-dims --num-cells 15
``` -->

Plots are saved in the plots directory by default, which can be changed with the --plot-file-path argument (Make sure to create the specified directory if it does not already exist)
```bash
julia --project cli.jl --show --initial-states blinker --plot-file-path plots2
```
<br/>

The plot at the top was created using the following command
```bash
julia --project cli.jl --show --initial-states single --plot-sse --num-steps 100 --num-cells 13 --plot-rounded --plot-cbe --file-formats svg
```
<br/>

## Work with the REPL
Julia uses a just-in-time compiler which takes extra time when the code is executed the first time. Subsequent executions will reuse the compiled functions and run a lot faster, even with different input parameters. However, when using the cli script, the compiled functions are lost between executions and have to be recompiled every time.

To prevent that, you might want to work from inside the julia REPL, especially if you plan to create many quick plots.
To do that, first open the julia REPL
```bash
julia
```
then include the instantiation file and use the project.
```julia
include("instantiate.jl")
using QuantumGameOfLife
```
<br/>

Afterwards, you can use the same command line options as with the cli by passing them to the start function.
```julia
QuantumGameOfLife.start("--show --initial-states blinker --file-formats pdf jpg --plot-sse --plot-rounded")
```

<br/>

---
> **NOTE**  
> This is a new and improved version of [BenjaminDecker/quantum-cellular-automaton](https://github.com/BenjaminDecker/quantum-cellular-automaton), re-written in the Julia language.
