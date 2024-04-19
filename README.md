# QuantumGameOfLife.jl

### A classical simulation of the Quantum Game of Life.

Simulate and create beautiful plots of quantum cellular automata, inspired by their classical counterpart.
The translation of classical rules into a quantum framework follows the method presented in [[1]](#1).



![](plots/plot.svg)

## Usage
List all available parameters with
```bash
$ julia --project cli.jl --help
```
<br/>

Create and show a plot
```bash
$ julia --project cli.jl --show
```
<br/>

Use different initial state vectors
```bash
$ julia --project cli.jl --show --initial-state triple_blinker
```
<br/>

Use a superposition of initial state vectors
```bash
$ julia --project cli.jl --show --initial-state blinker triple_blinker
```
<br/>

Plot additional measurement information
```bash
$ julia --project cli.jl --show --plot expect sse rounded cbe
```
<br/>

Use different QCA rules
```bash
$ julia --project cli.jl --show --distance 2 --activation-interval 2 3
```
<br/>

Write to different file formats
```bash
$ julia --project cli.jl --show --file-formats svg png pdf
```
<br/>

Plot the classical evolution and mps bond dimension
```bash
$ julia --project cli.jl --show --plot classical expect bond_dims
```
<br/>

Try the TDVP algorithm (This command was used to create the plot at the top)
```bash
$ julia --project cli.jl --show --initial-state single --algorithm tdvp --num-cells 33 --max-bond-dim 5 --num-steps 250 --sweeps-per-time-step 10 --plot classical expect sse rounded --step-size 0.4 --file-formats svg
```

Plots are saved in the plots directory by default, which can be changed with the --plotting-file-path argument. (Make sure to create the specified directory first if it does not already exist.)
```bash
$ julia --project cli.jl --show --plotting-file-path plots2
```
<br/>

## Work with the REPL
Julia uses a just-in-time compiler which takes extra time when code is executed the first time to compile functions before executing them. Subsequent executions will reuse the compiled functions and run a lot faster, even with different input parameters. However, when using the CLI script, the compiled functions are lost between executions and have to be recompiled every time.

To prevent that, you might want to work from inside the julia REPL, especially if you plan to run many quick simulations.

To do that, make sure that your working directory is the project root directory and open the julia REPL
```bash
$ cd QuantumGameOfLife.jl/
$ julia
```
then include the instantiation file and use the project.
```julia
julia> include("instantiate.jl")
julia> using QuantumGameOfLife
```
<br/>

Afterwards, you can use the same command line options as with the CLI by passing them to the start function. To see the effect, compare the runtimes of two executions of the same function.
```julia
julia> @time QuantumGameOfLife.start("--show --distance 2 --activation-interval 2 3")

julia> @time QuantumGameOfLife.start("--show --distance 2 --activation-interval 2 3")

julia> @time QuantumGameOfLife.start("--show --initial-state blinker --file-formats pdf jpg --plot expect sse rounded")
```

<br/>

## References
<a id="1">[1]</a> 
Ney, P. M., Notarnicola, S., Montangero, S., & Morigi, G. (2022). Entanglement in the quantum Game of Life. Physical Review A, 105(1), 012416, DOI: [10.1103/physreva.105.012416](http://dx.doi.org/10.1103/PhysRevA.105.012416)