# QuantumGameOfLife.jl

### A classical simulation of the Quantum Game of Life.

Simulate and create beautiful plots of quantum cellular automata, inspired by classical elementary cellular automata.
Given some rule from the Wolfram Code, the translation into a quantum framework follows the method presented in [[1]](#1).



![](plots/plot.svg)

## Setup
Before running the code for the first time, install all dependencies by running the instantiation script from the project directory

```bash
$ cd QuantumGameOfLife.jl/
$ julia instantiate.jl
```

The instantiation script is also called from the cli script to make sure everything is set up correctly every time.

## Usage
List all available parameters with
```bash
$ julia cli.jl --help
```
<br/>

Create and show a plot
```bash
$ julia cli.jl --show
```
<br/>

Use different initial state vectors
```bash
$ julia cli.jl --show --initial-state blinker triple_blinker
```
<br/>

Use a superposition of initial state vectors
```bash
$ julia cli.jl --show --initial-state blinker triple_blinker --superposition
```
<br/>

Plot additional measurement information
```bash
$ julia cli.jl --show --plot expect sse rounded cbe autocorrelation
```
<br/>

Write to different file formats
```bash
$ julia cli.jl --show --file-formats svg png pdf
```
<br/>

Try the TDVP algorithm (This command was used to create the plot at the top)
```bash
$ julia cli.jl --show --initial-state single --algorithm tdvp1 --num-cells 33 --max-bond-dim 5 --num-steps 350 --sweeps-per-time-step 10 --plot classical expect sse rounded autocorrelation cbe --step-size 0.4 --file-formats svg
```

Plots are saved in the plots directory by default, which can be changed with the --plotting-file-path argument. (Make sure to create the specified directory first if it does not already exist.)
```bash
$ julia cli.jl --show --plotting-file-path plots2
```
<br/>

### Use different QCA rules
Provide a Wolfram code via the rule parameter to create quantum analogues to classical elementary cellular automata and compare the results to the classical case.
Use a distance of 1 to get equivalents to elementary cellular automata where the rule labels the corresponding Wolfram code.
```bash
$ julia cli.jl --show --plot classical expect --distance 1 --rule 108
$ julia cli.jl --show --plot classical expect --distance 1 --rule 30
$ julia cli.jl --show --plot classical expect --distance 1 --rule 150
```
<br/>

Use higher distances and extend the Wolfram code to larger numbers
```bash
$ julia cli.jl --show --plot classical expect --distance 2 --rule 2266898040
```
<br/>

Use periodic boundary conditions
```bash
$ julia cli.jl --show --periodic-boundaries
```
<br/>

## Work with the REPL
Julia uses a just-in-time compiler which takes extra time when code is executed the first time to compile functions before executing them. Subsequent executions will reuse the compiled functions and run a lot faster, even with different input parameters. However, when using the CLI script, the compiled functions are lost between executions and have to be recompiled every time.

To prevent that, you might want to work from inside the julia REPL, especially if you plan to run many quick simulations.

To do so, make sure that your working directory is the project root directory and open the julia REPL
```bash
$ cd QuantumGameOfLife.jl/
$ julia
```
then, include the instantiation file and use the project.
```julia
julia> include("instantiate.jl")
julia> using QuantumGameOfLife
```
<br/>

Afterwards, you can use the same command line options as with the CLI by passing them to the start function. To see the effect, compare the runtimes of two consecutive executions of the same function.
```julia
julia> @time QuantumGameOfLife.start()

julia> @time QuantumGameOfLife.start()

julia> @time QuantumGameOfLife.start("--show --plot classical expect --distance 1 --rule 150")
```

<br/>

## References
<a id="1">[1]</a> 
Benjamin Decker. “Creating a Quantum Analogue to an Arbitrary Classical Elementary Cellular Automaton”. en. MA thesis. Technical University of Munich, 2024-10. URL: https://mediatum.ub.tum.de/1756463
