#!/usr/bin/env julia
# Command line entry point for QuantumGameOfLife.jl.
#
# Activates the project environment next to this script, makes sure all dependencies are
# instantiated and then starts the simulation with the given command line arguments.
#
# Usage: julia cli.jl --help
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using QuantumGameOfLife
QuantumGameOfLife.start()
