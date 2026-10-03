#!/usr/bin/env julia
# Command line entry point for QuantumGameOfLife.jl.
#
# Activates the project environment next to this script, makes sure all dependencies are
# instantiated and then starts the simulation from a TOML config file.
#
# Usage:
#   julia cli.jl run.toml     # start a run defined by a config file
#   julia cli.jl --help       # start a run with command line options instead; this also
#                             # writes the effective configuration to a timestamped config
#                             # file so the exact same run can be started later
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using QuantumGameOfLife
QuantumGameOfLife.start()
