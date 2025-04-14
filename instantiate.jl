using Pkg
Pkg.activate(".")
ENV["PYTHON"] = ""
Pkg.instantiate()
Pkg.build("PyCall")
using QuantumGameOfLife