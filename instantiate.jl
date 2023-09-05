ENV["PYTHON"] = ""
using Pkg
Pkg.activate(".")
Pkg.instantiate()