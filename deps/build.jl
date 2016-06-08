# deps.jl is created at the end of a successful build, so rm
# to ensure that failed builds are missing this file.
if isfile("deps.jl")
  rm("deps.jl")
end

Pkg.checkout("MPI", "master")
Pkg.build("MPI")
include("build_petscs.jl")
