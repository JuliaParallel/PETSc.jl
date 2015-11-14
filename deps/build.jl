# deps.jl is created at the end of a successful build, so rm
# to ensure that failed builds are missing this file.
if isfile("deps.jl")
  rm("deps.jl")
end

pkg_dict = Pkg.installed()  # get dictionary of installed package names to version numbers

if !haskey(pkg_dict, "MPI")
  Pkg.clone("MPI")
  Pkg.build("MPI")
end

include("build_petscs.jl")
