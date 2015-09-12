
pkg_dict = Pkg.installed()  # get dictionary of installed package names to version numbers

if !haskey(pkg_dict, "MPI")
  Pkg.clone("MPI")
  Pkg.build("MPI")
end

include("build_petscs.jl")
