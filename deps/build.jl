
pkg_dict = Pkg.installed()  # get dictionary of installed package names to version numbers

if !haskey(pkg_dict, "MP")
  Pkg.clone("MPI")
  Pkg.build("MPI")
end

if !haskey(ENV, "PETSC_DIR")  && !haskey(ENV, "PETSC_ARCH")
 run(`./install_petsc.sh`)
end 
