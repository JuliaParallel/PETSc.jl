if !haskey(ENV, "PETSC_DIR")  && !haskey(ENV, "PETSC_ARCH")
 run(`./install_petsc.sh`)
end 
