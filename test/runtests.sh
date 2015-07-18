#!/bin/bash

# run all test configuration
julia --check-bounds=yes -e 'Pkg.clone(pwd()); Pkg.build("PETSc")'

sum=0 
julia --check-bounds=yes ./test/runtests.jl
sum=$(expr $sum + $?)

julia --check-bounds=yes ./test/test_doublereal.jl
sum=$(expr $sum + $?)

exit $sum
