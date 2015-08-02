#!/bin/bash
# this script tests all 6 configurations of Petsc
# this script should be executed from one directory above the tests/ directory
# the build_dir is the location of the deps/ directory, where Petsc will
# be built


test_dir=`pwd`
build_dir="./deps"
echo $test_dir
# run all test configuration
sum=0 

cd $build_dir
./install_petsc.sh
cd $test_dir



# 32 bit integer, double precision real
echo "--- testing 32 bit integer, double precision, real---"
julia --check-bounds=yes ./test/test_doublereal32.jl
sum=$(expr $sum + $?)
echo $sum

julia --check-bounds=yes ./test/runtests.jl
sum=$(expr $sum + $?)
echo $sum

mpirun -np 2 julia --check-bounds=yes ./test/runtests.jl
sum=$(expr $sum + $?)
echo $sum


# 64 bit intger, double precision real
cd $build_dir
./install_petsc.sh --with-64-bit-indices
cd $test_dir

echo "--- testing 64 bit integer, double precision, real---"
julia --check-bounds=yes ./test/test_doublereal64.jl
sum=$(expr $sum + $?)
echo $sum


julia --check-bounds=yes ./test/runtests.jl
sum=$(expr $sum + $?)
echo $sum

mpirun -np 2 julia --check-bounds=yes ./test/runtests.jl
sum=$(expr $sum + $?)
echo $sum



# 32 bit intger, single precision real
cd $build_dir
./install_petsc.sh --with-precision=single
cd $test_dir

echo "--- testing 32 bit integer, single precision, real---"
julia --check-bounds=yes ./test/test_singlereal32.jl
sum=$(expr $sum + $?)
echo $sum


julia --check-bounds=yes ./test/runtests.jl
sum=$(expr $sum + $?)
echo $sum

mpirun -np 2 julia --check-bounds=yes ./test/runtests.jl
sum=$(expr $sum + $?)
echo $sum



# 64 bit intger, single precision real
cd $build_dir
./install_petsc.sh --with-precision=single  --with-64-bit-indices
cd $test_dir

echo "--- testing 64 bit integer, single precision, real---"
julia --check-bounds=yes ./test/test_singlereal64.jl
sum=$(expr $sum + $?)
echo $sum


julia --check-bounds=yes ./test/runtests.jl
sum=$(expr $sum + $?)
echo $sum

mpirun -np 2 julia --check-bounds=yes ./test/runtests.jl
sum=$(expr $sum + $?)
echo $sum



# 32 bit intger, double precision complex
cd $build_dir
./install_petsc.sh --with-scalar-type=complex
cd $test_dir

echo "--- testing 32 bit integer, double precision, complex---"
julia --check-bounds=yes ./test/test_doublecomplex32.jl
sum=$(expr $sum + $?)
echo $sum


julia --check-bounds=yes ./test/runtests.jl
sum=$(expr $sum + $?)
echo $sum

mpirun -np 2 julia --check-bounds=yes ./test/runtests.jl
sum=$(expr $sum + $?)
echo $sum



# 64 bit intger, double precision complex
cd $build_dir
./install_petsc.sh --with-scalar-type=complex --with-64-bit-indices
cd $test_dir

echo "--- testing 64 bit integer, double precision, complex---"
julia --check-bounds=yes ./test/test_doublecomplex64.jl
sum=$(expr $sum + $?)
echo $sum


julia --check-bounds=yes ./test/runtests.jl
sum=$(expr $sum + $?)
echo $sum

mpirun -np 2 julia --check-bounds=yes ./test/runtests.jl
sum=$(expr $sum + $?)
echo $sum


# Petsc does not support single precision complex
## 32 bit intger, single precision complex
#cd $build_dir
#./install_petsc.sh --with-scalar-type=complex --with-precision=single
#cd $test_dir
#
#echo "--- testing 32 bit integer, single precision, complex---"
#julia --check-bounds=yes ./test/test_singlecomplex32.jl
#sum=$(expr $sum + $?)
#
#
#julia --check-bounds=yes ./test/runtests.jl
#sum=$(expr $sum + $?)
#
#
#
## 64 bit intger, single precision complex
#cd $build_dir
#./install_petsc.sh --with-scalar-type=complex --with-64-bit-indices --with-precision=single
#cd $test_dir
#
#echo "--- testing 64 bit integer, single precision, complex---"
#julia --check-bounds=yes ./test/test_singlecomplex64.jl
#sum=$(expr $sum + $?)
#
#
#julia --check-bounds=yes ./test/runtests.jl
#sum=$(expr $sum + $?)
#

echo "exit status sum = $sum"
exit $sum
