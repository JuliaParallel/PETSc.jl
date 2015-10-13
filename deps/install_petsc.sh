#!/bin/bash
set -ex

#Use all available system cores
np=`getconf _NPROCESSORS_ONLN`

# all arguments are passed to the configure script

petsc_name=$1

# unset the environmental variables if they exist
unset PETSC_DIR
unset PETSC_ARCH

# some useful options
# --with-64-bit-indices=<true, false>
# --with-precision=<single,double>
# --wtih-scalar-type=<real, complex>
#cd ./$petsc_name
./configure $1 $2 $3 > fout
#./configure --with-debugging=0 COPTFLAGS='-O3 -march=native -mtrune=native' CXXOPTFLAGS='-O3 -march=native -mtune=native' FOPTFLAGS='-O3 -march=native -mtune=native' > fout

echo "finished configure"

# get PETSC_ARCH from the printout
PETSC_ARCH=$(cat ./fout | grep "PETSC_ARCH:" | awk '{print $2}')

# get the command printed out on the second to last line
cmd=$(tail -n 2 ./fout | head -n 1)
# execute the command

$cmd MAKE_NP=$np > fout2

echo "finished first command"

cmd2=$(tail -n 2 ./fout2 | head -n 1)
# execute the command

$cmd2 MAKE_NP=$np > fout3


echo "finished second command"


# cmd3 the slashes in cmd3a substitution are causing problems
# for julia, so we skip them
#cmd3=$(tail --lines=1 ./fout3)
#substr='<number of MPI processes you intend to use>'
#nprocs=np
#cmd3a="${cmd3/$substr/$nprocs}"
# execute the command

#$cmd3a | tee fout4


#echo "finished third command"

#export PETSC_DIR=`pwd`
#export PETSC_ARCH

petsc_dir=`pwd`
cd ..
echo "export PETSC_DIR=$petsc_dir" > use_petsc.sh
echo "export PETSC_ARCH=$PETSC_ARCH" >> use_petsc.sh

echo "$petsc_dir" > petsc_evars
echo "$PETSC_ARCH" >> petsc_evars

