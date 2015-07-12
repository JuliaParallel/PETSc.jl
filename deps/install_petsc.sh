#1/bin/sh
set -ex
petsc_name=petsc-3.6.0
fmt=.tar.gz
wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/$petsc_name$fmt
tar xfz ./$petsc_name$fmt

#export PETSC_DIR=`pwd`/$petsc_name
#export PETSC_ARCH=
unset PETSC_DIR
unset PETSC_ARCH

cd ./$petsc_name
./configure | tee fout

echo "finished configure"

# get PETSC_ARCH from the printout
PETSC_ARCH=$(cat ./fout | grep "PETSC_ARCH:" | awk '{print $2}')

# get the command printed out on the second to last line
cmd=$(tail --lines=2 ./fout | head --lines=1)
# execute the command

$cmd | tee fout2

echo "finished first command"

cmd2=$(tail --lines=2 ./fout2 | head --lines=1)
# execute the command

$cmd2 | tee fout3


echo "finished second command"


# cmd3 the slashes in cmd3a substitution are causing problems
# for julia, so we skip them
#cmd3=$(tail --lines=1 ./fout3)
#substr='<number of MPI processes you intend to use>'
#nprocs=4
#cmd3a="${cmd3/$substr/$nprocs}"
# execute the command

#$cmd3a | tee fout4


#echo "finished third command"

#export PETSC_DIR=`pwd`
#export PETSC_ARCH

petsc_dir=`pwd`
cd ..
echo "export PETSC_DIR=$petsc_dir" >> use_petsc.sh
echo "export PETSC_ARCH=$PETSC_ARCH" >> use_petsc.sh

echo "$petsc_dir" >> petsc_evars
echo "$PETSC_ARCH" >> petsc_evars




