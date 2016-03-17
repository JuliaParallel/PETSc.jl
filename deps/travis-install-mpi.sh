#!/bin/sh
# this conf file is taken from the mpi4py project
# http://mpi4py.scipy.org/

# because Travis' apt-get is broken, manually walk dependency tree
libnuma_pkg=libnuma1_2.0.9~rc5-1ubuntu2_amd64.deb
libnuma_url=http://mirrors.kernel.org/ubuntu/pool/main/n/numactl/

libhwlock_pkg=libhwloc5_1.8-1ubuntu1_amd64.deb
libhwlock_url=http://mirrors.kernel.org/ubuntu/pool/universe/h/hwloc/

hwlock_pkg=hwloc-nox_1.8-1ubuntu1_amd64.deb
hwlock_url=http://mirrors.kernel.org/ubuntu/pool/universe/h/hwloc/

libmpich_pkg=libmpich10_3.0.4-6ubuntu1_amd64.deb
libmpich_url=http://mirrors.kernel.org/ubuntu/pool/universe/m/mpich/

pkg_name=mpich_3.0.4-6ubuntu1_amd64.deb
url=http://mirrors.kernel.org/ubuntu/pool/universe/m/mpich/

mpichdev_pkg=libmpich-dev_3.0.4-6ubuntu1_amd64.deb
mpichdev_url=http://mirrors.kernel.org/ubuntu/pool/universe/m/mpich/

set -e
os=`uname`
case "$os" in
    Darwin)
        brew update
        brew upgrade cmake
        brew upgrade gcc
        brew install openmpi --build-from-source --verbose |
            sed -e 's/^.*$/./'
        ;;
    Linux)
        case $1 in
          mpich1) set -x;
            sudo apt-get install -q gfortran mpich-shmem-bin libmpich-shmem1.0-dev;;
          mpich2) set -x;
            sudo apt-get install -q gfortran mpich2 libmpich2-3 libmpich2-dev;;
          mpich3) set -x;
            sudo apt-get install -q gfortran libcr0 libcr-dev default-jdk;
            wget $libnuma_url$libnuma_pkg;
            sudo dpkg -i $libnuma_pkg;
            wget $libhwlock_url$libhwlock_pkg;
            sudo dpkg -i ./$libhwlock_pkg;
            wget $hwlock_url$hwlock_pkg;
            sudo dpkg -i ./$hwlock_pkg;
            wget $libmpich_url$libmpich_pkg;
            sudo dpkg -i ./$libmpich_pkg;
            wget $url$pkg_name;
            sudo dpkg -i ./$pkg_name;
            wget $mpichdev_url$mpichdev_pkg;
            sudo dpkg -i ./$mpichdev_pkg;;
          openmpi) set -x;
            sudo apt-get install -q gfortran openmpi-bin openmpi-common libopenmpi-dev;;
          *)
            echo "Unknown MPI implementation:" $1; exit 1;;
        esac
    ;;
    *)
        echo "Unknown operating system: $os"
        exit 1
    ;;
esac
