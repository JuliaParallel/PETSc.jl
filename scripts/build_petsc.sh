#!/usr/bin/env bash
# Build PETSc as a shared library for `PETSc.set_library!`, linked against the
# MPI library that MPI.jl uses, so the two never disagree on MPI.
#
# Usage:
#   scripts/build_petsc.sh PREFIX [-j N] [-- CONFIGURE_OPTIONS...]
#
#   PREFIX        install directory; the source is unpacked in PREFIX/src
#   -j N          parallel make jobs (default 2)
#   -- ...        extra options passed to PETSc's configure
#
# Environment:
#   PETSC_VERSION  PETSc release to build (default 3.25.4, the version the
#                  wrappers are generated for)
#   JULIA_PROJECT  the Julia environment whose MPI.jl is used (default: this
#                  repository)
#
# The last line printed is the path of the built library, for set_library!.
set -euo pipefail

# ============================================================================
#   Arguments
# ============================================================================

# Print the usage block above and stop.
usage() { sed -n '5,18p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }

[ $# -ge 1 ] || usage
prefix=$(realpath -m "$1"); shift
jobs=2
extra=()
while [ $# -gt 0 ]; do
    case "$1" in
        -j) jobs="$2"; shift 2 ;;
        --) shift; extra=("$@"); break ;;
        *) usage ;;
    esac
done

version="${PETSC_VERSION:-3.25.4}"
repo=$(cd "$(dirname "$0")/.." && pwd)
project="${JULIA_PROJECT:-$repo}"

# Run a Julia snippet in the chosen project and print its output.
jl() { julia --startup-file=no --project="$project" -e "$1"; }

# ============================================================================
#   Libraries from the Julia side
# ============================================================================

# MPI: the library MPI.jl loads. Its headers sit in ../include next to it.
libmpi=$(jl 'using MPI; print(MPI.API.libmpi)')
mpi_dir=$(dirname "$(dirname "$libmpi")")
[ -f "$mpi_dir/include/mpi.h" ] || { echo "no mpi.h in $mpi_dir/include (MPI library $libmpi)" >&2; exit 1; }

# Directories of the artifacts libmpi depends on (hwloc, …). The linker needs
# them through -rpath-link; at run time Julia has them loaded already.
link_dirs=$(jl 'using MPI, Libdl
    print(join(unique(dirname.(filter(contains("artifacts"), Libdl.dllist()))), ":"))')

# BLAS and LAPACK: Julia's libblastrampoline, through its ILP64 names
# (dgemm_64_, …), as in PETSc_jll. A BLAS of PETSc's own would export dgemm_,
# and inside Julia those calls bind to libblastrampoline's LP64 entry points,
# which have no library behind them.
lbt=$(jl 'using Libdl; print(realpath(only(filter(contains("libblastrampoline"), Libdl.dllist()))))')

# ============================================================================
#   Source
# ============================================================================

# Download and unpack the release once; later runs reuse it.
src="$prefix/src/petsc-$version"
if [ ! -d "$src" ]; then
    mkdir -p "$prefix/src"
    curl -fsSL "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-$version.tar.gz" |
        tar -xz -C "$prefix/src"
fi

# ============================================================================
#   Configure options
# ============================================================================

opts=(
    # Optimized shared library with Int64 indices, as PETSc.jl expects
    --prefix="$prefix"
    --with-shared-libraries=1
    --with-debugging=0
    --with-64-bit-indices=1
    --with-fc=0
    COPTFLAGS=-O2 CXXOPTFLAGS=-O2

    # MPI from MPI.jl
    --with-mpi-include="$mpi_dir/include"
    --with-mpi-lib="$libmpi"
    LDFLAGS="-Wl,-rpath-link,$link_dirs"

    # BLAS and LAPACK from libblastrampoline
    --with-blaslapack-lib="$lbt"
    --with-blaslapack-suffix=_64_
    --known-64-bit-blas-indices=1

    # Mesh generators for simplex meshes in the DMPlex tests
    --download-triangle=1
    --download-ctetgen=1

    # CPU only: configure enables CUDA or HIP on its own when it finds them
    --with-cuda=0
    --with-hip=0
)

# ============================================================================
#   Build and install
# ============================================================================

cd "$src"
./configure "${opts[@]}" "${extra[@]}"
make -j "$jobs" all
make install

echo "$prefix/lib/libpetsc.so"
