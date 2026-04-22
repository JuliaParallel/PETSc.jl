#!/bin/bash
# submit_scaling.sh
# Usage: ./submit_scaling.sh <ntasks> <Nx> <Ny> <Nz> <mg_levels> [ncoarse]
#
# Arguments:
#   ntasks     - total MPI tasks (default: 8)
#   Nx/Ny/Nz   - global grid size in each dimension (default: 129)
#   mg_levels  - number of multigrid levels (default: 4)
#   ncoarse    - number of ranks for coarse solve (default: 16)
#
# Coarse grid size guidance:
#   mg_levels=4, fine=129^3  -> coarse ~16^3 (~4096 DOF)
#   mg_levels=5, fine=257^3  -> coarse ~17^3 (~4913 DOF)
#
# Examples:
#   ./submit_scaling.sh 64  257 257 257 5 16
#   ./submit_scaling.sh 512 1025 1025 1025 6 16
#   ./submit_scaling.sh 4096 2049 2049 2049 7 16

NTASKS=${1:-8}
NX=${2:-129}
NY=${3:-129}
NZ=${4:-129}
NMGLEVELS=${5:-4}
NCOARSE=${6:-16}

if [ ${NCOARSE} -gt ${NTASKS} ]; then
  ec
  ho "ERROR: NCOARSE=${NCOARSE} > NTASKS=${NTASKS}"
  exit 1
fi

REDUCTION_FACTOR=$((NTASKS / NCOARSE))

NTASKS_PER_NODE=64
NNODES=$(( (NTASKS + NTASKS_PER_NODE - 1) / NTASKS_PER_NODE ))
ACTUAL_NTASKS_PER_NODE=$(( NTASKS < NTASKS_PER_NODE ? NTASKS : NTASKS_PER_NODE ))

echo "Submitting: ntasks=${NTASKS}, nodes=${NNODES}, ntasks-per-node=${ACTUAL_NTASKS_PER_NODE}"
echo "            Nx=${NX}, Ny=${NY}, Nz=${NZ}, mg_levels=${NMGLEVELS}"
echo "            ncoarse=${NCOARSE}, reduction_factor=${REDUCTION_FACTOR}"

sbatch \
  --ntasks=${NTASKS} \
  --nodes=${NNODES} \
  --ntasks-per-node=${ACTUAL_NTASKS_PER_NODE} \
  --output="scaling_%j_n${NTASKS}_nx${NX}_ny${NY}_nz${NZ}_mg${NMGLEVELS}_nc${NCOARSE}.out" \
  --export=ALL,NTASKS=${NTASKS},NX=${NX},NY=${NY},NZ=${NZ},NMGLEVELS=${NMGLEVELS},NCOARSE=${NCOARSE},REDUCTION_FACTOR=${REDUCTION_FACTOR} \
  job.sh