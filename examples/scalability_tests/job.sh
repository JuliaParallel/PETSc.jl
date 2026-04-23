#!/bin/bash -l
#SBATCH --job-name=scaling
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=200000M
#SBATCH --time=0-00:10:00
#SBATCH --account=??? # add your account here
# --ntasks, --nodes, --ntasks-per-node, --output passed by submit_scaling.sh

module --force purge
module load CrayEnv
module load craype
module load gcc-native/13.2
module load cray-mpich/8.1.32

export HOME=/users/kausbori
export JULIA_DEPOT_PATH=/users/kausbori/.julia
export MPITRAMPOLINE_LIB=/users/kausbori/mpiwrapper/lib64/libmpiwrapper.so
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_THREADS=1
export JULIA_CPU_TARGET="generic"

export MPICH_OFI_STARTUP_CONNECT=1
export MPICH_OFI_NIC_POLICY=NUMA
export FI_CXI_DEFAULT_CQ_SIZE=131072
export MPICH_OFI_STARTUP_CONNECT=0
export FI_CXI_DISABLE_CQ_HUGETLB=1
export MPICH_MAX_THREAD_SAFETY=single
export LD_LIBRARY_PATH=/opt/cray/pe/mpich/8.1.32/ofi/cray/17.0/lib:$LD_LIBRARY_PATH

JULIA=/users/kausbori/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/bin/julia

echo "=== Scaling run: ntasks=${NTASKS}, Nx=${NX}, Ny=${NY}, Nz=${NZ}, mg_levels=${NMGLEVELS}, ncoarse=${NCOARSE}, reduction_factor=${REDUCTION_FACTOR} ==="
echo "=== Running on ${SLURM_NNODES} nodes, ${SLURM_NTASKS} tasks, $(( SLURM_NTASKS / SLURM_NNODES )) tasks/node ==="

echo "=== MPI diagnostic ==="
srun -n 1 $JULIA --project=/users/kausbori/PETSc_jl_scalability -e '
using MPI
MPI.Init()
println("MPI library: ", MPI.libmpi)
println("MPI version: ", MPI.MPI_LIBRARY_VERSION_STRING)
println("hostname: ", gethostname())
'
echo "=== End MPI diagnostic ==="

echo "MPI run started at $(date)"
srun -n ${NTASKS} $JULIA --heap-size-hint=1G \
    --project=/users/kausbori/PETSc_jl_scalability ex45.jl \
  -Nx ${NX} \
  -Ny ${NY} \
  -Nz ${NZ} \
  -pc_mg_levels ${NMGLEVELS} \
  -ksp_type cg \
  -pc_type mg \
  -mg_levels_ksp_type chebyshev \
  -mg_levels_pc_type sor \
  -ksp_rtol 1e-8 \
  -mg_coarse_ksp_type preonly \
  -mg_coarse_pc_type telescope \
  -mg_coarse_telescope_reduction_factor ${REDUCTION_FACTOR} \
  -mg_coarse_telescope_pc_type lu \
  -mg_coarse_telescope_pc_factor_mat_solver_type superlu_dist \
  -pc_mg_log \
  -malloc_view \
  -memory_view \
  -ksp_converged_reason \
  -ksp_monitor_true_residual \
  -log_view \
  2>&1
echo "MPI run finished at $(date)"

echo "=== Memory usage ==="
sacct -j ${SLURM_JOB_ID} --format=JobID,MaxRSS,AveRSS,NNodes,NCPUs --units=G

echo "exit: $?"
echo "Job finished at $(date)"