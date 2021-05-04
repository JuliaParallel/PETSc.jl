using Test
using PETSc, MPI

if ~MPI.Initialized()
    MPI.Init()
end
PETSc.initialize()

#@testset "DMSTAG routines" begin

# Create 1D DMStag
dm = PETSc.DMStagCreate1d(MPI.COMM_SELF,PETSc.DM_BOUNDARY_NONE,20,2,2,PETSc.DMSTAG_STENCIL_BOX,2,[])

# Test get size
@test PETSc.DMStagGetGlobalSizes(dm) == (20,0,0)

# Test gettype
@test PETSc.gettype(dm) == "stag"               

# Boundary
@test PETSc.DMStagGetBoundaryTypes(dm)==PETSc.DM_BOUNDARY_NONE

# Corners
@test PETSc.DMStagGetCorners(dm) == (0,0,0)

# Destroy
PETSc.destroy(dm)

# Create new struc
dm = PETSc.DMStagCreate1d(MPI.COMM_SELF,PETSc.DM_BOUNDARY_PERIODIC,200,2,2,PETSc.DMSTAG_STENCIL_BOX,2,[])


#end

#PETSc.finalize()