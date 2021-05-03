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

# Destroy
PETSc.destroy(dm)

# Destroy
dm = PETSc.DMStagCreate1d(MPI.COMM_SELF,PETSc.DM_BOUNDARY_NONE,200,2,2,PETSc.DMSTAG_STENCIL_BOX,2,[])


#end

#PETSc.finalize()