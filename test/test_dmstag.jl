using Test
using PETSc, MPI

if ~MPI.Initialized()
    MPI.Init()
end
PETSc.initialize()

@testset "DMSTAG routines" begin

# Create 1D DMStag
dm = PETSc.DMStagCreate1d(MPI.COMM_SELF,PETSc.DM_BOUNDARY_NONE,20,2,2,PETSc.DMSTAG_STENCIL_BOX,2,[])


# Test that we can retrieve sizes again
M,N,P = PETSc.DMStagGetGlobalSizes(dm);
@test PETSc.DMStagGetGlobalSizes(dm) == (20,0,0)



end

#PETSc.finalize()