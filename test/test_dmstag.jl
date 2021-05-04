using Test
using PETSc, MPI

if ~MPI.Initialized()
    MPI.Init()
end
PETSc.initialize()

#@testset "DMSTAG routines" begin

# Create 1D DMStag
dm = PETSc.DMStagCreate1d(MPI.COMM_SELF,PETSc.DM_BOUNDARY_NONE,20,2,2,PETSc.DMSTAG_STENCIL_BOX,2)
PETSc.destroy(dm)

# Create 1D DMStag with array of local @ of points
dm = PETSc.DMStagCreate1d(MPI.COMM_SELF,PETSc.DM_BOUNDARY_NONE,20,2,2,PETSc.DMSTAG_STENCIL_BOX,2,[20])

# Test get size
@test PETSc.DMStagGetGlobalSizes(dm) == 20

# Test gettype
@test PETSc.gettype(dm) == "stag"               

# Boundary
@test PETSc.DMStagGetBoundaryTypes(dm)==PETSc.DM_BOUNDARY_NONE

# Corners
@test PETSc.DMStagGetCorners(dm) == (0,20,1)

# Destroy
err = PETSc.destroy(dm)

# Create new struct and pass keyword arguments
dm = PETSc.DMStagCreate1d(MPI.COMM_SELF,PETSc.DM_BOUNDARY_PERIODIC,200,2,2; stag_grid_x=199);
@test PETSc.DMStagGetGlobalSizes(dm) == 199

# Set coordinates
PETSc.DMStagSetUniformCoordinates(dm, 0, 10)

# retrieve coordinate and value slots
@test PETSc.DMStagGetProductCoordinateLocationSlot(dm, PETSc.DMSTAG_RIGHT) == 2
@test PETSc.DMStagGetLocationSlot(dm, PETSc.DMSTAG_RIGHT, 0) ==4


vec_test = PETSc.DMCreateGlobalVector(dm)


#r_pv = Ref{Ptr{$PetscScalar}}()

#end

#PETSc.finalize()



#x = randn(11)
#V = PETSc.VecSeq(x)
