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

# DOF
@test PETSc.DMStagGetDOF(dm) == (2,2)

# Destroy
PETSc.destroy(dm)

# Create new struct and pass keyword arguments
dm = PETSc.DMStagCreate1d(MPI.COMM_SELF,PETSc.DM_BOUNDARY_NONE,200,2,2; stag_grid_x=10);
@test PETSc.DMStagGetGlobalSizes(dm) == 10

dm_2D = PETSc.DMStagCreate2d(MPI.COMM_SELF,PETSc.DM_BOUNDARY_NONE,PETSc.DM_BOUNDARY_NONE,20,21,1,1,1,1,1,PETSc.DMSTAG_STENCIL_BOX,2)
@test PETSc.DMStagGetGlobalSizes(dm_2D) == (20, 21)

dm_3D = PETSc.DMStagCreate3d(MPI.COMM_SELF,PETSc.DM_BOUNDARY_NONE,PETSc.DM_BOUNDARY_NONE,PETSc.DM_BOUNDARY_NONE,20,21,22,1,1,1,2,2,2,2,PETSc.DMSTAG_STENCIL_BOX,1,[],[],[])
@test PETSc.DMStagGetGlobalSizes(dm_3D) == (20, 21, 22)

dmnew = PETSc.DMStagCreateCompatibleDMStag(dm_3D,1,1,2,2)
@test PETSc.DMStagGetGlobalSizes(dmnew) == (20, 21, 22)

# Set coordinates 
PETSc.DMStagSetUniformCoordinates(dm, 0, 10)

# retrieve DM with coordinates
#subDM = PETSc.DMProductGetDM(dm, 0)

# retrieve coordinate and value slots
#@test PETSc.DMStagGetProductCoordinateLocationSlot(dm, PETSc.DMSTAG_RIGHT) == 1
#@test PETSc.DMStagGetLocationSlot(dm, PETSc.DMSTAG_RIGHT, 0) ==1

# Create a global and local Vec from the DMStag
vec_test_global     = PETSc.DMCreateGlobalVector(dm)
vec_test            = PETSc.DMCreateLocalVector(dm)
vec_test_2D         = PETSc.DMCreateLocalVector(dm_2D)

# Simply extract an array from the local vector
#x = PETSc.unsafe_localarray(Float64, vec_test.ptr; read=true, write=false)

entriesPerElement = PETSc.DMStagGetEntriesPerElement(dm)

x,m = PETSc.DMStagGetGhostCorners(dm)


# testing how to set values in a local vector:
#
# Note; this test really belongs to a Vec test & should be pushed to a different test file
X = rand(10)
V = PETSc.VecSeq(X)

# create a local Julia array from the vector which we can modify (write=true)
x_local = PETSc.unsafe_localarray(Float64, V.ptr);  # create a local array

x_local[8:10] .= x_local[8:10]*2 .+ 100       # modify the julia array

finalize(x_local)                             # delete local array after local use

V   # this correctly shows the modified array values in the vector

# What I don't understand is that even in the case that we read the array
# as read-only, changing the values in the julia array modifies them in the PetscVec 
# (that seems to defy the purpose of having a read-only option)
#
# In practice this is likely not hugely important; we should simply keep in mind to not 
# change the values locally


# Test retrieving an array from the DMStag:
X = PETSc.DMStagVecGetArray(dm_2D,vec_test_2D);




#PETSc.DMStagVecGetValuesStencil(dm, vec_test_global.ptr, [pos], [12.0], PETSc.INSERT_VALUES)

vec_test.array[1:10] = 1:10

pos1 = PETSc.DMStagStencil_c(PETSc.DMSTAG_RIGHT,3,0,0,1)

#PETSc.DMStagVecGetValueStencil(dm, vec_test.ptr, pos1) # this sets a single value
PETSc.DMStagVecGetValueStencil(dm, vec_test.ptr, PETSc.DMStagStencil_c(PETSc.DMSTAG_RIGHT,3,0,0,1)) # this sets a single value

X_1D = PETSc.DMStagVecGetArray(dm,vec_test);




#PETSc.DMStagVecGetValuesStencil(dm, vec_test.ptr, [pos1; pos2])


#PETSc.DMStagGetProductCoordinateArrays(dm)
