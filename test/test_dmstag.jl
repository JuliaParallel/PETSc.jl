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
dm_1D = PETSc.DMStagCreate1d(MPI.COMM_SELF,PETSc.DM_BOUNDARY_NONE,200,2,2; stag_grid_x=10);
@test PETSc.DMStagGetGlobalSizes(dm_1D) == 10

dm_2D = PETSc.DMStagCreate2d(MPI.COMM_SELF,PETSc.DM_BOUNDARY_NONE,PETSc.DM_BOUNDARY_NONE,20,21,1,1,1,1,1,PETSc.DMSTAG_STENCIL_BOX,2)
@test PETSc.DMStagGetGlobalSizes(dm_2D) == (20, 21)

dm_3D = PETSc.DMStagCreate3d(MPI.COMM_SELF,PETSc.DM_BOUNDARY_NONE,PETSc.DM_BOUNDARY_NONE,PETSc.DM_BOUNDARY_NONE,20,21,22,1,1,1,2,2,2,2,PETSc.DMSTAG_STENCIL_BOX,1,[],[],[])
@test PETSc.DMStagGetGlobalSizes(dm_3D) == (20, 21, 22)

# copy struct
dmnew = PETSc.DMStagCreateCompatibleDMStag(dm_3D,1,1,2,2)
@test PETSc.DMStagGetGlobalSizes(dmnew) == (20, 21, 22)

# Set coordinates 
PETSc.DMStagSetUniformCoordinates(dm_1D, 0, 10)

# retrieve DM with coordinates
#subDM = PETSc.DMProductGetDM(dm, 0)

# retrieve coordinate and value slots
#@test PETSc.DMStagGetProductCoordinateLocationSlot(dm, PETSc.DMSTAG_RIGHT) == 1
#@test PETSc.DMStagGetLocationSlot(dm, PETSc.DMSTAG_RIGHT, 0) ==1

# Create a global and local Vec from the DMStag
vec_test_global     = PETSc.DMCreateGlobalVector(dm_1D)
vec_test            = PETSc.DMCreateLocalVector(dm_1D)
vec_test_2D         = PETSc.DMCreateLocalVector(dm_2D)

# Simply extract an array from the local vector
#x = PETSc.unsafe_localarray(Float64, vec_test.ptr; read=true, write=false)

entriesPerElement = PETSc.DMStagGetEntriesPerElement(dm_1D)

x,m = PETSc.DMStagGetGhostCorners(dm_1D)


# VEC test
# testing how to set values in a local vector:
#
# Note; this test really belongs to a Vec test & should be pushed to a different test file
v       =   rand(10)
v[10]   =   1;
V       =   PETSc.VecSeq(v)
@test V[10] == 1

# VEC test
# create a local Julia array from the vector which we can modify (write=true)
x_local =   PETSc.unsafe_localarray(Float64, V.ptr, write=true);    # create a local array from the vector
x_local[8:10] .= x_local[8:10]*2 .+ 100                             # modify the julia array
finalize(x_local)                                                   # delete local array after local use
@test v[10] == 102                                                  # check

# Note: What I don't understand is that even in the case that we read the array
# as read-only, changing the values in the julia array modifies them in the PetscVec 
# (that seems to defy the purpose of having a read-only option)
#
# In practice this is likely not hugely important; we should simply keep in mind to not 
# change the values locally

# Test retrieving an array from the DMStag:
X = PETSc.DMStagVecGetArray(dm_2D,vec_test_2D);
X[1,1,1] = 1;

# See if DMLocalToGlobal works
vec_test_global .= 0;
vec_test        .= 0;
vec_test[1:end] = 1:length(vec_test);
PETSc.DMLocalToGlobal(dm_1D, vec_test, PETSc.INSERT_VALUES, vec_test_global)
@test vec_test_global[20]==20

# test GlobalToLocal as well. 
# NOTE: as we currently only have VecSeq, parallel halos are not yet tested with this

# Test DMStagVecGetArray for a 1D case
vec_test.array[1:10] = 1:10
X_1D = PETSc.DMStagVecGetArray(dm_1D,vec_test);
@test X_1D[2,3] == 7.0

# Create two stencil locations
pos1 = PETSc.DMStagStencil(PETSc.DMSTAG_LEFT,1,0,0,1)
@test pos1.c == 1
pos2 = PETSc.DMStagStencil(PETSc.DMSTAG_RIGHT,4,0,0,0)
@test pos2.loc == PETSc.DMSTAG_RIGHT
@test pos2.i == 4

# Retrieve value from stencil
val = PETSc.DMStagVecGetValueStencil(dm_1D, vec_test, pos1) # this gets a single value
@test val==6

# Set single value in global vector using stencil
PETSc.DMStagVecSetValueStencil(dm_1D, vec_test_global, pos2, 2222.2, PETSc.INSERT_VALUES)
@test vec_test_global[21] == 2222.2

pos3 = PETSc.DMStagStencil_c(PETSc.DMSTAG_LEFT,1,0,0,1)

# NOTE: setting/getting multiple values is somehow not working for me. Can be called
#  by creating a wrapper
#val = PETSc.DMStagVecGetValuesStencil(dm_1D, vec_test, [pos3; pos3]) 


# Create matrix from dm object, Note: can only be viewed once it is assembled!
A = PETSc.DMCreateMatrix(dm_1D);  # 
@test size(A) == (42,42)
PETSc.assembled(A)

# set some values using normal indices:
A[1,1]= 1.0
A[1,10]= 1.0

# Set values using the DMStagStencil indices
PETSc.DMStagMatSetValueStencil(dm_1D, A, pos1, pos1, 11.1, PETSc.INSERT_VALUES)

# Assemble matrix
PETSc.assemble(A)
@test A[1,10] == 1.0 


# Reads a value from the matrix, using the stencil structure
@test PETSc.DMStagMatGetValueStencil(dm_1D, A, pos1, pos1)==11.1


#PETSc.DMStagVecGetValuesStencil(dm_1D, vec_test.ptr, [pos2]) # this sets a single valu

#PETSc.DMStagVecGetValuesStencil(dm_1D, vec_test.ptr, [pos1; pos2])


#end
