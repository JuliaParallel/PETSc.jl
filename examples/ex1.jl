# 1D staggered FD example
using PETSc, MPI, LinearAlgebra, SparseArrays# ForwardDiff

if ~MPI.Initialized()
    MPI.Init()
end

petsclib = PETSc.petsclibs[1];
PETSc.initialize(petsclib)
PetscScalar = PETSc.scalartype(petsclib)
PetscInt = PETSc.inttype(petsclib)
PetscReal = real(PetscScalar)

# ==========================================
# Main code 

comm = MPI.COMM_WORLD

# Compute initial solution
nx   =   10;
x0   =   0;
xend =   1;

# create dmstag for solution and setup
dm  = PETSc.DMStag(petsclib, comm,
                            (PETSc.DM_BOUNDARY_NONE,),
                            (nx,),
                            (1,1),
                            1,
                            PETSc.DMSTAG_STENCIL_BOX)

# create uniform coordinates
LibPETSc.DMStagSetUniformCoordinatesProduct(petsclib, dm, PetscReal(x0), PetscReal(xend), PetscReal(0.0), PetscReal(1.0), PetscReal(0.0), PetscReal(1.0))
#determine boundary type
bnd =   LibPETSc.DMStagGetBoundaryTypes(petsclib, dm); 

a = 1.0; b = 2.0; c = 1.0;
if bnd[1] == PETSc.DM_BOUNDARY_PERIODIC  # bnd is now a tuple, check first element
    b = a;
    c = 0.0;
end

#Compute reference solution on the grid, using direct array access
xa        = PETSc.DMGlobalVec(dm);
xa_Local  = PETSc.DMLocalVec(dm);
xa_array  = LibPETSc.DMStagVecGetArray(petsclib, dm, xa_Local);
dm_coord  = LibPETSc.DMGetCoordinateDM(petsclib, dm)
X_coord,_,_ = LibPETSc.DMStagGetProductCoordinateArrays(petsclib, dm)

corners    =    PETSc.getcorners(dm);

# Get the correct entries for each of our variables in local element-wise storage
iu = LibPETSc.DMStagGetLocationSlot(petsclib, dm, LibPETSc.DMSTAG_LEFT, PetscInt(0));
ip = LibPETSc.DMStagGetLocationSlot(petsclib, dm, LibPETSc.DMSTAG_ELEMENT, PetscInt(0));

# These lines segfault - unclear why
#ixu = LibPETSc.DMStagGetLocationSlot(petsclib, dm_coord, LibPETSc.DMSTAG_LEFT, PetscInt(0));
#ixp = LibPETSc.DMStagGetLocationSlot(petsclib, dm_coord, LibPETSc.DMSTAG_ELEMENT, PetscInt(0));


ixu = iu
ixp = ip

# analytical solution:
xa_array[1:end  ,iu+1] .= a  .+ (b .- a .- (c./2.0)) .* X_coord[1:end,ixu+1] .+ (c./2.0).*X_coord[1:end,ixu+1].*X_coord[1:end,ixu+1];
xa_array[1:end-1,ip+1] .= b .- a .- (c./2.0) .+ c .* X_coord[1:end-1,ixp+1];

LibPETSc.DMStagVecRestoreArray(petsclib, dm, xa_Local, xa_array)
LibPETSc.DMLocalToGlobalBegin(petsclib, dm, xa_Local, LibPETSc.INSERT_VALUES, xa)
LibPETSc.DMLocalToGlobalEnd(petsclib, dm, xa_Local, LibPETSc.INSERT_VALUES, xa)

dmForcing = LibPETSc.DMStagCreateCompatibleDMStag(petsclib, dm, PetscInt(1), PetscInt(0), PetscInt(0), PetscInt(0));
f         = PETSc.DMGlobalVec(dmForcing);
fLocal    = PETSc.DMLocalVec(dmForcing);
f        .= c;
fLocal   .= c;

A   = LibPETSc.DMCreateMatrix(petsclib, dm);
rhs = PETSc.DMGlobalVec(dm);


# Construct rhs vector
for e in corners.lower[1]-1:corners.upper[1]-1

    pos1 = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_ELEMENT, PetscInt(e), PetscInt(0), PetscInt(0), PetscInt(0));
    val1 = PetscScalar(0.0);
    LibPETSc.DMStagVecSetValuesStencil(petsclib, dm, rhs, PetscInt(1), [pos1], [val1], PETSc.INSERT_VALUES);

    pos2 = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_LEFT, PetscInt(e), PetscInt(0), PetscInt(0), PetscInt(0));
    if e == corners.lower[1]-1
        val2 = PetscScalar(a);
    else
        val2 = LibPETSc.DMStagVecGetValuesStencil(petsclib, dmForcing, fLocal, PetscInt(1), [pos2])[1];
    end
    LibPETSc.DMStagVecSetValuesStencil(petsclib, dm, rhs, PetscInt(1), [pos2], [val2], PETSc.INSERT_VALUES);

    if e == corners.upper[1]-1
        pos3 = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_RIGHT, PetscInt(e), PetscInt(0), PetscInt(0), PetscInt(0));
        val3 = PetscScalar(b);
        LibPETSc.DMStagVecSetValuesStencil(petsclib, dm, rhs, PetscInt(1), [pos3], [val3], PETSc.INSERT_VALUES);
    end
end

PETSc.assemble!(rhs)


# Cinstruct matrix A
for e in corners.lower[1]-1:corners.upper[1]-1
    # Note that PETSc ordering is zero-based 
    row = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_LEFT, PetscInt(e), PetscInt(0), PetscInt(0), PetscInt(0));
    if e == corners.lower[1]-1
        val1 = PetscScalar(1.0);
        LibPETSc.DMStagMatSetValuesStencil(petsclib, dm, A, PetscInt(1), [row], PetscInt(1), [row], [val1], PETSc.INSERT_VALUES);
    else
        col1 = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_ELEMENT, PetscInt(e), PetscInt(0), PetscInt(0), PetscInt(0));
        col2 = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_ELEMENT, PetscInt(e-1), PetscInt(0), PetscInt(0), PetscInt(0));

        h   =  X_coord[2,1]-X_coord[1,1];
        
        val1 = PetscScalar(1.0/h);
        val2 = PetscScalar(-1.0/h);
        val3 = PetscScalar(0.0);

        LibPETSc.DMStagMatSetValuesStencil(petsclib, dm, A, PetscInt(1), [row], PetscInt(1), [col1], [val1], PETSc.INSERT_VALUES);
        LibPETSc.DMStagMatSetValuesStencil(petsclib, dm, A, PetscInt(1), [row], PetscInt(1), [col2], [val2], PETSc.INSERT_VALUES);
        LibPETSc.DMStagMatSetValuesStencil(petsclib, dm, A, PetscInt(1), [row], PetscInt(1), [row] , [val3], PETSc.INSERT_VALUES);
    end
    if e == corners.upper[1]-1
        row2 = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_RIGHT, PetscInt(e), PetscInt(0), PetscInt(0), PetscInt(0));
        val4 = PetscScalar(1.0)
        LibPETSc.DMStagMatSetValuesStencil(petsclib, dm, A, PetscInt(1), [row2], PetscInt(1), [row2], [val4], PETSc.INSERT_VALUES);
    end

    row  = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_ELEMENT, PetscInt(e), PetscInt(0), PetscInt(0), PetscInt(0));
    col1 = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_RIGHT, PetscInt(e), PetscInt(0), PetscInt(0), PetscInt(0));
    col2 = LibPETSc.DMStagStencil(LibPETSc.DMSTAG_LEFT, PetscInt(e), PetscInt(0), PetscInt(0), PetscInt(0));

    #xu1 = LibPETSc.DMStagVecGetValuesStencil(petsclib, dm_coord, X_coord, PetscInt(1), [col1])[1];  # FIXME: Need to check if coordinate access is correct
    #xu2 = LibPETSc.DMStagVecGetValuesStencil(petsclib, dm_coord, X_coord, PetscInt(1), [col2])[1];  # FIXME: Need to check if coordinate access is correct
    xu1 = X_coord[2,1]
    xu2 = X_coord[1,1]
    h   = xu1-xu2;


    val1 = PetscScalar(-1.0/h);
    val2 = PetscScalar(1.0/h);
    val3 = PetscScalar(1.0);

    LibPETSc.DMStagMatSetValuesStencil(petsclib, dm, A, PetscInt(1), [row], PetscInt(1), [col1], [val1], PETSc.INSERT_VALUES);
    LibPETSc.DMStagMatSetValuesStencil(petsclib, dm, A, PetscInt(1), [row], PetscInt(1), [col2], [val2], PETSc.INSERT_VALUES);
    LibPETSc.DMStagMatSetValuesStencil(petsclib, dm, A, PetscInt(1), [row], PetscInt(1), [row] , [val3], PETSc.INSERT_VALUES);
end

PETSc.assemble!(A)



x   = PETSc.DMGlobalVec(dm);
ksp = PETSc.KSP(A);
PETSc.solve!(x,ksp,rhs);



xLocal    = PETSc.DMLocalVec(dm);
LibPETSc.DMGlobalToLocalBegin(petsclib, dm, x, LibPETSc.INSERT_VALUES, xLocal)
LibPETSc.DMGlobalToLocalEnd(petsclib, dm, x, LibPETSc.INSERT_VALUES, xLocal)

# get the local solution array:
xsolution_array  = LibPETSc.DMStagVecGetArray(petsclib, dm, xLocal);

# TODO: DMStagGetGhostArrayLocationSlot needs review - may not be available or may need different signature
#xu = LibPETSc.DMStagGetGhostArrayLocationSlot(petsclib, dm, xLocal, PETSc.DMSTAG_LEFT, PetscInt(0));
#xp = LibPETSc.DMStagGetGhostArrayLocationSlot(petsclib, dm, xLocal, PETSc.DMSTAG_ELEMENT, PetscInt(0));

iu = LibPETSc.DMStagGetLocationSlot(petsclib, dm, LibPETSc.DMSTAG_LEFT, PetscInt(0))
ip = LibPETSc.DMStagGetLocationSlot(petsclib, dm, LibPETSc.DMSTAG_ELEMENT, PetscInt(0))
#ip = LibPETSc.DMStagGetGhostArrayLocationSlot(petsclib, dm, xLocal, PETSc.DMSTAG_ELEMENT, PetscInt(0));

#print("u_array = ",xu,"\np_array = ",xp,"\n");

xa_norm    = LinearAlgebra.norm(xa);
error      = xa[:].-x[:]
error_norm = LinearAlgebra.norm(error);
errRel     = error_norm/xa_norm;

print("Error (abs): ",error_norm,"\nError (rel): ",errRel,"\n");

#PETSc.finalize(petsclib)
