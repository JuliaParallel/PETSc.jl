using PETSc, MPI, LinearAlgebra, SparseArrays, Plots, ForwardDiff

if ~MPI.Initialized()
    MPI.Init()
end

PETSc.initialize()

# ==========================================
# Main code 


# Compute initial solution
nx   =   3;
x0   =   0;
xend =   1;

# create dmstag for solution and setup
dm  = PETSc.DMStagCreate1d(MPI.COMM_SELF,PETSc.DM_BOUNDARY_NONE,nx,1,1,PETSc.DMSTAG_STENCIL_BOX,1);
# creat uniform coordinates
PETSc.DMStagSetUniformCoordinatesExplicit(dm, x0, xend);
#determine boundary type
bnd =   PETSc.DMStagGetBoundaryTypes(dm);

a = 1.0; b = 2.0; c = 1.0;
if bnd == PETSc.DM_BOUNDARY_PERIODIC
    b = a;
    c = 0.0;
end

#Compute reference solution on the grid, using direct array access
xa        = PETSc.DMCreateGlobalVector(dm);
xa_Local  = PETSc.DMCreateLocalVector(dm);
xa_array  = PETSc.DMStagVecGetArray(dm,xa_Local);
dm_coord  = PETSc.DMGetCoordinateDM(dm);
vec_coord = PETSc.DMGetCoordinatesLocal(dm);
X_coord = PETSc.DMStagVecGetArray(dm_coord, vec_coord);
start,n,nExtra = PETSc.DMStagGetCorners(dm);

# Get the correct entries for each of our variables in local element-wise storage
iu = PETSc.DMStagGetLocationSlot(dm, PETSc.DMSTAG_LEFT, 0);
ip = PETSc.DMStagGetLocationSlot(dm, PETSc.DMSTAG_ELEMENT, 0);
ixu = PETSc.DMStagGetLocationSlot(dm_coord, PETSc.DMSTAG_LEFT, 0);
ixp = PETSc.DMStagGetLocationSlot(dm_coord, PETSc.DMSTAG_ELEMENT, 0);

xa_array[1:end  ,iu+1] .= a  .+ (b .- a .- (c./2.0)) .* X_coord[1:end,ixu+1] .+ (c./2.0).*X_coord[1:end,ixu+1].*X_coord[1:end,ixu+1];
xa_array[1:end-1,ip+1] .= b .- a .- (c./2.0) .+ c .* X_coord[1:end-1,ixp+1];

PETSc.DMLocalToGlobal(dm, xa_Local, PETSc.INSERT_VALUES, xa)
dmForcing = PETSc.DMStagCreateCompatibleDMStag(dm,1,0);
f         = PETSc.DMCreateGlobalVector(dmForcing);
fLocal    = PETSc.DMCreateLocalVector(dmForcing);
f        .= c;
fLocal   .= c;

A   = PETSc.DMCreateMatrix(dm);
rhs = PETSc.DMCreateGlobalVector(dm);

for e in start:start+n-1

    pos1 = PETSc.DMStagStencil(PETSc.DMSTAG_ELEMENT,e,0,0,0);
    val1 = 0.0;
    PETSc.DMStagVecSetValueStencil(dm, rhs, pos1, val1, PETSc.INSERT_VALUES);

    pos2 = PETSc.DMStagStencil(PETSc.DMSTAG_LEFT,e,0,0,0);
    if e == start
        val2 = a;
    else
        val2 = PETSc.DMStagVecGetValueStencil(dmForcing, fLocal, pos2);
    end
    PETSc.DMStagVecSetValueStencil(dm, rhs, pos2, val2, PETSc.INSERT_VALUES);

    if e == start+n-1
        pos3 = PETSc.DMStagStencil(PETSc.DMSTAG_RIGHT,e,0,0,0);
        val3 = b;
        PETSc.DMStagVecSetValueStencil(dm, rhs, pos3, val3, PETSc.INSERT_VALUES);
    end
end

PETSc.assemble(rhs)

for e in start:start+n-1
    row = PETSc.DMStagStencil(PETSc.DMSTAG_LEFT,e,0,0,0);
    if e == start
        val1 = 1.0;
        PETSc.DMStagMatSetValueStencil(dm, A, row, row, val1, PETSc.INSERT_VALUES);
    else
        col1 = PETSc.DMStagStencil(PETSc.DMSTAG_ELEMENT,e,0,0,0);
        col2 = PETSc.DMStagStencil(PETSc.DMSTAG_ELEMENT,e-1,0,0,0);

        xp1 = PETSc.DMStagVecGetValueStencil(dm_coord, vec_coord, col1);
        xp2 = PETSc.DMStagVecGetValueStencil(dm_coord, vec_coord, col2);
        h   = xp1-xp2;
        #print("h = ", h,", xp1 = ",xp1,", xp2 = ",xp2,"\n")
        
        val1 = 1.0/h;
        val2 = -1.0/h;
        val3 = 0.0;

        PETSc.DMStagMatSetValueStencil(dm, A, row, col1, val1, PETSc.INSERT_VALUES);
        PETSc.DMStagMatSetValueStencil(dm, A, row, col2, val2, PETSc.INSERT_VALUES);
        PETSc.DMStagMatSetValueStencil(dm, A, row, row , val3, PETSc.INSERT_VALUES);
    end
    if e == start+n-1
        row2 = PETSc.DMStagStencil(PETSc.DMSTAG_RIGHT,e,0,0,0);
        val4 = 1.0
        PETSc.DMStagMatSetValueStencil(dm, A, row2, row2, val4, PETSc.INSERT_VALUES);
    end

    row  = PETSc.DMStagStencil(PETSc.DMSTAG_ELEMENT,e,0,0,0);
    col1 = PETSc.DMStagStencil(PETSc.DMSTAG_RIGHT,e,0,0,0);
    col2 = PETSc.DMStagStencil(PETSc.DMSTAG_LEFT,e,0,0,0);

    xu1 = PETSc.DMStagVecGetValueStencil(dm_coord, vec_coord, col1);
    xu2 = PETSc.DMStagVecGetValueStencil(dm_coord, vec_coord, col2);
    h   = xu1-xu2;

    val1 = -1.0/h;
    val2 = 1.0/h;
    val3 = 1.0;

    PETSc.DMStagMatSetValueStencil(dm, A, row, col1, val1, PETSc.INSERT_VALUES);
    PETSc.DMStagMatSetValueStencil(dm, A, row, col2, val2, PETSc.INSERT_VALUES);
    PETSc.DMStagMatSetValueStencil(dm, A, row, row , val3, PETSc.INSERT_VALUES);
end

PETSc.assemble(A)

x   = PETSc.DMCreateGlobalVector(dm);
ksp = PETSc.KSP(A);
PETSc.solve!(x,ksp,rhs);

xLocal    = PETSc.DMCreateLocalVector(dm);
PETSc.DMGlobalToLocal(dm,x, PETSc.INSERT_VALUES,xLocal)
xu         =   PETSc.DMStagGetGhostArrayLocationSlot(dm,xLocal, PETSc.DMSTAG_LEFT, 0);
xp         =   PETSc.DMStagGetGhostArrayLocationSlot(dm,xLocal, PETSc.DMSTAG_ELEMENT, 0);

print("u_array = ",xu,"\np_array = ",xp,"\n");

xa_norm    = LinearAlgebra.norm(xa);
error      = xa.-x
error_norm = LinearAlgebra.norm(error);
errRel     = error_norm/xa_norm;

print("Error (abs): ",error_norm,"\nError (rel): ",errRel,"\n");

#PETSc.finalize()
