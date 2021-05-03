using PETSc, MPI

if ~MPI.Initialized()
    MPI.Init()
end

DmSol = DMStagCreate1d(MPI.COMM_SELF,PETSc.DM_BOUNDARY_NONE,3,1,1,PETSc.DMSTAG_STENCIL_BOX,1,[])