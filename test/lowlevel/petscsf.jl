# Low-level PetscSF: graph setup and the hand-written communication wrappers
# (PetscSFBcast*/PetscSFReduce*/PetscSFFetchAndOp*, which getAPI.py does not export)
using Test
using PETSc, MPI
using PETSc: LibPETSc

MPI.Initialized() || MPI.Init()

@testset "Low-level PetscSF communication" begin
    petsclib = PETSc.getlib()
    PETSc.initialize(petsclib)
    PetscInt = petsclib.PetscInt

    # every leaf i reads root (n - 1 - i) on this rank: a permutation star forest
    n = 5
    sf = LibPETSc.PetscSFCreate(petsclib, MPI.COMM_SELF)
    ilocal = PetscInt.(0:n-1)
    iremote = [LibPETSc.PetscSFNode(0, n - 1 - i) for i in 0:n-1]
    LibPETSc.PetscSFSetGraph(petsclib, sf, n, n, ilocal, LibPETSc.PETSC_COPY_VALUES,
                             iremote, LibPETSc.PETSC_COPY_VALUES)
    LibPETSc.PetscSFSetUp(petsclib, sf)

    nroots, nleaves, il, ir = LibPETSc.PetscSFGetGraph(petsclib, sf)
    @test nroots == n && nleaves == n
    @test il === nothing || il == ilocal      # NULL ilocal means the leaves are contiguous [0, nleaves)
    @test [r.index for r in ir] == [n - 1 - i for i in 0:n-1]

    # broadcast roots -> leaves
    rootdata = Float64.(1:n)
    leafdata = zeros(Float64, n)
    LibPETSc.PetscSFBcastBegin(petsclib, sf, MPI.Datatype(Float64), rootdata, leafdata, MPI.REPLACE)
    LibPETSc.PetscSFBcastEnd(petsclib, sf, MPI.Datatype(Float64), rootdata, leafdata, MPI.REPLACE)
    @test leafdata == reverse(rootdata)

    # reduce leaves -> roots (sum)
    acc = zeros(Float64, n)
    LibPETSc.PetscSFReduceBegin(petsclib, sf, MPI.Datatype(Float64), leafdata, acc, MPI.SUM)
    LibPETSc.PetscSFReduceEnd(petsclib, sf, MPI.Datatype(Float64), leafdata, acc, MPI.SUM)
    @test acc == rootdata

    # fetch-and-op: leafupdate receives the old root value, roots are incremented
    roots = PetscInt.(10:10:10n)
    leaves = ones(PetscInt, n)
    update = zeros(PetscInt, n)
    LibPETSc.PetscSFFetchAndOpBegin(petsclib, sf, MPI.Datatype(PetscInt), roots, leaves, update, MPI.SUM)
    LibPETSc.PetscSFFetchAndOpEnd(petsclib, sf, MPI.Datatype(PetscInt), roots, leaves, update, MPI.SUM)
    @test update == reverse(PetscInt.(10:10:10n))
    @test roots == PetscInt.(10:10:10n) .+ 1

    @test (@inferred LibPETSc.PetscSFGetLeafRange(petsclib, sf)) == (0, n - 1)
    LibPETSc.PetscSFDestroy(petsclib, sf)
    PETSc.finalize(petsclib)
end
