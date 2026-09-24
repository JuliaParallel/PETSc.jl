# test/lifetimes.jl
# A Vec or Mat built on Julia arrays uses their memory without copying it, so the
# arrays must live as long as the PETSc object. Each test drops every Julia
# reference to the arrays, forces a collection and reuses the memory, then reads
# the PETSc object: a freed array shows up as changed values.

using Test
using PETSc
using MPI

MPI.Initialized() || MPI.Init()

# Collect twice and fill the freed memory, in the small pools as well as the
# large ones, with a value the tests never use
function scrub!()
    GC.gc(true)
    GC.gc(true)
    large = [fill(7.0, 1000) for _ in 1:2000]
    small = [fill(7.0, n) for _ in 1:20_000 for n in 1:8]
    ints = [fill(7, n) for _ in 1:20_000 for n in 1:8]
    return sum(first, large) + sum(first, small) + sum(first, ints)
end

@testset "lifetimes of wrapped arrays" begin
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
    comm = MPI.COMM_SELF
    PetscScalar = petsclib.PetscScalar
    PetscInt = petsclib.PetscInt
    LibPETSc = PETSc.LibPETSc

    @testset "PetscVec on a Julia array" begin
        vs = [PETSc.PetscVec(petsclib, fill(PetscScalar(1), 1000)) for _ in 1:20]
        ws = [PETSc.PetscVec(petsclib, comm, fill(PetscScalar(2), 1000)) for _ in 1:20]
        scrub!()
        @test all(v -> v[1] == 1 && v[1000] == 1, vs)
        @test all(w -> w[1] == 2 && w[1000] == 2, ws)
        foreach(PETSc.destroy!, vs)
        foreach(PETSc.destroy!, ws)
    end

    @testset "PetscMat on CSR arrays, outliving its handle inside a KSP" begin
        # 2I on 3×3, from arrays nothing else refers to once this returns
        function csr_solver()
            A = PETSc.PetscMat(petsclib, PetscInt[0, 1, 2, 3], PetscInt[0, 1, 2], fill(PetscScalar(2), 3))
            ksp = PETSc.KSP(A)
            PETSc.destroy!(A)      # the KSP still holds the matrix
            return ksp
        end
        ksp = csr_solver()
        scrub!()
        x = ksp \ PetscScalar[2, 4, 6]
        @test x ≈ PetscScalar[1, 2, 3]
        PETSc.destroy!(ksp)
    end

    @testset "temporaries survive a collection during the operation" begin
        S = PETSc.PetscMat(petsclib, PetscScalar[2 0 0; 0 2 0; 0 0 2])
        ksp = PETSc.KSP(S)
        p = PETSc.pc(ksp)
        PETSc.set_type!(p, :shell)
        PETSc.set_shell_apply!(p) do y, _p, x
            scrub!()               # while `ksp \ b` holds its copy of b
            PETSc.with_local_array!((y, x); write = (true, false)) do ya, xa
                ya .= xa
            end
        end
        @test ksp \ PetscScalar[2, 4, 6] ≈ PetscScalar[1, 2, 3]
        PETSc.destroy!(ksp)
        PETSc.destroy!(S)

        function f!(y, x)
            scrub!()               # while `M * x` holds its copy of x
            PETSc.with_local_array!((y, x); write = (true, false)) do ya, xa
                ya .= 2 .* xa
            end
        end
        M = PETSc.MatShell(petsclib, f!, comm, 3, 3)
        @test M * PetscScalar[1, 2, 3] == PetscScalar[2, 4, 6]
        PETSc.destroy!(M)
    end

    @testset "KSP from a SparseMatrixCSC keeps no extra reference to its matrix" begin
        refcount(ksp) = LibPETSc.PetscObjectGetReference(petsclib, first(LibPETSc.KSPGetOperators(petsclib, ksp)))
        S = PETSc.SparseArrays.sparse(PetscScalar[2 0; 0 2])
        ksp = PETSc.KSP(petsclib, comm, S)
        # the same solver built by hand, with the caller's reference released
        A = PETSc.PetscMat(petsclib, comm, S)
        ksp_ref = PETSc.KSP(A)
        PETSc.destroy!(A)
        @test refcount(ksp) == refcount(ksp_ref)
        @test ksp \ PetscScalar[2, 4] ≈ PetscScalar[1, 2]
        PETSc.destroy!(ksp)
        PETSc.destroy!(ksp_ref)
    end

    PETSc.finalize(petsclib)
end
