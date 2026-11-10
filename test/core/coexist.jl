# test/core/coexist.jl
# Libraries for different scalar types live side by side in one process: a
# Float32 and a ComplexF64 solver, each with a Julia callback, are set up
# together and used in turn. Every wrapper is generated once per library with
# that library's path in its ccall, and dispatch picks it from the object's type.

using Test
using PETSc, MPI, LinearAlgebra, SparseArrays

MPI.Initialized() || MPI.Init()

@testset "Float32 and ComplexF64 libraries side by side" begin
    library(T) = findfirst(lib -> lib.PetscScalar === T, PETSc.petsclibs)
    i_f, i_c = library(Float32), library(ComplexF64)
    if i_f === nothing || i_c === nothing
        # a custom library (set_library!) provides one build only
        @test_skip "needs both a Float32 and a ComplexF64 library"
    else
        lib_f, lib_c = PETSc.petsclibs[i_f], PETSc.petsclibs[i_c]
        PETSc.initialize(lib_f)
        PETSc.initialize(lib_c)
        comm = MPI.COMM_SELF

        n = 20
        S = spdiagm(-1 => fill(-1.0, n - 1), 0 => fill(2.0, n), 1 => fill(-1.0, n - 1))
        A_f = PETSc.PetscMat(lib_f, comm, SparseMatrixCSC{Float32}(S))
        A_c = PETSc.PetscMat(lib_c, comm, SparseMatrixCSC{ComplexF64}(S))
        @test typeof(A_f) != typeof(A_c)   # the library is part of the type

        ksp_f = PETSc.KSP(A_f; ksp_type = "cg", pc_type = "jacobi", ksp_rtol = 1e-6)
        ksp_c = PETSc.KSP(A_c; ksp_type = "gmres", pc_type = "jacobi", ksp_rtol = 1e-12)

        # the solves alternate between the two libraries
        b_f = ones(Float32, n)
        b_c = ComplexF64.(1:n) .+ im
        x_f = ksp_f \ b_f
        x_c = ksp_c \ b_c
        x_f2 = ksp_f \ b_f
        @test eltype(x_f) === Float32
        @test eltype(x_c) === ComplexF64
        @test S * x_f ≈ b_f rtol = 1e-4
        @test S * x_c ≈ b_c rtol = 1e-10
        @test x_f2 == x_f

        # a Julia closure as the operator of each library, also called in turn
        factor = 3
        shell_f = PETSc.MatShell((y, x) -> y .= factor .* x, lib_f, comm, n, n)
        shell_c = PETSc.MatShell((y, x) -> y .= (factor * im) .* x, lib_c, comm, n, n)
        z_f = rand(Float32, n)
        z_c = rand(ComplexF64, n)
        @test shell_f * z_f ≈ 3z_f
        @test shell_c * z_c ≈ 3im * z_c
        @test shell_f * z_f ≈ 3z_f

        foreach(PETSc.destroy!, (ksp_f, ksp_c, shell_f, shell_c, A_f, A_c))
        PETSc.finalize(lib_f)
        PETSc.finalize(lib_c)
    end
end
