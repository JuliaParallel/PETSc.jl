# test/solvers/callback_contract.jl
# What a callback hands to PETSc.jl and gets back: its return value is ignored,
# and the norms given to a convergence test have the library's PetscReal type.

using Test
using PETSc
using MPI
using Logging

MPI.Initialized() || MPI.Init()

# F(x) = x² − 4 on one unknown, root 2
function square_snes(petsclib, comm, residual_return, jacobian_return, test = nothing)
    PetscScalar = petsclib.PetscScalar
    PetscInt = petsclib.PetscInt
    snes = PETSc.SNES(petsclib, comm; pc_type = "lu", snes_rtol = 1e-6)
    r = PETSc.PetscVec(petsclib, zeros(PetscScalar, 1))
    PETSc.set_function!(snes, r) do fx, _snes, x
        fx[1] = x[1]^2 - 4
        return residual_return
    end
    J = PETSc.PetscMat(petsclib, zeros(PetscScalar, 1, 1))
    PETSc.set_snes_jacobian!(snes, J) do J, _snes, x
        J[1, 1] = 2x[1]
        PETSc.assemble!(J)
        return jacobian_return
    end
    isnothing(test) || PETSc.set_convergence_test!(test, snes)
    return snes, r, J
end

@testset "Callback contract" begin
    comm = MPI.COMM_SELF

    for petsclib in PETSc.petsclibs
        PETSc.initialize(petsclib)
        PetscScalar = petsclib.PetscScalar
        PetscReal = petsclib.PetscReal

        @testset "a callback's return value is ignored, $PetscScalar" begin
            for (rf, rj) in ((nothing, nothing), (0, 0), (1, 7), (-3, PETSc.LibPETSc.PetscErrorCode(76)))
                snes, r, J = square_snes(petsclib, comm, rf, rj)
                x = PETSc.PetscVec(petsclib, PetscScalar[3])
                logger = Test.TestLogger(min_level = Logging.Warn)
                with_logger(logger) do
                    PETSc.solve!(x, snes)
                end
                @test isempty(logger.logs)
                @test Int(PETSc.converged_reason(snes)) > 0
                @test x[1] ≈ 2 rtol = 1e-4
                foreach(PETSc.destroy!, (snes, x, r, J))
            end
        end

        @testset "convergence test norms are PetscReal, $PetscReal" begin
            types = Set{DataType}()
            test = (_snes, it, xnorm, gnorm, fnorm) -> begin
                push!(types, typeof(xnorm), typeof(gnorm), typeof(fnorm))
                it >= 1 && return PETSc.LibPETSc.SNES_CONVERGED_USER
                return PETSc.LibPETSc.SNES_CONVERGED_ITERATING
            end
            snes, r, J = square_snes(petsclib, comm, nothing, nothing, test)
            x = PETSc.PetscVec(petsclib, PetscScalar[3])
            PETSc.solve!(x, snes)
            @test types == Set([PetscReal])
            foreach(PETSc.destroy!, (snes, x, r, J))
        end
    end
end
