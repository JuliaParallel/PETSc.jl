# test/lowlevel/readers.jl
# LibPETSc readers on an object that has not solved anything return plain values,
# not handles or references.

using Test
using PETSc

@testset "LibPETSc readers before a solve" begin
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
    comm = LibPETSc.PETSC_COMM_SELF

    snes = LibPETSc.SNESCreate(petsclib, comm)
    LibPETSc.SNESSetFromOptions(petsclib, snes)
    @test LibPETSc.SNESGetConvergedReason(petsclib, snes) isa LibPETSc.SNESConvergedReason
    @test LibPETSc.SNESGetIterationNumber(petsclib, snes) == 0
    LibPETSc.SNESDestroy(petsclib, snes)

    # Tao objects work only where the Tao types survive a PETSc cycle (not on Windows)
    if PETSc.tao_usable_after_reinitialize()
        tao = LibPETSc.TaoCreate(petsclib, comm)
        @test LibPETSc.TaoGetConvergedReason(petsclib, tao) isa LibPETSc.TaoConvergedReason
        LibPETSc.TaoDestroy(petsclib, tao)
    end

    PETSc.finalize(petsclib)
end
