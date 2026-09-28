using PETSc

# Minimal SNES helper test
petsclib = PETSc.petsclibs[1]
PETSc.initialize(petsclib)

snes = LibPETSc.SNESCreate(petsclib, LibPETSc.PETSC_COMM_SELF)
LibPETSc.SNESSetFromOptions(petsclib, snes)

# Call helpers (should not MethodError)
reason = LibPETSc.SNESGetConvergedReason(petsclib, snes)
iters = LibPETSc.SNESGetIterationNumber(petsclib, snes)

@assert typeof(reason) <: Integer || typeof(reason) <: Enum
@assert typeof(iters) <: Integer

LibPETSc.SNESDestroy(petsclib, snes)
PETSc.finalize(petsclib)
