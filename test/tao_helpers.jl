using PETSc

# Minimal TAO helper test
petsclib = PETSc.petsclibs[1]
PETSc.initialize(petsclib)

if !PETSc.tao_usable_after_reinitialize()
    @info "Skipping Tao tests: PETSc 3.25.x loses the Tao types at PetscFinalize and the workaround is not possible with these binaries (Windows)"
else

tao = LibPETSc.TaoCreate(petsclib, LibPETSc.PETSC_COMM_SELF)
# Avoid calling TaoSetType string wrapper; keep minimal

# Call helper (should not MethodError)
reason = LibPETSc.TaoGetConvergedReason(petsclib, tao)
@assert typeof(reason) <: Integer || typeof(reason) <: Enum

LibPETSc.TaoDestroy(petsclib, tao)
PETSc.finalize(petsclib)
end
