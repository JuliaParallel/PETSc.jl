"""
   initialized(petsclib)

Check if `petsclib` is initialized

# External Links
$(_doc_external("Sys/PetscInitialized"))
"""
function Initialized end
@for_petsc function initialized(::$UnionPetscLib)
    r_flag = Ref{PetscBool}()
    @chk ccall(
        (:PetscInitialized, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscBool},),
        r_flag,
    )
    return r_flag[] == PETSC_TRUE
end

"""
   initialize([petsclib])

Initialized the `petsclib`, if no `petsclib` is given then all `PETSc.petsclibs`
will be initialized.

Additionally:
 - This will initialize MPI if it has not already been initialized.
 - It will disable the PETSc signal handler (via
   $(_petsc_link("Sys/PetscPopSignalHandler"))
 - Add an [`atexit`](https://docs.julialang.org/en/v1/base/base/#Base.atexit)
   hook to call [`PETSc.finalize`](@ref).

# External Links
$(_doc_external("Sys/PetscInitializeNoArguments"))
"""
function Initialize end

function initialize()
    map(initialize, petsclibs)
    return nothing
end

@for_petsc function initialize(::$UnionPetscLib)
    if !initialized($petsclib)
        MPI.Initialized() || MPI.Init()
        @chk ccall(
            (:PetscInitializeNoArguments, $petsc_library),
            PetscErrorCode,
            (),
        )

        # disable signal handler
        @chk ccall((:PetscPopSignalHandler, $petsc_library), PetscErrorCode, ())

        atexit(() -> finalize($petsclib))
    end
    return nothing
end

"""
   finalize(petsclib)

Finalize the `petsclib`, if no `petsclib` is given then all `PETSc.petsclibs`
will be finalized.

# External Links
$(_doc_external("Sys/PetscFinalize"))
"""
function finalize end

function finalize()
    map(finalize, petsclibs)
    return nothing
end

@for_petsc function finalize(::$UnionPetscLib)
    if !finalized($petsclib)
        @chk ccall((:PetscFinalize, $petsc_library), PetscErrorCode, ())
    end
    return nothing
end

"""
   finalized(petsclib)

Check if `petsclib` is finalized

# External Links
$(_doc_external("Sys/PetscFinalized"))
"""
function finalized end
@for_petsc function finalized(::$UnionPetscLib)
    r_flag = Ref{PetscBool}()
    @chk ccall(
        (:PetscFinalized, $petsc_library),
        PetscErrorCode,
        (Ptr{PetscBool},),
        r_flag,
    )
    return r_flag[] == PETSC_TRUE
end
