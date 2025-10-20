"""
   initialized(petsclib)

Check if `petsclib` is initialized

# External Links
$(_doc_external("Sys/PetscInitialized"))
"""
initialized(petsclib) = LibPETSc.PetscInitialized(petsclib)

"""
   initialize([petsclib])
?
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
function initialize()
    map(initialize, petsclibs)
    return nothing
end

function initialize(petsclib)
    if !initialized(petsclib)
        MPI.Initialized() || MPI.Init()
        petsclib.age += 1
        LibPETSc.PetscInitializeNoArguments(petsclib)

        # disable signal handler
        LibPETSc.PetscPopSignalHandler(petsclib)

        atexit(() -> finalize(petsclib))
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
function finalize()
    map(finalize, petsclibs)
    return nothing
end

function finalize(petsclib)
    if !finalized(petsclib)
        petsclib.age += 1
        LibPETSc.PetscFinalize(petsclib)
    end
    return nothing
end

"""
   finalized(petsclib)

Check if `petsclib` is finalized

# External Links
$(_doc_external("Sys/PetscFinalized"))
"""
finalized(petsclib) = LibPETSc.PetscFinalized(petsclib)


"""
    scalartype(petsclib::PetscLibType)

return the scalar type for the associated `petsclib`
"""
scalartype(::LibPETSc.PetscLibType{ST}) where {ST} = ST
scalartype(
    ::Type{PetscLib},
) where {PetscLib <: PetscLibType{ST}} where {ST} = ST


"""
    inttype(petsclib::PetscLibType)

return the int type for the associated `petsclib`
"""
inttype(::LibPETSc.PetscLibType{ST, IT}) where {ST, IT} = IT
inttype(
    ::Type{PetscLib},
) where {PetscLib <: PetscLibType{ST, IT}} where {ST, IT} = IT