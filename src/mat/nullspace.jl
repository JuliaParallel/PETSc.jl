"""
    mat_nullspace_create(petsclib, comm; has_const = true) -> MatNullSpace

Create a null space containing the constant vector (Neumann / saddle-point problems).
"""
function mat_nullspace_create(petsclib, comm; has_const::Bool = true)
    return LibPETSc.MatNullSpaceCreate(petsclib, comm,
        LibPETSc.PetscBool(has_const),
        petsclib.PetscInt(0),
        LibPETSc.PetscVec{typeof(petsclib)}[])
end

"""
    mat_nullspace_create(petsclib, comm, vecs) -> MatNullSpace

Create a null space spanned by the given vectors.  Each element of `vecs`
must have a `.ptr` field holding the underlying PETSc Vec handle (`CVec`).
The vectors should be orthonormal; call `VecNormalize` beforehand if needed.
"""
function mat_nullspace_create end

LibPETSc.@for_petsc function mat_nullspace_create(
    petsclib::$UnionPetscLib,
    comm,
    vecs,
)
    # the wrapper converts the handles itself; accept anything with a `.ptr`
    pvecs = [LibPETSc.PetscVec{$PetscLib}(v.ptr; own = false) for v in vecs]
    return LibPETSc.MatNullSpaceCreate(petsclib, comm,
        LibPETSc.PETSC_FALSE, $PetscInt(length(pvecs)), pvecs)
end

"""
    destroy!(nullsp::LibPETSc.MatNullSpace)
    destroy!(petsclib, nullsp::LibPETSc.MatNullSpace)

Destroy a `MatNullSpace` created by [`mat_nullspace_create`](@ref).

`MatNullSpace` is a bare `Ptr` and carries no library, so the one-argument form
uses the default library. Pass `petsclib` explicitly when several are in use.

# External Links
$(doc_external("Mat/MatNullSpaceDestroy"))
"""
function destroy!(petsclib::LibPETSc.PetscLibType, nullsp::LibPETSc.MatNullSpace)
    LibPETSc.MatNullSpaceDestroy(petsclib, nullsp)
    return nothing
end

destroy!(nullsp::LibPETSc.MatNullSpace) = destroy!(getlib(), nullsp)

"""
    set_nullspace!(mat, nullsp)

Attach `nullsp` to `mat` so the linear solver removes it each iteration.
"""
function set_nullspace!(mat, nullsp)
    LibPETSc.MatSetNullSpace(LibPETSc.getlib(typeof(mat).parameters[1]), mat, nullsp)
    return mat
end
