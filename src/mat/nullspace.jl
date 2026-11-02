"""
    mat_nullspace_create(petsclib, comm; has_const = true) -> MatNullSpace

Create a null space containing the constant vector (Neumann / saddle-point problems).
On a one-rank communicator the garbage collector destroys it if you do not.
"""
function mat_nullspace_create(petsclib, comm; has_const::Bool = true)
    nullsp = LibPETSc.MatNullSpaceCreate(petsclib, comm,
        LibPETSc.PetscBool(has_const),
        petsclib.PetscInt(0),
        LibPETSc.PetscVec{typeof(petsclib)}[])
    return own_serial!(nullsp, comm)
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
    nullsp = LibPETSc.MatNullSpaceCreate(petsclib, comm,
        LibPETSc.PETSC_FALSE, $PetscInt(length(pvecs)), pvecs)
    return own_serial!(nullsp, comm)
end

"""
    destroy!(nullsp::MatNullSpace)

Destroy a `MatNullSpace` created by [`mat_nullspace_create`](@ref). Does nothing
on a borrowed handle.

# External Links
$(doc_external("Mat/MatNullSpaceDestroy"))
"""
function destroy!(nullsp::LibPETSc.AbstractMatNullSpace{PetscLib}) where {PetscLib}
    return destroy_handle!(LibPETSc.MatNullSpaceDestroy, nullsp, PetscLib)
end

"""
    set_nullspace!(mat, nullsp)

Attach `nullsp` to `mat` so the linear solver removes it each iteration.
"""
function set_nullspace!(
    mat::AbstractPetscMat{PetscLib},
    nullsp::LibPETSc.AbstractMatNullSpace{PetscLib},
) where {PetscLib}
    LibPETSc.MatSetNullSpace(getlib(PetscLib), mat, nullsp)
    return mat
end
