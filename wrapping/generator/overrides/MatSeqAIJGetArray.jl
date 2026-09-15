# override for MatSeqAIJGetArray; C signature: MatSeqAIJGetArray(Mat A, PetscScalar* array[])
"""
    array::Vector{PetscScalar} = MatSeqAIJGetArray(petsclib::PetscLibType, A::AbstractPetscMat)

Returns a 1D Julia Array (Vector) providing direct read/write access to the internal
numerical values of a `MATSEQAIJ` (Sparse) matrix.

Note: This only provides access to the values array. It does not allow changing
the sparsity pattern (row pointers or column indices).

Not Collective

Input Parameter:
- `A` - a `MATSEQAIJ` matrix

Output Parameter:
- `array` - A `Vector` view of the non-zero values.

Level: intermediate

See also: `MatSeqAIJRestoreArray()`, `MatDenseGetArray()`

# External Links
$(_doc_external("Mat/MatSeqAIJGetArray"))
"""
function MatSeqAIJGetArray(petsclib::PetscLibType, A::AbstractPetscMat) end

@for_petsc function MatSeqAIJGetArray(petsclib::$UnionPetscLib, A::AbstractPetscMat)
    array_ = Ref{Ptr{$PetscScalar}}()

    # 1. Get the pointer to the numerical values
    @chk ccall(
               (:MatSeqAIJGetArray, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{Ptr{$PetscScalar}}),
               A, array_,
              )

    # 2. Get the matrix info to determine the number of nonzeros
    info_ref = Ref{LibPETSc.MatInfo}()
    LibPETSc.MatGetInfo(petsclib, A, LibPETSc.MAT_LOCAL, info_ref)

    # 3. Extract the nnz count using the [] syntax to unwrap the Ref
    nnz = Int(info_ref[].nz_used)

    return unsafe_wrap(Array, array_[], nnz; own = false)
end

