# override for SNESGetJacobianMat; C signature: SNESGetJacobianMat(<not in API snapshot>)
"""
    SNESGetJacobianMat(petsclib, snes) -> PetscMat

Return the assembled system (A) matrix from `snes`, passing `NULL` for the
preconditioner matrix, Jacobian function, and context.  Use this when only
the system matrix is needed (e.g. to attach a null space via `MatSetNullSpace`).
"""
function SNESGetJacobianMat(petsclib::PetscLibType, snes::AbstractSNES) end

@for_petsc function SNESGetJacobianMat(petsclib::$UnionPetscLib, snes::AbstractSNES)
    J_ref = Ref{CMat}(C_NULL)
    @chk ccall(
        (:SNESGetJacobian, $petsc_library), PetscErrorCode,
        (CSNES, Ptr{CMat}, Ptr{CMat}, Ptr{Cvoid}, Ptr{Ptr{Cvoid}}),
        snes, J_ref, C_NULL, C_NULL, C_NULL)
    return PetscMat{$PetscLib}(J_ref[])
end

