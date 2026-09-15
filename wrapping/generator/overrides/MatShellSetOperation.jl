# override for MatShellSetOperation; C signature: MatShellSetOperation(<not in API snapshot>)
"""
    MatShellSetOperation(petsclib::PetscLibType, mat::AbstractPetscMat, op::MatOperation, g::Ptr)

Allows user to set a matrix operation for a `MATSHELL` shell matrix.

Logically Collective

Input Parameters:
`mat` - the `MATSHELL` shell matrix
`op`  - the name of the operation
`g`   - a pointer to the function that provides the operation created with `@cfunction`

Level: advanced

-seealso: `Mat`, `MATSHELL`, `MatCreateShell()`, `MatShellGetContext()`, `MatShellGetOperation()`, `MatShellSetContext()`, `MatSetOperation()`, `MatShellSetManageScalingShifts()`, `MatShellSetMatProductOperation()`

# External Links
$(_doc_external("Mat/MatShellSetOperation"))
"""
function MatShellSetOperation(petsclib::PetscLibType, mat::AbstractPetscMat, op::MatOperation, g::Ptr) end

@for_petsc function MatShellSetOperation(petsclib::$UnionPetscLib, mat::AbstractPetscMat, op::MatOperation, g::Ptr)

    @chk ccall(
               (:MatShellSetOperation, $petsc_library),
               PetscErrorCode,
               (CMat, MatOperation, Ptr{Cvoid}),
               mat, op, g,
              )

	return nothing
end

