# Hand-written wrappers for functions/macros not described by getAPI.py (wrapping/overrides/)
# override for DMComputeL2Diff; C signature: DMComputeL2Diff(<not in API snapshot>)
"""
    DMComputeL2Diff(petsclib, dm, time, funcs, ctxs, X) -> PetscReal

Compute the L² norm of the difference between the global vector `X` and the
pointwise exact functions `funcs`.  `ctxs` is a `Vector{Ptr{Cvoid}}` of
context pointers (use `C_NULL` entries for no context).
"""
function DMComputeL2Diff(petsclib::PetscLibType, dm::AbstractPetscDM, time::Real,
                         funcs::Vector{Ptr{Cvoid}}, ctxs::Vector{Ptr{Cvoid}},
                         X::AbstractPetscVec) end

@for_petsc function DMComputeL2Diff(petsclib::$UnionPetscLib, dm::AbstractPetscDM,
                                    time::Real,
                                    funcs::Vector{Ptr{Cvoid}},
                                    ctxs::Vector{Ptr{Cvoid}},
                                    X::AbstractPetscVec)
    diff_ref = Ref{$PetscReal}(0)
    GC.@preserve funcs ctxs @chk ccall(
        (:DMComputeL2Diff, $petsc_library), PetscErrorCode,
        (CDM, $PetscReal, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}, CVec, Ptr{$PetscReal}),
        dm, $PetscReal(time), funcs, ctxs, X, diff_ref)
    return diff_ref[]
end

# override for DMProjectFunction; C signature: DMProjectFunction(<not in API snapshot>)
"""
    DMProjectFunction(petsclib, dm, time, funcs, ctxs, mode, X)

Project the pointwise functions `funcs` (one `Ptr{Cvoid}` per field, matching
`PetscSimplePointFn` signature) into the global vector `X`.  `ctxs` is a
matching `Vector{Ptr{Cvoid}}` of context pointers (use `C_NULL` entries for
no context).
"""
function DMProjectFunction(petsclib::PetscLibType, dm::AbstractPetscDM, time::Real,
                           funcs::Vector{Ptr{Cvoid}}, ctxs::Vector{Ptr{Cvoid}},
                           mode::InsertMode, X::AbstractPetscVec) end

@for_petsc function DMProjectFunction(petsclib::$UnionPetscLib, dm::AbstractPetscDM,
                                      time::Real,
                                      funcs::Vector{Ptr{Cvoid}},
                                      ctxs::Vector{Ptr{Cvoid}},
                                      mode::InsertMode, X::AbstractPetscVec)
    GC.@preserve funcs ctxs @chk ccall(
        (:DMProjectFunction, $petsc_library), PetscErrorCode,
        (CDM, $PetscReal, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}, InsertMode, CVec),
        dm, $PetscReal(time), funcs, ctxs, mode, X)
    return nothing
end

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

# override for PETSC_VIEWER_STDERR_SELF; C signature: PETSC_VIEWER_STDERR_SELF(<not in API snapshot>)
"""
	viewer::PetscViewer = PETSC_VIEWER_STDERR_SELF(petsclib::PetscLibType)

Get the default PETSc `STDERR` viewer for `MPI.COMM_SELF`

Not Collective

Input Parameter:
- `petsclib` - the PETSc library instance

Output Parameter:
- `viewer` - the viewer

Level: beginner

-seealso: `PETSC_VIEWER_STDERR_WORLD`, `PetscViewerASCIIGetStderr()`
"""
function PETSC_VIEWER_STDERR_SELF(petsclib::PetscLibType) end

@for_petsc function PETSC_VIEWER_STDERR_SELF(petsclib::$UnionPetscLib)
	viewer_ref = Ref{PetscViewer}()
	@chk ccall(
		(:PetscViewerASCIIGetStderr, $petsc_library),
		PetscErrorCode,
		(MPI_Comm, Ptr{PetscViewer}),
		MPI.COMM_SELF, viewer_ref,
	)
	return viewer_ref[]
end

# override for PETSC_VIEWER_STDERR_WORLD; C signature: PETSC_VIEWER_STDERR_WORLD(<not in API snapshot>)
"""
	viewer::PetscViewer = PETSC_VIEWER_STDERR_WORLD(petsclib::PetscLibType)

Get the default PETSc `STDERR` viewer for `MPI.COMM_WORLD`

Collective on `MPI.COMM_WORLD`

Input Parameter:
- `petsclib` - the PETSc library instance

Output Parameter:
- `viewer` - the viewer

Level: beginner

-seealso: `PETSC_VIEWER_STDERR_SELF`, `PetscViewerASCIIGetStderr()`
"""
function PETSC_VIEWER_STDERR_WORLD(petsclib::PetscLibType) end

@for_petsc function PETSC_VIEWER_STDERR_WORLD(petsclib::$UnionPetscLib)
	viewer_ref = Ref{PetscViewer}()
	@chk ccall(
		(:PetscViewerASCIIGetStderr, $petsc_library),
		PetscErrorCode,
		(MPI_Comm, Ptr{PetscViewer}),
		MPI.COMM_WORLD, viewer_ref,
	)
	return viewer_ref[]
end

# override for PETSC_VIEWER_STDOUT_SELF; C signature: PETSC_VIEWER_STDOUT_SELF(<not in API snapshot>)
"""
	viewer::PetscViewer = PETSC_VIEWER_STDOUT_SELF(petsclib::PetscLibType)

Get the default PETSc `STDOUT` viewer for `MPI.COMM_SELF`

Not Collective

Input Parameter:
- `petsclib` - the PETSc library instance

Output Parameter:
- `viewer` - the viewer

Level: beginner

-seealso: `PETSC_VIEWER_STDOUT_WORLD`, `PetscViewerASCIIGetStdout()`
"""
function PETSC_VIEWER_STDOUT_SELF(petsclib::PetscLibType) end

@for_petsc function PETSC_VIEWER_STDOUT_SELF(petsclib::$UnionPetscLib)
	viewer_ref = Ref{PetscViewer}()
	@chk ccall(
		(:PetscViewerASCIIGetStdout, $petsc_library),
		PetscErrorCode,
		(MPI_Comm, Ptr{PetscViewer}),
		MPI.COMM_SELF, viewer_ref,
	)
	return viewer_ref[]
end

# override for PETSC_VIEWER_STDOUT_WORLD; C signature: PETSC_VIEWER_STDOUT_WORLD(<not in API snapshot>)
"""
	viewer::PetscViewer = PETSC_VIEWER_STDOUT_WORLD(petsclib::PetscLibType)

Get the default PETSc `STDOUT` viewer for `MPI.COMM_WORLD`

Collective on `MPI.COMM_WORLD`

Input Parameter:
- `petsclib` - the PETSc library instance

Output Parameter:
- `viewer` - the viewer

Level: beginner

-seealso: `PETSC_VIEWER_STDOUT_SELF`, `PetscViewerASCIIGetStdout()`
"""
function PETSC_VIEWER_STDOUT_WORLD(petsclib::PetscLibType) end

@for_petsc function PETSC_VIEWER_STDOUT_WORLD(petsclib::$UnionPetscLib)
	viewer_ref = Ref{PetscViewer}()
	@chk ccall(
		(:PetscViewerASCIIGetStdout, $petsc_library),
		PetscErrorCode,
		(MPI_Comm, Ptr{PetscViewer}),
		MPI.COMM_WORLD, viewer_ref,
	)
	return viewer_ref[]
end

# override for SNESGetJacobianMat; C signature: SNESGetJacobianMat(<not in API snapshot>)
"""
    SNESGetJacobianMat(petsclib, snes) -> PetscMat

Return the assembled system (A) matrix from `snes`, passing `NULL` for the
preconditioner matrix, Jacobian function, and context.  Use this when only
the system matrix is needed (e.g. to attach a null space via `MatSetNullSpace`).
"""
function SNESGetJacobianMat(petsclib::PetscLibType, snes::AbstractPetscSNES) end

@for_petsc function SNESGetJacobianMat(petsclib::$UnionPetscLib, snes::AbstractPetscSNES)
    J_ref = Ref{CMat}(C_NULL)
    @chk ccall(
        (:SNESGetJacobian, $petsc_library), PetscErrorCode,
        (CSNES, Ptr{CMat}, Ptr{CMat}, Ptr{Cvoid}, Ptr{Ptr{Cvoid}}),
        snes, J_ref, C_NULL, C_NULL, C_NULL)
    return PetscMat{$PetscLib}(J_ref[])
end

