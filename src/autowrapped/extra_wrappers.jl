# Hand-written wrappers for functions/macros not described by getAPI.py (wrapping/overrides/)
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

