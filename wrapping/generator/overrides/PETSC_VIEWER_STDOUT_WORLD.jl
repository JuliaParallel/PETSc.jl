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

