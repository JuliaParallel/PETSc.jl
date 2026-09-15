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

