"""
    save_vtk(vec::AbstractPetscVec, filename)
    save_vtk(vecs, filename)

Write the global vector `vec` to a VTK unstructured-grid file (`.vtu`),
readable by ParaView and VisIt.  Works for both serial and MPI runs; PETSc
gathers all ranks into a single file.  The filename must end in `.vtu`.
The mesh geometry is taken from the `DM` attached to `vec` by PETSc.

The second form takes any iterable of `AbstractPetscVec` and writes each as a
separate point-data array in one file.  The field name shown in ParaView/VisIt
comes from the name set on the vector via [`set_name!`](@ref).

Neither `petsclib` nor `comm` is an argument: the vector carries the library as
a type parameter and PETSc is asked for its communicator (§8).

`vtk_save!` and `vtk_save_fields!` were the two v0.4 names for these two
methods, and both shim the old argument list.

# External Links
$(doc_external("Viewer/PetscViewerVTKOpen"))
"""
function save_vtk end

LibPETSc.@for_petsc function save_vtk(
    vec::AbstractPetscVec{$PetscLib},
    filename::AbstractString,
)
    petsclib = getlib($PetscLib)
    comm = LibPETSc.PetscObjectGetComm(petsclib, vec)
    vtk = LibPETSc.PetscViewerVTKOpen(petsclib, comm, String(filename), LibPETSc.FILE_MODE_WRITE)
    LibPETSc.VecView(petsclib, vec, vtk)
    vtk_ref = Ref{LibPETSc.PetscViewer}(vtk)
    LibPETSc.PetscViewerDestroy(petsclib, vtk_ref)
    return nothing
end

# The multi-field method of `save_vtk`, documented with the single-vector one
# above. The communicator comes from the first vector, so the iterable must not
# be empty.
function save_vtk(
    vecs,   # iterable of AbstractPetscVec
    filename::AbstractString,
)
    isempty(vecs) && throw(
        ArgumentError("save_vtk needs at least one vector to take the communicator from"),
    )
    v1 = first(vecs)
    petsclib = petsclib_of(v1)
    comm = LibPETSc.PetscObjectGetComm(petsclib, v1)
    vtk = LibPETSc.PetscViewerVTKOpen(petsclib, comm, String(filename), LibPETSc.FILE_MODE_WRITE)
    for vec in vecs
        LibPETSc.VecView(petsclib, vec, vtk)
    end
    vtk_ref = Ref{LibPETSc.PetscViewer}(vtk)
    LibPETSc.PetscViewerDestroy(petsclib, vtk_ref)
    return nothing
end
