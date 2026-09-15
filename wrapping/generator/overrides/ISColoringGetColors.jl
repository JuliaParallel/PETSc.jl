# override for ISColoringGetColors; C signature: ISColoringGetColors(ISColoring iscoloring, PetscInt* n, PetscInt* nc, ISColoringValue** colors)
"""
	n::PetscInt, nc::PetscInt, colors::Vector{UInt16} = ISColoringGetColors(petsclib::PetscLibType, iscoloring::ISColoring)
Returns an array with the color for each local node

Not Collective

Input Parameter:
- `iscoloring` - the coloring context

Output Parameters:
- `n`      - number of nodes (DOFs)
- `nc`     - number of colors
- `colors` - copy of the color array (one `UInt16` entry per DOF)

Level: advanced

-seealso: `ISColoring`, `ISColoringValue`, `ISColoringRestoreIS()`, `ISColoringView()`, `ISColoringGetIS()`

# External Links
$(_doc_external("Vec/ISColoringGetColors"))
"""
function ISColoringGetColors(petsclib::PetscLibType, iscoloring::ISColoring) end

@for_petsc function ISColoringGetColors(petsclib::$UnionPetscLib, iscoloring::ISColoring)
    n_  = Ref{$PetscInt}()
    nc_ = Ref{$PetscInt}()
    # PETSc returns a pointer into its own storage; we must not free it.
    colors_ptr_ = Ref{ISColoringValue}(C_NULL)

    @chk ccall(
               (:ISColoringGetColors, $petsc_library),
               PetscErrorCode,
               (ISColoring, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{ISColoringValue}),
               iscoloring, n_, nc_, colors_ptr_,
              )

    n  = n_[]
    nc = nc_[]
    # ISColoringValue is a PETSc opaque pointer alias for `unsigned short *`.
    # Reinterpret as Ptr{UInt16} and copy to a Julia-owned Vector so the
    # caller can safely use the data after ISColoringDestroy.
    colors_raw = Ptr{UInt16}(UInt(colors_ptr_[]))
    colors = copy(unsafe_wrap(Vector{UInt16}, colors_raw, Int(n); own = false))

    return n, nc, colors
end

