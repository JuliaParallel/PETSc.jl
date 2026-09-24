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

See also: `PETSC_VIEWER_STDERR_WORLD`, `PetscViewerASCIIGetStderr()`
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

See also: `PETSC_VIEWER_STDERR_SELF`, `PetscViewerASCIIGetStderr()`
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

See also: `PETSC_VIEWER_STDOUT_WORLD`, `PetscViewerASCIIGetStdout()`
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

See also: `PETSC_VIEWER_STDOUT_SELF`, `PetscViewerASCIIGetStdout()`
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

# override for PetscSFBcastBegin; C signature: PetscSFBcastBegin(PetscSF sf, MPI_Datatype unit, const void *rootdata, void *leafdata, MPI_Op op)
"""
    PetscSFBcastBegin(petsclib::PetscLibType, sf::PetscSF, unit::MPI.Datatype, rootdata, leafdata, op::MPI.Op)

Begin broadcasting root values to leaves: for each leaf, `leafdata[leaf] = op(leafdata[leaf], rootdata[root])`.

`unit` is the MPI datatype of one entry (e.g. `MPI.Datatype(Float64)` or `MPI.Datatype(petsclib.PetscInt)`),
`rootdata` and `leafdata` are `Array`s (or raw pointers) of that type with the SF's root and leaf counts.
`op` is usually `MPI.REPLACE`; use `MPI.SUM` etc. to combine with the existing leaf values. Both arrays must stay alive until `PetscSFBcastEnd` returns.

See also: `PetscSF`, `PetscSFBcastEnd()`, `PetscSFSetGraph()`

# External Links
$(_doc_external("Vec/PetscSFBcastBegin"))
"""
function PetscSFBcastBegin(petsclib::PetscLibType, sf::PetscSF, unit::MPI.Datatype, rootdata::Union{Ptr, AbstractArray}, leafdata::Union{Ptr, AbstractArray}, op::MPI.Op) end

@for_petsc function PetscSFBcastBegin(petsclib::$UnionPetscLib, sf::PetscSF, unit::MPI.Datatype, rootdata::Union{Ptr, AbstractArray}, leafdata::Union{Ptr, AbstractArray}, op::MPI.Op)
    GC.@preserve rootdata leafdata begin
        @chk ccall(
            (:PetscSFBcastBegin, $petsc_library), PetscErrorCode,
            (PetscSF, MPI_Datatype, Ptr{Cvoid}, Ptr{Cvoid}, MPI_Op),
            sf, unit, rootdata isa Ptr ? rootdata : pointer(rootdata), leafdata isa Ptr ? leafdata : pointer(leafdata), op)
    end
    return nothing
end
# override for PetscSFBcastEnd; C signature: PetscSFBcastEnd(PetscSF sf, MPI_Datatype unit, const void *rootdata, void *leafdata, MPI_Op op)
"""
    PetscSFBcastEnd(petsclib::PetscLibType, sf::PetscSF, unit::MPI.Datatype, rootdata, leafdata, op::MPI.Op)

End a broadcast started with `PetscSFBcastBegin`; must be called with the same arguments.

`unit` is the MPI datatype of one entry (e.g. `MPI.Datatype(Float64)` or `MPI.Datatype(petsclib.PetscInt)`),
`rootdata` and `leafdata` are `Array`s (or raw pointers) of that type with the SF's root and leaf counts.


See also: `PetscSF`, `PetscSFBcastBegin()`, `PetscSFSetGraph()`

# External Links
$(_doc_external("Vec/PetscSFBcastEnd"))
"""
function PetscSFBcastEnd(petsclib::PetscLibType, sf::PetscSF, unit::MPI.Datatype, rootdata::Union{Ptr, AbstractArray}, leafdata::Union{Ptr, AbstractArray}, op::MPI.Op) end

@for_petsc function PetscSFBcastEnd(petsclib::$UnionPetscLib, sf::PetscSF, unit::MPI.Datatype, rootdata::Union{Ptr, AbstractArray}, leafdata::Union{Ptr, AbstractArray}, op::MPI.Op)
    GC.@preserve rootdata leafdata begin
        @chk ccall(
            (:PetscSFBcastEnd, $petsc_library), PetscErrorCode,
            (PetscSF, MPI_Datatype, Ptr{Cvoid}, Ptr{Cvoid}, MPI_Op),
            sf, unit, rootdata isa Ptr ? rootdata : pointer(rootdata), leafdata isa Ptr ? leafdata : pointer(leafdata), op)
    end
    return nothing
end
# override for PetscSFFetchAndOpBegin; C signature: PetscSFFetchAndOpBegin(PetscSF sf, MPI_Datatype unit, void *rootdata, const void *leafdata, void *leafupdate, MPI_Op op)
"""
    PetscSFFetchAndOpBegin(petsclib::PetscLibType, sf::PetscSF, unit::MPI.Datatype, rootdata, leafdata, leafupdate, op::MPI.Op)

Begin a fetch-and-op: every leaf fetches the current root value into `leafupdate` and then applies `rootdata[root] = op(rootdata[root], leafdata[leaf])`, atomically per root.

`unit` is the MPI datatype of one entry (e.g. `MPI.Datatype(Float64)`); `rootdata` has the SF's root count,
`leafdata` and `leafupdate` its leaf count. The arrays must stay alive until `PetscSFFetchAndOpEnd` returns.

See also: `PetscSF`, `PetscSFFetchAndOpEnd()`, `PetscSFReduceBegin()`

# External Links
$(_doc_external("Vec/PetscSFFetchAndOpBegin"))
"""
function PetscSFFetchAndOpBegin(petsclib::PetscLibType, sf::PetscSF, unit::MPI.Datatype, rootdata::Union{Ptr, AbstractArray}, leafdata::Union{Ptr, AbstractArray}, leafupdate::Union{Ptr, AbstractArray}, op::MPI.Op) end

@for_petsc function PetscSFFetchAndOpBegin(petsclib::$UnionPetscLib, sf::PetscSF, unit::MPI.Datatype, rootdata::Union{Ptr, AbstractArray}, leafdata::Union{Ptr, AbstractArray}, leafupdate::Union{Ptr, AbstractArray}, op::MPI.Op)
    GC.@preserve rootdata leafdata leafupdate begin
        @chk ccall(
            (:PetscSFFetchAndOpBegin, $petsc_library), PetscErrorCode,
            (PetscSF, MPI_Datatype, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, MPI_Op),
            sf, unit, rootdata isa Ptr ? rootdata : pointer(rootdata), leafdata isa Ptr ? leafdata : pointer(leafdata),
            leafupdate isa Ptr ? leafupdate : pointer(leafupdate), op)
    end
    return nothing
end
# override for PetscSFFetchAndOpEnd; C signature: PetscSFFetchAndOpEnd(PetscSF sf, MPI_Datatype unit, void *rootdata, const void *leafdata, void *leafupdate, MPI_Op op)
"""
    PetscSFFetchAndOpEnd(petsclib::PetscLibType, sf::PetscSF, unit::MPI.Datatype, rootdata, leafdata, leafupdate, op::MPI.Op)

End a fetch-and-op started with `PetscSFFetchAndOpBegin`; must be called with the same arguments.

`unit` is the MPI datatype of one entry (e.g. `MPI.Datatype(Float64)`); `rootdata` has the SF's root count,
`leafdata` and `leafupdate` its leaf count. The arrays must stay alive until `PetscSFFetchAndOpEnd` returns.

See also: `PetscSF`, `PetscSFFetchAndOpBegin()`, `PetscSFReduceBegin()`

# External Links
$(_doc_external("Vec/PetscSFFetchAndOpEnd"))
"""
function PetscSFFetchAndOpEnd(petsclib::PetscLibType, sf::PetscSF, unit::MPI.Datatype, rootdata::Union{Ptr, AbstractArray}, leafdata::Union{Ptr, AbstractArray}, leafupdate::Union{Ptr, AbstractArray}, op::MPI.Op) end

@for_petsc function PetscSFFetchAndOpEnd(petsclib::$UnionPetscLib, sf::PetscSF, unit::MPI.Datatype, rootdata::Union{Ptr, AbstractArray}, leafdata::Union{Ptr, AbstractArray}, leafupdate::Union{Ptr, AbstractArray}, op::MPI.Op)
    GC.@preserve rootdata leafdata leafupdate begin
        @chk ccall(
            (:PetscSFFetchAndOpEnd, $petsc_library), PetscErrorCode,
            (PetscSF, MPI_Datatype, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, MPI_Op),
            sf, unit, rootdata isa Ptr ? rootdata : pointer(rootdata), leafdata isa Ptr ? leafdata : pointer(leafdata),
            leafupdate isa Ptr ? leafupdate : pointer(leafupdate), op)
    end
    return nothing
end
# override for PetscSFReduceBegin; C signature: PetscSFReduceBegin(PetscSF sf, MPI_Datatype unit, const void *leafdata, void *rootdata, MPI_Op op)
"""
    PetscSFReduceBegin(petsclib::PetscLibType, sf::PetscSF, unit::MPI.Datatype, leafdata, rootdata, op::MPI.Op)

Begin reducing leaf values into their roots: `rootdata[root] = op(rootdata[root], leafdata[leaf])` over all leaves of a root.

`unit` is the MPI datatype of one entry (e.g. `MPI.Datatype(Float64)` or `MPI.Datatype(petsclib.PetscInt)`),
`leafdata` and `rootdata` are `Array`s (or raw pointers) of that type with the SF's root and leaf counts.
Typical `op`s are `MPI.SUM`, `MPI.MAX`, `MPI.MIN` and `MPI.REPLACE`. Both arrays must stay alive until `PetscSFReduceEnd` returns.

See also: `PetscSF`, `PetscSFReduceEnd()`, `PetscSFSetGraph()`

# External Links
$(_doc_external("Vec/PetscSFReduceBegin"))
"""
function PetscSFReduceBegin(petsclib::PetscLibType, sf::PetscSF, unit::MPI.Datatype, leafdata::Union{Ptr, AbstractArray}, rootdata::Union{Ptr, AbstractArray}, op::MPI.Op) end

@for_petsc function PetscSFReduceBegin(petsclib::$UnionPetscLib, sf::PetscSF, unit::MPI.Datatype, leafdata::Union{Ptr, AbstractArray}, rootdata::Union{Ptr, AbstractArray}, op::MPI.Op)
    GC.@preserve leafdata rootdata begin
        @chk ccall(
            (:PetscSFReduceBegin, $petsc_library), PetscErrorCode,
            (PetscSF, MPI_Datatype, Ptr{Cvoid}, Ptr{Cvoid}, MPI_Op),
            sf, unit, leafdata isa Ptr ? leafdata : pointer(leafdata), rootdata isa Ptr ? rootdata : pointer(rootdata), op)
    end
    return nothing
end
# override for PetscSFReduceEnd; C signature: PetscSFReduceEnd(PetscSF sf, MPI_Datatype unit, const void *leafdata, void *rootdata, MPI_Op op)
"""
    PetscSFReduceEnd(petsclib::PetscLibType, sf::PetscSF, unit::MPI.Datatype, leafdata, rootdata, op::MPI.Op)

End a reduction started with `PetscSFReduceBegin`; must be called with the same arguments.

`unit` is the MPI datatype of one entry (e.g. `MPI.Datatype(Float64)` or `MPI.Datatype(petsclib.PetscInt)`),
`leafdata` and `rootdata` are `Array`s (or raw pointers) of that type with the SF's root and leaf counts.


See also: `PetscSF`, `PetscSFReduceBegin()`, `PetscSFSetGraph()`

# External Links
$(_doc_external("Vec/PetscSFReduceEnd"))
"""
function PetscSFReduceEnd(petsclib::PetscLibType, sf::PetscSF, unit::MPI.Datatype, leafdata::Union{Ptr, AbstractArray}, rootdata::Union{Ptr, AbstractArray}, op::MPI.Op) end

@for_petsc function PetscSFReduceEnd(petsclib::$UnionPetscLib, sf::PetscSF, unit::MPI.Datatype, leafdata::Union{Ptr, AbstractArray}, rootdata::Union{Ptr, AbstractArray}, op::MPI.Op)
    GC.@preserve leafdata rootdata begin
        @chk ccall(
            (:PetscSFReduceEnd, $petsc_library), PetscErrorCode,
            (PetscSF, MPI_Datatype, Ptr{Cvoid}, Ptr{Cvoid}, MPI_Op),
            sf, unit, leafdata isa Ptr ? leafdata : pointer(leafdata), rootdata isa Ptr ? rootdata : pointer(rootdata), op)
    end
    return nothing
end
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
    return PetscMat(J_ref[], petsclib; own = false)
end

