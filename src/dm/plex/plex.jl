# ── DMPlex constructors ──────────────────────────────────────────────────────

"""
    DMPlex(
        petsclib::PetscLib,
        comm::MPI.Comm;
        set_from_options = true,
        setup            = false,
        prefix           = "",
        options...,
    )

Create an empty `DMPLEX` mesh.  Use the keyword `options...` (or PETSc command-line
options) such as `dm_plex_box_faces`, `dm_plex_dim`, `dm_plex_simplex`,
`dm_plex_filename`, etc. to populate the mesh — this is the most flexible entry
point and matches the `DMCreate` + `DMSetType(DMPLEX)` + `DMSetFromOptions` flow
used in PETSc's tutorials (e.g. `snes/ex12.c`).

For an explicit box mesh, use the
[`DMPlex`](@ref) method taking `(petsclib, comm, dim, simplex, faces; ...)`,
below.

# External Links
$(doc_external("DMPlex/DMPlexCreate"))
$(doc_external("DM/DMSetType"))
$(doc_external("DM/DMSetFromOptions"))
"""
function DMPlex(
    petsclib::PetscLib,
    comm::MPI.Comm;
    set_from_options::Bool = true,
    setup::Bool = false,
    prefix::AbstractString = "",
    options...,
) where {PetscLib}
    check_initialized(getlib(PetscLib))

    dm = LibPETSc.DMCreate(petsclib, comm)
    LibPETSc.DMSetType(petsclib, dm, "plex")
    dm = DMPlex{PetscLib}(dm.ptr, getlib(PetscLib).age, true)
    return finish_dm!(dm, comm, options; prefix, set_from_options, setup)
end

"""
    DMPlex(
        petsclib::PetscLib,
        comm::MPI.Comm,
        dim::Integer,
        simplex::Bool,
        faces::AbstractVector{<:Integer};
        lower         = ntuple(_ -> 0.0, dim),
        upper         = ntuple(_ -> 1.0, dim),
        periodicity   = ntuple(_ -> DM_BOUNDARY_NONE, dim),
        interpolate   = true,
        localize_height = 0,
        sparse_localize = false,
        set_from_options = true,
        setup           = false,
        prefix          = "",
        options...,
    )

Create a structured-grid box `DMPLEX` mesh in `dim` dimensions with `faces[d]`
cells along axis `d`.  When `simplex == true` each box cell is split into
simplices (triangles in 2D, tetrahedra in 3D); otherwise tensor-product cells
are used.

# External Links
$(doc_external("DMPlex/DMPlexCreateBoxMesh"))
"""
function DMPlex(
    petsclib::PetscLib,
    comm::MPI.Comm,
    dim::Integer,
    simplex::Bool,
    faces::AbstractVector{<:Integer};
    lower            = ntuple(_ -> 0.0, dim),
    upper            = ntuple(_ -> 1.0, dim),
    periodicity      = ntuple(_ -> DM_BOUNDARY_NONE, dim),
    interpolate      = true,
    localize_height  = 0,
    sparse_localize  = false,
    set_from_options::Bool = true,
    setup::Bool      = false,
    prefix::AbstractString = "",
    options...,
) where {PetscLib}
    check_initialized(getlib(PetscLib))
    for (name, v) in
        (("faces", faces), ("lower", lower), ("upper", upper), ("periodicity", periodicity))
        length(v) == dim || throw(
            DimensionMismatch("$name has length $(length(v)), expected dim = $dim"),
        )
    end

    PetscInt  = inttype(PetscLib)
    PetscReal = real(scalartype(PetscLib))

    faces_v       = collect(PetscInt.(faces))
    lower_v       = collect(PetscReal.(lower))
    upper_v       = collect(PetscReal.(upper))
    periodicity_v = collect(periodicity)

    dm = LibPETSc.DMPlexCreateBoxMesh(
        petsclib, comm,
        PetscInt(dim),
        LibPETSc.PetscBool(simplex),
        faces_v, lower_v, upper_v, periodicity_v,
        LibPETSc.PetscBool(interpolate),
        PetscInt(localize_height),
        LibPETSc.PetscBool(sparse_localize),
    )
    dm = DMPlex{PetscLib}(dm.ptr, getlib(PetscLib).age, true)
    return finish_dm!(dm, comm, options; prefix, set_from_options, setup)
end


# ── Convenience helpers ──────────────────────────────────────────────────────

"""
    issimplex(dm::DMPlex) -> Bool

Return `true` when the cells of `dm` (assumed to be a `DMPLEX`) are simplices.

# External Links
$(doc_external("DMPlex/DMPlexIsSimplex"))
"""
function issimplex(dm::DMPlex{PetscLib}) where {PetscLib}
    return Bool(LibPETSc.DMPlexIsSimplex(getlib(PetscLib), dm))
end

"""
    distribute(dm::DMPlex; overlap = 0) -> DMPlex

Distribute `dm` across its communicator, with `overlap` layers of cells shared
between neighboring ranks, and return the distributed mesh as a new, owned DM.
`dm` itself is left unchanged and still has to be destroyed. When PETSc doesn't
redistribute (on one rank) the returned DM has a null pointer, and `destroy!` on
it does nothing.

# External Links
$(doc_external("DMPlex/DMPlexDistribute"))
"""
function distribute(
    dm::DMPlex{PetscLib};
    overlap::Integer = 0,
) where {PetscLib}
    petsclib = getlib(PetscLib)
    PetscInt = inttype(PetscLib)
    _, dm_par = LibPETSc.DMPlexDistribute(petsclib, dm, PetscInt(overlap))
    # The distributed DM is a new object handed to the caller, not a borrowed
    # handle, so it is narrowed as owned.
    return narrow(dm_par; own = true)
end


# ── Wire DMPlex FEM residual/Jacobian assembly into SNES ─────────────────────

"""
    set_snes_local_fem!(dm; use_obj = false, ctx = C_NULL)

Tell SNES to use `DMPlex`'s built-in FEM residual / Jacobian assembly on `dm`.

`petsclib` is not an argument: the `dm` carries it as a type parameter (§8).
The v0.4 name `plex_set_snes_local_fem!` shims the old form.

# External Links
$(doc_external("DMPlex/DMPlexSetSNESLocalFEM"))
"""
function set_snes_local_fem! end

LibPETSc.@for_petsc function set_snes_local_fem!(
    dm::AbstractPetscDM{$PetscLib};
    use_obj::Bool = false,
    ctx::Ptr{Cvoid} = C_NULL,
)
    petsclib = getlib($PetscLib)
    LibPETSc.DMPlexSetSNESLocalFEM(petsclib, dm, LibPETSc.PetscBool(use_obj), ctx)
    return dm
end
