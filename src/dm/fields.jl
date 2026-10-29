"""
    ds(dm::AbstractPetscDM) -> PetscDS

Return the `PetscDS` (discrete system) attached to `dm`.

$(doc_borrowed())

# External Links
$(doc_external("DM/DMGetDS"))
"""
function ds end

LibPETSc.@for_petsc function ds(dm::AbstractPetscDM{$PetscLib})
    return PetscDS{$PetscLib}(LibPETSc.DMGetDS(getlib($PetscLib), dm))
end

"""
    create_ds!(dm::AbstractPetscDM)

Build the `PetscDS` for `dm` from the currently-attached fields.

# External Links
$(doc_external("DM/DMCreateDS"))
"""
function create_ds!(dm::AbstractPetscDM{PetscLib}) where {PetscLib}
    LibPETSc.DMCreateDS(getlib(PetscLib), dm)
    return dm
end

"""
    set_field!(dm::AbstractPetscDM, field::Integer, fe; label = C_NULL)

Attach a `PetscFE` (or other discretization object) as the `field`-th field of
`dm` (0-based field index, matching PETSc).

# External Links
$(doc_external("DM/DMSetField"))
"""
function set_field! end

LibPETSc.@for_petsc function set_field!(
    dm::AbstractPetscDM{$PetscLib},
    field::Integer,
    fe::Ptr;
    label::Ptr{Cvoid} = C_NULL,
)
    LibPETSc.DMSetField(getlib($PetscLib), dm, $PetscInt(field), label, Ptr{Cvoid}(fe))
    return dm
end

"""
    add_boundary!(dm, bctype, name, label, values, field, comps,
                  bcfunc_ptr, bcfunc_t_ptr = C_NULL, ctx = C_NULL)

Attach a boundary condition to `dm`, and return `dm`. `LibPETSc.DMAddBoundary`
also returns the number PETSc gives the boundary.  `label` is a `Ptr{Cvoid}` from
[`label`](@ref).  `bcfunc_ptr` / `bcfunc_t_ptr` are `@cfunction(...)` results
matching PETSc's boundary callback signature.

`petsclib` is not an argument: the `dm` carries it as a type parameter (§8).
The v0.4 form with `petsclib` first still works in v0.5 and warns.

The `values` are label values PETSc issued and `field` is the number
[`set_field!`](@ref) gave the field, so both keep PETSc's numbering (§12.1).

# External Links
$(doc_external("DM/DMAddBoundary"))
"""
function add_boundary! end

LibPETSc.@for_petsc function add_boundary!(
    dm::AbstractPetscDM{$PetscLib},
    bctype::LibPETSc.DMBoundaryConditionType,
    name::AbstractString,
    label::Ptr{Cvoid},
    values::AbstractVector{<:Integer},
    field::Integer,
    comps::AbstractVector{<:Integer},
    bcfunc_ptr::Ptr{Cvoid},
    bcfunc_t_ptr::Ptr{Cvoid} = C_NULL,
    ctx::Ptr{Cvoid}          = C_NULL,
)
    petsclib = getlib($PetscLib)
    values_v = collect($PetscInt.(values))
    comps_v  = collect($PetscInt.(comps))
    LibPETSc.DMAddBoundary(
        petsclib, dm, bctype, String(name), label,
        $PetscInt(length(values_v)), values_v,
        $PetscInt(field),
        $PetscInt(length(comps_v)), comps_v,
        bcfunc_ptr, bcfunc_t_ptr, ctx,
    )
    return dm
end


# ── DMProjectFunction / DMComputeL2Diff ─────────────────────────────────────

"""
    project_function!(X, dm, time, funcs; ctxs = nothing, mode = INSERT_VALUES)

Project a collection of pointwise functions into the global vector `X` and
return `X`.

The written vector comes first and `petsclib` is gone (§8); v0.4 took `X` last
and `petsclib` first. A shim forwards the v0.4 name `dm_project_function!`
from the old order (§17.2).
`funcs` holds one `@cfunction` pointer per field, matching PETSc's
`PetscSimplePointFn` signature. `ctxs` is a matching vector of context pointers,
or `nothing` to use `C_NULL` for every field.

`mode` is an `InsertMode`: `INSERT_VALUES` sets the free DOFs only,
`INSERT_ALL_VALUES` the constrained ones too.

# External Links
$(doc_external("DM/DMProjectFunction"))
"""
function project_function!(
    X::AbstractPetscVec{PetscLib},
    dm::AbstractPetscDM{PetscLib},
    time::Real,
    funcs::AbstractVector;
    ctxs::Union{AbstractVector, Nothing} = nothing,
    mode::InsertMode = INSERT_VALUES,
) where {PetscLib}
    fptrs = Ptr{Cvoid}[Ptr{Cvoid}(f) for f in funcs]
    cptrs = ctxs === nothing ?
        fill(C_NULL, length(fptrs)) :
        Ptr{Cvoid}[Ptr{Cvoid}(c) for c in ctxs]
    length(cptrs) == length(fptrs) || throw(
        DimensionMismatch("got $(length(fptrs)) funcs but $(length(cptrs)) ctxs"),
    )
    LibPETSc.DMProjectFunction(getlib(PetscLib), dm, time, fptrs, cptrs, mode, X)
    return X
end

"""
    l2diff(dm, time, funcs, ctxs, X) -> PetscReal

Compute the L² error between the global vector `X` and the pointwise exact
functions `funcs` (one per field, matching `PetscSimplePointFn` signature).
`ctxs` is a matching vector of context pointers, or `nothing` to use `C_NULL`
for every field.  Requires that [`set_exact_solution!`](@ref) has been called.

`petsclib` is not an argument: the `dm` carries it as a type parameter (§8).
The v0.4 name `dm_compute_l2diff` shims the old form.

# External Links
$(doc_external("DM/DMComputeL2Diff"))
"""
function l2diff end

LibPETSc.@for_petsc function l2diff(
    dm::AbstractPetscDM{$PetscLib},
    time::Real,
    funcs::AbstractVector{Ptr{Cvoid}},
    ctxs::AbstractVector{Ptr{Cvoid}},
    X::AbstractPetscVec{$PetscLib},
)
    petsclib = getlib($PetscLib)
    n = length(funcs)
    length(ctxs) == n || throw(
        DimensionMismatch("got $n funcs but $(length(ctxs)) ctxs"),
    )
    funcs_v = collect(funcs)
    ctxs_v  = collect(ctxs)
    return LibPETSc.DMComputeL2Diff(petsclib, dm, time, funcs_v, ctxs_v, X)
end

function l2diff(
    dm::AbstractPetscDM,
    time::Real,
    funcs::AbstractVector,
    ctxs::Union{AbstractVector, Nothing},
    X::AbstractPetscVec,
)
    fptrs = Ptr{Cvoid}[Ptr{Cvoid}(f) for f in funcs]
    cptrs = ctxs === nothing ?
        fill(C_NULL, length(fptrs)) :
        Ptr{Cvoid}[Ptr{Cvoid}(c) for c in ctxs]
    return l2diff(dm, time, fptrs, cptrs, X)
end

# `dm_create_global_vec`/`dm_create_local_vec` here and `DMGlobalVec`/`DMLocalVec` in
# dm.jl were the same two calls under two names. Both pairs collapse onto
# `global_vec`/`local_vec`, defined once in src/dm.jl.

# `global_to_local!(dm, gvec, lvec; mode)` lived here as a second shape of the
# same call. §8 collapses both onto `global_to_local!(lvec, dm, gvec; mode)`,
# defined once in src/dm.jl.

"""
    set_auxiliary_vec!(dm, aux_local)

Attach a local auxiliary vector `aux_local` to `dm` (global label / value 0 / part 0).
The auxiliary field values are forwarded to all pointwise functions as the `a` argument.
"""
function set_auxiliary_vec!(dm::AbstractPetscDM{PetscLib}, aux_local) where {PetscLib}
    petsclib = getlib(PetscLib)
    LibPETSc.DMSetAuxiliaryVec(petsclib, dm,
        Ptr{Cvoid}(C_NULL),
        PetscLib.PetscInt(0), PetscLib.PetscInt(0),
        aux_local)
    return dm
end

"""
    project_field!(X, dm, time, U, funcs; mode = INSERT_VALUES)

Project a function of the fields in the input vector `U` into the FE space of `dm`,
writing the result into `X`, and return `X`.

The written vector comes first and `petsclib` is gone (§8); v0.4 took `X` last
and `petsclib` first. A shim forwards the v0.4 name `dm_project_field!` from
the old order (§17.2).  `funcs` is a vector of `Ptr{Cvoid}` function pointers
(one per field in `dm`), each with the `PetscPointFn` / `@residual_fn` signature.
`U` must be associated with a DM that shares the same mesh as `dm` (e.g. obtained via
[`clone`](@ref)).

Use this to compute derived quantities (e.g. stress from displacement gradient) and
project them onto a new FE field for visualisation. `mode` is an `InsertMode`, as
for [`project_function!`](@ref).

# External Links
$(doc_external("DM/DMProjectField"))
"""
function project_field!(
    X::AbstractPetscVec{PetscLib},
    dm::AbstractPetscDM{PetscLib},
    time::Real,
    U::AbstractPetscVec{PetscLib},
    funcs;
    mode::InsertMode = INSERT_VALUES,
) where {PetscLib}
    fptrs = Ptr{Cvoid}[Ptr{Cvoid}(f) for f in funcs]
    GC.@preserve fptrs LibPETSc.DMProjectField(
        getlib(PetscLib), dm, PetscLib.PetscReal(time), U, pointer(fptrs), mode, X,
    )
    return X
end

"""
    add_natural_boundary!(dm, ds, name, label, label_value, field, f0_ptr, f1_ptr = C_NULL)

Register a natural (Neumann) boundary condition on `dm` following the PETSc 3.22+
pattern used in the C tutorials, and return `dm`.  Equivalent to:
```c
DMAddBoundary(dm, DM_BC_NATURAL, name, label, 1, &val, field, 0, NULL, NULL, NULL, NULL, &bd);
PetscDSGetBoundary(ds, bd, &wf, NULL, ...);
PetscWeakFormSetIndexBdResidual(wf, label, val, field, 0, 0, f0, 0, NULL);
```

`f0_ptr` is a `@bd_fn`-generated C-callable function pointer implementing
the boundary integrand (`PetscBdPointFn` signature, with the outward normal `n[]`
between `x[]` and `numConstants`).

`petsclib` is not an argument: the `dm` carries it as a type parameter (§8).
The v0.4 form with `petsclib` first still works in v0.5 and warns.

# External Links
$(doc_external("DM/DMAddBoundary"))
$(doc_external("Dm/PetscWeakFormSetIndexBdResidual"))
"""
function add_natural_boundary! end

LibPETSc.@for_petsc function add_natural_boundary!(
    dm::AbstractPetscDM{$PetscLib},
    ds::PetscDS{$PetscLib},
    name::AbstractString,
    label::Ptr{Cvoid},
    label_value::Integer,
    field::Integer,
    f0_ptr::Ptr{Cvoid},
    f1_ptr::Ptr{Cvoid} = C_NULL,
)
    petsclib = getlib($PetscLib)
    # Step 1: register the boundary (NULL function — function set via weak form below)
    bd = LibPETSc.DMAddBoundary(petsclib, dm, LibPETSc.DM_BC_NATURAL, String(name),
        label, $PetscInt(1), $PetscInt[label_value],
        $PetscInt(field), $PetscInt(0), $PetscInt[],
        C_NULL, C_NULL, C_NULL)
    # Step 2: get the per-boundary PetscWeakForm
    wf = LibPETSc.PetscDSGetBoundary(petsclib, ds.ptr, bd).wf
    # Step 3: install the boundary integrand on the per-boundary weak form
    LibPETSc.PetscWeakFormSetIndexBdResidual(petsclib, wf, label,
        $PetscInt(label_value), $PetscInt(field),
        $PetscInt(0), $PetscInt(0), f0_ptr, $PetscInt(0), f1_ptr)
    return dm
end

"""
    copy_disc!(src::AbstractPetscDM, dst::AbstractPetscDM)

Copy the discretisation (fields, `PetscDS`, BCs) from `src` to `dst`.  Useful when
propagating FEM setup to coarser levels of a multigrid hierarchy.

# External Links
$(doc_external("DM/DMCopyDisc"))
"""
function copy_disc! end

LibPETSc.@for_petsc function copy_disc!(
    src::AbstractPetscDM{$PetscLib},
    dst::AbstractPetscDM{$PetscLib},
)
    LibPETSc.DMCopyDisc(getlib($PetscLib), src, dst)
    return dst
end
