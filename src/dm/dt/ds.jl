abstract type AbstractPetscDS{T} end

mutable struct PetscDS{PetscLib} <: AbstractPetscDS{PetscLib}
    ptr::LibPETSc.PetscDS   # Ptr{_n_PetscDS}
end

Base.unsafe_convert(::Type{LibPETSc.PetscDS}, v::PetscDS) = v.ptr


# ── PetscDS point-function setters ──────────────────────────────────────────

"""
    set_residual!(ds, field, f0_ptr, f1_ptr)

Attach the pointwise residual point-function pair to field `field` (0-based) of
the `PetscDS` `ds`.  `f0_ptr` and `f1_ptr` must be `@cfunction(...)`-style
function pointers (`Ptr{Cvoid}`) matching PETSc's `PetscPointFn` signature.

# External Links
$(doc_external("Dm/PetscDSSetResidual"))
"""
function set_residual!(
    ds::PetscDS{PetscLib},
    field::Integer,
    f0_ptr::Ptr{Cvoid},
    f1_ptr::Ptr{Cvoid},
) where {PetscLib}
    petsclib = getlib(PetscLib)
    LibPETSc.PetscDSSetResidual(petsclib, ds.ptr, PetscLib.PetscInt(field), f0_ptr, f1_ptr)
    return ds
end

"""
    set_jacobian!(ds, fieldI, fieldJ, g0_ptr, g1_ptr, g2_ptr, g3_ptr)

Attach the pointwise Jacobian point-function quadruple for the block
`(fieldI, fieldJ)` of the `PetscDS` `ds`.  Pass `C_NULL` for any unused term.

# External Links
$(doc_external("Dm/PetscDSSetJacobian"))
"""
function set_jacobian!(
    ds::PetscDS{PetscLib},
    fieldI::Integer,
    fieldJ::Integer,
    g0_ptr::Ptr{Cvoid},
    g1_ptr::Ptr{Cvoid},
    g2_ptr::Ptr{Cvoid},
    g3_ptr::Ptr{Cvoid},
) where {PetscLib}
    petsclib = getlib(PetscLib)
    LibPETSc.PetscDSSetJacobian(petsclib, ds.ptr, PetscLib.PetscInt(fieldI), PetscLib.PetscInt(fieldJ),
                                g0_ptr, g1_ptr, g2_ptr, g3_ptr)
    return ds
end

"""
    set_jacobian_preconditioner!(ds, field_i, field_j, g0, g1, g2, g3)

Set the Jacobian preconditioner terms for the weak form between fields `field_i` and `field_j`.
Each `gN` argument is either a `Ptr{Cvoid}` function pointer or `C_NULL`.

This is the preconditioner equivalent of [`set_jacobian!`](@ref).

# External Links
$(doc_external("Dm/PetscDSSetJacobianPreconditioner"))
"""
function set_jacobian_preconditioner!(
    ds::PetscDS{PetscLib},
    field_i::Integer, field_j::Integer,
    g0::Ptr{Cvoid}, g1::Ptr{Cvoid}, g2::Ptr{Cvoid}, g3::Ptr{Cvoid},
) where {PetscLib}
    LibPETSc.PetscDSSetJacobianPreconditioner(
        getlib(PetscLib), ds.ptr,
        PetscLib.PetscInt(field_i), PetscLib.PetscInt(field_j),
        g0, g1, g2, g3,
    )
    return ds
end

"""
    set_exact_solution!(ds, field, sol_ptr, ctx = C_NULL)

Register the pointwise exact-solution function `sol_ptr` for field `field`
(0-based) of the `PetscDS` `ds`.  The function pointer must match PETSc's
`PetscSimplePointFn` / `PetscErrorCode f(dim, time, x, Nc, u, ctx)` signature.

This is required for `DMComputeL2Diff` and `DMComputeExactSolution` to work.

# External Links
$(doc_external("Dm/PetscDSSetExactSolution"))
"""
function set_exact_solution!(
    ds::PetscDS{PetscLib},
    field::Integer,
    sol_ptr::Ptr{Cvoid},
    ctx::Ptr{Cvoid} = C_NULL,
) where {PetscLib}
    petsclib = getlib(PetscLib)
    LibPETSc.PetscDSSetExactSolution(petsclib, ds.ptr, PetscLib.PetscInt(field), sol_ptr, ctx)
    return ds
end

"""
    set_constants!(ds::PetscDS, constants)

Set the named constants (accessible as `constants[i]` in pointwise functions) on the
`PetscDS` `ds`.  `constants` must be an `AbstractVector` whose elements are convertible
to `PetscScalar`.

# External Links
$(doc_external("Dm/PetscDSSetConstants"))
"""
function set_constants!(ds::PetscDS{PetscLib}, constants::AbstractVector) where {PetscLib}
    petsclib = getlib(PetscLib)
    cst_v = collect(PetscLib.PetscScalar, constants)
    LibPETSc.PetscDSSetConstants(petsclib, ds.ptr, PetscLib.PetscInt(length(cst_v)), cst_v)
    return ds
end
