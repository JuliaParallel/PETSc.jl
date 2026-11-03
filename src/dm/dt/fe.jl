"""
    fe_create_default(petsclib, comm, dim, Nc, simplex; degree = 1, prefix = "", qorder = -1)

Create a default `PetscFE` (Lagrange) for a single field with `Nc` components on
either a simplex (`simplex == true`) or tensor cell.

`degree` is the polynomial degree (1 = P1/Q1, 2 = P2/Q2, etc.).  It is injected
into the options database under `prefix` as `-\${prefix}petscspace_degree` before
the call to `PetscFECreateDefault`, so any options explicitly set by the user
will take precedence.

# External Links
$(doc_external("FE/PetscFECreateDefault"))
"""
function fe_create_default(
    petsclib::PetscLib,
    comm::MPI.Comm,
    dim::Integer,
    Nc::Integer,
    simplex::Bool;
    degree::Integer  = 1,
    prefix::AbstractString = "",
    qorder::Integer = -1,
) where {PetscLib}
    PetscInt = inttype(PetscLib)
    # PetscFECreateDefault reads the degree from the options database only
    entries = ["-$(prefix)petscspace_degree" => string(degree)]
    fe = with_layered_options(PetscLib, entries) do
        LibPETSc.PetscFECreateDefault(
            petsclib, comm,
            PetscInt(dim), PetscInt(Nc),
            LibPETSc.PetscBool(simplex),
            String(prefix),
            PetscInt(qorder),
        )
    end
    return own_serial!(fe, comm)
end

"""
    fe_create_lagrange(petsclib, comm, dim, Nc, simplex, degree; qorder = -1)

Create a `PetscFE` Lagrange element of polynomial degree `degree` for a single
field with `Nc` components on either a simplex (`simplex == true`) or tensor cell.
Unlike `fe_create_default`, the degree is specified explicitly and does not depend
on the options database.

# External Links
$(doc_external("FE/PetscFECreateLagrange"))
"""
function fe_create_lagrange(
    petsclib::PetscLib,
    comm::MPI.Comm,
    dim::Integer,
    Nc::Integer,
    simplex::Bool,
    degree::Integer;
    qorder::Integer = -1,
) where {PetscLib}
    PetscInt = inttype(PetscLib)
    fe = LibPETSc.PetscFECreateLagrange(
        petsclib, comm,
        PetscInt(dim), PetscInt(Nc),
        LibPETSc.PetscBool(simplex),
        PetscInt(degree),
        PetscInt(qorder),
    )
    return own_serial!(fe, comm)
end

"""
    copy_quadrature!(dst_fe, src_fe)

Copy the quadrature rule from `src_fe` to `dst_fe` so both fields share
the same integration points, and return `dst_fe`.

# External Links
$(doc_external("FE/PetscFECopyQuadrature"))
"""
function copy_quadrature!(
    dst_fe::LibPETSc.AbstractPetscFE{PetscLib},
    src_fe::LibPETSc.AbstractPetscFE{PetscLib},
) where {PetscLib}
    LibPETSc.PetscFECopyQuadrature(getlib(PetscLib), src_fe, dst_fe)
    return dst_fe
end

"""
    compose_constant_nullspace!(obj)

Create a trivial (constant) `MatNullSpace` and attach it to `obj`, a `PetscFE` or
an index set, under the key `"nullspace"` via `PetscObjectCompose`.  This signals to
the fieldsplit preconditioner that the field or split has a constant null space.
Returns `obj`.
"""
function compose_constant_nullspace!(
    obj::Union{LibPETSc.AbstractPetscFE{PetscLib}, LibPETSc.AbstractIS{PetscLib}},
) where {PetscLib}
    petsclib = getlib(PetscLib)
    nsp = mat_nullspace_create(petsclib, LibPETSc.PetscObjectGetComm(petsclib, obj); has_const = true)
    try
        LibPETSc.PetscObjectCompose(petsclib,
            convert(Ptr{Cvoid}, obj), "nullspace", convert(Ptr{Cvoid}, nsp))
    finally
        destroy!(nsp)   # `obj` keeps its own reference
    end
    return obj
end

"""
    destroy!(fe::PetscFE)

Destroy a `PetscFE` created by [`fe_create_default`](@ref) or
[`fe_create_lagrange`](@ref). Does nothing on a borrowed handle.

# External Links
$(doc_external("FE/PetscFEDestroy"))
"""
function destroy!(fe::LibPETSc.AbstractPetscFE{PetscLib}) where {PetscLib}
    return destroy_handle!(LibPETSc.PetscFEDestroy, fe, PetscLib)
end
