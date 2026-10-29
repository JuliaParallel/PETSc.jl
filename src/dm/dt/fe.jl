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
    # Push the polynomial degree into the options DB so PetscFECreateDefault
    # picks up the correct space even when no command-line option is given.
    opt_key = isempty(prefix) ? "petscspace_degree" : "$(prefix)petscspace_degree"
    opts = PetscOptions(petsclib; Symbol(opt_key) => degree)
    return with_options(opts) do
        LibPETSc.PetscFECreateDefault(
            petsclib, comm,
            PetscInt(dim), PetscInt(Nc),
            LibPETSc.PetscBool(simplex),
            String(prefix),
            PetscInt(qorder),
        )
    end
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
    return LibPETSc.PetscFECreateLagrange(
        petsclib, comm,
        PetscInt(dim), PetscInt(Nc),
        LibPETSc.PetscBool(simplex),
        PetscInt(degree),
        PetscInt(qorder),
    )
end

"""
    copy_quadrature!(petsclib, src_fe, dst_fe)

Copy the quadrature rule from `src_fe` to `dst_fe` so both fields share
the same integration points.
"""
function copy_quadrature!(petsclib, src_fe, dst_fe)
    LibPETSc.PetscFECopyQuadrature(petsclib, src_fe, dst_fe)
    return dst_fe
end

"""
    compose_constant_nullspace!(petsclib, comm, fe)

Create a trivial (constant) `MatNullSpace` and attach it to the FE object `fe`
under the key `"nullspace"` via `PetscObjectCompose`.  This signals to the
fieldsplit preconditioner that the field has a constant null space.
"""
function compose_constant_nullspace! end

LibPETSc.@for_petsc function compose_constant_nullspace!(
    petsclib::$UnionPetscLib,
    comm,
    fe,
)
    nsp = mat_nullspace_create(petsclib, comm; has_const = true)
    LibPETSc.PetscObjectCompose(petsclib,
        convert(Ptr{Cvoid}, fe), "nullspace", convert(Ptr{Cvoid}, nsp))
    destroy!(petsclib, nsp)
    return fe
end
