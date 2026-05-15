import .LibPETSc: AbstractPetscDM, PetscDM, CDM, getlib

abstract type AbstractPetscDS{T} end

mutable struct PetscDS{PetscLib} <: AbstractPetscDS{PetscLib}
    ptr::LibPETSc.PetscDS   # Ptr{_n_PetscDS}
end

Base.unsafe_convert(::Type{LibPETSc.PetscDS}, v::PetscDS) = v.ptr

# ── @cfunction helpers for PETSc FEM pointwise callbacks ─────────────────────
#
# These macros generate a thin C-callable wrapper that unwraps all C pointer
# arguments into Julia Vectors, then calls the user's pure-Julia function.
#
# PetscInt, PetscReal, PetscScalar must be concrete types in the calling scope
# (e.g. `const PetscInt = petsclib.PetscInt`).
#
# Usage:
#   function my_fn(t, x, u, ctx);  u[1] = x[1]^2 + x[2]^2;  end
#   const my_fn_ptr = PETSc.@petsc_simple_fn(my_fn)
#
#   function f0(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
#               aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, out)
#       out .= 4
#   end
#   const f0_ptr = PETSc.@petsc_residual_fn(f0, Nf)        # outsz = Nf
#   const f1_ptr = PETSc.@petsc_residual_fn(f1, dim_*Nf)   # outsz = dim*Nf
#   const g3_ptr = PETSc.@petsc_jacobian_fn(g3, dim_*dim_) # outsz = dim² (Nf=1)

# Helper: substitute gensymmed param names into an output-size expression.
_petsc_subst(ex, m) =
    ex isa Symbol ? get(m, ex, ex) :
    ex isa Expr   ? Expr(ex.head, map(a -> _petsc_subst(a, m), ex.args)...) :
    ex

"""
    @petsc_simple_fn(f)

Generate a C-callable wrapper for `f(t, x, u, ctx)` and return its C pointer,
matching the `PetscSimplePointFn` signature.  Arguments seen by `f`:
  - `t`   — time (`PetscReal`)
  - `x`   — coordinates (`Vector{PetscReal}`, length `dim`)
  - `u`   — output field values (`Vector{PetscScalar}`, length `Nc`)
  - `ctx` — user context (`Ptr{Cvoid}`)

Used by `dm_project_function!`, `add_boundary!`, `set_exact_solution!`.
`PetscInt`, `PetscReal`, `PetscScalar` must be in scope.
"""
macro petsc_simple_fn(f)
    cfn = Symbol(f, :__cfn__)
    gdim = gensym("dim"); gNc = gensym("Nc")
    PI = Core.eval(__module__, :PetscInt)
    PS = Core.eval(__module__, :PetscScalar)
    PR = Core.eval(__module__, :PetscReal)
    cfn_ptr_expr = :(Base.@cfunction($cfn, $PI,
        ($PI, $PR, Ptr{$PR}, $PI, Ptr{$PS}, Ptr{Cvoid})))
    quote
        function $(esc(cfn))(
            $gdim::$PI, t::$PR, xp::Ptr{$PR},
            $gNc::$PI, up::Ptr{$PS}, ctx::Ptr{Cvoid},
        )::$PI
            x = unsafe_wrap(Vector{$PR}, xp, $gdim)
            u = unsafe_wrap(Vector{$PS}, up, $gNc)
            $(esc(f))(t, x, u, ctx)
            return $PI(0)
        end
        @eval $cfn_ptr_expr
    end
end

"""
    @petsc_residual_fn(f, outsz)

Generate a C-callable wrapper for the pure-Julia residual function `f` and
return its C pointer, matching the `PetscPointFn` signature.  All pointer
arguments are unwrapped into Julia Vectors before calling `f`.

`outsz` is an expression giving the output array length in terms of `dim_`,
`Nf`, `NfAux`, or `numConstants` (e.g. `Nf` for `f0`-type functions,
`dim_*Nf` for `f1`-type functions).

User function signature:
```
f(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
  aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, out)
```
`PetscInt`, `PetscReal`, `PetscScalar` must be in scope.
"""
macro petsc_residual_fn(f, outsz)
    cfn = Symbol(f, :__cfn__)
    gdim = gensym("dim"); gNf = gensym("Nf"); gNfAux = gensym("NfAux"); gnC = gensym("nC")
    m = Dict{Symbol,Any}(:dim_ => gdim, :Nf => gNf, :NfAux => gNfAux, :numConstants => gnC)
    outsz_s = _petsc_subst(outsz, m)
    PI = Core.eval(__module__, :PetscInt)
    PS = Core.eval(__module__, :PetscScalar)
    PR = Core.eval(__module__, :PetscReal)
    cfn_ptr_expr = :(Base.@cfunction($cfn, Cvoid,
        ($PI, $PI, $PI, Ptr{$PI}, Ptr{$PI},
         Ptr{$PS}, Ptr{$PS}, Ptr{$PS}, Ptr{$PI}, Ptr{$PI},
         Ptr{$PS}, Ptr{$PS}, Ptr{$PS},
         $PR, Ptr{$PR}, $PI, Ptr{$PS}, Ptr{$PS})))
    quote
        function $(esc(cfn))(
            $gdim::$PI, $gNf::$PI, $gNfAux::$PI,
            uOff_p::Ptr{$PI}, uOff_xp::Ptr{$PI},
            u_p::Ptr{$PS}, u_tp::Ptr{$PS}, u_xp::Ptr{$PS},
            aOff_p::Ptr{$PI}, aOff_xp::Ptr{$PI},
            a_p::Ptr{$PS}, a_tp::Ptr{$PS}, a_xp::Ptr{$PS},
            t::$PR, x_p::Ptr{$PR}, $gnC::$PI, cst_p::Ptr{$PS},
            out_p::Ptr{$PS},
        )::Cvoid
            uOff   = unsafe_wrap(Vector{$PI}, uOff_p,  $gNf + 1)
            uOff_x = unsafe_wrap(Vector{$PI}, uOff_xp, $gNf + 1)
            u      = unsafe_wrap(Vector{$PS}, u_p,     uOff[end])
            u_t    = unsafe_wrap(Vector{$PS}, u_tp,    uOff[end])
            u_x    = unsafe_wrap(Vector{$PS}, u_xp,    uOff_x[end])
            aOff   = $gNfAux > 0 ? unsafe_wrap(Vector{$PI}, aOff_p,  $gNfAux + 1) : $PI[]
            aOff_x = $gNfAux > 0 ? unsafe_wrap(Vector{$PI}, aOff_xp, $gNfAux + 1) : $PI[]
            a      = $gNfAux > 0 ? unsafe_wrap(Vector{$PS}, a_p,  aOff[end])   : $PS[]
            a_t    = $gNfAux > 0 ? unsafe_wrap(Vector{$PS}, a_tp, aOff[end])   : $PS[]
            a_x    = $gNfAux > 0 ? unsafe_wrap(Vector{$PS}, a_xp, aOff_x[end]) : $PS[]
            x      = unsafe_wrap(Vector{$PR}, x_p,     $gdim)
            constants = $gnC > 0 ? unsafe_wrap(Vector{$PS}, cst_p, $gnC) : $PS[]
            out    = unsafe_wrap(Vector{$PS}, out_p,   $outsz_s)
            $(esc(f))($gdim, $gNf, $gNfAux, uOff, uOff_x, u, u_t, u_x,
                      aOff, aOff_x, a, a_t, a_x, t, x, $gnC, constants, out)
            return
        end
        @eval $cfn_ptr_expr
    end
end

"""
    @petsc_jacobian_fn(f, outsz)

Generate a C-callable wrapper for the pure-Julia Jacobian function `f` and
return its C pointer, matching the `PetscPointJacFn` signature.  Identical to
`@petsc_residual_fn` but with an extra `u_tShift` scalar between `t` and `x`.

User function signature:
```
f(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
  aOff, aOff_x, a, a_t, a_x, t, u_tShift, x, numConstants, constants, out)
```
`PetscInt`, `PetscReal`, `PetscScalar` must be in scope.
"""
macro petsc_jacobian_fn(f, outsz)
    cfn = Symbol(f, :__cfn__)
    gdim = gensym("dim"); gNf = gensym("Nf"); gNfAux = gensym("NfAux"); gnC = gensym("nC")
    m = Dict{Symbol,Any}(:dim_ => gdim, :Nf => gNf, :NfAux => gNfAux, :numConstants => gnC)
    outsz_s = _petsc_subst(outsz, m)
    PI = Core.eval(__module__, :PetscInt)
    PS = Core.eval(__module__, :PetscScalar)
    PR = Core.eval(__module__, :PetscReal)
    cfn_ptr_expr = :(Base.@cfunction($cfn, Cvoid,
        ($PI, $PI, $PI, Ptr{$PI}, Ptr{$PI},
         Ptr{$PS}, Ptr{$PS}, Ptr{$PS}, Ptr{$PI}, Ptr{$PI},
         Ptr{$PS}, Ptr{$PS}, Ptr{$PS},
         $PR, $PR, Ptr{$PR}, $PI, Ptr{$PS}, Ptr{$PS})))
    quote
        function $(esc(cfn))(
            $gdim::$PI, $gNf::$PI, $gNfAux::$PI,
            uOff_p::Ptr{$PI}, uOff_xp::Ptr{$PI},
            u_p::Ptr{$PS}, u_tp::Ptr{$PS}, u_xp::Ptr{$PS},
            aOff_p::Ptr{$PI}, aOff_xp::Ptr{$PI},
            a_p::Ptr{$PS}, a_tp::Ptr{$PS}, a_xp::Ptr{$PS},
            t::$PR, utShift::$PR, x_p::Ptr{$PR},
            $gnC::$PI, cst_p::Ptr{$PS},
            out_p::Ptr{$PS},
        )::Cvoid
            uOff   = unsafe_wrap(Vector{$PI}, uOff_p,  $gNf + 1)
            uOff_x = unsafe_wrap(Vector{$PI}, uOff_xp, $gNf + 1)
            u      = unsafe_wrap(Vector{$PS}, u_p,     uOff[end])
            u_t    = unsafe_wrap(Vector{$PS}, u_tp,    uOff[end])
            u_x    = unsafe_wrap(Vector{$PS}, u_xp,    uOff_x[end])
            aOff   = $gNfAux > 0 ? unsafe_wrap(Vector{$PI}, aOff_p,  $gNfAux + 1) : $PI[]
            aOff_x = $gNfAux > 0 ? unsafe_wrap(Vector{$PI}, aOff_xp, $gNfAux + 1) : $PI[]
            a      = $gNfAux > 0 ? unsafe_wrap(Vector{$PS}, a_p,  aOff[end])   : $PS[]
            a_t    = $gNfAux > 0 ? unsafe_wrap(Vector{$PS}, a_tp, aOff[end])   : $PS[]
            a_x    = $gNfAux > 0 ? unsafe_wrap(Vector{$PS}, a_xp, aOff_x[end]) : $PS[]
            x      = unsafe_wrap(Vector{$PR}, x_p,     $gdim)
            constants = $gnC > 0 ? unsafe_wrap(Vector{$PS}, cst_p, $gnC) : $PS[]
            out    = unsafe_wrap(Vector{$PS}, out_p,   $outsz_s)
            $(esc(f))($gdim, $gNf, $gNfAux, uOff, uOff_x, u, u_t, u_x,
                      aOff, aOff_x, a, a_t, a_x, t, utShift, x, $gnC, constants, out)
            return
        end
        @eval $cfn_ptr_expr
    end
end

"""
    @petsc_bd_fn(f, outsz)

Generate a C-callable wrapper for the pure-Julia boundary pointwise function `f`
and return its C pointer, matching the `PetscBdPointFn` signature.  Identical to
`@petsc_residual_fn` but with an extra outward-normal vector `n` between `x` and
`numConstants`.

User function signature:
```
f(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
  aOff, aOff_x, a, a_t, a_x, t, x, n, numConstants, constants, out)
```
`PetscInt`, `PetscReal`, `PetscScalar` must be in scope.
"""
macro petsc_bd_fn(f, outsz)
    cfn = Symbol(f, :__cfn__)
    gdim = gensym("dim"); gNf = gensym("Nf"); gNfAux = gensym("NfAux"); gnC = gensym("nC")
    m = Dict{Symbol,Any}(:dim_ => gdim, :Nf => gNf, :NfAux => gNfAux, :numConstants => gnC)
    outsz_s = _petsc_subst(outsz, m)
    PI = Core.eval(__module__, :PetscInt)
    PS = Core.eval(__module__, :PetscScalar)
    PR = Core.eval(__module__, :PetscReal)
    cfn_ptr_expr = :(Base.@cfunction($cfn, Cvoid,
        ($PI, $PI, $PI, Ptr{$PI}, Ptr{$PI},
         Ptr{$PS}, Ptr{$PS}, Ptr{$PS}, Ptr{$PI}, Ptr{$PI},
         Ptr{$PS}, Ptr{$PS}, Ptr{$PS},
         $PR, Ptr{$PR}, Ptr{$PR}, $PI, Ptr{$PS}, Ptr{$PS})))
    quote
        function $(esc(cfn))(
            $gdim::$PI, $gNf::$PI, $gNfAux::$PI,
            uOff_p::Ptr{$PI}, uOff_xp::Ptr{$PI},
            u_p::Ptr{$PS}, u_tp::Ptr{$PS}, u_xp::Ptr{$PS},
            aOff_p::Ptr{$PI}, aOff_xp::Ptr{$PI},
            a_p::Ptr{$PS}, a_tp::Ptr{$PS}, a_xp::Ptr{$PS},
            t::$PR, x_p::Ptr{$PR}, n_p::Ptr{$PR},
            $gnC::$PI, cst_p::Ptr{$PS},
            out_p::Ptr{$PS},
        )::Cvoid
            uOff   = unsafe_wrap(Vector{$PI}, uOff_p,  $gNf + 1)
            uOff_x = unsafe_wrap(Vector{$PI}, uOff_xp, $gNf + 1)
            u      = unsafe_wrap(Vector{$PS}, u_p,     uOff[end])
            u_t    = unsafe_wrap(Vector{$PS}, u_tp,    uOff[end])
            u_x    = unsafe_wrap(Vector{$PS}, u_xp,    uOff_x[end])
            aOff   = $gNfAux > 0 ? unsafe_wrap(Vector{$PI}, aOff_p,  $gNfAux + 1) : $PI[]
            aOff_x = $gNfAux > 0 ? unsafe_wrap(Vector{$PI}, aOff_xp, $gNfAux + 1) : $PI[]
            a      = $gNfAux > 0 ? unsafe_wrap(Vector{$PS}, a_p,  aOff[end])   : $PS[]
            a_t    = $gNfAux > 0 ? unsafe_wrap(Vector{$PS}, a_tp, aOff[end])   : $PS[]
            a_x    = $gNfAux > 0 ? unsafe_wrap(Vector{$PS}, a_xp, aOff_x[end]) : $PS[]
            x      = unsafe_wrap(Vector{$PR}, x_p, $gdim)
            n      = unsafe_wrap(Vector{$PR}, n_p, $gdim)
            constants = $gnC > 0 ? unsafe_wrap(Vector{$PS}, cst_p, $gnC) : $PS[]
            out    = unsafe_wrap(Vector{$PS}, out_p, $outsz_s)
            $(esc(f))($gdim, $gNf, $gNfAux, uOff, uOff_x, u, u_t, u_x,
                      aOff, aOff_x, a, a_t, a_x, t, x, n, $gnC, constants, out)
            return
        end
        @eval $cfn_ptr_expr
    end
end

# ── DMPlex constructors ──────────────────────────────────────────────────────

"""
    DMPlex(
        petsclib::PetscLib,
        comm::MPI.Comm;
        setfromoptions = true,
        dmsetup        = false,
        prefix         = "",
        options...,
    )

Create an empty `DMPLEX` mesh.  Use the keyword `options...` (or PETSc command-line
options) such as `dm_plex_box_faces`, `dm_plex_dim`, `dm_plex_simplex`,
`dm_plex_filename`, etc. to populate the mesh — this is the most flexible entry
point and matches the `DMCreate` + `DMSetType(DMPLEX)` + `DMSetFromOptions` flow
used in PETSc's tutorials (e.g. `snes/ex12.c`).

For an explicit box mesh, use the
[`DMPlex(petsclib, comm, dim, simplex, faces; ...)`](@ref) method below.

# External Links
$(_doc_external("DMPlex/DMPlexCreate"))
$(_doc_external("DM/DMSetType"))
$(_doc_external("DM/DMSetFromOptions"))
"""
function DMPlex(
    petsclib::PetscLib,
    comm::MPI.Comm;
    setfromoptions = true,
    dmsetup        = false,
    prefix         = "",
    options...,
) where {PetscLib}
    @assert initialized(getlib(PetscLib))

    dm = LibPETSc.DMCreate(petsclib, comm)
    LibPETSc.DMSetType(petsclib, dm, "plex")

    if !isempty(prefix)
        LibPETSc.DMSetOptionsPrefix(petsclib, dm, prefix)
    end

    if setfromoptions
        opts = Options(petsclib; options...)
        push!(opts)
        LibPETSc.DMSetFromOptions(petsclib, dm)
        pop!(opts)
    end

    if dmsetup
        setup!(dm)
    end

    if MPI.Comm_size(comm) == 1
        finalizer(destroy, dm)
    end
    return dm
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
        setfromoptions  = true,
        dmsetup         = false,
        prefix          = "",
        options...,
    )

Create a structured-grid box `DMPLEX` mesh in `dim` dimensions with `faces[d]`
cells along axis `d`.  When `simplex == true` each box cell is split into
simplices (triangles in 2D, tetrahedra in 3D); otherwise tensor-product cells
are used.

# External Links
$(_doc_external("DMPlex/DMPlexCreateBoxMesh"))
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
    setfromoptions   = true,
    dmsetup          = false,
    prefix           = "",
    options...,
) where {PetscLib}
    @assert initialized(getlib(PetscLib))
    @assert length(faces) == dim
    @assert length(lower) == dim
    @assert length(upper) == dim
    @assert length(periodicity) == dim

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

    if !isempty(prefix)
        LibPETSc.DMSetOptionsPrefix(petsclib, dm, prefix)
    end

    if setfromoptions
        opts = Options(petsclib; options...)
        push!(opts)
        LibPETSc.DMSetFromOptions(petsclib, dm)
        pop!(opts)
    end

    if dmsetup
        setup!(dm)
    end

    if MPI.Comm_size(comm) == 1
        finalizer(destroy, dm)
    end
    return dm
end


# ── Convenience helpers ──────────────────────────────────────────────────────

"""
    isplexsimplex(dm::AbstractPetscDM) -> Bool

Return `true` when the cells of `dm` (assumed to be a `DMPLEX`) are simplices.

# External Links
$(_doc_external("DMPlex/DMPlexIsSimplex"))
"""
function isplexsimplex(dm::AbstractPetscDM{PetscLib}) where {PetscLib}
    return Bool(LibPETSc.DMPlexIsSimplex(getlib(PetscLib), dm))
end

"""
    plexdistribute!(dm::AbstractPetscDM; overlap = 0) -> Union{Nothing, AbstractPetscDM}

Distribute the (serial) `DMPLEX` `dm` across the communicator with the given
point-overlap.  Returns the new distributed `DM` on the owning communicator, or
`nothing` if no redistribution occurred (e.g. serial run).

The original `dm` is *not* destroyed — the caller is responsible for that.

# External Links
$(_doc_external("DMPlex/DMPlexDistribute"))
"""
function plexdistribute!(
    dm::AbstractPetscDM{PetscLib};
    overlap::Integer = 0,
) where {PetscLib}
    PetscInt = inttype(PetscLib)
    # DMPlexDistribute returns (sf, dmParallel); here we only need dmParallel.
    _, dm_par = LibPETSc.DMPlexDistribute(getlib(PetscLib), dm, PetscInt(overlap))
    return dm_par
end


# ── PetscDS / PetscFE convenience ────────────────────────────────────────────

"""
    getds(dm::AbstractPetscDM) -> PetscDS

Return the `PetscDS` (discrete system) attached to `dm`.

# External Links
$(_doc_external("DM/DMGetDS"))
"""
function getds end

LibPETSc.@for_petsc function getds(dm::AbstractPetscDM{$PetscLib})
    return PetscDS{$PetscLib}(LibPETSc.DMGetDS(getlib($PetscLib), dm))
end

"""
    createds!(dm::AbstractPetscDM)

Build the `PetscDS` for `dm` from the currently-attached fields.

# External Links
$(_doc_external("DM/DMCreateDS"))
"""
function createds!(dm::AbstractPetscDM{PetscLib}) where {PetscLib}
    LibPETSc.DMCreateDS(getlib(PetscLib), dm)
    return nothing
end

"""
    setfield!(dm::AbstractPetscDM, field::Integer, fe; label = C_NULL)

Attach a `PetscFE` (or other discretization object) as the `field`-th field of
`dm` (0-based field index, matching PETSc).

# External Links
$(_doc_external("DM/DMSetField"))
"""
function setfield! end

LibPETSc.@for_petsc function setfield!(
    dm::AbstractPetscDM{$PetscLib},
    field::Integer,
    fe::Ptr;
    label::Ptr{Cvoid} = C_NULL,
)
    LibPETSc.DMSetField(getlib($PetscLib), dm, $PetscInt(field), label, Ptr{Cvoid}(fe))
    return nothing
end

"""
    fe_create_default(petsclib, comm, dim, Nc, simplex; degree = 1, prefix = "", qorder = -1)

Create a default `PetscFE` (Lagrange) for a single field with `Nc` components on
either a simplex (`simplex == true`) or tensor cell.

`degree` is the polynomial degree (1 = P1/Q1, 2 = P2/Q2, etc.).  It is injected
into the options database under `prefix` as `-\${prefix}petscspace_degree` before
the call to `PetscFECreateDefault`, so any options explicitly set by the user
will take precedence.

# External Links
$(_doc_external("FE/PetscFECreateDefault"))
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
    opts = Options(petsclib; Symbol(opt_key) => degree)
    push!(opts)
    fe = LibPETSc.PetscFECreateDefault(
        petsclib, comm,
        PetscInt(dim), PetscInt(Nc),
        LibPETSc.PetscBool(simplex),
        String(prefix),
        PetscInt(qorder),
    )
    pop!(opts)
    return fe
end

"""
    fe_create_lagrange(petsclib, comm, dim, Nc, simplex, degree; qorder = -1)

Create a `PetscFE` Lagrange element of polynomial degree `degree` for a single
field with `Nc` components on either a simplex (`simplex == true`) or tensor cell.
Unlike `fe_create_default`, the degree is specified explicitly and does not depend
on the options database.

# External Links
$(_doc_external("FE/PetscFECreateLagrange"))
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


# ── Label / Boundary helpers ─────────────────────────────────────────────────

"""
    getlabel(dm::AbstractPetscDM, name::AbstractString) -> Ptr{Cvoid}

Look up a `DMLabel` on `dm` by name and return its raw pointer (`Ptr{Cvoid}`).
Returns `C_NULL` if no label with that name exists.

# External Links
$(_doc_external("DM/DMGetLabel"))
"""
function getlabel end

LibPETSc.@for_petsc function getlabel(
    dm::AbstractPetscDM{$PetscLib},
    name::AbstractString,
)
    return LibPETSc.DMGetLabel(getlib($PetscLib), dm, String(name))
end

"""
    add_boundary!(petsclib, dm, bctype, name, label, values, field, comps,
                  bcfunc_ptr, bcfunc_t_ptr = C_NULL, ctx = C_NULL) -> bd::PetscInt

Attach a boundary condition to `dm`.  `label` is a `Ptr{Cvoid}` from
[`getlabel`](@ref).  `bcfunc_ptr` / `bcfunc_t_ptr` are `@cfunction(...)` results
matching PETSc's boundary callback signature.

# External Links
$(_doc_external("DM/DMAddBoundary"))
"""
function add_boundary! end

LibPETSc.@for_petsc function add_boundary!(
    petsclib::$UnionPetscLib,
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
    values_v = collect($PetscInt.(values))
    comps_v  = collect($PetscInt.(comps))
    return LibPETSc.DMAddBoundary(
        petsclib, dm,
        bctype, String(name), label,
        $PetscInt(length(values_v)), values_v,
        $PetscInt(field),
        $PetscInt(length(comps_v)), comps_v,
        bcfunc_ptr, bcfunc_t_ptr, ctx,
    )
end


# ── Wire DMPlex FEM residual/Jacobian assembly into SNES ─────────────────────

"""
    plex_set_snes_local_fem!(petsclib, dm; use_obj = false, ctx = C_NULL)

Tell SNES to use `DMPlex`'s built-in FEM residual / Jacobian assembly on `dm`.

# External Links
$(_doc_external("DMPlex/DMPlexSetSNESLocalFEM"))
"""
function plex_set_snes_local_fem! end

LibPETSc.@for_petsc function plex_set_snes_local_fem!(
    petsclib::$UnionPetscLib,
    dm::AbstractPetscDM{$PetscLib};
    use_obj::Bool = false,
    ctx::Ptr{Cvoid} = C_NULL,
)
    LibPETSc.DMPlexSetSNESLocalFEM(petsclib, dm, LibPETSc.PetscBool(use_obj), ctx)
    return nothing
end


# ── PetscDS point-function setters ──────────────────────────────────────────

"""
    set_residual!(ds, field, f0_ptr, f1_ptr)

Attach the pointwise residual point-function pair to field `field` (0-based) of
the `PetscDS` `ds`.  `f0_ptr` and `f1_ptr` must be `@cfunction(...)`-style
function pointers (`Ptr{Cvoid}`) matching PETSc's `PetscPointFn` signature.

# External Links
$(_doc_external("Dm/PetscDSSetResidual"))
"""
function set_residual!(
    ds::PetscDS{PetscLib},
    field::Integer,
    f0_ptr::Ptr{Cvoid},
    f1_ptr::Ptr{Cvoid},
) where {PetscLib}
    petsclib = getlib(PetscLib)
    LibPETSc.PetscDSSetResidual(petsclib, ds.ptr, PetscLib.PetscInt(field), f0_ptr, f1_ptr)
    return nothing
end

"""
    set_jacobian!(ds, fieldI, fieldJ, g0_ptr, g1_ptr, g2_ptr, g3_ptr)

Attach the pointwise Jacobian point-function quadruple for the block
`(fieldI, fieldJ)` of the `PetscDS` `ds`.  Pass `C_NULL` for any unused term.

# External Links
$(_doc_external("Dm/PetscDSSetJacobian"))
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
    return nothing
end

"""
    set_exact_solution!(ds, field, sol_ptr, ctx = C_NULL)

Register the pointwise exact-solution function `sol_ptr` for field `field`
(0-based) of the `PetscDS` `ds`.  The function pointer must match PETSc's
`PetscSimplePointFn` / `PetscErrorCode f(dim, time, x, Nc, u, ctx)` signature.

This is required for `DMComputeL2Diff` and `DMComputeExactSolution` to work.

# External Links
$(_doc_external("Dm/PetscDSSetExactSolution"))
"""
function set_exact_solution!(
    ds::PetscDS{PetscLib},
    field::Integer,
    sol_ptr::Ptr{Cvoid},
    ctx::Ptr{Cvoid} = C_NULL,
) where {PetscLib}
    petsclib = getlib(PetscLib)
    LibPETSc.PetscDSSetExactSolution(petsclib, ds.ptr, PetscLib.PetscInt(field), sol_ptr, ctx)
    return nothing
end


# ── DMProjectFunction / DMComputeL2Diff ─────────────────────────────────────

"""
    dm_project_function!(petsclib, dm, time, funcs, ctxs, mode, X)

Project a collection of pointwise functions into the global vector `X`.
`funcs` is a vector of `@cfunction` pointers (one per field) matching
PETSc's `PetscSimplePointFn` signature.  `ctxs` is a matching vector of
context pointers, or `nothing` to use `C_NULL` for every field.

`mode` is an `InsertMode` enum value, typically `INSERT_ALL_VALUES` (sets
constrained DOFs too) or `INSERT_VALUES` (free DOFs only).

# External Links
$(_doc_external("DM/DMProjectFunction"))
"""
function dm_project_function! end

LibPETSc.@for_petsc function dm_project_function!(
    petsclib::$UnionPetscLib,
    dm::AbstractPetscDM{$PetscLib},
    time::Real,
    funcs::AbstractVector{Ptr{Cvoid}},
    ctxs::AbstractVector{Ptr{Cvoid}},
    mode::LibPETSc.InsertMode,
    X::AbstractPetscVec{$PetscLib},
)
    n = length(funcs)
    @assert length(ctxs) == n "funcs and ctxs must have the same length"
    funcs_v = collect(funcs)
    ctxs_v  = collect(ctxs)
    GC.@preserve funcs_v ctxs_v LibPETSc.@chk ccall(
        (:DMProjectFunction, $petsc_library),
        LibPETSc.PetscErrorCode,
        (LibPETSc.CDM, $PetscReal, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}},
         LibPETSc.InsertMode, LibPETSc.CVec),
        dm, $PetscReal(time), funcs_v, ctxs_v, mode, X,
    )
    return nothing
end

function dm_project_function!(
    petsclib, dm::AbstractPetscDM,
    time::Real,
    funcs::AbstractVector,
    ctxs::Union{AbstractVector, Nothing},
    mode::LibPETSc.InsertMode,
    X,
)
    fptrs = Ptr{Cvoid}[Ptr{Cvoid}(f) for f in funcs]
    cptrs = ctxs === nothing ?
        fill(C_NULL, length(fptrs)) :
        Ptr{Cvoid}[Ptr{Cvoid}(c) for c in ctxs]
    dm_project_function!(petsclib, dm, time, fptrs, cptrs, mode, X)
end

"""
    dm_compute_l2diff(petsclib, dm, time, funcs, ctxs, X) -> PetscReal

Compute the L² error between the global vector `X` and the pointwise exact
functions `funcs` (one per field, matching `PetscSimplePointFn` signature).
`ctxs` is a matching vector of context pointers, or `nothing` to use `C_NULL`
for every field.  Requires that [`set_exact_solution!`](@ref) has been called.

# External Links
$(_doc_external("DM/DMComputeL2Diff"))
"""
function dm_compute_l2diff end

LibPETSc.@for_petsc function dm_compute_l2diff(
    petsclib::$UnionPetscLib,
    dm::AbstractPetscDM{$PetscLib},
    time::Real,
    funcs::AbstractVector{Ptr{Cvoid}},
    ctxs::AbstractVector{Ptr{Cvoid}},
    X::AbstractPetscVec{$PetscLib},
)
    n = length(funcs)
    @assert length(ctxs) == n "funcs and ctxs must have the same length"
    funcs_v  = collect(funcs)
    ctxs_v   = collect(ctxs)
    diff_ref = Ref{$PetscReal}(0)
    GC.@preserve funcs_v ctxs_v LibPETSc.@chk ccall(
        (:DMComputeL2Diff, $petsc_library),
        LibPETSc.PetscErrorCode,
        (LibPETSc.CDM, $PetscReal, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}},
         LibPETSc.CVec, Ptr{$PetscReal}),
        dm, $PetscReal(time), funcs_v, ctxs_v, X, diff_ref,
    )
    return diff_ref[]
end

function dm_compute_l2diff(
    petsclib, dm::AbstractPetscDM,
    time::Real,
    funcs::AbstractVector,
    ctxs::Union{AbstractVector, Nothing},
    X,
)
    fptrs = Ptr{Cvoid}[Ptr{Cvoid}(f) for f in funcs]
    cptrs = ctxs === nothing ?
        fill(C_NULL, length(fptrs)) :
        Ptr{Cvoid}[Ptr{Cvoid}(c) for c in ctxs]
    dm_compute_l2diff(petsclib, dm, time, fptrs, cptrs, X)
end

# ── Auxiliary-field helpers ───────────────────────────────────────────────────

"""
    dmclone(dm::AbstractPetscDM) -> AbstractPetscDM

Return a new DM that is a clone of `dm` (same topology, no fields or DS).
"""
function dmclone(dm::AbstractPetscDM{PetscLib}) where {PetscLib}
    newdm = LibPETSc.PetscDM(getlib(PetscLib))
    LibPETSc.DMClone(getlib(PetscLib), dm, newdm)
    return newdm
end

"""
    dm_create_global_vec(dm::AbstractPetscDM) -> AbstractPetscVec

Allocate a global vector matching the layout of `dm`.
"""
function dm_create_global_vec(dm::AbstractPetscDM{PetscLib}) where {PetscLib}
    return LibPETSc.DMCreateGlobalVector(getlib(PetscLib), dm)
end

"""
    dm_create_local_vec(dm::AbstractPetscDM) -> AbstractPetscVec

Allocate a local (ghosted) vector matching the layout of `dm`.
"""
function dm_create_local_vec(dm::AbstractPetscDM{PetscLib}) where {PetscLib}
    return LibPETSc.DMCreateLocalVector(getlib(PetscLib), dm)
end

"""
    dm_global_to_local!(dm, global_vec, local_vec; mode = INSERT_VALUES)

Scatter `global_vec` into `local_vec` (including ghost values).
"""
function dm_global_to_local!(
    dm::AbstractPetscDM{PetscLib},
    gvec,
    lvec;
    mode::LibPETSc.InsertMode = LibPETSc.INSERT_VALUES,
) where {PetscLib}
    LibPETSc.DMGlobalToLocal(getlib(PetscLib), dm, gvec, mode, lvec)
    return nothing
end

"""
    dm_set_auxiliary_vec!(dm, aux_local)

Attach a local auxiliary vector `aux_local` to `dm` (global label / value 0 / part 0).
The auxiliary field values are forwarded to all pointwise functions as the `a` argument.
"""
function dm_set_auxiliary_vec!(dm::AbstractPetscDM{PetscLib}, aux_local) where {PetscLib}
    petsclib = getlib(PetscLib)
    LibPETSc.DMSetAuxiliaryVec(petsclib, dm,
        Ptr{Cvoid}(C_NULL),
        PetscLib.PetscInt(0), PetscLib.PetscInt(0),
        aux_local)
    return nothing
end

"""
    mat_null_space_create(petsclib, comm; has_const = true) -> MatNullSpace

Create a null space containing the constant vector (Neumann problems).
"""
function mat_null_space_create(petsclib, comm; has_const::Bool = true)
    return LibPETSc.MatNullSpaceCreate(petsclib, comm,
        LibPETSc.PetscBool(has_const),
        petsclib.PetscInt(0),
        LibPETSc.PetscVec[])
end

"""
    mat_set_null_space!(mat, nullsp)

Attach `nullsp` to `mat` so the linear solver removes it each iteration.
"""
function mat_set_null_space!(mat, nullsp)
    LibPETSc.MatSetNullSpace(LibPETSc.getlib(typeof(mat).parameters[1]), mat, nullsp)
    return nothing
end
