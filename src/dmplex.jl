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
    petsclib = getlib(PetscLib)
    PetscInt = inttype(PetscLib)
    dm_par = LibPETSc.PetscDM(petsclib)
    LibPETSc.DMPlexDistribute(petsclib, dm, PetscInt(overlap),
                              Ptr{LibPETSc.PetscSF}(C_NULL), dm_par)
    return dm_par
end


# ── PetscDS / PetscFE convenience ────────────────────────────────────────────

"""
    petsc_setname!(petsclib, obj, name)

Set the name of any PETSc object (DM, Vec, FE, …) to `name`.
`obj` can be any pointer type that is convertible to `Ptr{Cvoid}`.

Thin wrapper around `PetscObjectSetName` that avoids needing an explicit
`convert(Ptr{Cvoid}, obj)` at the call site.

# External Links
$(_doc_external("Sys/PetscObjectSetName"))
"""
function petsc_setname!(petsclib::LibPETSc.PetscLibType, obj, name::AbstractString)
    LibPETSc.PetscObjectSetName(petsclib, convert(Ptr{Cvoid}, obj), String(name))
    return nothing
end

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
        petsclib, dm, bctype, String(name), label,
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
    set_jacobian_preconditioner!(ds, field_i, field_j, g0, g1, g2, g3)

Set the Jacobian preconditioner terms for the weak form between fields `field_i` and `field_j`.
Each `gN` argument is either a `Ptr{Cvoid}` function pointer or `C_NULL`.

This is the preconditioner equivalent of [`set_jacobian!`](@ref).

# External Links
$(_doc_external("Dm/PetscDSSetJacobianPreconditioner"))
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
    LibPETSc.DMProjectFunction(petsclib, dm, time, funcs_v, ctxs_v, mode, X)
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
    funcs_v = collect(funcs)
    ctxs_v  = collect(ctxs)
    return LibPETSc.DMComputeL2Diff(petsclib, dm, time, funcs_v, ctxs_v, X)
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
    dm_coarsen_hook_add!(dm, coarsenhook, restricthook = C_NULL)

Register a C callback invoked each time `dm` is coarsened (e.g. during FAS
hierarchy setup).  The coarsenhook signature is:
```
hook(fine::Ptr{Cvoid}, coarse::Ptr{Cvoid}, ctx::Ptr{Cvoid}) → Cint
```
Use `@cfunction(f, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))` to create
the pointer.  `restricthook` (optional) is called on each nonlinear solve
restriction step.
"""
function dm_coarsen_hook_add! end

# Convenience: infer petsclib from the DM type parameter
function dm_coarsen_hook_add!(
    dm::AbstractPetscDM{PL},
    coarsenhook::Ptr{Cvoid},
    restricthook::Ptr{Cvoid} = C_NULL,
) where {PL}
    LibPETSc.DMCoarsenHookAdd(getlib(PL), dm, coarsenhook, restricthook)
end

"""
    fe_copy_quadrature!(petsclib, src_fe, dst_fe)

Copy the quadrature rule from `src_fe` to `dst_fe` so both fields share
the same integration points.
"""
function fe_copy_quadrature!(petsclib, src_fe, dst_fe)
    LibPETSc.PetscFECopyQuadrature(petsclib, src_fe, dst_fe)
    return nothing
end

"""
    create_split_boundary_labels!(dm)

Split the single `"marker"` label on a 2D box mesh into four per-wall labels
(`"markerBottom"`, `"markerRight"`, `"markerTop"`, `"markerLeft"`, marker IDs 1–4)
so that corner DOFs can belong to multiple wall labels simultaneously.
"""
function create_split_boundary_labels! end

LibPETSc.@for_petsc function create_split_boundary_labels!(
    dm::AbstractPetscDM{$PetscLib},
)
    petsclib = getlib($PetscLib)
    names = ("markerBottom", "markerRight", "markerTop", "markerLeft")
    for (name, id) in zip(names, $PetscInt[1, 2, 3, 4])
        LibPETSc.DMCreateLabel(petsclib, dm, name)
        is = LibPETSc.IS{$PetscLib}()
        LibPETSc.DMGetStratumIS(petsclib, dm, "marker", id, is)
        is.ptr == C_NULL && continue
        label = LibPETSc.DMGetLabel(petsclib, dm, name)
        LibPETSc.DMLabelInsertIS(petsclib, label, is, $PetscInt(1))
        LibPETSc.ISDestroy(petsclib, is)
    end
    return nothing
end

"""
    mat_null_space_create(petsclib, comm; has_const = true) -> MatNullSpace

Create a null space containing the constant vector (Neumann / saddle-point problems).
"""
function mat_null_space_create(petsclib, comm; has_const::Bool = true)
    return LibPETSc.MatNullSpaceCreate(petsclib, comm,
        LibPETSc.PetscBool(has_const),
        petsclib.PetscInt(0),
        LibPETSc.CVec[])
end

"""
    mat_null_space_create(petsclib, comm, vecs) -> MatNullSpace

Create a null space spanned by the given vectors.  Each element of `vecs`
must have a `.ptr` field holding the underlying PETSc Vec handle (`CVec`).
The vectors should be orthonormal; call `VecNormalize` beforehand if needed.
"""
function mat_null_space_create end

LibPETSc.@for_petsc function mat_null_space_create(
    petsclib::$UnionPetscLib,
    comm,
    vecs,
)
    cvecs = LibPETSc.CVec[v.ptr for v in vecs]
    return LibPETSc.MatNullSpaceCreate(petsclib, comm,
        LibPETSc.PETSC_FALSE, $PetscInt(length(cvecs)), cvecs)
end

"""
    mat_null_space_destroy!(petsclib, nullsp)

Destroy a `MatNullSpace` created by `mat_null_space_create`.
"""
function mat_null_space_destroy!(petsclib, nullsp::LibPETSc.MatNullSpace)
    LibPETSc.MatNullSpaceDestroy(petsclib, nullsp)
    return nothing
end

"""
    mat_set_null_space!(mat, nullsp)

Attach `nullsp` to `mat` so the linear solver removes it each iteration.
"""
function mat_set_null_space!(mat, nullsp)
    LibPETSc.MatSetNullSpace(LibPETSc.getlib(typeof(mat).parameters[1]), mat, nullsp)
    return nothing
end

"""
    fe_compose_constant_null_space!(petsclib, comm, fe)

Create a trivial (constant) `MatNullSpace` and attach it to the FE object `fe`
under the key `"nullspace"` via `PetscObjectCompose`.  This signals to the
fieldsplit preconditioner that the field has a constant null space.
"""
function fe_compose_constant_null_space! end

LibPETSc.@for_petsc function fe_compose_constant_null_space!(
    petsclib::$UnionPetscLib,
    comm,
    fe,
)
    nsp = mat_null_space_create(petsclib, comm; has_const = true)
    LibPETSc.PetscObjectCompose(petsclib,
        convert(Ptr{Cvoid}, fe), "nullspace", convert(Ptr{Cvoid}, nsp))
    mat_null_space_destroy!(petsclib, nsp)
    return nothing
end

"""
    snes_set_jacobian_null_space!(snes, nullsp)

Retrieve the assembled Jacobian matrix from `snes` (after `SNESSetUp`) and
attach `nullsp` to it.  Must be called after `SNESSetUp` and before `SNESSolve`.
"""
function snes_set_jacobian_null_space! end

LibPETSc.@for_petsc function snes_set_jacobian_null_space!(
    snes::LibPETSc.PetscSNES{$PetscLib},
    nullsp::LibPETSc.MatNullSpace,
)
    petsclib = getlib($PetscLib)
    J = LibPETSc.SNESGetJacobianMat(petsclib, snes)
    LibPETSc.MatSetNullSpace(petsclib, J, nullsp)
    return nothing
end

"""
    vtk_save!(petsclib, comm, filename, vec)

Write the global vector `vec` to a VTK unstructured-grid file (`.vtu`),
readable by ParaView and VisIt.  Works for both serial and MPI runs; PETSc
gathers all ranks into a single file.  The filename must end in `.vtu`.
The mesh geometry is taken from the `DM` attached to `vec` by PETSc.
"""
function vtk_save! end

LibPETSc.@for_petsc function vtk_save!(
    petsclib::$UnionPetscLib,
    comm::MPI.Comm,
    filename::AbstractString,
    vec::AbstractPetscVec{$PetscLib},
)
    vtk = LibPETSc.PetscViewerVTKOpen(petsclib, comm, String(filename), LibPETSc.FILE_MODE_WRITE)
    LibPETSc.VecView(petsclib, vec, vtk)
    vtk_ref = Ref{LibPETSc.PetscViewer}(vtk)
    LibPETSc.PetscViewerDestroy(petsclib, vtk_ref)
    return nothing
end

"""
    dm_project_field!(petsclib, dm, time, U, funcs, mode, X)

Project a function of the fields in the input vector `U` into the FE space of `dm`,
writing the result into `X`.  `funcs` is a vector of `Ptr{Cvoid}` function pointers
(one per field in `dm`), each with the `PetscPointFn` / `@petsc_residual_fn` signature.
`U` must be associated with a DM that shares the same mesh as `dm` (e.g. obtained via
[`dmclone`](@ref)).

Use this to compute derived quantities (e.g. stress from displacement gradient) and
project them onto a new FE field for visualisation.

# External Links
$(_doc_external("DM/DMProjectField"))
"""
function dm_project_field! end

LibPETSc.@for_petsc function dm_project_field!(
    petsclib::$UnionPetscLib,
    dm::AbstractPetscDM{$PetscLib},
    time::Real,
    U::AbstractPetscVec{$PetscLib},
    funcs::AbstractVector{Ptr{Cvoid}},
    mode::LibPETSc.InsertMode,
    X::AbstractPetscVec{$PetscLib},
)
    funcs_v = collect(Ptr{Cvoid}, funcs)
    GC.@preserve funcs_v LibPETSc.DMProjectField(petsclib, dm, $PetscReal(time), U, pointer(funcs_v), mode, X)
    return nothing
end

function dm_project_field!(
    petsclib, dm::AbstractPetscDM,
    time::Real, U, funcs, mode, X,
)
    fptrs = Ptr{Cvoid}[Ptr{Cvoid}(f) for f in funcs]
    dm_project_field!(petsclib, dm, time, U, fptrs, mode, X)
end

"""
    vtk_save_fields!(petsclib, comm, filename, vecs)

Write multiple global vectors to a single VTK file (`.vtu`).  `vecs` is any
iterable of `AbstractPetscVec`; each is written as a separate point-data array.
The field name shown in ParaView/VisIt comes from the name set on the vector via
`PetscObjectSetName`.

This is the multi-field equivalent of [`vtk_save!`](@ref).
"""
function vtk_save_fields! end

LibPETSc.@for_petsc function vtk_save_fields!(
    petsclib::$UnionPetscLib,
    comm::MPI.Comm,
    filename::AbstractString,
    vecs,   # iterable of AbstractPetscVec
)
    vtk = LibPETSc.PetscViewerVTKOpen(petsclib, comm, String(filename), LibPETSc.FILE_MODE_WRITE)
    for vec in vecs
        LibPETSc.VecView(petsclib, vec, vtk)
    end
    vtk_ref = Ref{LibPETSc.PetscViewer}(vtk)
    LibPETSc.PetscViewerDestroy(petsclib, vtk_ref)
    return nothing
end

"""
    set_constants!(ds::PetscDS, constants)

Set the named constants (accessible as `constants[i]` in pointwise functions) on the
`PetscDS` `ds`.  `constants` must be an `AbstractVector` whose elements are convertible
to `PetscScalar`.

# External Links
$(_doc_external("Dm/PetscDSSetConstants"))
"""
function set_constants!(ds::PetscDS{PetscLib}, constants::AbstractVector) where {PetscLib}
    petsclib = getlib(PetscLib)
    cst_v = collect(PetscLib.PetscScalar, constants)
    LibPETSc.PetscDSSetConstants(petsclib, ds.ptr, PetscLib.PetscInt(length(cst_v)), cst_v)
    return nothing
end

"""
    add_natural_boundary!(petsclib, dm, ds, name, label, label_value, field, f0_ptr, f1_ptr = C_NULL) -> bd::PetscInt

Register a natural (Neumann) boundary condition on `dm` following the PETSc 3.22+
pattern used in the C tutorials.  Equivalent to:
```c
DMAddBoundary(dm, DM_BC_NATURAL, name, label, 1, &val, field, 0, NULL, NULL, NULL, NULL, &bd);
PetscDSGetBoundary(ds, bd, &wf, NULL, ...);
PetscWeakFormSetIndexBdResidual(wf, label, val, field, 0, 0, f0, 0, NULL);
```

`f0_ptr` is a `@petsc_bd_fn`-generated C-callable function pointer implementing
the boundary integrand (`PetscBdPointFn` signature, with the outward normal `n[]`
between `x[]` and `numConstants`).

# External Links
$(_doc_external("DM/DMAddBoundary"))
$(_doc_external("Dm/PetscWeakFormSetIndexBdResidual"))
"""
function add_natural_boundary! end

LibPETSc.@for_petsc function add_natural_boundary!(
    petsclib::$UnionPetscLib,
    dm::AbstractPetscDM{$PetscLib},
    ds::PetscDS{$PetscLib},
    name::AbstractString,
    label::Ptr{Cvoid},
    label_value::Integer,
    field::Integer,
    f0_ptr::Ptr{Cvoid},
    f1_ptr::Ptr{Cvoid} = C_NULL,
)
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
    return bd
end

"""
    dm_copy_disc!(src::AbstractPetscDM, dst::AbstractPetscDM)

Copy the discretisation (fields, `PetscDS`, BCs) from `src` to `dst`.  Useful when
propagating FEM setup to coarser levels of a multigrid hierarchy.

# External Links
$(_doc_external("DM/DMCopyDisc"))
"""
function dm_copy_disc! end

LibPETSc.@for_petsc function dm_copy_disc!(
    src::AbstractPetscDM{$PetscLib},
    dst::AbstractPetscDM{$PetscLib},
)
    LibPETSc.DMCopyDisc(getlib($PetscLib), src, dst)
    return nothing
end

"""
    dm_get_coarse(dm::AbstractPetscDM) -> AbstractPetscDM

Return the coarse `DM` from which `dm` was obtained by refinement (e.g. via
`-dm_refine_hierarchy`).  The returned DM is a borrowed reference owned by PETSc;
do **not** call `destroy` on it.  Check `convert(Ptr{Cvoid}, cdm) == C_NULL` to
detect when there is no coarser level.

# External Links
$(_doc_external("DM/DMGetCoarseDM"))
"""
function dm_get_coarse end

LibPETSc.@for_petsc function dm_get_coarse(dm::AbstractPetscDM{$PetscLib})
    petsclib = getlib($PetscLib)
    cdm = LibPETSc.PetscDM(petsclib)
    LibPETSc.DMGetCoarseDM(petsclib, dm, cdm)
    return cdm
end

"""
    vtk_merge_tensor!(fname, names...)

Post-process a VTK `.vtu` file written by PETSc's VTK viewer to merge 9 separate
scalar `DataArray`s named `name.0`…`name.8` (produced when a field has more than 3
components) into a single `NumberOfComponents="9"` DataArray, for each name given.
Also sets `Tensors="name1 name2 …"` on the `<PointData>` tag so ParaView treats the
arrays as 3×3 tensors.  The file is rewritten in-place.

PETSc's VTK writer always splits fields with more than 3 components into separate
scalar arrays; this function reassembles them for proper tensor visualisation.
"""
function vtk_merge_tensor!(fname::AbstractString, names::AbstractString...)
    for name in names
        _vtk_merge_one_tensor!(fname, name)
    end
    return nothing
end

function _vtk_merge_one_tensor!(fname::AbstractString, name::AbstractString)
    isfile(fname) || return
    raw = read(fname)

    needle = b"<AppendedData"
    app_i  = let found = nothing
        for i in 1:length(raw)-length(needle)
            if @views raw[i:i+length(needle)-1] == needle
                found = i; break
            end
        end
        found
    end
    app_i === nothing && return
    gt  = findnext(isequal(UInt8('>')), raw, app_i + length(needle))
    gt  === nothing && return
    sep = findnext(isequal(UInt8('_')), raw, gt + 1)
    sep === nothing && return

    xml = String(raw[1:sep-1])
    bin = raw[sep+1:end]

    hdr_size = contains(xml, "header_type=\"UInt64\"") ? 8 : 4
    read_hdr(off) = hdr_size == 8 ?
        Int(read(IOBuffer(bin[off+1:off+8]), UInt64)) :
        Int(read(IOBuffer(bin[off+1:off+4]), UInt32))

    pat = Regex("""<DataArray\\b[^>]*\\bName=\"$(name)\\.(\\d+)\"[^>]*/?>""")
    ms  = collect(eachmatch(pat, xml))
    length(ms) == 9 || return

    entries = [(comp    = parse(Int, m[1]),
                xml_beg = m.offset,
                xml_end = m.offset + length(m.match) - 1,
                bin_off = parse(Int, match(r"offset=\"(\d+)\"", m.match)[1]))
               for m in ms]

    xml_first = minimum(e.xml_beg for e in entries)
    xml_last  = maximum(e.xml_end for e in entries)

    sorted   = sort(entries; by = e -> e.comp)
    offsets  = [e.bin_off for e in sorted]

    val_size = contains(ms[1].match, "Float64") ? 8 : 4
    RT       = val_size == 8 ? Float64 : Float32

    function read_scalar(off)
        nbytes = read_hdr(off)
        nvals  = nbytes ÷ val_size
        io     = IOBuffer(bin[off+hdr_size+1:off+hdr_size+nbytes])
        [read(io, RT) for _ in 1:nvals]
    end

    scalars = [read_scalar(off) for off in offsets]
    npoints = length(scalars[1])

    merged = Vector{RT}(undef, npoints * 9)
    for i in 0:npoints-1, k in 0:8
        merged[i*9+k+1] = scalars[k+1][i+1]
    end

    merged_nbytes = UInt64(npoints * 9 * val_size)
    merged_block  = let buf = IOBuffer()
        hdr_size == 8 ? write(buf, merged_nbytes) : write(buf, UInt32(merged_nbytes))
        for v in merged; write(buf, v); end
        take!(buf)
    end

    last_off    = offsets[9]
    old_bin_end = last_off + hdr_size + npoints * val_size
    Δ           = length(merged_block) - (old_bin_end - offsets[1])

    new_bin  = vcat(bin[1:offsets[1]], merged_block, bin[old_bin_end+1:end])

    type_str = RT == Float64 ? "Float64" : "Float32"
    new_da   = "<DataArray Name=\"$name\" NumberOfComponents=\"9\"" *
               " type=\"$type_str\" format=\"appended\" offset=\"$(offsets[1])\"/>"

    post = xml[xml_last+1:end]
    post_fixed = replace(post, r"offset=\"(\d+)\"" => function(s)
        m2  = match(r"\"(\d+)\"", s)
        off = parse(Int, m2[1])
        off > last_off ? "offset=\"$(off + Δ)\"" : s
    end)

    pre = xml[1:xml_first-1]
    # Julia's replace passes SubString (not RegexMatch) to the function, so we
    # must call match() explicitly to extract capture groups.
    pre = replace(pre, r"<PointData\b([^>]*)>" =>
        function(s)
            if contains(s, "Tensors=")
                m2 = match(r"Tensors=\"([^\"]*)\"", String(s))
                replace(String(s), m2.match => "Tensors=\"$(m2[1]) $name\"")
            else
                replace(String(s), "<PointData" => "<PointData Tensors=\"$name\"")
            end
        end)

    open(fname, "w") do io
        write(io, pre)
        write(io, new_da)
        write(io, post_fixed)
        write(io, UInt8('_'))
        write(io, new_bin)
    end
    return nothing
end
