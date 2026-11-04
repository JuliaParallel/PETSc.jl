# ── @cfunction helpers for PETSc FEM pointwise callbacks ─────────────────────
#
# These macros generate a thin C-callable wrapper that hands every C pointer
# argument to the user's pure-Julia function as a `PtrVector`, then calls it.
#
# The first argument names the library (a PetscLibType value or type). It is
# evaluated when the macro expands, so it must be a global defined before the
# macro call; the C signature uses that library's PetscInt, PetscReal and
# PetscScalar.
#
# Usage:
#   const petsclib = PETSc.getlib(; PetscScalar = Float64)
#   function my_fn(t, x, u, ctx);  u[1] = x[1]^2 + x[2]^2;  end
#   const my_fn_ptr = PETSc.@simple_fn(petsclib, my_fn)
#
#   function f0(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
#               aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, out)
#       out .= 4
#   end
#   const f0_ptr = PETSc.@residual_fn(petsclib, f0, Nf)        # outsz = Nf
#   const f1_ptr = PETSc.@residual_fn(petsclib, f1, dim_*Nf)   # outsz = dim*Nf
#   const g3_ptr = PETSc.@jacobian_fn(petsclib, g3, dim_*dim_) # outsz = dim² (Nf=1)
#
# The wrapper's C pointer comes from an `@eval` at run time: within a single
# top-level block the `@cfunction` is compiled before the wrapper exists.

"""
    PtrVector{T}(ptr, n)

An `AbstractVector{T}` over `n` values PETSc owns at `ptr`, as the pointwise
macros pass each array argument. It is valid only during the call. Indexing is
bounds-checked unless the caller uses `@inbounds`, and `pointer` and `ccall` see
`ptr`. It does not allocate, which matters at one call per quadrature point.
"""
struct PtrVector{T} <: AbstractVector{T}
    ptr::Ptr{T}
    n::Int
    PtrVector{T}(ptr, n) where {T} = new{T}(ptr, n)
end
PtrVector(ptr::Ptr{T}, n) where {T} = PtrVector{T}(ptr, n)

Base.size(v::PtrVector) = (v.n,)
Base.IndexStyle(::Type{<:PtrVector}) = IndexLinear()
Base.@propagate_inbounds function Base.getindex(v::PtrVector, i::Int)
    @boundscheck checkbounds(v, i)
    return unsafe_load(v.ptr, i)
end
Base.@propagate_inbounds function Base.setindex!(v::PtrVector, x, i::Int)
    @boundscheck checkbounds(v, i)
    unsafe_store!(v.ptr, x, i)
    return v
end
Base.pointer(v::PtrVector) = v.ptr
Base.unsafe_convert(::Type{Ptr{T}}, v::PtrVector{T}) where {T} = v.ptr

# Look up the PetscInt, PetscScalar and PetscReal of the library `lib` names
# in module `mod`.
function pointwise_types(mod::Module, lib)
    L = Core.eval(mod, lib)
    L isa LibPETSc.PetscLibType || (L isa Type && L <: LibPETSc.PetscLibType && isconcretetype(L)) ||
        throw(ArgumentError("the first argument of a pointwise macro must name a PETSc library or its type, got `$lib` = $L"))
    return LibPETSc.inttype(L), LibPETSc.scalartype(L), LibPETSc.realtype(L)
end

# Helper: substitute gensymmed param names into an output-size expression.
petsc_subst(ex, m) =
    ex isa Symbol ? get(m, ex, ex) :
    ex isa Expr   ? Expr(ex.head, map(a -> petsc_subst(a, m), ex.args)...) :
    ex

"""
    @simple_fn(petsclib, f)

Generate a C-callable wrapper for `f(t, x, u, ctx)` and return its C pointer,
matching the `PetscSimplePointFn` signature.  Arguments seen by `f`:
  - `t`   — time (`PetscReal`)
  - `x`   — coordinates (a `PtrVector{PetscReal}` of length `dim`)
  - `u`   — output field values (a `PtrVector{PetscScalar}` of length `Nc`)
  - `ctx` — user context (`Ptr{Cvoid}`)

Used by `project_function!`, `add_boundary!`, `set_exact_solution!`.
`petsclib` is a PETSc library or its type, and must be a global defined before
the macro call: it sets the C types of the arguments.
"""
macro simple_fn(lib, f)
    return simple_fn_expr(f, pointwise_types(__module__, lib)...)
end

function simple_fn_expr(f, PI, PS, PR)
    run_callback = GlobalRef(@__MODULE__, :run_callback)
    PV = GlobalRef(@__MODULE__, :PtrVector)
    cfn = Symbol(f, :__cfn__)
    gdim = gensym("dim"); gNc = gensym("Nc")
    cfn_ptr_expr = :(Base.@cfunction($cfn, Cint,
        ($PI, $PR, Ptr{$PR}, $PI, Ptr{$PS}, Ptr{Cvoid})))
    quote
        function $(esc(cfn))(
            $gdim::$PI, t::$PR, xp::Ptr{$PR},
            $gNc::$PI, up::Ptr{$PS}, ctx::Ptr{Cvoid},
        )::Cint
            x = $PV(xp, $gdim)
            u = $PV(up, $gNc)
            return $run_callback(() -> $(esc(f))(t, x, u, ctx), "pointwise function")
        end
        @eval $cfn_ptr_expr
    end
end

"""
    @residual_fn(petsclib, f, outsz)

Generate a C-callable wrapper for the pure-Julia residual function `f` and
return its C pointer, matching the `PetscPointFn` signature.  Each pointer
argument reaches `f` as a `PETSc.PtrVector`, an `AbstractVector` over
PETSc's memory, so annotate the arguments as
`AbstractVector` or not at all.

`outsz` is an expression giving the output array length in terms of `dim_`,
`Nf`, `NfAux`, or `numConstants` (e.g. `Nf` for `f0`-type functions,
`dim_*Nf` for `f1`-type functions).

User function signature:
```
f(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
  aOff, aOff_x, a, a_t, a_x, t, x, numConstants, constants, out)
```
`petsclib` is as in `@simple_fn`.
"""
macro residual_fn(lib, f, outsz)
    return residual_fn_expr(f, outsz, pointwise_types(__module__, lib)...)
end

function residual_fn_expr(f, outsz, PI, PS, PR)
    run_callback = GlobalRef(@__MODULE__, :run_callback)
    PV = GlobalRef(@__MODULE__, :PtrVector)
    cfn = Symbol(f, :__cfn__)
    gdim = gensym("dim"); gNf = gensym("Nf"); gNfAux = gensym("NfAux"); gnC = gensym("nC")
    m = Dict{Symbol,Any}(:dim_ => gdim, :Nf => gNf, :NfAux => gNfAux, :numConstants => gnC)
    outsz_s = petsc_subst(outsz, m)
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
            uOff   = $PV(uOff_p,  $gNf + 1)
            uOff_x = $PV(uOff_xp, $gNf + 1)
            u      = $PV(u_p,     uOff[end])
            u_t    = $PV(u_tp,    uOff[end])
            u_x    = $PV(u_xp,    uOff_x[end])
            # with no auxiliary fields PETSc may pass null pointers: length 0
            naux   = $gNfAux > 0 ? $gNfAux + 1 : 0
            aOff   = $PV(aOff_p,  naux)
            aOff_x = $PV(aOff_xp, naux)
            a      = $PV(a_p,  naux > 0 ? aOff[end]   : 0)
            a_t    = $PV(a_tp, naux > 0 ? aOff[end]   : 0)
            a_x    = $PV(a_xp, naux > 0 ? aOff_x[end] : 0)
            x      = $PV(x_p,     $gdim)
            constants = $PV(cst_p, $gnC)
            out    = $PV(out_p,   $outsz_s)
            $run_callback(() -> $(esc(f))($gdim, $gNf, $gNfAux, uOff, uOff_x, u, u_t, u_x,
                      aOff, aOff_x, a, a_t, a_x, t, x, $gnC, constants, out),
                "pointwise function")
            return
        end
        @eval $cfn_ptr_expr
    end
end

"""
    @jacobian_fn(petsclib, f, outsz)

Generate a C-callable wrapper for the pure-Julia Jacobian function `f` and
return its C pointer, matching the `PetscPointJacFn` signature.  Identical to
`@residual_fn` but with an extra `u_tShift` scalar between `t` and `x`.

User function signature:
```
f(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
  aOff, aOff_x, a, a_t, a_x, t, u_tShift, x, numConstants, constants, out)
```
`petsclib` is as in `@simple_fn`.
"""
macro jacobian_fn(lib, f, outsz)
    return jacobian_fn_expr(f, outsz, pointwise_types(__module__, lib)...)
end

function jacobian_fn_expr(f, outsz, PI, PS, PR)
    run_callback = GlobalRef(@__MODULE__, :run_callback)
    PV = GlobalRef(@__MODULE__, :PtrVector)
    cfn = Symbol(f, :__cfn__)
    gdim = gensym("dim"); gNf = gensym("Nf"); gNfAux = gensym("NfAux"); gnC = gensym("nC")
    m = Dict{Symbol,Any}(:dim_ => gdim, :Nf => gNf, :NfAux => gNfAux, :numConstants => gnC)
    outsz_s = petsc_subst(outsz, m)
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
            uOff   = $PV(uOff_p,  $gNf + 1)
            uOff_x = $PV(uOff_xp, $gNf + 1)
            u      = $PV(u_p,     uOff[end])
            u_t    = $PV(u_tp,    uOff[end])
            u_x    = $PV(u_xp,    uOff_x[end])
            # with no auxiliary fields PETSc may pass null pointers: length 0
            naux   = $gNfAux > 0 ? $gNfAux + 1 : 0
            aOff   = $PV(aOff_p,  naux)
            aOff_x = $PV(aOff_xp, naux)
            a      = $PV(a_p,  naux > 0 ? aOff[end]   : 0)
            a_t    = $PV(a_tp, naux > 0 ? aOff[end]   : 0)
            a_x    = $PV(a_xp, naux > 0 ? aOff_x[end] : 0)
            x      = $PV(x_p,     $gdim)
            constants = $PV(cst_p, $gnC)
            out    = $PV(out_p,   $outsz_s)
            $run_callback(() -> $(esc(f))($gdim, $gNf, $gNfAux, uOff, uOff_x, u, u_t, u_x,
                      aOff, aOff_x, a, a_t, a_x, t, utShift, x, $gnC, constants, out),
                "pointwise function")
            return
        end
        @eval $cfn_ptr_expr
    end
end

"""
    @bd_fn(petsclib, f, outsz)

Generate a C-callable wrapper for the pure-Julia boundary pointwise function `f`
and return its C pointer, matching the `PetscBdPointFn` signature.  Identical to
`@residual_fn` but with an extra outward-normal vector `n` between `x` and
`numConstants`.

User function signature:
```
f(dim_, Nf, NfAux, uOff, uOff_x, u, u_t, u_x,
  aOff, aOff_x, a, a_t, a_x, t, x, n, numConstants, constants, out)
```
`petsclib` is as in `@simple_fn`.
"""
macro bd_fn(lib, f, outsz)
    return bd_fn_expr(f, outsz, pointwise_types(__module__, lib)...)
end

function bd_fn_expr(f, outsz, PI, PS, PR)
    run_callback = GlobalRef(@__MODULE__, :run_callback)
    PV = GlobalRef(@__MODULE__, :PtrVector)
    cfn = Symbol(f, :__cfn__)
    gdim = gensym("dim"); gNf = gensym("Nf"); gNfAux = gensym("NfAux"); gnC = gensym("nC")
    m = Dict{Symbol,Any}(:dim_ => gdim, :Nf => gNf, :NfAux => gNfAux, :numConstants => gnC)
    outsz_s = petsc_subst(outsz, m)
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
            uOff   = $PV(uOff_p,  $gNf + 1)
            uOff_x = $PV(uOff_xp, $gNf + 1)
            u      = $PV(u_p,     uOff[end])
            u_t    = $PV(u_tp,    uOff[end])
            u_x    = $PV(u_xp,    uOff_x[end])
            # with no auxiliary fields PETSc may pass null pointers: length 0
            naux   = $gNfAux > 0 ? $gNfAux + 1 : 0
            aOff   = $PV(aOff_p,  naux)
            aOff_x = $PV(aOff_xp, naux)
            a      = $PV(a_p,  naux > 0 ? aOff[end]   : 0)
            a_t    = $PV(a_tp, naux > 0 ? aOff[end]   : 0)
            a_x    = $PV(a_xp, naux > 0 ? aOff_x[end] : 0)
            x      = $PV(x_p, $gdim)
            n      = $PV(n_p, $gdim)
            constants = $PV(cst_p, $gnC)
            out    = $PV(out_p, $outsz_s)
            $run_callback(() -> $(esc(f))($gdim, $gNf, $gNfAux, uOff, uOff_x, u, u_t, u_x,
                      aOff, aOff_x, a, a_t, a_x, t, x, n, $gnC, constants, out),
                "pointwise function")
            return
        end
        @eval $cfn_ptr_expr
    end
end
