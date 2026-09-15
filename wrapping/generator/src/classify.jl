# Argument classification and the per-argument code fragments

struct FArg
    name::String          # Julia-side name
    name_ccall::String    # expression passed to ccall
    typename::String      # Julia type in the signature / docs (concrete handle names)
    ccall_str::String     # type in the ccall tuple
    output::Bool
    init::String
    extract::String
    isarray::Bool
    stars::Int
    isfunction::Bool
end

const SIMPLE_TYPES = ("PetscScalar", "PetscBool", "PetscReal", "PetscComplex", "PetscInt")

"""
Is the argument an output? The manual page's `Input Parameters` / `Output Parameters` sections
decide when they mention the argument; otherwise a pointer is guessed to be an output, with the
original generator's heuristics restricting that guess.
"""
function is_output(r::Rules, fn::String, name::String, typename::String, stars::Int, output_vars, input_vars)
    name in input_vars && return false
    name in output_vars && return stars >= 1
    output = stars == 1
    if output
        simple = typename in SIMPLE_TYPES || typename in r.enum_types || typename in r.string_types ||
                 typename == "MPI_Comm"
        if (!occursin("Create", fn) && !occursin("Duplicate", fn) && !occursin("Type", fn) && !simple) ||
           occursin("Restore", fn) || occursin("Copy", fn)
            output = false
        end
    end
    return output
end

is_simple(r::Rules, t) = t in SIMPLE_TYPES || t in r.enum_types || t == "MPI_Comm" || t in ("Cint", "Csize_t", "Cdouble", "Cfloat", "Bool", "Cchar", "Int32", "PetscInt64", "PetscInt32", "PetscMPIInt", "PetscCount", "PetscLogDouble", "PetscObjectId", "PetscObjectState", "PetscClassId")

function init_extract(r::Rules, typename::String, name::String, isarray::Bool, isoutput::Bool, stars::Int)
    init, extract, name_ccall = "", "", name
    if !isarray && isoutput && is_handle(r, typename)
        h = r.handles[typename]
        name_ccall = "$(name)_"
        init = "$name_ccall = Ref{$(h.c)}()"
        extract = "$name = $(h.julia)($(name_ccall)[], petsclib)"
    elseif !isarray && !isoutput && stars == 1 && is_handle(r, typename)
        name_ccall = "$(name)_"
        init = "$name_ccall = Ref($(name).ptr)"
        extract = "$(name).ptr = C_NULL"
    elseif isarray && isoutput && stars > 0
        name_ccall = "$(name)_"
        init = "$name_ccall = Ref{" * "Ptr{"^stars * typename * "}"^stars * "}()"
        extract = "$name = unsafe_wrap(Array, $name_ccall[], VecGetLocalSize(petsclib, x); own = false)"
    elseif isarray && isoutput && stars == 0
        init = "$name = Vector{$typename}(undef, ni);  # CHECK SIZE!!"
    elseif isarray && !isoutput && stars == 1
        name_ccall = "$(name)_"
        init = "$name_ccall = Ref(pointer($name))"
    elseif !isarray && isoutput && typename in r.string_types
        name_ccall = "$(name)_"
        init = "$name_ccall = Ref{$typename}()"
        extract = "$name = unsafe_string($(name_ccall)[])"
    elseif !isarray && isoutput
        name_ccall = "$(name)_"
        init = "$name_ccall = Ref{$typename}()"
        extract = "$name = $(name_ccall)[]"
    end
    return init, extract, name_ccall
end

is_destroy(fn::Fn) = endswith(fn.name, "Destroy")

function classify(r::Rules, fn::Fn, a::Arg, input_vars, output_vars)
    stars = a.stars
    typename = map_type(r, a.typename)
    name = rename_arg(r, a.name)
    isarray = a.array
    isoutput = is_output(r, fn.name, name, typename, stars, output_vars, input_vars)
    ov = get(get(r.args, fn.name, Dict{String,Dict{String,Any}}()), name, Dict{String,Any}())
    if haskey(ov, "direction")
        isoutput = ov["direction"] == "out"
    end
    # --- function pointers and void pointers -----------------------------------------
    # A function pointer is passed as a Ptr{Cvoid} (from @cfunction); `Fn **out` returns one.
    isfnptr = a.isfunction || typename == "Ptr{Cvoid}" || (endswith(typename, "Fn") && !isarray)
    if typename == "Ptr{Cvoid}"                      # inline `void (*f)(...)`: name may sit in the C type
        isempty(name) && (name = rename_arg(r, String(match(r"\(\**(\w+)", a.typename).captures[1])))
        return FArg(name, name, typename, typename, false, "", "", false, stars, true)
    end
    if isfnptr || typename == "Cvoid"
        if stars >= 2 && !isarray
            return FArg(name, "$(name)_", "Ptr{Cvoid}", "Ptr{Ptr{Cvoid}}", true,
                        "$(name)_ = Ref{Ptr{Cvoid}}()", "$name = $(name)_[]", false, stars, isfnptr)
        elseif !isarray
            return FArg(name, name, "Ptr{Cvoid}", "Ptr{Cvoid}", false, "", "", false, stars, isfnptr)
        end
    end
    # caller-allocated output arrays are only allocated when their length is known
    if isarray && isoutput && stars == 0 && !haskey(ov, "len") && !any(x.name == "ni" for x in fn.args)
        isoutput = false
    end
    # `char **name` output: a C string PETSc owns
    if typename == "Cchar" && stars == 2 && !isarray && isoutput
        return FArg(name, "$(name)_", "String", "Ptr{Ptr{Cchar}}", true,
                    "$(name)_ = Ref{Ptr{Cchar}}()", "$name = unsafe_string($(name)_[])", false, stars, false)
    end
    # `T **out` output that is not marked as an array: hand back the raw pointer
    if stars == 2 && !isarray && isoutput && !is_handle(r, typename)
        return FArg(name, "$(name)_", "Ptr{$typename}", "Ptr{Ptr{$typename}}", true,
                    "$(name)_ = Ref{Ptr{$typename}}()", "$name = $(name)_[]", false, stars, false)
    end
    # scalar passed by reference as an input (`PetscInt *n` documented as input, Restore functions)
    if stars == 1 && !isarray && !isoutput && is_simple(r, typename) && (name in input_vars || occursin("Restore", fn.name))
        return FArg(name, "$(name)_", typename, "Ptr{$typename}", false,
                    "$(name)_ = Ref{$typename}($name)", "", false, stars, false)
    end
    typename_ccall = is_handle(r, typename) ? r.handles[typename].c : typename
    init, extract, name_ccall = init_extract(r, typename, name, isarray, isoutput, stars)
    ccall_str = "Ptr{"^stars * typename_ccall * "}"^stars
    if isarray
        ccall_str = "Ptr{$ccall_str}"
        if !isoutput
            typename = typename == "Cchar" ? "String" : "Vector{$typename}"
        else
            typename = "Vector{$typename}"
        end
    end
    # --- generic rules beyond the original heuristics ---------------------------------
    if !isarray && stars == 1 && is_handle(r, typename) && !isoutput && !is_destroy(fn)
        # a handle written back into a caller-supplied object (Get*, Duplicate into existing, ...)
        extract = "$(name).ptr = $(name_ccall)[]"
    end
    byref = get(ov, "byref", false) || (is_destroy(fn) && !isarray && stars == 1 && !is_handle(r, typename) &&
                                         !isoutput && typename != "Ptr{Cvoid}" && !a.isfunction &&
                                         !endswith(typename, "Fn") && typename != "Cvoid")
    if byref
        # opaque handle passed by reference: accept the handle or a Ref to it
        typename = "Union{$typename, Ref{$typename}}"
        name_ccall = "$(name)_"
        init = "$name_ccall = $name isa Base.RefValue ? $name : Ref{$(a.typename == typename ? typename : map_type(r, a.typename))}($name)"
        extract = ""
    end
    if isarray && isoutput && stars > 0 && (haskey(ov, "size") || haskey(ov, "prelude") || get(ov, "nullinit", false))
        base = "Ptr{"^stars * (is_handle(r, typename) ? r.handles[typename].c : replace(typename, "Vector{" => "", "}" => "")) * "}"^stars
        init = get(ov, "nullinit", false) ? "$name_ccall = Ref{$base}(C_NULL)" : "$name_ccall = Ref{$base}()"
        sz = get(ov, "size", "VecGetLocalSize(petsclib, x)")
        pre = get(ov, "prelude", "")
        extract = (isempty(pre) ? "" : pre * "\n\t") * "$name = unsafe_wrap(Array, $name_ccall[], $sz; own = false)"
    end
    if isarray && isoutput && stars == 0 && haskey(ov, "len")
        init = "$name = Vector{$(replace(typename, "Vector{" => "", "}" => ""))}(undef, $(ov["len"]))"
    end
    if get(ov, "nullable", false) && !isoutput
        typename = startswith(typename, "Union{") ? typename : "Union{Ptr, $typename}"
    end
    FArg(name, name_ccall, typename, ccall_str, isoutput, init, extract, isarray, stars, a.isfunction)
end

classify_all(r::Rules, fn::Fn, input_vars, output_vars) =
    FArg[classify(r, fn, a, input_vars, output_vars) for a in fn.args]
