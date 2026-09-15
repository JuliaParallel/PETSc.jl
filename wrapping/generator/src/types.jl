# Rules loading and C -> Julia type mapping

struct Handle
    c_name::String
    julia::String
    abstract::String
    c::String
end

struct Rules
    replace_order::Vector{String}
    replace_map::Dict{String,String}
    fix_substring::Bool
    rename_args::Dict{String,String}
    handles::Dict{String,Handle}          # keyed by the Julia struct name (after mapping)
    handle_c_names::Dict{String,Handle}   # keyed by the C name
    predeclared::Set{String}
    dispatch_types::Vector{String}
    senum_overrides::Dict{String,String}
    exclude::Dict{String,String}
    files::Vector{Dict{String,Any}}
    args::Dict{String,Dict{String,Dict{String,Any}}}   # fn => arg => overrides
    enum_types::Set{String}       # C enums (filled from the API snapshot)
    string_types::Set{String}     # PETSc string enums such as KSPType (filled from the API snapshot)
    struct_types::Set{String}     # C structs passed by value/pointer (filled from the API snapshot)
end

function load_rules(dir::AbstractString)
    t = TOML.parsefile(joinpath(dir, "types.toml"))
    f = TOML.parsefile(joinpath(dir, "files.toml"))
    # mined rules first, hand-written rules on top (argument tables merge, hand keys win)
    a = Dict{String,Any}()
    for f in ("args_mined.toml", "args.toml")
        isfile(joinpath(dir, f)) || continue
        for (fn, tab) in TOML.parsefile(joinpath(dir, f))
            dst = get!(a, fn, Dict{String,Any}())
            for (an, ov) in tab
                merge!(get!(dst, an, Dict{String,Any}()), ov)
            end
        end
    end
    rep = t["replace"]
    order = String.(rep["order"])
    rmap = Dict(k => String(rep[k]) for k in order)
    handles = Dict{String,Handle}()
    hc = Dict{String,Handle}()
    for h in t["handles"]
        H = Handle(h["c_name"], h["julia"], h["abstract"], h["c"])
        handles[H.julia] = H
        hc[H.c_name] = H
    end
    args = Dict{String,Dict{String,Dict{String,Any}}}()
    for (fn, tab) in a
        args[fn] = Dict(String(an) => Dict{String,Any}(ov) for (an, ov) in tab)
    end
    Rules(order, rmap, get(t, "fix_substring_replacements", false),
        Dict(String(k) => String(v) for (k, v) in get(t, "rename_args", Dict())),
        handles, hc, Set(String.(get(t["predeclared"], "names", String[]))),
        String.(t["dispatch_types"]),
        Dict(String(k) => String(v) for (k, v) in get(t, "senum_overrides", Dict())),
        Dict(String(k) => String(v) for (k, v) in get(f, "exclude", Dict())),
        Vector{Dict{String,Any}}(f["file"]), args, Set{String}(), Set{String}(), Set{String}())
end

"""C type name -> Julia type name (the original generator's `replace_types`)."""
function map_type(r::Rules, typename::AbstractString)
    s = replace(replace(String(typename), "std::" => ""), "::" => "_")   # std::size_t -> size_t, moab::Range -> moab_Range
    occursin("(", s) && return "Ptr{Cvoid}"          # inline function pointer `void (*f0)(...)`
    for k in r.replace_order
        pat = r.fix_substring ? Regex("\\b" * k * "\\b") : k
        s = replace(s, pat => r.replace_map[k])
    end
    for (cn, h) in r.handle_c_names
        cn == h.julia && continue
        s = replace(s, Regex("\\b" * cn * "\\b") => h.julia)
    end
    return s
end

rename_arg(r::Rules, name::AbstractString) = get(r.rename_args, name, String(name))

"""Add `\$` in front of the library-dependent scalar/int types (`replace_dispatch_types`)."""
function dispatch(r::Rules, s::AbstractString)
    out = String(s)
    for t in r.dispatch_types
        out = replace(out, t => "\$" * t)
    end
    return out
end

is_handle(r::Rules, jtype::AbstractString) = haskey(r.handles, jtype)

"""Type annotation of an input argument: abstract supertype for handles (PR #263)."""
function abstract_arg_type(r::Rules, typename::AbstractString)
    m = match(r"^AbstractArray\{(\w+)\}$", typename)
    m !== nothing && return String(typename)
    m = match(r"^Union\{Ptr, (.*)\}$", typename)
    m !== nothing && return "Union{Ptr, $(abstract_arg_type(r, m.captures[1]))}"
    m = match(r"^Vector\{(\w+)\}$", typename)
    if m !== nothing
        inner = m.captures[1]
        return is_handle(r, inner) ? "Vector{<:$(r.handles[inner].abstract)}" : String(typename)
    end
    return is_handle(r, typename) ? r.handles[typename].abstract : String(typename)
end
