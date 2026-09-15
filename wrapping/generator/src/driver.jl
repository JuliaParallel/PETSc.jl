# Putting it together: which functions go where, and writing the output directory

struct Plan
    file::String
    functions::Vector{String}
end

function plan_files(api::API, r::Rules)
    plans = Plan[]
    seen = Set{String}()
    for f in r.files
        names = String[]
        if get(f, "standalone", false)
            append!(names, api.standalone)
        end
        for c in get(f, "classes", String[])
            haskey(api.classes, c) || (@warn "class $c not in API snapshot"; continue)
            append!(names, api.classes[c])
        end
        filter!(n -> !occursin("_", n) && !haskey(r.exclude, n) && !(n in seen), names)
        union!(seen, names)
        push!(plans, Plan(f["name"], sort!(names)))
    end
    unmapped = setdiff(Set(keys(api.classes)), Set(c for f in r.files for c in get(f, "classes", String[])))
    isempty(unmapped) || @warn "classes not assigned to any file" unmapped = sort!(collect(unmapped))
    return plans
end

"""Split a Julia argument list on commas that are not inside braces/parentheses."""
function split_toplevel(str::AbstractString)
    parts = String[]; depth = 0; buf = IOBuffer()
    for c in str
        (c == '{' || c == '(') && (depth += 1)
        (c == '}' || c == ')') && (depth -= 1)
        if c == ',' && depth == 0
            push!(parts, String(take!(buf)))
        else
            write(buf, c)
        end
    end
    push!(parts, String(take!(buf)))
    return parts
end

"""Type names an argument list refers to (as the original `extract_function_typeargs_from_class`)."""
function referenced_types(args::Vector{FArg})
    ts = String[]
    for a in args
        for m in eachmatch(r"[A-Za-z_]\w*", a.typename)
            t = m.match
            t in ("Union", "Ptr", "Ref", "Vector") && continue
            push!(ts, t)
        end
    end
    return ts
end

function load_overrides(dir::AbstractString)
    ov = Dict{String,String}()
    isdir(dir) || return ov
    for f in readdir(dir)
        endswith(f, ".jl") || continue
        ov[f[1:end-3]] = read(joinpath(dir, f), String)
    end
    return ov
end

"""
    generate(; api_json, petsc_dir, outdir, wrapping_dir)

Generates the complete `src/autowrapped` directory into `outdir`.
"""
function generate(; api_json::AbstractString, petsc_dir::AbstractString, outdir::AbstractString,
                  wrapping_dir::AbstractString = dirname(@__DIR__), verbose::Bool = true)
    t0 = time()
    api = load_api(api_json)
    r = load_rules(joinpath(wrapping_dir, "rules"))
    union!(r.enum_types, keys(api.enums))
    union!(r.string_types, keys(api.senums))
    union!(r.struct_types, keys(api.structs))
    overrides = load_overrides(joinpath(wrapping_dir, "overrides"))
    verbose && println("API $(api.version): $(length(api.functions)) functions; building docstring index ...")
    docs = build_docindex(petsc_dir)
    verbose && println("  indexed $(length(docs.blocks)) manual pages in $(round(time()-t0, digits=1)) s")
    mkpath(outdir)

    prologue = read(joinpath(wrapping_dir, "prologue.jl"), String)
    structs = read(joinpath(wrapping_dir, "structs.jl"), String)
    petscbool = read(joinpath(wrapping_dir, "petscbool.jl"), String)
    # typedefs the prologue already defines (PetscInt, PetscScalar, ... are library dependent) or hand-defined
    typedef_exclude = Set(["PetscGeom", "PetscInt32", "PetscBool", "PetscComplex"]) ∪ defined_names(prologue)
    known = defined_names(prologue) ∪ defined_names(structs) ∪ Set(["PetscBool", "PETSC_TRUE", "PETSC_FALSE"])
    union!(known, keys(api.enums), keys(api.senums), setdiff(keys(api.typedefs), typedef_exclude), keys(api.structs))
    union!(known, r.predeclared)
    for h in values(r.handles)
        push!(known, h.julia); push!(known, h.abstract); push!(known, h.c)
    end
    isknown(t) = t in known || t == "String" || isdefined(Base, Symbol(t)) || isdefined(Core, Symbol(t)) ||
                 startswith(t, "Libc.") || startswith(t, "\$")

    td_names, td_values, td_skipped = plan_typedefs(api, r, typedef_exclude, isknown)
    setdiff!(known, td_skipped)          # so they get declared as opaque types below

    plans = plan_files(api, r)
    opaque = Set{String}()
    nfun = 0
    used_overrides = Set{String}()
    for p in plans
        open(joinpath(outdir, p.file), "w") do io
            for name in p.functions
                fn = api.functions[name]
                if haskey(overrides, name)
                    print(io, overrides[name])
                    endswith(overrides[name], "\n") || println(io)
                    push!(used_overrides, name)
                    nfun += 1
                    continue
                end
                input_vars, output_vars, doc_lines = function_docs(docs, name)
                args = classify_all(r, fn, input_vars, output_vars)
                for t in referenced_types(args)
                    isknown(t) || push!(opaque, t)
                end
                render_function(io, r, fn, args, doc_lines, manual_section(docs, name, fn.mansec))
                nfun += 1
            end
        end
    end
    # hand-written extra functions (not in the API) live in overrides/ too
    extras = sort!(collect(setdiff(keys(overrides), used_overrides)))
    open(joinpath(outdir, "extra_wrappers.jl"), "w") do io
        println(io, "# Hand-written wrappers for functions/macros not described by getAPI.py (wrapping/overrides/)")
        for name in extras
            print(io, overrides[name])
            endswith(overrides[name], "\n") || println(io)
        end
    end

    open(joinpath(outdir, "enums_wrappers.jl"), "w") do io
        write_enums(io, api, Set(["KSPConvergedReason", "PetscMemType"]))
    end
    open(joinpath(outdir, "senums_wrappers.jl"), "w") do io
        write_senums(io, api, r)
    end
    open(joinpath(outdir, "typedefs_wrappers.jl"), "w") do io
        write_typedefs(io, td_names, td_values, [petscbool])
    end
    write(joinpath(outdir, "struct_wrappers.jl"), structs)
    write(joinpath(outdir, "petscarray.jl"), read(joinpath(wrapping_dir, "petscarray.jl"), String))
    open(joinpath(outdir, "opaque_types.jl"), "w") do io
        write_opaque_types(io, collect(opaque))
    end
    open(joinpath(outdir, "petsc_wrappers_version.jl"), "w") do io
        write_version(io, api.version)
    end
    includes = Tuple{String,Bool}[("petscarray.jl", true), ("enums_wrappers.jl", true), ("senums_wrappers.jl", true),
                                  ("typedefs_wrappers.jl", true), ("struct_wrappers.jl", true), ("opaque_types.jl", true)]
    for f in r.files
        push!(includes, (f["name"], get(f, "include", true)))
    end
    push!(includes, ("extra_wrappers.jl", true))
    open(joinpath(outdir, "petsc_library.jl"), "w") do io
        write_library_file(io, prologue, includes)
    end
    # static sanity check: every type name used in a signature or ccall tuple must be defined somewhere
    defined = copy(known)
    union!(defined, opaque)
    for f in readdir(outdir)
        endswith(f, ".jl") || continue
        f in ("petsc_library.jl",) && continue
        union!(defined, defined_names(read(joinpath(outdir, f), String)))
    end
    missing = Set{String}()
    for p in plans
        txt = read(joinpath(outdir, p.file), String)
        for m in eachmatch(r"(?m)^@for_petsc function \w+\((.*)\)\s*$", txt)
            for arg in split_toplevel(m.captures[1])
                occursin("::", arg) || continue
                for id in eachmatch(r"[A-Za-z_]\w*", last(split(arg, "::"; limit = 2)))
                    n = id.match
                    (isknown(n) || n in defined || n in ("Union", "Ptr", "Ref", "Vector", "Nothing", "AbstractString", "AbstractArray", "UnionPetscLib", "MPI")) || push!(missing, n)
                end
            end
        end
        for m in eachmatch(r"(?m)^\s+\(([^:].*)\),\s*$", txt), id in eachmatch(r"[A-Za-z_]\w*", m.captures[1])
            n = id.match
            (isknown(n) || n in defined || n in ("Ptr", "Ref", "Nothing", "MPI")) || push!(missing, n)
        end
    end
    isempty(missing) || @warn "type names used but not defined anywhere in the output" missing = sort!(collect(missing))
    verbose && println("wrote $nfun functions in $(length(plans)) files, $(length(opaque)) opaque types, $(length(extras)) extra wrappers to $outdir in $(round(time()-t0, digits=1)) s")
    return nothing
end
