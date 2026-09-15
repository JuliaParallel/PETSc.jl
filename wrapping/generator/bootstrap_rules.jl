# Mine the hand-edited golden wrappers for per-argument rules and write rules/args.toml.
#
#   julia --project=wrapping/generator wrapping/generator/bootstrap_rules.jl GOLDEN_DIR
#
# Recognised patterns in a golden @for_petsc body:
#   X_ = Ref{Ptr{T}}(C_NULL)                          -> nullinit
#   X = unsafe_wrap(Array, X_[], EXPR; own = false)    -> size = EXPR (unless the placeholder)
#   assignment lines between the ccall and the wrap    -> prelude
#   X = Vector{T}(undef, EXPR)                         -> len = EXPR (unless "ni")
#   name::Union{X, Ref{X}} in the signature            -> byref (non-Destroy functions only)
#   name::Union{..., Ptr} / Union{Ptr, ...}            -> nullable
include(joinpath(@__DIR__, "src", "blocks.jl"))

function mine(golden::AbstractString)
    blocks = load_blocks(golden)
    rules = Dict{String,Dict{String,Dict{String,Any}}}()
    add!(fn, arg, k, v) = (get!(get!(rules, fn, Dict{String,Dict{String,Any}}()), arg, Dict{String,Any}())[k] = v)
    for (key, b) in blocks
        startswith(key, "fn:") || continue
        fn = key[4:end]
        occursin("#", fn) && continue
        isfile(joinpath(@__DIR__, "overrides", fn * ".jl")) && continue   # verbatim override, no rules needed
        _, code = split_doc(b.text)
        sm = match(r"(?m)^@for_petsc function \w+\(petsclib::\$UnionPetscLib,?", code)
        sm === nothing && continue
        # argument list up to the matching parenthesis (signatures may span several lines)
        rest = code[sm.offset + length(sm.match):end]
        depth = 1; stop = 0
        for (k, c) in enumerate(rest)
            c == '(' && (depth += 1); c == ')' && (depth -= 1)
            if depth == 0; stop = k - 1; break; end
        end
        arglist = replace(rest[1:stop], "\n" => " ")
        args = String[]; depth = 0; cur = IOBuffer()
        for c in arglist
            if c == '{' || c == '('; depth += 1 elseif c == '}' || c == ')'; depth -= 1 end
            if c == ',' && depth == 0
                push!(args, String(take!(cur)))
            else
                write(cur, c)
            end
        end
        push!(args, String(take!(cur)))
        for a in args
            m = match(r"^\s*(\w+)::(.+?)\s*$", a)
            m === nothing && continue
            name, typ = m.captures
            mb = match(r"^Union\{(\w+), Ref\{\1\}\}$", typ)
            if mb !== nothing
                # Union{X, Ref{X}} inputs of Get functions are superseded by returning the value
            elseif occursin(r"Union\{.*\bPtr\b.*\}", typ)
                add!(fn, name, "nullable", true)
            end
        end
        # a handle output the heuristics would treat as an input, hand-changed to return a new object
        for m in eachmatch(r"^\s*(\w+)\s*=\s*(PetscVec|PetscMat|PetscDM|PetscKSP|PetscSNES|IS|TS|AO|Tao|PF|PetscOptions)(?:\{PetscLib\})?\((?:\1_\[\], petsclib|C_NULL, \w+\.age)\)"m, code)
            if !occursin("Create", fn) && !occursin("Duplicate", fn)
                add!(fn, m.captures[1], "direction", "out")
            end
        end
        lines = [strip(l) for l in split(code, '\n')]
        ccall_end = findfirst(l -> l == ")", lines)
        for (i, l) in enumerate(lines)
            m = match(r"^(\w+) = unsafe_wrap\(Array, (\w+)_\[\], (.+); own = false\)", l)
            if m !== nothing && m.captures[3] != "VecGetLocalSize(petsclib, x)" && m.captures[1] == m.captures[2]
                add!(fn, m.captures[1], "size", String(m.captures[3]))
                if ccall_end !== nothing
                    pre = String[]
                    for j in ccall_end+1:i-1
                        pl = lines[j]
                        (isempty(pl) || startswith(pl, "#")) && continue
                        occursin(r"^[\w, ]+ = ", pl) && !occursin("unsafe_wrap", pl) && push!(pre, pl)
                    end
                    isempty(pre) || add!(fn, m.captures[1], "prelude", join(pre, "\n\t"))
                end
            end
            m = match(r"^(\w+)_ = Ref\{Ptr\{.*\}\}\(C_NULL\)", l)
            m === nothing || add!(fn, m.captures[1], "nullinit", true)
            m = match(r"^(\w+) = Vector\{[^}]*\}\(undef, (.+)\)\s*;?\s*$", l)
            if m !== nothing && strip(m.captures[2]) != "ni" && !occursin("CHECK SIZE", l)
                add!(fn, m.captures[1], "len", String(strip(m.captures[2])))
            end
        end
    end
    # fallback independent of block splitting: nullable inputs from the untyped stubs
    for f in readdir(golden)
        endswith(f, "_wrappers.jl") || continue
        for m in eachmatch(r"(?m)^function (\w+)\(petsclib::PetscLibType,?(.*)\) end\s*$", read(joinpath(golden, f), String))
            fn = m.captures[1]
            isfile(joinpath(@__DIR__, "overrides", fn * ".jl")) && continue
            for a in split_toplevel_args(m.captures[2])
                am = match(r"^\s*(\w+)::(.+?)\s*$", a)
                am === nothing && continue
                name, typ = am.captures
                if occursin(r"Union\{.*\bPtr\b.*\}", typ) && match(r"^Union\{(\w+), Ref\{\1\}\}$", typ) === nothing
                    add!(fn, name, "nullable", true)
                end
            end
        end
    end
    return rules
end

function split_toplevel_args(str)
    parts = String[]; depth = 0; buf = IOBuffer()
    for c in str
        (c == '{' || c == '(') && (depth += 1); (c == '}' || c == ')') && (depth -= 1)
        if c == ',' && depth == 0; push!(parts, String(take!(buf))) else write(buf, c) end
    end
    push!(parts, String(take!(buf)))
    return parts
end

function write_toml(rules, path)
    open(path, "w") do io
        println(io, "# Per-function argument overrides MINED from the hand-edited wrappers by bootstrap_rules.jl.")
        println(io, "# Regenerated wholesale by that script: do not edit by hand, use rules/args.toml instead.\n")
        for fn in sort(collect(keys(rules)))
            for arg in sort(collect(keys(rules[fn])))
                println(io, "[$fn.$arg]")
                for (k, v) in sort(collect(rules[fn][arg]); by = first)
                    if v isa Bool
                        println(io, "$k = $v")
                    else
                        println(io, "$k = ", repr(String(v)))
                    end
                end
                println(io)
            end
        end
    end
end

rules = mine(ARGS[1])
write_toml(rules, joinpath(@__DIR__, "rules", "args_mined.toml"))
n = sum(length(v) for v in values(rules))
println("wrote $n argument rules for $(length(rules)) functions to rules/args_mined.toml")
