# Splitting autowrapped files into named blocks (used by golden_diff.jl and bootstrap_rules.jl)

struct Block
    file::String
    text::String
end

function split_wrappers(file::AbstractString, text::AbstractString, blocks::Dict{String,Block})
    lines = split(text, '\n'; keepempty=true)
    starts = Int[]
    for i in 1:length(lines)-1
        if strip(lines[i]) == "\"\"\"" && occursin(r"^\s+\S.*\(petsclib", lines[i+1])
            push!(starts, i)
        end
    end
    header_end = isempty(starts) ? length(lines) : starts[1] - 1
    header = join(lines[1:header_end], '\n')
    for m in eachmatch(r"(?m)^(?:const|mutable struct)\s+(\w+)", header)
        n = m.captures[1]
        startswith(n, "_n_") && continue
        blocks["type:$n"] = Block(file, m.match)
    end
    push!(starts, length(lines) + 1)
    for k in 1:length(starts)-1
        chunk = join(lines[starts[k]:starts[k+1]-1], '\n')
        m = match(r"(?m)^@for_petsc function (\w+)\(", chunk)
        m === nothing && (m = match(r"(?m)^function (\w+)\(", chunk))
        m === nothing && (m = match(r"(?m)^(?:LibPETSc\.)?@for_petsc function (?:LibPETSc\.)?(\w+)\(", chunk))
        name = m === nothing ? "chunk:$file:$k" : m.captures[1]
        key = "fn:$name"
        haskey(blocks, key) && (key = "fn:$name#$(file)")
        blocks[key] = Block(file, chunk)
    end
end

function split_simple(file, text, blocks, re, prefix)
    for m in eachmatch(re, text)
        blocks["$prefix:$(m.captures[1])"] = Block(file, m.match)
    end
end

function load_blocks(dir::AbstractString)
    blocks = Dict{String,Block}()
    for f in sort(readdir(dir))
        endswith(f, ".jl") || continue
        text = read(joinpath(dir, f), String)
        if f == "enums_wrappers.jl"
            split_simple(f, text, blocks, r"(?ms)^@enum (\w+).*?^end", "enum")
        elseif f == "senums_wrappers.jl"
            split_simple(f, text, blocks, r"(?m)^(\w+)=.*$", "senum")
        elseif f == "typedefs_wrappers.jl"
            split_simple(f, text, blocks, r"(?m)^const (\w+) = .*$", "typedef")
            m = match(r"(?ms)^primitive type (PetscBool).*?^Base\.iszero.*?$", text)
            m === nothing || (blocks["block:PetscBool"] = Block(f, m.match))
        elseif f == "struct_wrappers.jl"
            split_simple(f, text, blocks, r"(?ms)^(?:mutable )?struct (\w+).*?^end", "struct")
        elseif f == "opaque_types.jl"
            for m in eachmatch(r"(?m)^(?:const|mutable struct)\s+(\w+)", text)
                n = m.captures[1]
                startswith(n, "_n_") && continue
                blocks["type:$n"] = Block(f, m.match)
            end
        elseif f in ("petsc_library.jl", "petscarray.jl", "petsc_wrappers_version.jl")
            blocks["file:$f"] = Block(f, text)
        else
            split_wrappers(f, text, blocks)
        end
    end
    return blocks
end

"""Lines of a block with trailing whitespace removed, blank/comment/doc-link lines dropped."""
function code_lines(s::AbstractString)
    out = String[]
    for l in split(s, '\n')
        t = rstrip(l)
        isempty(t) && continue
        startswith(t, "\$(_doc_external(") && continue
        startswith(strip(t), "#") && continue
        push!(out, String(strip(t)))
    end
    return out
end
normalize(s) = join(code_lines(s), '\n')

"""Docstring region (between the opening and closing triple quotes) and the rest of a block."""
function split_doc(s::AbstractString)
    m = match(r"(?s)^\s*\"\"\"\n(.*?)\n\"\"\"\n(.*)$", s)
    m === nothing && return "", String(s)
    return String(m.captures[1]), String(m.captures[2])
end

"""Coarse classification of why a golden block and a new block differ."""
function categorize(g::AbstractString, n::AbstractString)
    gd, gc = split_doc(g); nd, nc = split_doc(n)
    gl = code_lines(gc); nl = code_lines(nc)
    if gl == nl
        return "docs-only"
    end
    go = setdiff(gl, nl); no = setdiff(nl, gl)   # lines only on one side
    j(v) = join(v, "\n")
    G = j(go); N = j(no)
    if all(occursin(".ptr = C_NULL", l) for l in go) && all(occursin(r"\.ptr = \w+_\[\]", l) for l in no) && !isempty(go)
        return "writeback-fix"
    end
    if occursin(r"Ptr\{(Ptr\{)?(IS|PetscVec|PetscMat|PetscDM|PetscKSP|PetscSNES|TS|AO|Tao|PF)\}", G) && occursin(r"Ptr\{(Ptr\{)?C(IS|Vec|Mat|DM|KSP|SNES|TS|AO|Tao|PF)\}", N)
        return "handle-ccall-fix"
    end
    if occursin("::Ptr{Cvoid}", G) && occursin(r"::\w*Fn\b", N)
        return "fnptr-typed"
    end
    if occursin(r"Union\{\w+, Ref\{\w+\}\}", G) && !occursin("Union{", N)
        return "byref-missing"
    end
    if occursin(r"Union\{\w+, Ref\{\w+\}\}", N) && !occursin("Union{", G)
        return "byref-new"
    end
    if occursin("unsafe_wrap(Array", G) && occursin("VecGetLocalSize(petsclib, x)", N)
        return "size-rule"
    end
    if occursin(r"Vector\{[^}]*\}\(undef, ", G) && occursin("(undef, ni)", N)
        return "len-rule"
    end
    if occursin("Union{Ptr", G) && !occursin("Union{Ptr", N)
        return "nullable-rule"
    end
    if occursin(r"\bRef\{\w+\}\(\)", G) && occursin(r"\.ptr = \w+_\[\]", N) && occursin(r"return \w", G)
        return "returns-new-object"
    end
    if occursin("::Vector{Cchar}", G) && occursin("::String", N)
        return "string-arg-fix"
    end
    # golden applied substring renames (local->loc, global->glob, end->end_, function->fnc) inside identifiers
    ren(l) = replace(l, "function" => "fnc", "end" => "end_", "global" => "glob", "local" => "loc")
    if [ren(l) for l in nl] == gl
        return "arg-rename-fix"
    end
    if occursin(r"Ref\(pointer\(", N) && occursin(r"_ = Ref\{Ptr\{", G) && occursin("Restore", N)
        return "restore-input-fix"
    end
    if occursin(r"= Ref\(\w+\.ptr\)", N) && !occursin(r"= Ref\(\w+\.ptr\)", G) && occursin(r"\.ptr = \w+_\[\]", N)
        return "handle-out-byvalue-fix"
    end
    if occursin("::Cvoid", G) && occursin("::Ptr{Cvoid}", N)
        return "voidptr-fix"
    end
    if occursin(r"::\w+Fn\b", G) && occursin("::Ptr{Cvoid}", N)
        return "fnptr-ptrcvoid-fix"
    end
    if occursin("unsafe_string(", G) && occursin(r"^\s*\w+ = \w+_\[\]$"m, N) && !occursin("unsafe_string(", N)
        return "enum-out-fix"
    end
    if (occursin("Union{Ptr{", G) || occursin("Union{Ptr,", G)) && occursin("Union{Ptr, ", N)
        return "nullable-form"
    end
    if occursin(r"Ref\{(IS|TS|AO|PF|Tao|PetscVec|PetscMat|PetscDM|PetscKSP|PetscSNES)\}\(\)", G) && occursin(r"Ref\{C\w+\}\(\)", N)
        return "handle-return-fix"
    end
    if occursin("(undef, ni)", G) && !occursin("(undef, ni)", N)
        return "ni-alloc-removed"
    end
    istuple(l) = occursin(r"^\s*\(.*\),\s*$", l)
    gt = [l for l in gl if istuple(l)]; nt = [l for l in nl if istuple(l)]
    if !isempty(gt) && !isempty(nt) && count("Ptr{", nt[1]) > count("Ptr{", gt[1]) && !any(occursin("Ref", l) for l in go)
        return "byvalue-ptr-fix"
    end
    outnames(d) = begin
        h = match(r"^\s*(.*?)\s*=\s*\w+\(petsclib", d)
        h === nothing ? String[] : [String(first(split(strip(x), "::"))) for x in split(h.captures[1], ",")]
    end
    gouts = outnames(gd); nouts = outnames(nd)
    if gouts != nouts
        if length(nouts) > length(gouts)
            extra = setdiff(nouts, gouts)
            # the new outputs are enums (or other plain by-value types) that the old heuristics refused
            hdr = match(r"^\s*(.*?)\s*=\s*\w+\(petsclib", nd)
            types = Dict(String(first(split(strip(x), "::"))) => String(last(split(strip(x), "::"))) for x in split(hdr.captures[1], ","))
            if all(!occursin(r"Vector|Ptr|Petsc(Int|Real|Scalar|Bool)|^(IS|TS|AO|PF|Tao|PetscVec|PetscMat|PetscDM|PetscKSP|PetscSNES)$", get(types, e, "?")) for e in extra)
                return "enum-direction-fix"
            end
            return "direction-new-output"
        end
        return "direction-golden-output"
    end
    return "other"
end
