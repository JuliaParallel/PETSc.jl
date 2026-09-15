# Docstring index: one pass over PETSc's src/**/*.c collecting the /*@ ... @*/ manual-page
# blocks, then the same cleaning rules the original find_doc_strings.jl applied (kept
# faithful so the generated docstrings are reproduced byte for byte).

struct DocIndex
    blocks::Dict{String,Vector{Tuple{String,Vector{String}}}}   # fn name => [(file, raw lines)] in walk order
    definers::Dict{String,Vector{String}}                       # fn name => files containing "PetscErrorCode NAME("
    petsc_dir::String
    mansec_cache::Dict{String,String}                           # directory => manual page section
end

const _EXCLUDED_DIRS = ("benchmarks", "tutorials", "tests", "petsc4py")

function _c_files(petsc_dir)
    files = String[]
    for (root, dirs, fs) in walkdir(joinpath(petsc_dir, "src"))
        for f in fs
            endswith(f, ".c") && push!(files, joinpath(root, f))
        end
    end
    filter!(f -> !any(occursin(d, f) for d in _EXCLUDED_DIRS), files)
    return files
end

const _DEF_RE = r"PetscErrorCode (\w+)\("

function build_docindex(petsc_dir::AbstractString)
    blocks = Dict{String,Vector{Tuple{String,Vector{String}}}}()
    definers = Dict{String,Vector{String}}()
    for file in _c_files(petsc_dir)
        txt = read(file, String)
        occursin("/*@", txt) || occursin("PetscErrorCode ", txt) || continue
        for m in eachmatch(_DEF_RE, txt)
            push!(get!(definers, m.captures[1], String[]), file)
        end
        lines = split(txt, '\n')
        i = 1
        while i <= length(lines)
            if occursin("/*@", strip(lines[i]))
                # first line after the marker names the function
                i += 1
                i > length(lines) && break
                first = lines[i]
                toks = split(strip(first))
                if !isempty(toks)
                    name = String(toks[1])
                    raw = String[]
                    # a block ends at "@*/"; a plain "*/" or the next "/*@" also ends it (PETSc typos)
                    while i <= length(lines) && !occursin("@*/", lines[i]) && strip(lines[i]) != "*/" &&
                          !occursin("/*@", lines[i])
                        push!(raw, String(lines[i]))
                        i += 1
                    end
                    if i <= length(lines) && occursin("/*@", lines[i])
                        i -= 1   # let the outer loop see the new block marker
                    end
                    push!(get!(blocks, name, Vector{Tuple{String,Vector{String}}}()), (file, raw))
                end
            end
            i += 1
        end
    end
    for v in values(definers)
        unique!(v)
    end
    DocIndex(blocks, definers, String(petsc_dir), Dict{String,String}())
end

# ---- faithful port of find_doc_strings.jl -------------------------------------------

function _add_backticks(line::AbstractString)
    parts = split(strip(line), "- ")
    if length(parts) > 1
        var_name = parts[2]
        n = length(var_name)
        parts[2] = rpad("`" * strip(var_name) * "`", n + 2)
        return join(parts .* "- ")[1:end-2]
    end
    return String(line)
end

# processing applied by the original read_c_function_docs to one raw block
function _process_block(raw::Vector{String})
    out = String[]
    for l in raw
        line = strip(l)
        if !isempty(line) && (startswith(line, ".") || startswith(line, "+"))
            line = "-" * line[2:end]
        end
        line = strip(line)
        if startswith(line, '-') && length(findall('-', line)) > 1
            line = _add_backticks(line)
        end
        push!(out, String(line))
    end
    return out
end

function _finish_block(comment::Vector{String})
    if !isempty(comment)
        comment[1] = strip(comment[1])
        parts = split(comment[1], "-")
        comment[1] = length(parts) >= 2 ? String(strip(parts[2])) : ""
    end
    comment = replace.(comment, "\$" => "")
    comment = replace.(comment, "[](ch_stag)," => "")
    comment = replace.(comment, "[](ch_dmbase)," => "")
    comment = replace.(comment, "-seealso:  " => "See also: \n=== \n")
    comment = replace.(comment, "seealso:  " => "See also: \n=== \n")
    return String.(comment)
end

# original: iterate defining files in walk order, return the first file's non-empty block
function raw_docs(idx::DocIndex, name::String)
    files = get(idx.definers, name, String[])
    blks = get(idx.blocks, name, Tuple{String,Vector{String}}[])
    for f in files
        lines = String[]
        for (bf, raw) in blks
            bf == f && append!(lines, _process_block(raw))
        end
        isempty(lines) || return _finish_block(lines)
    end
    return nothing
end

function _get_last_line(start, comment)
    for l in start:length(comment)
        isempty(strip(comment[l])) && return l - 1
    end
    return length(comment)
end

function _extract_variable(str::AbstractString)
    occursin('-', str) || return ""
    parts = split(str, '-')
    length(parts) < 2 && return ""
    var_part = replace(strip(parts[2]), "`" => "")
    return String(first(split(var_part, ' ')))
end

function extract_input_output_vars(comment::Vector{String})
    input_vars, output_vars = String[], String[]
    for (key, acc) in (("Input Parameter", input_vars), ("Output Parameter", output_vars))
        any(startswith.(comment, key)) || continue
        start = first(findall(startswith.(strip.(comment), key))) + 1
        stop = _get_last_line(start, comment)
        for l in start:stop
            v = _extract_variable(comment[l])
            isempty(v) || push!(acc, v)
        end
    end
    return input_vars, output_vars
end

function _remove_notes!(comment::Vector{String}, keyword::String)
    note_start = findfirst(c -> occursin(keyword, c), comment)
    note_start === nothing && return comment
    l = findfirst(c -> occursin("-seealso:", c), comment)
    l === nothing && return comment
    note_end = l - 1
    note_start <= note_end && deleteat!(comment, note_start:note_end)
    return comment
end

function remove_notes(comment::Vector{String})
    for k in ("Note:", "Notes:", "Developer Note:", "Fortran Notes:", "Example Usage:", "-vb", "Example Usage\\:")
        _remove_notes!(comment, k)
    end
    for i in eachindex(comment)
        comment[i] = replace(comment[i], "\\:" => ":", "``" => "`")
    end
    return comment
end

"""
    input_vars, output_vars, doc_lines = function_docs(idx, name)

`doc_lines` is `nothing` when no manual page was found.
"""
function function_docs(idx::DocIndex, name::String)
    comment = raw_docs(idx, name)
    comment === nothing && return String[], String[], nothing
    input_vars, output_vars = extract_input_output_vars(comment)
    return input_vars, output_vars, remove_notes(copy(comment))
end

const _FALLBACK_MANSEC = Dict("vec" => "Vec", "mat" => "Mat", "dm" => "DM", "ksp" => "KSP", "snes" => "SNES",
    "ts" => "TS", "tao" => "Tao", "sys" => "Sys", "viewer" => "Viewer", "ml" => "PC")

"""
    manual_section(idx, name, mansec)

Directory of the function's manual page on petsc.org (`DMPlex`, `PC`, ...): the `SUBMANSEC`, or
else the `MANSEC`, of the makefile in the directory of the defining C file, searching upwards.
"""
function manual_section(idx::DocIndex, name::String, mansec::String)
    files = get(idx.definers, name, String[])
    for f in files
        dir = dirname(f)
        while startswith(dir, idx.petsc_dir) && length(dir) > length(idx.petsc_dir)
            if haskey(idx.mansec_cache, dir)
                sec = idx.mansec_cache[dir]
                isempty(sec) || return sec
            else
                mk = joinpath(dir, "makefile")
                sec = ""
                if isfile(mk)
                    txt = read(mk, String)
                    m = match(r"(?m)^\s*SUBMANSEC\s*=\s*(\S+)", txt)
                    m === nothing && (m = match(r"(?m)^\s*MANSEC\s*=\s*(\S+)", txt))
                    m === nothing || (sec = String(m.captures[1]))
                end
                idx.mansec_cache[dir] = sec
                isempty(sec) || return sec
            end
            dir = dirname(dir)
        end
    end
    return get(_FALLBACK_MANSEC, lowercase(mansec), _titlecase_fallback(mansec))
end
_titlecase_fallback(s) = isempty(s) ? s : uppercase(s[1:1]) * lowercase(s[2:end])
