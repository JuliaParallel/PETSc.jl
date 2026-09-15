# Compare every `LibPETSc.X(...)` call in src/ (high-level), ext/, test/ and examples/ with the
# generated stub of X: reports calls whose positional argument count does not match, and
# references to names that no longer exist.
#
#   julia wrapping/generator/check_callers.jl [AUTOWRAPPED_DIR]
root = dirname(dirname(@__DIR__))
wrapdir = length(ARGS) >= 1 ? ARGS[1] : joinpath(root, "src", "autowrapped")

# arity of generated stubs: function X(petsclib::PetscLibType, a::T, b::U) end
stubs = Dict{String,Int}()
defined = Set{String}()
for f in readdir(wrapdir)
    endswith(f, ".jl") || continue
    txt = read(joinpath(wrapdir, f), String)
    for m in eachmatch(r"(?m)^function (\w+)\(petsclib::PetscLibType(.*)\)(?: end)?\s*$", txt)
        args = strip(m.captures[2])
        n = isempty(args) ? 0 : count(",", replace(args, r"\{[^{}]*(\{[^{}]*\})*[^{}]*\}" => "")) # commas outside braces
        stubs[m.captures[1]] = n
    end
    for m in eachmatch(r"(?m)^(?:const|mutable struct|struct|abstract type|primitive type|@enum|@for_petsc function|function)\s+(\w+)", txt)
        push!(defined, m.captures[1])
    end
    for m in eachmatch(r"(?m)^(\w+)\s*=\s*", txt); push!(defined, m.captures[1]); end
    for m in eachmatch(r"(?m)^\s+(\w+) = \d+\s*$", txt); push!(defined, m.captures[1]); end   # enum members
end

function split_args(s)
    parts = String[]; depth = 0; buf = IOBuffer()
    for c in s
        (c in "({[") && (depth += 1); (c in ")}]") && (depth -= 1)
        if c == ',' && depth == 0; push!(parts, String(take!(buf))) else write(buf, c) end
    end
    push!(parts, String(take!(buf)))
    filter!(p -> !isempty(strip(p)), parts)
end

problems = String[]
for dir in ("src", "ext", "test", "examples")
    d = joinpath(root, dir)
    isdir(d) || continue
    for (r, _, fs) in walkdir(d), f in fs
        endswith(f, ".jl") || continue
        occursin("autowrapped", r) && continue
        path = joinpath(r, f)
        lines = split(read(path, String), '\n')
        for (i, l) in enumerate(lines)
            for m in eachmatch(r"LibPETSc\.(\w+)(\()?", l)
                name = m.captures[1]
                if !(name in defined) && !(name in ("PetscLibType", "petsclibs", "PETSC_COMM_SELF", "@for_petsc", "@chk", "getlib", "PetscError", "PetscErrorCode", "UnionPetscLibType", "scalartype", "inttype", "realtype", "MPI_Comm", "libs", "PetscInt", "PetscScalar", "PetscReal", "PetscBool", "PetscLib", "PETSC_DECIDE", "PETSC_DETERMINE"))
                    push!(problems, "$path:$i: LibPETSc.$name is not defined by the generated files")
                end
                (m.captures[2] === nothing || !haskey(stubs, name)) && continue
                # collect the call's argument text (may span lines)
                start = m.offset + length(m.match)
                buf = IOBuffer(); depth = 1; j = i; pos = start
                while depth > 0 && j <= length(lines)
                    s = lines[j]
                    while pos <= lastindex(s)
                        c = s[pos]
                        c == '(' && (depth += 1); c == ')' && (depth -= 1)
                        depth == 0 && break
                        write(buf, c); pos = nextind(s, pos)
                    end
                    depth == 0 && break
                    write(buf, ' '); j += 1; pos = 1
                end
                args = split_args(String(take!(buf)))
                nargs = length(args) - 1   # minus petsclib
                any(startswith(strip(a), "...") || endswith(strip(a), "...") for a in args) && continue
                any(occursin("=", a) && !occursin("==", a) for a in args) && continue   # keyword args
                if nargs != stubs[name]
                    push!(problems, "$path:$i: LibPETSc.$name called with $nargs args, generated stub takes $(stubs[name])")
                end
            end
        end
    end
end
foreach(println, problems)
println("\n$(length(problems)) problems")
