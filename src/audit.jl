"""
    audit_petsc_file(path::AbstractString)

Scan a PETSc.jl source file and report PETSc object creations and destroys.

This utility reads `path`, looks for common PETSc object creation patterns
(`KSP(...)`, `SNES(...)`, `DMDA(...)`, `DMStag(...)`, as well as
`VecCreate...`, `VecDuplicate...`, `Vec...WithArray(...)` and corresponding
`Mat...` creators via either `PETSc.` or `LibPETSc.`), plus common DM/Vec/Mat
allocators like `DMGlobalVec`, `DMLocalVec`, `DMGetCoordinateDM`,
`DMStagCreateCompatibleDMStag`, and `DMCreateMatrix`; and for explicit
`destroy(x)` calls (optionally prefixed with `PETSc.`).

It prints three sections:
 - CREATION statements: line numbers and variables assigned to created objects
 - DESTROY calls: line numbers and variables destroyed
 - POSSIBLY UNDESTROYED objects: variables that appear in creations but have no
     matching destroy call
    - FINALIZE calls: presence of `PETSc.finalize(petsclib)`; suggests adding it
        if missing

Notes:
 - Heuristic only: matches common constructors and `LibPETSc` creation routines;
     it does not follow control flow or scopes.
 - Creations without assignment cannot be cross-checked against destroys.
 - Objects produced via helpers (e.g., `similar(...)`) are not detected unless
     they use known creation patterns.

Returns `nothing`.
"""
function audit_petsc_file(path::AbstractString)
    # Read full content and strip block comments (#= ... =#) before line parsing
    content = read(path, String)
    content = replace(content, r"(?s)#=.*?=#" => "")
    lines = split(content, '\n')

    # Patterns for creation calls
    type_creators = r"^(?:\s*)(?:([A-Za-z_]\w*)\s*=\s*)?(?:(?:PETSc|LibPETSc)\.)?(KSP|SNES|DMDA|DMStag)\s*\("
    # Match common PETSc Vec/Mat creators, including suffixed forms like
    # VecCreateSeq, VecCreateMPI, VecDuplicateVecs, MatCreateSeqAIJ, etc.
    vec_create = r"^(?:\s*)(?:([A-Za-z_]\w*)\s*=\s*)?(?:(?:PETSc|LibPETSc)\.)?Vec(?:Create\w*|Duplicate\w*|Load\w*|\w*WithArray)\s*\("
    mat_create = r"^(?:\s*)(?:([A-Za-z_]\w*)\s*=\s*)?(?:(?:PETSc|LibPETSc)\.)?Mat(?:Create\w*|Duplicate\w*|Load\w*|\w*With\w*Arrays)\s*\("
    # Additional allocators that should be destroyed
    dm_vec_alloc = r"^(?:\s*)(?:([A-Za-z_]\w*)\s*=\s*)?(?:(?:PETSc|LibPETSc)\.)?(DMGlobalVec|DMLocalVec)\s*\("
    dm_alloc_dm = r"^(?:\s*)(?:([A-Za-z_]\w*)\s*=\s*)?(?:(?:PETSc|LibPETSc)\.)?(DMGetCoordinateDM|DMStagCreateCompatibleDMStag)\s*\("
    dm_alloc_mat = r"^(?:\s*)(?:([A-Za-z_]\w*)\s*=\s*)?(?:(?:PETSc|LibPETSc)\.)?(DMCreateMatrix)\s*\("

    # Some PETSc calls create multiple objects via tuple assignment
    # e.g. `x,b = LibPETSc.MatCreateVecs(petsclib, A)`
    mat_create_vecs = r"^(?:\s*)(.+?)\s*=\s*(?:(?:PETSc|LibPETSc)\.)?MatCreateVecs\s*\("

    # Pattern for destroy calls
    destroy_pat = r"^(?:\s*)(?:(?:PETSc)\.)?destroy\s*\(\s*([A-Za-z_]\w*)"
    # Pattern for finalize calls
    finalize_pat = r"^(?:\s*)(?:(?:PETSc)\.)?finalize\s*\(\s*([A-Za-z_]\w*)"

    created = Vector{Tuple{Int,Union{Nothing,String},String}}()
    destroyed = Vector{Tuple{Int,String}}()
    finalized = Vector{Tuple{Int,String}}()

    for (i, ln_raw) in enumerate(lines)
        # Skip empty and full-line comments
        s = strip(ln_raw)
        if isempty(s) || startswith(s, "#")
            continue
        end
        ln_proc = ln_raw

        if occursin(type_creators, ln_proc)
            m = match(type_creators, ln_proc)
            var = m.captures[1]
            func = m.captures[2]
            push!(created, (i, var, func))
            continue
        end
        if occursin(vec_create, ln_proc)
            m = match(vec_create, ln_proc)
            var = m.captures[1]
            push!(created, (i, var, "Vec"))
            continue
        end
        if occursin(mat_create, ln_proc)
            m = match(mat_create, ln_proc)
            var = m.captures[1]
            push!(created, (i, var, "Mat"))
            continue
        end
        if occursin(dm_vec_alloc, ln_proc)
            m = match(dm_vec_alloc, ln_proc)
            var = m.captures[1]
            push!(created, (i, var, "Vec"))
            continue
        end
        if occursin(dm_alloc_dm, ln_proc)
            m = match(dm_alloc_dm, ln_proc)
            var = m.captures[1]
            push!(created, (i, var, "DM"))
            continue
        end
        if occursin(dm_alloc_mat, ln_proc)
            m = match(dm_alloc_mat, ln_proc)
            var = m.captures[1]
            push!(created, (i, var, "Mat"))
            continue
        end

        if occursin(mat_create_vecs, ln_proc)
            m = match(mat_create_vecs, ln_proc)
            lhs = m.captures[1]
            # Strip parentheses and split on commas: "(x, b)" or "x,b" etc.
            lhs = replace(lhs, "(" => "", ")" => "")
            vars = [strip(v) for v in split(lhs, ',') if !isempty(strip(v))]
            if isempty(vars)
                push!(created, (i, nothing, "Vec"))
            else
                for v in vars
                    if occursin(r"^[A-Za-z_]\w*$", v)
                        push!(created, (i, v, "Vec"))
                    end
                end
            end
            continue
        end

        if occursin(destroy_pat, ln_proc)
            m = match(destroy_pat, ln_proc)
            var = m.captures[1]
            push!(destroyed, (i, var))
            continue
        end
        if occursin(finalize_pat, ln_proc)
            m = match(finalize_pat, ln_proc)
            var = m.captures[1]
            push!(finalized, (i, var))
            continue
        end
    end

    println("CREATION statements:")
    for (lnum, var, func) in created
        if isnothing(var)
            println("  line $(lnum): $(func)(...) (no assignment)")
        else
            println("  line $(lnum): $(var) = $(func)(...)")
        end
    end

    println("DESTROY calls:")
    for (lnum, var) in destroyed
        println("  line $(lnum): destroy($(var))")
    end

    created_vars = Set{String}()
    for (_, var, _) in created
        if !isnothing(var)
            push!(created_vars, var::String)
        end
    end
    destroyed_vars = Set{String}(map(x -> x[2], destroyed))

    missing = setdiff(created_vars, destroyed_vars)
    if !isempty(missing)
        println("POSSIBLY UNDESTROYED objects:")
        for v in sort(collect(missing))
            println("  $(v)")
        end
    else
        println("No obvious leaks: all assigned creations have a destroy call.")
    end

    if isempty(finalized)
        println("FINALIZE calls: none detected")
        println("Suggestion: add PETSc.finalize(petsclib) at the end of the routine.")
    else
        println("FINALIZE calls:")
        for (lnum, var) in finalized
            println("  line $(lnum): finalize($(var))")
        end
    end

    return nothing
end
