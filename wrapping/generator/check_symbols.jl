# Every wrapped function must be a symbol of the PETSc library, otherwise its ccall fails at run
# time with "could not load symbol". Functions that are `static inline` in a PETSc header or C
# macros show up in getAPI.py's output but have no symbol: put them in `[exclude]` of
# rules/files.toml. Functions of optional packages (CUDA, Kokkos, MOAB, ...) are also reported;
# they stay wrapped because other PETSc_jll builds may have them.
#
#   julia --project=. wrapping/generator/check_symbols.jl [autowrapped-dir]
#
# Runs in the package environment (needs PETSc_jll), not in wrapping/generator's.
using PETSc, MPI, Libdl

dir = isempty(ARGS) ? joinpath(@__DIR__, "..", "..", "src", "autowrapped") : ARGS[1]
petsclib = PETSc.getlib()
h = Libdl.dlopen(petsclib.petsc_library)
missing = String[]
for f in sort!(readdir(dir; join = true))
    endswith(f, "_wrappers.jl") && basename(f) != "extra_wrappers.jl" || continue
    for m in eachmatch(r"^function (\w+)\(petsclib::PetscLibType"m, read(f, String))
        Libdl.dlsym_e(h, m.captures[1]) == C_NULL && push!(missing, "$(basename(f)): $(m.captures[1])")
    end
end
println(length(missing), " wrapped functions have no symbol in ", basename(string(petsclib.petsc_library)))
foreach(println, missing)
