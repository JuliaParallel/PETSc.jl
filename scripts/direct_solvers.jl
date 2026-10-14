#!/usr/bin/env julia
#
# Which direct solvers in the PETSc_jll build give the right answer, in serial and in parallel.
#
#   julia --project=. scripts/direct_solvers.jl            # table of every case
#   julia --project=. scripts/direct_solvers.jl --check    # also fail on a case that differs from EXPECTED
#
# The solver processes use the active project, so `--project=env` tests the JLLs pinned in `env`.
#
# Every solver runs in a process of its own under `mpiexec` (1 and 2 ranks), so that one that
# crashes is reported as a crash instead of ending the run. Each process solves two systems, an
# SPD one and a saddle point, both with `ksp_type preonly` (the factorization alone) and with
# `ksp_type gmres` around it, and compares with Julia's own sparse LU. The table this prints is
# the one in docs/src/man/direct_solvers.md.

using MPI, PETSc, SparseArrays, LinearAlgebra, Printf

# The systems
# ----------------------------------------------------------------------------

# −∇·(η∇u) on an n×n grid, Dirichlet, η = contrast in a central square, scaled by 1/h²
function laplacian(n, contrast)
    η = ones(n, n)
    c, r = n ÷ 2, max(1, n ÷ 8)
    η[(c - r):(c + r), (c - r):(c + r)] .= contrast
    idx = LinearIndices((n, n))
    I, J, V = Int[], Int[], Float64[]
    for j in 1:n, i in 1:n
        k, d = idx[i, j], 0.0
        for (di, dj) in ((1, 0), (-1, 0), (0, 1), (0, -1))
            ii, jj = i + di, j + dj
            w = 2 / (1 / η[i, j] + 1 / η[clamp(ii, 1, n), clamp(jj, 1, n)])   # harmonic mean
            d += w
            (1 <= ii <= n && 1 <= jj <= n) && (push!(I, k); push!(J, idx[ii, jj]); push!(V, -w))
        end
        push!(I, k); push!(J, k); push!(V, d)
    end
    return sparse(I, J, V, n^2, n^2) * n^2
end

spd_system(n) = laplacian(n, 1e4)

# [A Bᵀ; B −εI]: two velocity components with a viscosity contrast of 1e6, a divergence B and a
# small pressure penalty ε, the shape of a penalized Stokes problem. The pressure block is the
# one that makes a factorization without row pivoting lose accuracy.
function saddle_system(n; contrast = 1e6, ε = 1e-12)
    A = blockdiag(laplacian(n, contrast), laplacian(n, contrast))
    Dx = spdiagm(n, n, 0 => -ones(n), 1 => ones(n - 1)) * n
    B = [kron(sparse(1.0I, n, n), Dx) kron(Dx, sparse(1.0I, n, n))]
    return [A B'; B -ε * sparse(1.0I, n^2, n^2)]
end

const PROBLEMS = (:spd, :saddle, :saddle_mild)
const MODES = (:preonly, :gmres)

# PETSc options per solver, and the problems each applies to
const SOLVERS = [
    ("petsc",        (pc_type = "lu", pc_factor_mat_solver_type = "petsc"),        PROBLEMS),
    ("mumps",        (pc_type = "lu", pc_factor_mat_solver_type = "mumps"),        PROBLEMS),
    ("superlu_dist", (pc_type = "lu", pc_factor_mat_solver_type = "superlu_dist"), PROBLEMS),
    ("umfpack",      (pc_type = "lu", pc_factor_mat_solver_type = "umfpack"),      PROBLEMS),
    ("klu",          (pc_type = "lu", pc_factor_mat_solver_type = "klu"),          PROBLEMS),
    ("cholmod",      (pc_type = "cholesky", pc_factor_mat_solver_type = "cholmod"), (:spd,)),
    ("hypre",        (pc_type = "hypre", pc_hypre_type = "boomeramg"),             (:spd,)),   # gmres only: one V-cycle is no solve
]

# One process: every problem and mode of one solver
# ----------------------------------------------------------------------------

# A distributed AIJ matrix holding K, each rank filling the rows it owns
function distributed(petsclib, comm, K)
    LibPETSc = PETSc.LibPETSc
    PetscInt = petsclib.PetscInt
    N = size(K, 1)
    nz = PetscInt(maximum(diff(SparseMatrixCSC(K').colptr)))
    A = LibPETSc.MatCreateAIJ(petsclib, comm, PetscInt(LibPETSc.PETSC_DECIDE), PetscInt(LibPETSc.PETSC_DECIDE),
                              PetscInt(N), PetscInt(N), nz, C_NULL, nz, C_NULL)
    rstart, rend = LibPETSc.MatGetOwnershipRange(petsclib, A)
    Kt = SparseMatrixCSC(K')                 # column j of Kt is row j of K
    for i in (rstart + 1):rend, p in nzrange(Kt, i)
        A[i, rowvals(Kt)[p]] = nonzeros(Kt)[p]
    end
    PETSc.assemble!(A)
    return A, (rstart + 1):rend
end

function run_case(name)
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
    _, options, problems = SOLVERS[findfirst(s -> s[1] == name, SOLVERS)]
    for problem in problems, mode in MODES
        name == "hypre" && mode === :preonly && continue
        K = problem === :spd ? spd_system(24) :
            problem === :saddle ? saddle_system(24) : saddle_system(32; contrast = 1e4, ε = 1e-6)
        b = ones(size(K, 1))
        reference = K \ b
        A, rows = distributed(petsclib, comm, K)
        ksp = PETSc.KSP(A; ksp_type = String(mode), ksp_rtol = 1e-10, ksp_max_it = 50, options...)
        x, bv = PETSc.LibPETSc.MatCreateVecs(petsclib, A)
        PETSc.with_local_array!(bv) do a
            a .= b[rows]
        end
        status, error = try
            PETSc.solve!(x, ksp, bv; check = true)
            local_x = PETSc.with_local_array!(copy, x; write = false)
            local_err = sum(abs2, local_x .- reference[rows])
            err = sqrt(MPI.Allreduce(local_err, +, comm)) / norm(reference)
            (err < 1e-8 ? "ok" : "wrong"), err
        catch e
            (e isa PETSc.ConvergenceError ? "diverged" : "error"), NaN
        end
        rank == 0 && println("RESULT $name $problem $mode $status $(@sprintf("%.1e", error))")
        foreach(PETSc.destroy!, (ksp, x, bv, A))
    end
    PETSc.finalize(petsclib)
    return nothing
end

# The driver
# ----------------------------------------------------------------------------

# What PETSc_jll 3.25.4 gives on Linux x86_64 with MPICH, once MUMPS_jll and MPICH_jll have the
# fix for parallel MUMPS (JuliaPackaging/Yggdrasil#15135; before it, MUMPS on 2 ranks hangs).
# `ok` matches Julia's LU to 1e-8, `wrong` reports success with a larger error, `diverged` or
# `error` reports the failure, `crash` ends the process and `hang` never finishes. `--check`
# compares against this.
const EXPECTED = let expected = Dict{Tuple{String, Symbol, Symbol, Int}, String}()
    serial_only = ("petsc", "umfpack", "klu", "cholmod")
    for (name, _, problems) in SOLVERS, problem in problems, mode in MODES, nranks in (1, 2)
        name == "hypre" && mode === :preonly && continue
        expected[(name, problem, mode, nranks)] = nranks == 2 && name in serial_only ? "error" : "ok"
    end
    # PETSc's own LU doesn't pivot: 1e-8 on the stiff saddle point, where the others give 1e-13
    expected[("petsc", :saddle, :preonly, 1)] = expected[("petsc", :saddle, :gmres, 1)] = "wrong"
    for nranks in (1, 2)
        expected[("superlu_dist", :saddle_mild, :preonly, nranks)] = "wrong"
        # Not checked: within a factor of 2 of the 1e-8 limit, or different from run to run
        delete!(expected, ("superlu_dist", :saddle, :preonly, nranks))
        delete!(expected, ("superlu_dist", :saddle, :gmres, nranks))
        delete!(expected, ("superlu_dist", :saddle_mild, :gmres, nranks))
    end
    expected
end

# Seconds one solver's process may take; it needs about 15 when nothing goes wrong
const TIMEOUT = 300

function run_all(; check::Bool)
    script = @__FILE__
    results = Dict{Tuple{String, Symbol, Symbol, Int}, Tuple{String, String}}()
    for nranks in (1, 2), (name, _, problems) in SOLVERS
        cmd = `$(MPI.mpiexec()) -n $nranks $(Base.julia_cmd()) --project=$(Base.active_project()) $script --case $name`
        out = IOBuffer()
        process = run(pipeline(ignorestatus(cmd); stdout = out, stderr = devnull); wait = false)
        hung = timedwait(() -> process_exited(process), TIMEOUT) === :timed_out
        hung && (kill(process); wait(process))
        seen = Set{Tuple{Symbol, Symbol}}()
        for line in eachline(IOBuffer(take!(out)))
            startswith(line, "RESULT ") || continue
            _, _, problem, mode, status, err = split(line)
            results[(name, Symbol(problem), Symbol(mode), nranks)] = (String(status), String(err))
            push!(seen, (Symbol(problem), Symbol(mode)))
        end
        for problem in problems, mode in MODES
            name == "hypre" && mode === :preonly && continue
            (problem, mode) in seen && continue
            results[(name, problem, mode, nranks)] = (hung ? "hang" : "crash", "")
        end
    end
    print_table(results)
    if check
        differs = [k for (k, v) in EXPECTED if get(results, k, ("missing", ""))[1] != v]
        isempty(differs) || error("cases that differ from EXPECTED: $(join(differs, ", "))")
    end
    return results
end

function print_table(results)
    println("| solver | ranks | SPD, preonly | SPD, gmres | saddle, preonly | saddle, gmres | mild saddle, preonly | mild saddle, gmres |")
    println("|---|---|---|---|---|---|---|---|")
    for (name, _, _) in SOLVERS, nranks in (1, 2)
        cells = String[]
        for problem in PROBLEMS, mode in MODES
            status, err = get(results, (name, problem, mode, nranks), ("n/a", ""))
            push!(cells, status in ("ok", "wrong") ? "$status ($err)" : status)
        end
        println("| $name | $nranks | ", join(cells, " | "), " |")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    i = findfirst(==("--case"), ARGS)
    if i === nothing
        run_all(; check = "--check" in ARGS)
    else
        run_case(ARGS[i + 1])
    end
end
