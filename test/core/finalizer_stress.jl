# test/core/finalizer_stress.jl
# Run by finalizer_queue.jl in a subprocess with several threads and
# PETSC_OPTIONS="-malloc_debug -malloc_test". One task makes every PETSc call: it
# drops serial Vecs, whose finalizers run on whichever thread collects, and solves.
# The other task makes no PETSc call; it only allocates and collects. Arguments:
# the number of seconds to run.

using PETSc, MPI, SparseArrays, LinearAlgebra
MPI.Init()
petsclib = PETSc.petsclibs[1]
PETSc.initialize(petsclib)
n = 2000
seconds = parse(Float64, ARGS[1])
A = spdiagm(-1 => fill(-1.0, n - 1), 0 => fill(2.0, n), 1 => fill(-1.0, n - 1))
ksp = PETSc.KSP(petsclib, MPI.COMM_SELF, A; ksp_type = "cg", pc_type = "jacobi", ksp_rtol = 1e-8)
b = PETSc.PetscVec(petsclib, ones(n))
x = PETSc.PetscVec(petsclib, zeros(n))
M = PETSc.PetscMat(petsclib, A)
stop = Threads.Atomic{Bool}(false)
solves = Threads.Atomic{Int}(0)
petsc_task = Threads.@spawn begin
    while !stop[]
        for _ in 1:50
            PETSc.PetscVec(petsclib, rand(50))
        end
        M * b
        PETSc.solve!(x, ksp, b)
        Threads.atomic_add!(solves, 1)
        yield()
    end
end
gc_task = Threads.@spawn begin
    while !stop[]
        junk = [rand(1000) for _ in 1:100]
        GC.gc(false)
        yield()
    end
end
sleep(seconds)
stop[] = true
wait(petsc_task)
wait(gc_task)
residual = norm(A * x[:] - ones(n)) / sqrt(n)
println("solves ", solves[], " residual ", residual)
residual < 1e-10 || error("wrong solution")
foreach(PETSc.destroy!, (ksp, b, x, M))
PETSc.finalize(petsclib)
