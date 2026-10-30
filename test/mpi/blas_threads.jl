# test/mpi/blas_threads.jl
# Run on several ranks of one node: `initialize` gives each rank one BLAS
# thread, unless `-blas_num_threads` or the environment sizes the pool.

using Test
using MPI
MPI.Initialized() || MPI.Init()
using PETSc
using LinearAlgebra: BLAS

@testset "BLAS threads under MPI" begin
    petsclib = PETSc.petsclibs[1]
    PETSc.isinitialized(petsclib) && PETSc.finalize(petsclib)
    environment_sets_it = any(v -> haskey(ENV, v), PETSc.BLAS_THREAD_VARIABLES)
    node = MPI.Comm_split_type(MPI.COMM_WORLD, MPI.COMM_TYPE_SHARED, MPI.Comm_rank(MPI.COMM_WORLD))
    shared = MPI.Comm_size(node) > 1
    MPI.free(node)

    BLAS.set_num_threads(4)
    PETSc.initialize(petsclib)
    @test BLAS.get_num_threads() == (shared && !environment_sets_it ? 1 : 4)
    # the user can still resize the pool after initialize
    BLAS.set_num_threads(3)
    @test BLAS.get_num_threads() == 3
    PETSc.finalize(petsclib)

    # the option wins over the default
    PETSc.initialize(petsclib; options = ["-blas_num_threads", "2"])
    @test BLAS.get_num_threads() == 2
    PETSc.finalize(petsclib)

    # `initialize` is collective, so ranks with different environments must not block
    # each other: only rank 0 has OMP_NUM_THREADS set here
    saved = Dict(v => get(ENV, v, nothing) for v in PETSc.BLAS_THREAD_VARIABLES)
    try
        foreach(v -> delete!(ENV, v), PETSc.BLAS_THREAD_VARIABLES)
        MPI.Comm_rank(MPI.COMM_WORLD) == 0 && (ENV["OMP_NUM_THREADS"] = "2")
        BLAS.set_num_threads(4)
        PETSc.initialize(petsclib)
        @test PETSc.isinitialized(petsclib)
        mine_sets_it = haskey(ENV, "OMP_NUM_THREADS")
        @test BLAS.get_num_threads() == (shared && !mine_sets_it ? 1 : 4)
        PETSc.finalize(petsclib)
    finally
        for (v, value) in saved
            value === nothing ? delete!(ENV, v) : (ENV[v] = value)
        end
    end
end
