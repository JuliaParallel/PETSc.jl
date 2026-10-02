# Performance

PETSc.jl aims for a fixed cost per call: the work a high-level call does on the Julia side does not grow with the size of the vector, matrix or grid it acts on. `test/core/performance.jl` holds that as a contract, checking every call below at two problem sizes. If a call in this page starts allocating more at the larger size, that is a bug worth reporting.

## Read and write arrays in a block

[`with_local_array!`](@ref) and [`with_field_views!`](@ref) check the arrays out once, hand them to your function and give them back, so the cost is a few allocations per vector no matter how long it is. Entry-by-entry access through `v[i]` goes through `VecGetValues`, one C call per entry, and only reaches the entries this rank owns:

```julia
# one checkout, then plain Julia on the array
PETSc.with_local_array!(x, y; write = (true, false)) do a, b
    @inbounds for i in eachindex(a)
        a[i] = 2b[i]
    end
end

# one C call per entry, and wrong on several ranks
for i in 1:length(x)
    x[i] = 2y[i]
end
```

Broadcasting does the checkout for you, in both directions. `x .= 2 .* y` writes the entries this rank owns. `w = 2 .* x` gives a plain `Vector`, and since that only means something where one rank holds the whole vector, it throws on a distributed one: broadcast into a `PetscVec` of the same layout instead.

Every view a block hands you is indexed the way the DM numbers its points, so a view and [`stencil`](@ref) address the same entry, and the arrays stay concretely typed inside the block.

## BLAS threads

Every `PETSc_jll` build calls BLAS through libblastrampoline, which means PETSc's dense kernels run in Julia's own BLAS pool. PETSc cannot size that pool, so `-mat_mumps_*` and any other dense work inherits whatever Julia set.

When several MPI ranks share a node, [`initialize`](@ref) sets the pool to one thread per rank, unless `-blas_num_threads`, `OPENBLAS_NUM_THREADS` or `OMP_NUM_THREADS` says otherwise. Without that, each rank runs a pool of busy-waiting threads competing for the same cores: a 4-rank DMStag Stokes solve measured 25 times slower. A serial run, or one rank per node, keeps Julia's default.

To choose for yourself:

```julia
PETSc.initialize(petsclib; options = ["-blas_num_threads", "2"])
```

## Keep PETSc objects concretely typed

A wrapper carries its library in its type, so a high-level call resolves the right method at compile time. Storing objects in a field typed `Any`, or in a `Ref{Any}`, throws that away and every call through it becomes a dynamic dispatch:

```julia
mutable struct Ctx
    dm::Any            # every call on ctx.dm dispatches at run time
end

mutable struct Ctx{D}
    dm::D              # resolved at compile time
end
```

Where a context has to stay untyped, put a function barrier between it and the loop: read the fields once into local variables and pass them to a second function that does the work.

## Precompile before an MPI launch

Precompile serially, once, before running on several ranks:

```bash
julia --project=. -e 'using PETSc'
mpiexec -n 4 julia --project=. run.jl
```

Ranks that start against a stale cache take turns: one precompiles while the others wait on its lock, so every rank pays for the compile before the run starts.

## Profiling a solve

`-log_view` gives PETSc's own breakdown, which attributes time to `KSPSolve`, `MatMult`, `PCApply` and the rest, and counts the calls:

```julia
PETSc.initialize(petsclib; log_view = true)
```

It is the right first measurement for a slow solve, because it separates time spent in PETSc from time spent in your callbacks. Julia's `@profile` then covers the callbacks themselves.
