# Fixture for test_audit.jl. Never executed, only parsed.
#
# Every construct here has a known answer. Only `never_freed` and `quoted_v`
# leak; everything else is created and released, or is not a tracked creation.

using PETSc

# released inside a loop
for i in 1:3
    loop_v = LibPETSc.VecCreateSeq(petsclib, comm, 10)
    LibPETSc.VecDestroy(petsclib, loop_v)
end

# released in a finally block
function in_try()
    try_v = LibPETSc.VecCreateSeq(petsclib, comm, 10)
    try
        use(try_v)
    finally
        PETSc.destroy!(try_v)
    end
end

# created on both branches, released once
if cond
    branch_v = LibPETSc.VecCreateSeq(petsclib, comm, 10)
else
    branch_v = LibPETSc.VecCreateSeq(petsclib, comm, 20)
end
PETSc.destroy!(branch_v)

# creation wrapped in a macro
@time macro_v = LibPETSc.VecCreateSeq(petsclib, comm, 10)
PETSc.destroy!(macro_v)

# released by broadcasting over a literal
bcast_a = LibPETSc.VecCreateSeq(petsclib, comm, 10)
bcast_b = LibPETSc.VecCreateSeq(petsclib, comm, 10)
PETSc.destroy!.([bcast_a, bcast_b])

# released through a splat
splat_v = LibPETSc.VecCreateSeq(petsclib, comm, 10)
PETSc.destroy!(splat_v...)

# release call carrying a keyword argument
kwarg_v = LibPETSc.VecCreateSeq(petsclib, comm, 10)
PETSc.destroy!(kwarg_v; force = true)

# handed to the garbage collector, in both finalizer spellings
fin_v = LibPETSc.VecCreateSeq(petsclib, comm, 10)
finalizer(destroy, fin_v)
fin_m = PETSc.MatSeqAIJ(petsclib, 10, 10, 3)
finalizer(m -> (destroy(m); data), fin_m)

# borrowed pointer: the caller does not own it
borrowed = PETSc.VecPtr(petsclib, some_ptr, false)

# let block, global binding, and a non-ASCII name
let
    let_v = LibPETSc.VecCreateSeq(petsclib, comm, 10)
    PETSc.destroy!(let_v)
end
global global_v = LibPETSc.VecCreateSeq(petsclib, comm, 10)
PETSc.destroy!(global_v)
Δv = LibPETSc.VecCreateSeq(petsclib, comm, 10)
PETSc.destroy!(Δv)

#=
block_v = LibPETSc.VecCreateSeq(petsclib, comm, 10)
=#

# the release below is quoted, so it releases nothing here
macro cleanup()
    return quote
        PETSc.destroy!(quoted_v)
    end
end
quoted_v = LibPETSc.VecCreateSeq(petsclib, comm, 10)

# the control: nothing releases this
never_freed = LibPETSc.VecCreateSeq(petsclib, comm, 10)
