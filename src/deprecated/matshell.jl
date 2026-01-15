"""
    MatShell(
        petsclib::PetscLib,
        obj::OType,
        comm::MPI.Comm,
        local_rows,
        local_cols,
        global_rows = LibPETSc.PETSC_DECIDE,
        global_cols = LibPETSc.PETSC_DECIDE,
    )

Create a `global_rows X global_cols` PETSc shell matrix object wrapping `obj`
with local size `local_rows X local_cols`.

The `obj` will be registered as an `MATOP_MULT` function and if if `obj` is a
`Function`, then the multiply action `obj(y,x)`; otherwise it calls `mul!(y,
obj, x)`.

if `comm == MPI.COMM_SELF` then the garbage connector can finalize the object,
otherwise the user is responsible for calling [`destroy`](@ref).

# External Links
$(_doc_external("Mat/MatCreateShell"))
$(_doc_external("Mat/MatShellSetOperation"))
$(_doc_external("Mat/MATOP_MULT"))
"""
mutable struct MatShell{PetscLib, PetscScalar, OType} <:
               AbstractMat{PetscLib, PetscScalar}
    ptr::CMat
    obj::OType
    age::Int
end

struct MatOp{PetscLib, PetscInt, Op} end

function (::MatOp{PetscLib, PetscInt, LibPETSc.MATOP_MULT})(
    M::CMat,
    cx::CVec,
    cy::CVec,
)::PetscInt where {PetscLib, PetscInt}
    r_ctx = Ref{Ptr{Cvoid}}()
    LibPETSc.MatShellGetContext(PetscLib, M, r_ctx)
    ptr = r_ctx[]
    mat = unsafe_pointer_to_objref(ptr)

    PetscScalar = PetscLib.PetscScalar
    x = VecPtr(PetscLib, cx, false)
    y = VecPtr(PetscLib, cy, false)

    _mul!(y, mat, x)

    return PetscInt(0)
end

function _mul!(
    y,
    mat::MatShell{PetscLib, PetscScalar, F},
    x,
) where {PetscLib, PetscScalar, F <: Function}
    mat.obj(y, x)
end

function _mul!(y, mat::MatShell, x)
    LinearAlgebra.mul!(y, mat.obj, x)
end

# We have to use the macro here because of the @cfunction
LibPETSc.@for_petsc function MatShell(
    petsclib::$PetscLib,
    obj::OType,
    comm::MPI.Comm,
    local_rows,
    local_cols,
    global_rows = LibPETSc.PETSC_DECIDE,
    global_cols = LibPETSc.PETSC_DECIDE,
) where {OType}
    mat = MatShell{$PetscLib, $PetscScalar, OType}(C_NULL, obj, petsclib.age)

    # we use the MatShell object itself
    ctx = pointer_from_objref(mat)

    LibPETSc.MatCreateShell(
        petsclib,
        comm,
        local_rows,
        local_cols,
        global_rows,
        global_cols,
        pointer_from_objref(mat),
        mat,
    )

    mulptr = @cfunction(
        MatOp{$PetscLib, $PetscInt, LibPETSc.MATOP_MULT}(),
        $PetscInt,
        (CMat, CVec, CVec)
    )
    LibPETSc.MatShellSetOperation(petsclib, mat, LibPETSc.MATOP_MULT, mulptr)

    if MPI.Comm_size(comm) == 1
        finalizer(destroy, mat)
    end

    return mat
end
