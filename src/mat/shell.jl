struct MatOp{PetscLib, Op} end

function (::MatOp{PetscLib, LibPETSc.MATOP_MULT})(
            M::CMat,
            cx::CVec,
            cy::CVec,
        ) where {PetscLib}
    return run_callback("MatShell multiply") do
        state = unsafe_pointer_to_objref(LibPETSc.MatShellGetContext(PetscLib, M))::MatShellState
        petsclib = getlib(PetscLib)
        _mul!(PetscVec(cy, petsclib; own = false), state.obj, PetscVec(cx, petsclib; own = false))
    end
end

function (::MatOp{PetscLib, LibPETSc.MATOP_MULT_TRANSPOSE})(
            M::CMat,
            cx::CVec,
            cy::CVec,
        ) where {PetscLib}
    return run_callback("MatShell mult_transpose") do
        state = unsafe_pointer_to_objref(LibPETSc.MatShellGetContext(PetscLib, M))::MatShellState
        petsclib = getlib(PetscLib)
        state.mult_transpose(PetscVec(cy, petsclib; own = false), PetscVec(cx, petsclib; own = false))
    end
end

function (::MatOp{PetscLib, LibPETSc.MATOP_GET_DIAGONAL})(
            M::CMat,
            cd::CVec,
        ) where {PetscLib}
    return run_callback("MatShell diagonal") do
        state = unsafe_pointer_to_objref(LibPETSc.MatShellGetContext(PetscLib, M))::MatShellState
        state.diagonal(PetscVec(cd, getlib(PetscLib); own = false))
    end
end

# MatMultAdd(A, x, y, z) is z = y + Ax; z may be y
function (::MatOp{PetscLib, LibPETSc.MATOP_MULT_ADD})(
            M::CMat,
            cx::CVec,
            cy::CVec,
            cz::CVec,
        ) where {PetscLib}
    return run_callback("MatShell mult_add") do
        state = unsafe_pointer_to_objref(LibPETSc.MatShellGetContext(PetscLib, M))::MatShellState
        petsclib = getlib(PetscLib)
        state.mult_add(
            PetscVec(cz, petsclib; own = false),
            PetscVec(cx, petsclib; own = false),
            PetscVec(cy, petsclib; own = false),
        )
    end
end




# The macro is needed because of the @cfunction


"""
    MatShell(obj, petsclib, comm::MPI.Comm, local_rows, local_cols;
             global_rows = PETSC_DECIDE, global_cols = PETSC_DECIDE,
             mult_transpose = nothing, diagonal = nothing, mult_add = nothing)

Create a `global_rows X global_cols` PETSc shell matrix object wrapping `obj`
with local size `local_rows X local_cols`.

`obj` is the multiply: a `Function` is called as `obj(y, x)` to fill `y`,
anything else goes through `mul!(y, obj, x)`. It comes first so `do` works:

```julia
A = MatShell(petsclib, comm, n, n) do y, x
    y .= 2 .* x
end
```

The keywords give the shell more operations, each a function that fills its
first argument, as `obj` does:

- `mult_transpose(y, x)`: `y = Aᵀx`. Then `transpose(A) * x` and `A' * x` work,
  the second through the conjugate of `x` in a complex build.
- `diagonal(d)`: `d` gets the diagonal of `A`. Then `pc_type = jacobi` works.
- `mult_add(z, x, y)`: `z = y + Ax`, where `z` may be `y`. Without it PETSc
  multiplies and adds. PETSc calls it as given, so a `MatScale` or `MatShift`
  on the shell does not reach it.

An operation left out throws a `PetscError` when PETSc needs it.

if `comm == MPI.COMM_SELF` then the garbage connector can finalize the object,
otherwise the user is responsible for calling [`destroy!`](@ref).

$(doc_callback())

# External Links
$(doc_external("Mat/MatCreateShell"))
$(doc_external("Mat/MatShellSetOperation"))
$(doc_external("Mat/MATOP_MULT"))
$(doc_external("Mat/MatMultTranspose"))
$(doc_external("Mat/MatGetDiagonal"))
$(doc_external("Mat/MatMultAdd"))
"""
mutable struct MatShell{PetscLib, OType} <: AbstractPetscMat{PetscLib}
    ptr::CMat
    obj::OType
    age
    own::Bool
end

# The Julia side of a MatShell, kept with the PETSc object.
# It holds the operators, not the wrapper, so the wrapper can still be collected
# and its finalizer run.
mutable struct MatShellState <: ObjectState
    obj::Any
    mult_transpose::Any
    diagonal::Any
    mult_add::Any
    alive::Bool
end
MatShellState() = MatShellState(nothing, nothing, nothing, nothing, true)
state_type(::Type{<:MatShell}) = MatShellState

function MatShell(
    obj,
    petsclib::PetscLibType,
    comm::MPI.Comm,
    local_rows::Integer,
    local_cols::Integer;
    global_rows::Integer = LibPETSc.PETSC_DECIDE,
    global_cols::Integer = LibPETSc.PETSC_DECIDE,
    mult_transpose = nothing,
    diagonal = nothing,
    mult_add = nothing,
)
    return create_matshell(petsclib, obj, comm, local_rows, local_cols, global_rows, global_cols,
                           mult_transpose, diagonal, mult_add)
end

# The per-library body of `MatShell`. The public method above takes `petsclib`
# as the abstract type, so the old argument order's shim can be told apart
# from it without an ambiguity.
LibPETSc.@for_petsc function create_matshell(
    petsclib::$PetscLib,
    obj::OType,
    comm::MPI.Comm,
    local_rows,
    local_cols,
    global_rows,
    global_cols,
    mult_transpose,
    diagonal,
    mult_add,
) where {OType}
    mat = MatShell{$PetscLib, OType}(C_NULL, obj, petsclib.age, true)

#=
    ccall(
        (:MatCreateShell, $petsc_library),
        LibPETSc.PetscErrorCode,
        (
            LibPETSc.MPI_Comm,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            $PetscInt,
            Ptr{Cvoid},
            Ptr{CMat},
        ),
        comm,
        local_rows,
        local_cols,
        global_rows,
        global_cols,
        pointer_from_objref(mat),
        A_,
    )
=#
    A_ = Ref{CMat}()
    ccall(
               (:MatCreateShell, $petsc_library),
               LibPETSc.PetscErrorCode,
               (LibPETSc.MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{Cvoid}, Ptr{CMat}),
               comm, local_rows, local_cols, global_rows, global_cols, C_NULL, A_,
              )

    mat.ptr = A_[]
    state = object_state!(mat)
    state.obj = obj
    LibPETSc.MatShellSetContext(petsclib, mat, pointer_from_objref(state))

    #=
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


    mat = LibPETSc.MatCreateShell(
        petsclib,
        comm,
        local_rows,
        local_cols,
        global_rows,
        global_cols,
        pointer_from_objref(mat),
    )
    =#
  
    mulptr = @cfunction(
        MatOp{$PetscLib, LibPETSc.MATOP_MULT}(),
        LibPETSc.PetscErrorCode,
        (CMat, CVec, CVec)
    )

    LibPETSc.MatShellSetOperation(petsclib, mat, LibPETSc.MATOP_MULT, mulptr)

    # The optional operations, each registered only when given
    state.mult_transpose = mult_transpose
    state.diagonal = diagonal
    state.mult_add = mult_add
    if !isnothing(mult_transpose)
        ptr = @cfunction(MatOp{$PetscLib, LibPETSc.MATOP_MULT_TRANSPOSE}(),
                         LibPETSc.PetscErrorCode, (CMat, CVec, CVec))
        LibPETSc.MatShellSetOperation(petsclib, mat, LibPETSc.MATOP_MULT_TRANSPOSE, ptr)
    end
    if !isnothing(diagonal)
        ptr = @cfunction(MatOp{$PetscLib, LibPETSc.MATOP_GET_DIAGONAL}(),
                         LibPETSc.PetscErrorCode, (CMat, CVec))
        LibPETSc.MatShellSetOperation(petsclib, mat, LibPETSc.MATOP_GET_DIAGONAL, ptr)
    end
    if !isnothing(mult_add)
        ptr = @cfunction(MatOp{$PetscLib, LibPETSc.MATOP_MULT_ADD}(),
                         LibPETSc.PetscErrorCode, (CMat, CVec, CVec, CVec))
        LibPETSc.MatShellSetOperation(petsclib, mat, LibPETSc.MATOP_MULT_ADD, ptr)
    end

    # The garbage collector can only destroy it when no other rank takes part
    if MPI.Comm_size(comm) == 1
        finalizer(destroy!, mat)
    end

    return mat
end

# The operator of a MatShell applied to `x`: a function is called as `obj(y, x)`,
# anything else goes through `mul!(y, obj, x)`
_mul!(y, obj::Function, x) = obj(y, x)
_mul!(y, obj, x) = LinearAlgebra.mul!(y, obj, x)


# Matrix-vector multiplication for MatShell with Julia arrays
function Base.:*(M::MatShell{PetscLib}, x::AbstractVector) where {PetscLib}
    PetscScalar = scalartype(PetscLib)
    PetscInt = inttype(PetscLib)
    
    # Get matrix dimensions
    m, n = size(M)
    
    # Create PETSc vectors wrapping the Julia arrays
    petsclib = getlib(PetscLib)
    # PETSc works on this copy of `x` in place, so it must outlive the product
    x_copy = PetscScalar.(x)
    return GC.@preserve x_copy begin
        petsc_x = seq_vec_with_array(petsclib, MPI.COMM_SELF, PetscInt(1), PetscInt(n), x_copy)
        petsc_y = LibPETSc.VecCreateSeq(petsclib, MPI.COMM_SELF, PetscInt(m))
        mul!(petsc_y, M, petsc_x)
        result = petsc_y[:]
        destroy!(petsc_x)
        destroy!(petsc_y)
        result
    end
end
