# ============================================================================
#   Operations on stencils
# ============================================================================

"""
    set_values!(J::AbstractPetscMat, dm::DMStag, rows, cols, vals; mode = INSERT_VALUES)

Write the dense block `vals` into `J` at the stencils `rows` × `cols` and return
`J`. `vals` is row-major, with `length(rows) * length(cols)` entries. Under
`ADD_VALUES`, entries whose row and column repeat are summed. `rows`, `cols` and
`vals` can be any `AbstractVector`. A `Vector` or a prefix view of one,
`view(buf, 1:n)`, is passed without a copy, so a scratch buffer can be reused for
blocks of any size; any other vector is copied first.

# External Links
$(doc_external("DMStag/DMStagMatSetValuesStencil"))
"""
function set_values!(
    J::AbstractPetscMat{PetscLib},
    dm::DMStag{PetscLib},
    rows::AbstractVector{LibPETSc.DMStagStencil},
    cols::AbstractVector{LibPETSc.DMStagStencil},
    vals::AbstractVector;
    mode::InsertMode = INSERT_VALUES,
) where {PetscLib}
    length(vals) == length(rows) * length(cols) || throw(
        DimensionMismatch(
            "a $(length(rows))x$(length(cols)) block needs " *
            "$(length(rows) * length(cols)) values, got $(length(vals))",
        ),
    )
    T = inttype(PetscLib)
    LibPETSc.DMStagMatSetValuesStencil(
        getlib(PetscLib), dm, J,
        T(length(rows)), stencil_vector(rows),
        T(length(cols)), stencil_vector(cols),
        value_vector(PetscLib, vals), mode,
    )
    return J
end

"""
    set_values!(v::AbstractPetscVec, dm::DMStag, positions, vals; mode = INSERT_VALUES)

Write `vals` into `v` at the stencils `positions` and return `v`. Both can be any
`AbstractVector` of the same length. A `Vector` or a prefix view of one,
`view(buf, 1:n)`, is passed without a copy; any other vector is copied first.

# External Links
$(doc_external("DMStag/DMStagVecSetValuesStencil"))
"""
function set_values!(
    v::AbstractPetscVec{PetscLib},
    dm::DMStag{PetscLib},
    positions::AbstractVector{LibPETSc.DMStagStencil},
    vals::AbstractVector;
    mode::InsertMode = INSERT_VALUES,
) where {PetscLib}
    length(vals) == length(positions) || throw(
        DimensionMismatch("$(length(positions)) positions, but $(length(vals)) values"),
    )
    LibPETSc.DMStagVecSetValuesStencil(
        getlib(PetscLib), dm, v, inttype(PetscLib)(length(positions)),
        stencil_vector(positions), value_vector(PetscLib, vals), mode,
    )
    return v
end

"""
    zero_rows_local!(J::AbstractPetscMat, dm::DMStag, rows; diag = 1, x = nothing, b = nothing)

[`zero_rows_local!`](@ref) with the rows given as stencils of `dm`, which `J` was
created from. Collective: a rank that owns none of the rows passes an empty vector.

# External Links
$(doc_external("DMStag/DMStagStencilToIndexLocal"))
"""
function zero_rows_local!(
    J::AbstractPetscMat{PetscLib},
    dm::DMStag{PetscLib, N},
    rows::AbstractVector{LibPETSc.DMStagStencil};
    diag = 1,
    x = nothing,
    b = nothing,
) where {PetscLib, N}
    T = inttype(PetscLib)
    rows_0b = LibPETSc.DMStagStencilToIndexLocal(
        getlib(PetscLib), dm, T(N), T(length(rows)), stencil_vector(rows),
    )
    return zero_rows_local!(J, rows_0b; diag, x, b)
end

# The C calls take a `Vector` and a count, so a `Vector` passes as it is, and so does
# the parent of a prefix view `view(buf, 1:n)`, the count being `n`. Any other vector
# is copied into one.
const PrefixView{T} = SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int}}, true}

stencil_vector(v::Vector{LibPETSc.DMStagStencil}) = v
stencil_vector(v::PrefixView{LibPETSc.DMStagStencil}) =
    first(parentindices(v)[1]) == 1 ? parent(v) : collect(v)
stencil_vector(v::AbstractVector{LibPETSc.DMStagStencil}) = collect(v)

value_vector(::Type{PetscLib}, v::Vector) where {PetscLib} =
    eltype(v) === scalartype(PetscLib) ? v : Vector{scalartype(PetscLib)}(v)
value_vector(::Type{PetscLib}, v::PrefixView) where {PetscLib} =
    eltype(v) === scalartype(PetscLib) && first(parentindices(v)[1]) == 1 ? parent(v) :
    Vector{scalartype(PetscLib)}(v)
value_vector(::Type{PetscLib}, v::AbstractVector) where {PetscLib} =
    Vector{scalartype(PetscLib)}(v)
