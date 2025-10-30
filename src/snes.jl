import .LibPETSc: AbstractPetscSNES, CSNES, PetscSNES

# Custom display for REPL
function Base.show(io::IO, v::AbstractPetscSNES{PetscLib}) where {PetscLib}
    if v.ptr == C_NULL
        print(io, "PETSc SNES (null pointer)")
        return
    else
        print(io, "PETSc SNES object")
    end
    return nothing
end