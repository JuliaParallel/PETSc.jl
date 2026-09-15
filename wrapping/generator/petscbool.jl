primitive type PetscBool 32 end

const PETSC_FALSE = Base.bitcast(PetscBool, Int32(0))
const PETSC_TRUE = Base.bitcast(PetscBool, Int32(1))

@inline _petscbool_bits(x::PetscBool) = Base.bitcast(Int32, x)

PetscBool(x::Bool) = convert(PetscBool, x)

Base.convert(::Type{PetscBool}, x::Bool) = x ? PETSC_TRUE : PETSC_FALSE
Base.convert(::Type{PetscBool}, x::PetscBool) = x
Base.convert(::Type{PetscBool}, x::Integer) = x == 0 ? PETSC_FALSE : PETSC_TRUE
Base.Bool(x::PetscBool) = _petscbool_bits(x) != 0
Base.convert(::Type{Bool}, x::PetscBool) = Base.Bool(x)
Base.convert(::Type{Int32}, x::PetscBool) = _petscbool_bits(x)
Base.cconvert(::Type{PetscBool}, x::Bool) = convert(PetscBool, x)
Base.cconvert(::Type{PetscBool}, x::PetscBool) = x
Base.unsafe_convert(::Type{PetscBool}, x::PetscBool) = x
Base.promote_rule(::Type{PetscBool}, ::Type{Bool}) = Bool
Base.promote_rule(::Type{Bool}, ::Type{PetscBool}) = Bool
Base.getindex(r::Base.RefValue{PetscBool}) = Bool(getfield(r, :x))

Base.:(==)(x::PetscBool, y::PetscBool) = _petscbool_bits(x) == _petscbool_bits(y)
Base.:(==)(x::PetscBool, y::Bool) = Bool(x) == y
Base.:(==)(x::Bool, y::PetscBool) = x == Bool(y)
Base.:(!)(x::PetscBool) = !Bool(x)
Base.hash(x::PetscBool, h::UInt) = hash(Bool(x), h)
Base.show(io::IO, x::PetscBool) = show(io, Bool(x))
Base.zero(::Type{PetscBool}) = PETSC_FALSE
Base.one(::Type{PetscBool}) = PETSC_TRUE
Base.iszero(x::PetscBool) = !Bool(x)
