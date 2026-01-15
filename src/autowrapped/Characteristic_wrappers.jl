# autodefined type arguments for class ------
mutable struct _n_Characteristic end
const Characteristic = Ptr{_n_Characteristic}

# -------------------------------------------------------
"""
	CharacteristicFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the `Characteristics` package. It is
called from `PetscFinalize()`.

Level: developer

-seealso: [](ch_ts), `PetscFinalize()`, `CharacteristicInitializePackage()`

# External Links
$(_doc_external("Ts/CharacteristicFinalizePackage"))
"""
function CharacteristicFinalizePackage(petsclib::PetscLibType) end

@for_petsc function CharacteristicFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:CharacteristicFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	CharacteristicInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the Characteristic package. It is called
from PetscDLLibraryRegister() when using dynamic libraries, and on the first call to CharacteristicCreate()
when using static libraries.

Level: developer

-seealso: [](ch_ts), `PetscInitialize()`, `CharacteristicFinalizePackage()`

# External Links
$(_doc_external("Ts/CharacteristicInitializePackage"))
"""
function CharacteristicInitializePackage(petsclib::PetscLibType) end

@for_petsc function CharacteristicInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:CharacteristicInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	CharacteristicDestroy(petsclib::PetscLibType,c::Characteristic) 
Destroys a `Characteristic` context created with `CharacteristicCreate()`

Collective

Input Parameter:
- `c` - the `Characteristic` context

Level: beginner

-seealso: `Characteristic`, `CharacteristicCreate()`

# External Links
$(_doc_external("Ts/CharacteristicDestroy"))
"""
function CharacteristicDestroy(petsclib::PetscLibType, c::Characteristic) end

@for_petsc function CharacteristicDestroy(petsclib::$UnionPetscLib, c::Characteristic )

    @chk ccall(
               (:CharacteristicDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{Characteristic},),
               c,
              )


	return nothing
end 

"""
	c::Characteristic = CharacteristicCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Creates a `Characteristic` context for use with the Method of Characteristics

Collective

Input Parameter:
- `comm` - MPI communicator

Output Parameter:
- `c` - the `Characteristic` context

Level: beginner

-seealso: `Characteristic`, `CharacteristicDestroy()`

# External Links
$(_doc_external("Ts/CharacteristicCreate"))
"""
function CharacteristicCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function CharacteristicCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	c_ = Ref{Characteristic}()

    @chk ccall(
               (:CharacteristicCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Characteristic}),
               comm, c_,
              )

	c = c_[]

	return c
end 

"""
	CharacteristicSetType(petsclib::PetscLibType,c::Characteristic, type::CharacteristicType) 
Builds Characteristic for a particular solver.

Logically Collective

Input Parameters:
- `c`    - the method of characteristics context
- `type` - a known method

Options Database Key:
- `-characteristic_type <method>` - Sets the method; use -help for a list
of available methods

Level: intermediate

-seealso: [](ch_ts), `CharacteristicType`

# External Links
$(_doc_external("Ts/CharacteristicSetType"))
"""
function CharacteristicSetType(petsclib::PetscLibType, c::Characteristic, type::CharacteristicType) end

@for_petsc function CharacteristicSetType(petsclib::$UnionPetscLib, c::Characteristic, type::CharacteristicType )

    @chk ccall(
               (:CharacteristicSetType, $petsc_library),
               PetscErrorCode,
               (Characteristic, CharacteristicType),
               c, type,
              )


	return nothing
end 

"""
	CharacteristicSetUp(petsclib::PetscLibType,c::Characteristic) 
Sets up the internal data structures for the
later use of a `Charactoristic` .

Collective

Input Parameter:
- `c` - context obtained from CharacteristicCreate()

Level: developer

-seealso: [](ch_ts), `Characteristic`, `CharacteristicCreate()`, `CharacteristicSolve()`, `CharacteristicDestroy()`

# External Links
$(_doc_external("Ts/CharacteristicSetUp"))
"""
function CharacteristicSetUp(petsclib::PetscLibType, c::Characteristic) end

@for_petsc function CharacteristicSetUp(petsclib::$UnionPetscLib, c::Characteristic )

    @chk ccall(
               (:CharacteristicSetUp, $petsc_library),
               PetscErrorCode,
               (Characteristic,),
               c,
              )


	return nothing
end 

"""
	CharacteristicRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Adds an approarch to the method of characteristics package.

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - name of a new approach
- `function` - routine to create method context

Level: advanced

-seealso: [](ch_ts), `CharacteristicRegisterAll()`, `CharacteristicRegisterDestroy()`

# External Links
$(_doc_external("Ts/CharacteristicRegister"))
"""
function CharacteristicRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function CharacteristicRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:CharacteristicRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	CharacteristicSetVelocityInterpolation(petsclib::PetscLibType,c::Characteristic, da::PetscDM, v::PetscVec, vOld::PetscVec, numComponents::PetscInt, components::Vector{PetscInt}, interp::external, ctx::Cvoid) 

# External Links
$(_doc_external("Ts/CharacteristicSetVelocityInterpolation"))
"""
function CharacteristicSetVelocityInterpolation(petsclib::PetscLibType, c::Characteristic, da::PetscDM, v::PetscVec, vOld::PetscVec, numComponents::PetscInt, components::Vector{PetscInt}, interp::external, ctx::Cvoid) end

@for_petsc function CharacteristicSetVelocityInterpolation(petsclib::$UnionPetscLib, c::Characteristic, da::PetscDM, v::PetscVec, vOld::PetscVec, numComponents::$PetscInt, components::Vector{$PetscInt}, interp::external, ctx::Cvoid )

    @chk ccall(
               (:CharacteristicSetVelocityInterpolation, $petsc_library),
               PetscErrorCode,
               (Characteristic, CDM, CVec, CVec, $PetscInt, Ptr{$PetscInt}, external, Ptr{Cvoid}),
               c, da, v, vOld, numComponents, components, interp, ctx,
              )


	return nothing
end 

"""
	CharacteristicSetVelocityInterpolationLocal(petsclib::PetscLibType,c::Characteristic, da::PetscDM, v::PetscVec, vOld::PetscVec, numComponents::PetscInt, components::Vector{PetscInt}, interp::external, ctx::Cvoid) 

# External Links
$(_doc_external("Ts/CharacteristicSetVelocityInterpolationLocal"))
"""
function CharacteristicSetVelocityInterpolationLocal(petsclib::PetscLibType, c::Characteristic, da::PetscDM, v::PetscVec, vOld::PetscVec, numComponents::PetscInt, components::Vector{PetscInt}, interp::external, ctx::Cvoid) end

@for_petsc function CharacteristicSetVelocityInterpolationLocal(petsclib::$UnionPetscLib, c::Characteristic, da::PetscDM, v::PetscVec, vOld::PetscVec, numComponents::$PetscInt, components::Vector{$PetscInt}, interp::external, ctx::Cvoid )

    @chk ccall(
               (:CharacteristicSetVelocityInterpolationLocal, $petsc_library),
               PetscErrorCode,
               (Characteristic, CDM, CVec, CVec, $PetscInt, Ptr{$PetscInt}, external, Ptr{Cvoid}),
               c, da, v, vOld, numComponents, components, interp, ctx,
              )


	return nothing
end 

"""
	CharacteristicSetFieldInterpolation(petsclib::PetscLibType,c::Characteristic, da::PetscDM, v::PetscVec, numComponents::PetscInt, components::Vector{PetscInt}, interp::external, ctx::Cvoid) 

# External Links
$(_doc_external("Ts/CharacteristicSetFieldInterpolation"))
"""
function CharacteristicSetFieldInterpolation(petsclib::PetscLibType, c::Characteristic, da::PetscDM, v::PetscVec, numComponents::PetscInt, components::Vector{PetscInt}, interp::external, ctx::Cvoid) end

@for_petsc function CharacteristicSetFieldInterpolation(petsclib::$UnionPetscLib, c::Characteristic, da::PetscDM, v::PetscVec, numComponents::$PetscInt, components::Vector{$PetscInt}, interp::external, ctx::Cvoid )

    @chk ccall(
               (:CharacteristicSetFieldInterpolation, $petsc_library),
               PetscErrorCode,
               (Characteristic, CDM, CVec, $PetscInt, Ptr{$PetscInt}, external, Ptr{Cvoid}),
               c, da, v, numComponents, components, interp, ctx,
              )


	return nothing
end 

"""
	CharacteristicSetFieldInterpolationLocal(petsclib::PetscLibType,c::Characteristic, da::PetscDM, v::PetscVec, numComponents::PetscInt, components::Vector{PetscInt}, interp::external, ctx::Cvoid) 

# External Links
$(_doc_external("Ts/CharacteristicSetFieldInterpolationLocal"))
"""
function CharacteristicSetFieldInterpolationLocal(petsclib::PetscLibType, c::Characteristic, da::PetscDM, v::PetscVec, numComponents::PetscInt, components::Vector{PetscInt}, interp::external, ctx::Cvoid) end

@for_petsc function CharacteristicSetFieldInterpolationLocal(petsclib::$UnionPetscLib, c::Characteristic, da::PetscDM, v::PetscVec, numComponents::$PetscInt, components::Vector{$PetscInt}, interp::external, ctx::Cvoid )

    @chk ccall(
               (:CharacteristicSetFieldInterpolationLocal, $petsc_library),
               PetscErrorCode,
               (Characteristic, CDM, CVec, $PetscInt, Ptr{$PetscInt}, external, Ptr{Cvoid}),
               c, da, v, numComponents, components, interp, ctx,
              )


	return nothing
end 

"""
	CharacteristicSolve(petsclib::PetscLibType,c::Characteristic, dt::PetscReal, solution::PetscVec) 
Apply the Method of Characteristics solver

Collective

Input Parameters:
- `c`        - context obtained from `CharacteristicCreate()`
- `dt`       - the time-step
- `solution` - vector holding the solution

Level: developer

-seealso: [](ch_ts), `Characteristic`, `CharacteristicCreate()`, `CharacteristicDestroy()`

# External Links
$(_doc_external("Ts/CharacteristicSolve"))
"""
function CharacteristicSolve(petsclib::PetscLibType, c::Characteristic, dt::PetscReal, solution::PetscVec) end

@for_petsc function CharacteristicSolve(petsclib::$UnionPetscLib, c::Characteristic, dt::$PetscReal, solution::PetscVec )

    @chk ccall(
               (:CharacteristicSolve, $petsc_library),
               PetscErrorCode,
               (Characteristic, $PetscReal, CVec),
               c, dt, solution,
              )


	return nothing
end 

