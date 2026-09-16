"""
	c::Characteristic = CharacteristicCreate(petsclib::PetscLibType, comm::MPI_Comm) 
Creates a `Characteristic` context for use with the Method of Characteristics

Collective

Input Parameter:
- `comm` - MPI communicator

Output Parameter:
- `c` - the `Characteristic` context

Level: beginner

See also: `Characteristic`, `CharacteristicDestroy()`

# External Links
$(_doc_external("Characteristic/CharacteristicCreate"))
"""
function CharacteristicCreate(petsclib::PetscLibType, comm::MPI_Comm)
    error("CharacteristicCreate: no generated method for these argument types")
end

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
	CharacteristicDestroy(petsclib::PetscLibType, c::Union{Characteristic, Ref{Characteristic}}) 
Destroys a `Characteristic` context created with `CharacteristicCreate()`

Collective

Input Parameter:
- `c` - the `Characteristic` context

Level: beginner

See also: `Characteristic`, `CharacteristicCreate()`

# External Links
$(_doc_external("Characteristic/CharacteristicDestroy"))
"""
function CharacteristicDestroy(petsclib::PetscLibType, c::Union{Characteristic, Ref{Characteristic}})
    error("CharacteristicDestroy: no generated method for these argument types")
end

@for_petsc function CharacteristicDestroy(petsclib::$UnionPetscLib, c::Union{Characteristic, Ref{Characteristic}} )
	c_ = c isa Base.RefValue ? c : Ref{Characteristic}(c)

    @chk ccall(
               (:CharacteristicDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{Characteristic},),
               c_,
              )


	return nothing
end 

"""
	CharacteristicFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the `Characteristics` package. It is
called from `PetscFinalize()`.

Level: developer

See also: `PetscFinalize()`, `CharacteristicInitializePackage()`

# External Links
$(_doc_external("Characteristic/CharacteristicFinalizePackage"))
"""
function CharacteristicFinalizePackage(petsclib::PetscLibType)
    error("CharacteristicFinalizePackage: no generated method for these argument types")
end

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

See also: `PetscInitialize()`, `CharacteristicFinalizePackage()`

# External Links
$(_doc_external("Sys/CharacteristicInitializePackage"))
"""
function CharacteristicInitializePackage(petsclib::PetscLibType)
    error("CharacteristicInitializePackage: no generated method for these argument types")
end

@for_petsc function CharacteristicInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:CharacteristicInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	CharacteristicRegister(petsclib::PetscLibType, sname::String, fnc::external) 
Adds an approarch to the method of characteristics package.

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - name of a new approach
- `function` - routine to create method context

Level: advanced

See also: `CharacteristicRegisterAll()`, `CharacteristicRegisterDestroy()`

# External Links
$(_doc_external("Characteristic/CharacteristicRegister"))
"""
function CharacteristicRegister(petsclib::PetscLibType, sname::String, fnc::external)
    error("CharacteristicRegister: no generated method for these argument types")
end

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
	CharacteristicSetFieldInterpolation(petsclib::PetscLibType, c::Characteristic, da::AbstractPetscDM, v::AbstractPetscVec, numComponents::PetscInt, components::Vector{PetscInt}, interp::external, ctx::Ptr{Cvoid}) 
Sets the routine used to interpolate the field being advected at the foot of a characteristic

Not Collective

Input Parameters:
- `c`             - the `Characteristic` context
- `da`            - the `DM` describing the layout of the field vector
- `v`             - the field vector to be advected
- `numComponents` - the number of field components to interpolate
- `components`    - the indices of the field components in `v`
- `interp`        - the interpolation routine, called with the global vector
- `ctx`           - context passed to the interpolation routine

Calling sequence of `interp`:
- `v`             - the field `Vec` from which to interpolate
- `interpIndices` - the coordinates at which to interpolate
- `numComponents` - the number of components to interpolate
- `components`    - the indices of the components in `v`
- `values`        - the interpolated values, one per component per point
- `ctx`           - the application context

Level: developer

See also: `Characteristic`, `CharacteristicSetFieldInterpolationLocal()`, `CharacteristicSetVelocityInterpolation()`

# External Links
$(_doc_external("Characteristic/CharacteristicSetFieldInterpolation"))
"""
function CharacteristicSetFieldInterpolation(petsclib::PetscLibType, c::Characteristic, da::AbstractPetscDM, v::AbstractPetscVec, numComponents::Integer, components::AbstractVector{<:Number}, interp::external, ctx::Ptr{Cvoid})
    error("CharacteristicSetFieldInterpolation: no generated method for these argument types")
end

@for_petsc function CharacteristicSetFieldInterpolation(petsclib::$UnionPetscLib, c::Characteristic, da::AbstractPetscDM, v::AbstractPetscVec, numComponents::$PetscInt, components::Vector{$PetscInt}, interp::external, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:CharacteristicSetFieldInterpolation, $petsc_library),
               PetscErrorCode,
               (Characteristic, CDM, CVec, $PetscInt, Ptr{$PetscInt}, external, Ptr{Cvoid}),
               c, da, v, numComponents, components, interp, ctx,
              )


	return nothing
end 

"""
	CharacteristicSetFieldInterpolationLocal(petsclib::PetscLibType, c::Characteristic, da::AbstractPetscDM, v::AbstractPetscVec, numComponents::PetscInt, components::Vector{PetscInt}, interp::external, ctx::Ptr{Cvoid}) 
Sets the routine used to interpolate the field being advected at the foot of a characteristic using a locally-accessible array

Not Collective

Input Parameters:
- `c`             - the `Characteristic` context
- `da`            - the `DM` describing the layout of the field vector
- `v`             - the field vector to be advected
- `numComponents` - the number of field components to interpolate
- `components`    - the indices of the field components in `v`
- `interp`        - the interpolation routine, called with a local array pointer rather than a `Vec`
- `ctx`           - context passed to the interpolation routine

Calling sequence of `interp`:
- `array`         - the locally-accessible array of the field vector obtained from the `DM`
- `interpIndices` - the coordinates at which to interpolate
- `numComponents` - the number of components to interpolate
- `components`    - the indices of the components in the array
- `values`        - the interpolated values, one per component per point
- `ctx`           - the application context

Level: developer

See also: `Characteristic`, `CharacteristicSetFieldInterpolation()`, `CharacteristicSetVelocityInterpolationLocal()`

# External Links
$(_doc_external("Characteristic/CharacteristicSetFieldInterpolationLocal"))
"""
function CharacteristicSetFieldInterpolationLocal(petsclib::PetscLibType, c::Characteristic, da::AbstractPetscDM, v::AbstractPetscVec, numComponents::Integer, components::AbstractVector{<:Number}, interp::external, ctx::Ptr{Cvoid})
    error("CharacteristicSetFieldInterpolationLocal: no generated method for these argument types")
end

@for_petsc function CharacteristicSetFieldInterpolationLocal(petsclib::$UnionPetscLib, c::Characteristic, da::AbstractPetscDM, v::AbstractPetscVec, numComponents::$PetscInt, components::Vector{$PetscInt}, interp::external, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:CharacteristicSetFieldInterpolationLocal, $petsc_library),
               PetscErrorCode,
               (Characteristic, CDM, CVec, $PetscInt, Ptr{$PetscInt}, external, Ptr{Cvoid}),
               c, da, v, numComponents, components, interp, ctx,
              )


	return nothing
end 

"""
	CharacteristicSetType(petsclib::PetscLibType, c::Characteristic, type::String) 
Builds Characteristic for a particular solver.

Logically Collective

Input Parameters:
- `c`    - the method of characteristics context
- `type` - a known method

Options Database Key:
- `-characteristic_type method` - Sets the method; use -help for a list
of available methods

Level: intermediate

See also: `CharacteristicType`

# External Links
$(_doc_external("Characteristic/CharacteristicSetType"))
"""
function CharacteristicSetType(petsclib::PetscLibType, c::Characteristic, type::String)
    error("CharacteristicSetType: no generated method for these argument types")
end

@for_petsc function CharacteristicSetType(petsclib::$UnionPetscLib, c::Characteristic, type::String )

    @chk ccall(
               (:CharacteristicSetType, $petsc_library),
               PetscErrorCode,
               (Characteristic, CharacteristicType),
               c, type,
              )


	return nothing
end 

"""
	CharacteristicSetUp(petsclib::PetscLibType, c::Characteristic) 
Sets up the internal data structures for the
later use of a `Charactoristic` .

Collective

Input Parameter:
- `c` - context obtained from CharacteristicCreate()

Level: developer

See also: `Characteristic`, `CharacteristicCreate()`, `CharacteristicSolve()`, `CharacteristicDestroy()`

# External Links
$(_doc_external("Characteristic/CharacteristicSetUp"))
"""
function CharacteristicSetUp(petsclib::PetscLibType, c::Characteristic)
    error("CharacteristicSetUp: no generated method for these argument types")
end

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
	CharacteristicSetVelocityInterpolation(petsclib::PetscLibType, c::Characteristic, da::AbstractPetscDM, v::AbstractPetscVec, vOld::AbstractPetscVec, numComponents::PetscInt, components::Vector{PetscInt}, interp::external, ctx::Ptr{Cvoid}) 
Sets the routine used to interpolate the velocity field at points along a characteristic

Not Collective

Input Parameters:
- `c`             - the `Characteristic` context
- `da`            - the `DM` describing the layout of the velocity vectors
- `v`             - the current velocity vector
- `vOld`          - the previous-time-step velocity vector
- `numComponents` - the number of velocity components to interpolate
- `components`    - the indices of the velocity components in `v` and `vOld`
- `interp`        - the interpolation routine, called with the global vector
- `ctx`           - context passed to the interpolation routine

Calling sequence of `interp`:
- `v`             - the velocity `Vec` from which to interpolate
- `interpIndices` - the coordinates at which to interpolate
- `numComponents` - the number of components to interpolate
- `components`    - the indices of the components in `v`
- `values`        - the interpolated values, one per component per point
- `ctx`           - the application context

Level: developer

See also: `Characteristic`, `CharacteristicSetVelocityInterpolationLocal()`, `CharacteristicSetFieldInterpolation()`

# External Links
$(_doc_external("Characteristic/CharacteristicSetVelocityInterpolation"))
"""
function CharacteristicSetVelocityInterpolation(petsclib::PetscLibType, c::Characteristic, da::AbstractPetscDM, v::AbstractPetscVec, vOld::AbstractPetscVec, numComponents::Integer, components::AbstractVector{<:Number}, interp::external, ctx::Ptr{Cvoid})
    error("CharacteristicSetVelocityInterpolation: no generated method for these argument types")
end

@for_petsc function CharacteristicSetVelocityInterpolation(petsclib::$UnionPetscLib, c::Characteristic, da::AbstractPetscDM, v::AbstractPetscVec, vOld::AbstractPetscVec, numComponents::$PetscInt, components::Vector{$PetscInt}, interp::external, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:CharacteristicSetVelocityInterpolation, $petsc_library),
               PetscErrorCode,
               (Characteristic, CDM, CVec, CVec, $PetscInt, Ptr{$PetscInt}, external, Ptr{Cvoid}),
               c, da, v, vOld, numComponents, components, interp, ctx,
              )


	return nothing
end 

"""
	CharacteristicSetVelocityInterpolationLocal(petsclib::PetscLibType, c::Characteristic, da::AbstractPetscDM, v::AbstractPetscVec, vOld::AbstractPetscVec, numComponents::PetscInt, components::Vector{PetscInt}, interp::external, ctx::Ptr{Cvoid}) 
Sets the routine used to interpolate the velocity field along a characteristic using a locally-accessible array

Not Collective

Input Parameters:
- `c`             - the `Characteristic` context
- `da`            - the `DM` describing the layout of the velocity vectors
- `v`             - the current velocity vector
- `vOld`          - the previous-time-step velocity vector
- `numComponents` - the number of velocity components to interpolate
- `components`    - the indices of the velocity components in `v` and `vOld`
- `interp`        - the interpolation routine, called with a local array pointer rather than a `Vec`
- `ctx`           - context passed to the interpolation routine

Calling sequence of `interp`:
- `array`         - the locally-accessible array of the velocity vector obtained from the `DM`
- `interpIndices` - the coordinates at which to interpolate
- `numComponents` - the number of components to interpolate
- `components`    - the indices of the components in the array
- `values`        - the interpolated values, one per component per point
- `ctx`           - the application context

Level: developer

See also: `Characteristic`, `CharacteristicSetVelocityInterpolation()`, `CharacteristicSetFieldInterpolationLocal()`

# External Links
$(_doc_external("Characteristic/CharacteristicSetVelocityInterpolationLocal"))
"""
function CharacteristicSetVelocityInterpolationLocal(petsclib::PetscLibType, c::Characteristic, da::AbstractPetscDM, v::AbstractPetscVec, vOld::AbstractPetscVec, numComponents::Integer, components::AbstractVector{<:Number}, interp::external, ctx::Ptr{Cvoid})
    error("CharacteristicSetVelocityInterpolationLocal: no generated method for these argument types")
end

@for_petsc function CharacteristicSetVelocityInterpolationLocal(petsclib::$UnionPetscLib, c::Characteristic, da::AbstractPetscDM, v::AbstractPetscVec, vOld::AbstractPetscVec, numComponents::$PetscInt, components::Vector{$PetscInt}, interp::external, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:CharacteristicSetVelocityInterpolationLocal, $petsc_library),
               PetscErrorCode,
               (Characteristic, CDM, CVec, CVec, $PetscInt, Ptr{$PetscInt}, external, Ptr{Cvoid}),
               c, da, v, vOld, numComponents, components, interp, ctx,
              )


	return nothing
end 

"""
	CharacteristicSolve(petsclib::PetscLibType, c::Characteristic, dt::PetscReal, solution::AbstractPetscVec) 
Apply the Method of Characteristics solver

Collective

Input Parameters:
- `c`        - context obtained from `CharacteristicCreate()`
- `dt`       - the time-step
- `solution` - vector holding the solution

Level: developer

See also: `Characteristic`, `CharacteristicCreate()`, `CharacteristicDestroy()`

# External Links
$(_doc_external("Characteristic/CharacteristicSolve"))
"""
function CharacteristicSolve(petsclib::PetscLibType, c::Characteristic, dt::Real, solution::AbstractPetscVec)
    error("CharacteristicSolve: no generated method for these argument types")
end

@for_petsc function CharacteristicSolve(petsclib::$UnionPetscLib, c::Characteristic, dt::$PetscReal, solution::AbstractPetscVec )

    @chk ccall(
               (:CharacteristicSolve, $petsc_library),
               PetscErrorCode,
               (Characteristic, $PetscReal, CVec),
               c, dt, solution,
              )


	return nothing
end 

