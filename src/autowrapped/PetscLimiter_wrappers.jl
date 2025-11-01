# autodefined type arguments for class ------
mutable struct _n_PetscLimiter end
const PetscLimiter = Ptr{_n_PetscLimiter}

# -------------------------------------------------------
"""
	PetscLimiterRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Adds a new `PetscLimiter` implementation

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - The name of a new user-defined creation routine
- `function` - The creation routine

-seealso: `PetscLimiter`, `PetscLimiterType`, `PetscLimiterRegisterAll()`, `PetscLimiterRegisterDestroy()`

# External Links
$(_doc_external("DM/PetscLimiterRegister"))
"""
function PetscLimiterRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function PetscLimiterRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:PetscLimiterRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	PetscLimiterSetType(petsclib::PetscLibType,lim::PetscLimiter, name::PetscLimiterType) 
Builds a `PetscLimiter` for a given `PetscLimiterType`

Collective

Input Parameters:
- `lim`  - The `PetscLimiter` object
- `name` - The kind of limiter

Options Database Key:
- `-petsclimiter_type <type>` - Sets the PetscLimiter type; use -help for a list of available types

Level: intermediate

-seealso: `PetscLimiter`, `PetscLimiterType`, `PetscLimiterGetType()`, `PetscLimiterCreate()`

# External Links
$(_doc_external("DM/PetscLimiterSetType"))
"""
function PetscLimiterSetType(petsclib::PetscLibType, lim::PetscLimiter, name::PetscLimiterType) end

@for_petsc function PetscLimiterSetType(petsclib::$UnionPetscLib, lim::PetscLimiter, name::PetscLimiterType )

    @chk ccall(
               (:PetscLimiterSetType, $petsc_library),
               PetscErrorCode,
               (PetscLimiter, PetscLimiterType),
               lim, name,
              )


	return nothing
end 

"""
	name::PetscLimiterType = PetscLimiterGetType(petsclib::PetscLibType,lim::PetscLimiter) 
Gets the `PetscLimiterType` name (as a string) from the `PetscLimiter`.

Not Collective

Input Parameter:
- `lim` - The `PetscLimiter`

Output Parameter:
- `name` - The `PetscLimiterType`

Level: intermediate

-seealso: `PetscLimiter`, `PetscLimiterType`, `PetscLimiterSetType()`, `PetscLimiterCreate()`

# External Links
$(_doc_external("DM/PetscLimiterGetType"))
"""
function PetscLimiterGetType(petsclib::PetscLibType, lim::PetscLimiter) end

@for_petsc function PetscLimiterGetType(petsclib::$UnionPetscLib, lim::PetscLimiter )
	name_ = Ref{PetscLimiterType}()

    @chk ccall(
               (:PetscLimiterGetType, $petsc_library),
               PetscErrorCode,
               (PetscLimiter, Ptr{PetscLimiterType}),
               lim, name_,
              )

	name = unsafe_string(name_[])

	return name
end 

"""
	PetscLimiterViewFromOptions(petsclib::PetscLibType,A::PetscLimiter, obj::PetscObject, name::String) 
View a `PetscLimiter` based on values in the options database

Collective

Input Parameters:
- `A`    - the `PetscLimiter` object to view
- `obj`  - Optional object that provides the options prefix to use
- `name` - command line option name

Level: intermediate

-seealso: `PetscLimiter`, `PetscLimiterView()`, `PetscObjectViewFromOptions()`, `PetscLimiterCreate()`

# External Links
$(_doc_external("DM/PetscLimiterViewFromOptions"))
"""
function PetscLimiterViewFromOptions(petsclib::PetscLibType, A::PetscLimiter, obj::PetscObject, name::String) end

@for_petsc function PetscLimiterViewFromOptions(petsclib::$UnionPetscLib, A::PetscLimiter, obj::PetscObject, name::String )

    @chk ccall(
               (:PetscLimiterViewFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscLimiter, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

"""
	PetscLimiterView(petsclib::PetscLibType,lim::PetscLimiter, v::PetscViewer) 
Views a `PetscLimiter`

Collective

Input Parameters:
- `lim` - the `PetscLimiter` object to view
- `v`   - the viewer

Level: beginner

-seealso: `PetscLimiter`, `PetscViewer`, `PetscLimiterDestroy()`, `PetscLimiterViewFromOptions()`

# External Links
$(_doc_external("DM/PetscLimiterView"))
"""
function PetscLimiterView(petsclib::PetscLibType, lim::PetscLimiter, v::PetscViewer) end

@for_petsc function PetscLimiterView(petsclib::$UnionPetscLib, lim::PetscLimiter, v::PetscViewer )

    @chk ccall(
               (:PetscLimiterView, $petsc_library),
               PetscErrorCode,
               (PetscLimiter, PetscViewer),
               lim, v,
              )


	return nothing
end 

"""
	PetscLimiterSetFromOptions(petsclib::PetscLibType,lim::PetscLimiter) 
sets parameters in a `PetscLimiter` from the options database

Collective

Input Parameter:
- `lim` - the `PetscLimiter` object to set options for

Level: intermediate

-seealso: `PetscLimiter`, `PetscLimiterView()`

# External Links
$(_doc_external("DM/PetscLimiterSetFromOptions"))
"""
function PetscLimiterSetFromOptions(petsclib::PetscLibType, lim::PetscLimiter) end

@for_petsc function PetscLimiterSetFromOptions(petsclib::$UnionPetscLib, lim::PetscLimiter )

    @chk ccall(
               (:PetscLimiterSetFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscLimiter,),
               lim,
              )


	return nothing
end 

"""
	PetscLimiterSetUp(petsclib::PetscLibType,lim::PetscLimiter) 
Construct data structures for the `PetscLimiter`

Collective

Input Parameter:
- `lim` - the `PetscLimiter` object to setup

Level: intermediate

-seealso: `PetscLimiter`, `PetscLimiterView()`, `PetscLimiterDestroy()`

# External Links
$(_doc_external("DM/PetscLimiterSetUp"))
"""
function PetscLimiterSetUp(petsclib::PetscLibType, lim::PetscLimiter) end

@for_petsc function PetscLimiterSetUp(petsclib::$UnionPetscLib, lim::PetscLimiter )

    @chk ccall(
               (:PetscLimiterSetUp, $petsc_library),
               PetscErrorCode,
               (PetscLimiter,),
               lim,
              )


	return nothing
end 

"""
	PetscLimiterDestroy(petsclib::PetscLibType,lim::PetscLimiter) 
Destroys a `PetscLimiter` object

Collective

Input Parameter:
- `lim` - the `PetscLimiter` object to destroy

Level: beginner

-seealso: `PetscLimiter`, `PetscLimiterView()`

# External Links
$(_doc_external("DM/PetscLimiterDestroy"))
"""
function PetscLimiterDestroy(petsclib::PetscLibType, lim::PetscLimiter) end

@for_petsc function PetscLimiterDestroy(petsclib::$UnionPetscLib, lim::PetscLimiter )

    @chk ccall(
               (:PetscLimiterDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscLimiter},),
               lim,
              )


	return nothing
end 

"""
	lim::PetscLimiter = PetscLimiterCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Creates an empty `PetscLimiter` object. The type can then be set with `PetscLimiterSetType()`.

Collective

Input Parameter:
- `comm` - The communicator for the `PetscLimiter` object

Output Parameter:
- `lim` - The `PetscLimiter` object

Level: beginner

-seealso: `PetscLimiter`, `PetscLimiterType`, `PetscLimiterSetType()`, `PETSCLIMITERSIN`

# External Links
$(_doc_external("DM/PetscLimiterCreate"))
"""
function PetscLimiterCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function PetscLimiterCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	lim_ = Ref{PetscLimiter}()

    @chk ccall(
               (:PetscLimiterCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscLimiter}),
               comm, lim_,
              )

	lim = lim_[]

	return lim
end 

"""
	phi::PetscReal = PetscLimiterLimit(petsclib::PetscLibType,lim::PetscLimiter, flim::PetscReal) 
Limit the flux

Input Parameters:
- `lim`  - The `PetscLimiter`
- `flim` - The input field

Output Parameter:
- `phi` - The limited field

Level: beginner

-seealso: `PetscLimiter`, `PetscLimiterType`, `PetscLimiterSetType()`, `PetscLimiterCreate()`

# External Links
$(_doc_external("DM/PetscLimiterLimit"))
"""
function PetscLimiterLimit(petsclib::PetscLibType, lim::PetscLimiter, flim::PetscReal) end

@for_petsc function PetscLimiterLimit(petsclib::$UnionPetscLib, lim::PetscLimiter, flim::$PetscReal )
	phi_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscLimiterLimit, $petsc_library),
               PetscErrorCode,
               (PetscLimiter, $PetscReal, Ptr{$PetscReal}),
               lim, flim, phi_,
              )

	phi = phi_[]

	return phi
end 

