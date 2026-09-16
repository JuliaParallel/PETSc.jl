"""
	lim::PetscLimiter = PetscLimiterCreate(petsclib::PetscLibType, comm::MPI_Comm) 
Creates an empty `PetscLimiter` object. The type can then be set with `PetscLimiterSetType()`.

Collective

Input Parameter:
- `comm` - The communicator for the `PetscLimiter` object

Output Parameter:
- `lim` - The `PetscLimiter` object

Level: beginner

See also: `PetscLimiter`, `PetscLimiterType`, `PetscLimiterSetType()`, `PETSCLIMITERSIN`

# External Links
$(_doc_external("FV/PetscLimiterCreate"))
"""
function PetscLimiterCreate(petsclib::PetscLibType, comm::MPI_Comm)
    error("PetscLimiterCreate: no generated method for these argument types")
end

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
	PetscLimiterDestroy(petsclib::PetscLibType, lim::Union{PetscLimiter, Ref{PetscLimiter}}) 
Destroys a `PetscLimiter` object

Collective

Input Parameter:
- `lim` - the `PetscLimiter` object to destroy

Level: beginner

See also: `PetscLimiter`, `PetscLimiterView()`

# External Links
$(_doc_external("FV/PetscLimiterDestroy"))
"""
function PetscLimiterDestroy(petsclib::PetscLibType, lim::Union{PetscLimiter, Ref{PetscLimiter}})
    error("PetscLimiterDestroy: no generated method for these argument types")
end

@for_petsc function PetscLimiterDestroy(petsclib::$UnionPetscLib, lim::Union{PetscLimiter, Ref{PetscLimiter}} )
	lim_ = lim isa Base.RefValue ? lim : Ref{PetscLimiter}(lim)

    @chk ccall(
               (:PetscLimiterDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscLimiter},),
               lim_,
              )


	return nothing
end 

"""
	name::PetscLimiterType = PetscLimiterGetType(petsclib::PetscLibType, lim::PetscLimiter) 
Gets the `PetscLimiterType` name (as a string) from the `PetscLimiter`.

Not Collective

Input Parameter:
- `lim` - The `PetscLimiter`

Output Parameter:
- `name` - The `PetscLimiterType`

Level: intermediate

See also: `PetscLimiter`, `PetscLimiterType`, `PetscLimiterSetType()`, `PetscLimiterCreate()`

# External Links
$(_doc_external("FV/PetscLimiterGetType"))
"""
function PetscLimiterGetType(petsclib::PetscLibType, lim::PetscLimiter)
    error("PetscLimiterGetType: no generated method for these argument types")
end

@for_petsc function PetscLimiterGetType(petsclib::$UnionPetscLib, lim::PetscLimiter )
	name_ = Ref{PetscLimiterType}()

    @chk ccall(
               (:PetscLimiterGetType, $petsc_library),
               PetscErrorCode,
               (PetscLimiter, Ptr{PetscLimiterType}),
               lim, name_,
              )

	name = name_[] == C_NULL ? "" : unsafe_string(name_[])

	return name
end 

"""
	phi::PetscReal = PetscLimiterLimit(petsclib::PetscLibType, lim::PetscLimiter, flim::PetscReal) 
Limit the flux

Input Parameters:
- `lim`  - The `PetscLimiter`
- `flim` - The input field

Output Parameter:
- `phi` - The limited field

Level: beginner

See also: `PetscLimiter`, `PetscLimiterType`, `PetscLimiterSetType()`, `PetscLimiterCreate()`

# External Links
$(_doc_external("FV/PetscLimiterLimit"))
"""
function PetscLimiterLimit(petsclib::PetscLibType, lim::PetscLimiter, flim::Real)
    error("PetscLimiterLimit: no generated method for these argument types")
end

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

"""
	PetscLimiterRegister(petsclib::PetscLibType, sname::String, fnc::external) 
Adds a new `PetscLimiter` implementation

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - The name of a new user-defined creation routine
- `function` - The creation routine

See also: `PetscLimiter`, `PetscLimiterType`, `PetscLimiterRegisterAll()`, `PetscLimiterRegisterDestroy()`

# External Links
$(_doc_external("FV/PetscLimiterRegister"))
"""
function PetscLimiterRegister(petsclib::PetscLibType, sname::String, fnc::external)
    error("PetscLimiterRegister: no generated method for these argument types")
end

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
	PetscLimiterSetFromOptions(petsclib::PetscLibType, lim::PetscLimiter) 
sets parameters in a `PetscLimiter` from the options database

Collective

Input Parameter:
- `lim` - the `PetscLimiter` object to set options for

Level: intermediate

See also: `PetscLimiter`, `PetscLimiterView()`

# External Links
$(_doc_external("FV/PetscLimiterSetFromOptions"))
"""
function PetscLimiterSetFromOptions(petsclib::PetscLibType, lim::PetscLimiter)
    error("PetscLimiterSetFromOptions: no generated method for these argument types")
end

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
	PetscLimiterSetType(petsclib::PetscLibType, lim::PetscLimiter, name::PetscLimiterType) 
Builds a `PetscLimiter` for a given `PetscLimiterType`

Collective

Input Parameters:
- `lim`  - The `PetscLimiter` object
- `name` - The kind of limiter

Options Database Key:
- `-petsclimiter_type type` - Sets the PetscLimiter type; use -help for a list of available types

Level: intermediate

See also: `PetscLimiter`, `PetscLimiterType`, `PetscLimiterGetType()`, `PetscLimiterCreate()`

# External Links
$(_doc_external("FV/PetscLimiterSetType"))
"""
function PetscLimiterSetType(petsclib::PetscLibType, lim::PetscLimiter, name::PetscLimiterType)
    error("PetscLimiterSetType: no generated method for these argument types")
end

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
	PetscLimiterSetUp(petsclib::PetscLibType, lim::PetscLimiter) 
Construct data structures for the `PetscLimiter`

Collective

Input Parameter:
- `lim` - the `PetscLimiter` object to setup

Level: intermediate

See also: `PetscLimiter`, `PetscLimiterView()`, `PetscLimiterDestroy()`

# External Links
$(_doc_external("FV/PetscLimiterSetUp"))
"""
function PetscLimiterSetUp(petsclib::PetscLibType, lim::PetscLimiter)
    error("PetscLimiterSetUp: no generated method for these argument types")
end

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
	PetscLimiterView(petsclib::PetscLibType, lim::PetscLimiter, v::PetscViewer) 
Views a `PetscLimiter`

Collective

Input Parameters:
- `lim` - the `PetscLimiter` object to view
- `v`   - the viewer

Level: beginner

See also: `PetscLimiter`, `PetscViewer`, `PetscLimiterDestroy()`, `PetscLimiterViewFromOptions()`

# External Links
$(_doc_external("FV/PetscLimiterView"))
"""
function PetscLimiterView(petsclib::PetscLibType, lim::PetscLimiter, v::PetscViewer)
    error("PetscLimiterView: no generated method for these argument types")
end

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
	PetscLimiterViewFromOptions(petsclib::PetscLibType, A::PetscLimiter, obj, name::String) 
View a `PetscLimiter` based on values in the options database

Collective

Input Parameters:
- `A`    - the `PetscLimiter` object to view
- `obj`  - Optional object that provides the options prefix to use
- `name` - command line option name

Options Database Key:
- `-name [viewertype][:...]` - option name and values. See `PetscObjectViewFromOptions()` for the possible arguments

Level: intermediate

See also: `PetscLimiter`, `PetscLimiterView()`, `PetscObjectViewFromOptions()`, `PetscLimiterCreate()`

# External Links
$(_doc_external("FV/PetscLimiterViewFromOptions"))
"""
function PetscLimiterViewFromOptions(petsclib::PetscLibType, A::PetscLimiter, obj, name::String)
    error("PetscLimiterViewFromOptions: no generated method for these argument types")
end

@for_petsc function PetscLimiterViewFromOptions(petsclib::$UnionPetscLib, A::PetscLimiter, obj, name::String )

    @chk ccall(
               (:PetscLimiterViewFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscLimiter, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

