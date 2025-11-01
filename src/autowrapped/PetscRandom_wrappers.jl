"""
	PetscRandomDestroy(petsclib::PetscLibType,r::PetscRandom) 
Destroys a `PetscRandom` object that was created by `PetscRandomCreate()`.

Collective

Input Parameter:
- `r` - the random number generator object

Level: intermediate

-seealso: `PetscRandom`, `PetscRandomGetValue()`, `PetscRandomCreate()`, `VecSetRandom()`

# External Links
$(_doc_external("Sys/PetscRandomDestroy"))
"""
function PetscRandomDestroy(petsclib::PetscLibType, r::PetscRandom) end

@for_petsc function PetscRandomDestroy(petsclib::$UnionPetscLib, r::PetscRandom )

    @chk ccall(
               (:PetscRandomDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscRandom},),
               r,
              )


	return nothing
end 

"""
	PetscRandomGetSeed(petsclib::PetscLibType,r::PetscRandom, seed::PetscInt64) 
Gets the random seed.

Not collective

Input Parameter:
- `r` - The random number generator context

Output Parameter:
- `seed` - The random seed

Level: intermediate

-seealso: `PetscRandom`, `PetscRandomCreate()`, `PetscRandomSetSeed()`, `PetscRandomSeed()`

# External Links
$(_doc_external("Sys/PetscRandomGetSeed"))
"""
function PetscRandomGetSeed(petsclib::PetscLibType, r::PetscRandom, seed::PetscInt64) end

@for_petsc function PetscRandomGetSeed(petsclib::$UnionPetscLib, r::PetscRandom, seed::$PetscInt64 )

    @chk ccall(
               (:PetscRandomGetSeed, $petsc_library),
               PetscErrorCode,
               (PetscRandom, Ptr{$PetscInt64}),
               r, seed,
              )


	return nothing
end 

"""
	PetscRandomSetSeed(petsclib::PetscLibType,r::PetscRandom, seed::PetscInt64) 
Sets the random seed. You MUST call `PetscRandomSeed()` after this call to have the new seed used.

Not collective

Input Parameters:
- `r`    - The random number generator context
- `seed` - The random seed

Level: intermediate

-seealso: `PetscRandom`, `PetscRandomCreate()`, `PetscRandomGetSeed()`, `PetscRandomSeed()`

# External Links
$(_doc_external("Sys/PetscRandomSetSeed"))
"""
function PetscRandomSetSeed(petsclib::PetscLibType, r::PetscRandom, seed::PetscInt64) end

@for_petsc function PetscRandomSetSeed(petsclib::$UnionPetscLib, r::PetscRandom, seed::$PetscInt64 )

    @chk ccall(
               (:PetscRandomSetSeed, $petsc_library),
               PetscErrorCode,
               (PetscRandom, $PetscInt64),
               r, seed,
              )


	return nothing
end 

"""
	PetscRandomSetFromOptions(petsclib::PetscLibType,rnd::PetscRandom) 
Configures the random number generator from the options database.

Collective

Input Parameter:
- `rnd` - The random number generator context

Options Database Keys:
- `-random_seed <integer>`    - provide a seed to the random number generator
- `-random_no_imaginary_part` - makes the imaginary part of the random number zero, this is useful when you want the
same code to produce the same result when run with real numbers or complex numbers for regression testing purposes

Level: beginner

-seealso: `PetscRandom`, `PetscRandomCreate()`, `PetscRandomSetType()`

# External Links
$(_doc_external("Sys/PetscRandomSetFromOptions"))
"""
function PetscRandomSetFromOptions(petsclib::PetscLibType, rnd::PetscRandom) end

@for_petsc function PetscRandomSetFromOptions(petsclib::$UnionPetscLib, rnd::PetscRandom )

    @chk ccall(
               (:PetscRandomSetFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscRandom,),
               rnd,
              )


	return nothing
end 

"""
	PetscRandomSetOptionsPrefix(petsclib::PetscLibType,r::PetscRandom, prefix::String) 
Sets the prefix used for searching for all
`PetscRandom` options in the database.

Logically Collective

Input Parameters:
- `r`      - the random number generator context
- `prefix` - the prefix to prepend to all option names

Level: advanced

-seealso: `PetscRandom`, `PetscRandomSetFromOptions()`

# External Links
$(_doc_external("Sys/PetscRandomSetOptionsPrefix"))
"""
function PetscRandomSetOptionsPrefix(petsclib::PetscLibType, r::PetscRandom, prefix::String) end

@for_petsc function PetscRandomSetOptionsPrefix(petsclib::$UnionPetscLib, r::PetscRandom, prefix::String )

    @chk ccall(
               (:PetscRandomSetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (PetscRandom, Ptr{Cchar}),
               r, prefix,
              )


	return nothing
end 

"""
	PetscRandomViewFromOptions(petsclib::PetscLibType,A::PetscRandom, obj::PetscObject, name::String) 
View a `PetscRandom` object based on the options database

Collective

Input Parameters:
- `A`    - the random number generator context
- `obj`  - Optional object
- `name` - command line option

Level: intermediate

-seealso: `PetscRandom`, `PetscRandomView`, `PetscObjectViewFromOptions()`, `PetscRandomCreate()`

# External Links
$(_doc_external("Sys/PetscRandomViewFromOptions"))
"""
function PetscRandomViewFromOptions(petsclib::PetscLibType, A::PetscRandom, obj::PetscObject, name::String) end

@for_petsc function PetscRandomViewFromOptions(petsclib::$UnionPetscLib, A::PetscRandom, obj::PetscObject, name::String )

    @chk ccall(
               (:PetscRandomViewFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscRandom, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

"""
	PetscRandomView(petsclib::PetscLibType,rnd::PetscRandom, viewer::PetscViewer) 
Views a random number generator object.

Collective

Input Parameters:
- `rnd`    - The random number generator context
- `viewer` - an optional visualization context

Level: beginner

-seealso: `PetscRandom`, `PetscRealView()`, `PetscScalarView()`, `PetscIntView()`

# External Links
$(_doc_external("Sys/PetscRandomView"))
"""
function PetscRandomView(petsclib::PetscLibType, rnd::PetscRandom, viewer::PetscViewer) end

@for_petsc function PetscRandomView(petsclib::$UnionPetscLib, rnd::PetscRandom, viewer::PetscViewer )

    @chk ccall(
               (:PetscRandomView, $petsc_library),
               PetscErrorCode,
               (PetscRandom, PetscViewer),
               rnd, viewer,
              )


	return nothing
end 

"""
	r::PetscRandom = PetscRandomCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Creates an object for generating random numbers,
and initializes the random-number generator.

Collective

Input Parameter:
- `comm` - MPI communicator

Output Parameter:
- `r` - the random number generator object

Level: intermediate

-seealso: `PetscRandomSetType()`, `PetscRandomGetValue()`, `PetscRandomGetValueReal()`, `PetscRandomSetInterval()`,
`PetscRandomDestroy()`, `VecSetRandom()`, `PetscRandomType`, `PetscRandom`

# External Links
$(_doc_external("Sys/PetscRandomCreate"))
"""
function PetscRandomCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function PetscRandomCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	r_ = Ref{PetscRandom}()

    @chk ccall(
               (:PetscRandomCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscRandom}),
               comm, r_,
              )

	r = r_[]

	return r
end 

"""
	PetscRandomSeed(petsclib::PetscLibType,r::PetscRandom) 
Seed the random number generator.

Not collective

Input Parameter:
- `r` - The random number generator context

Level: intermediate

-seealso: `PetscRandomCreate()`, `PetscRandomGetSeed()`, `PetscRandomSetSeed()`

# External Links
$(_doc_external("Sys/PetscRandomSeed"))
"""
function PetscRandomSeed(petsclib::PetscLibType, r::PetscRandom) end

@for_petsc function PetscRandomSeed(petsclib::$UnionPetscLib, r::PetscRandom )

    @chk ccall(
               (:PetscRandomSeed, $petsc_library),
               PetscErrorCode,
               (PetscRandom,),
               r,
              )


	return nothing
end 

"""
	PetscRandomSetType(petsclib::PetscLibType,rnd::PetscRandom, type::PetscRandomType) 
Builds a context for generating a particular type of random numbers.

Collective

Input Parameters:
- `rnd`  - The random number generator context
- `type` - The name of the random type

Options Database Key:
- `-random_type <type>` - Sets the random type; use -help for a list
of available types

Level: intermediate

-seealso: `PetscRandom`, `PetscRandomType`, `PetscRandomGetType()`, `PetscRandomCreate()`

# External Links
$(_doc_external("Sys/PetscRandomSetType"))
"""
function PetscRandomSetType(petsclib::PetscLibType, rnd::PetscRandom, type::PetscRandomType) end

@for_petsc function PetscRandomSetType(petsclib::$UnionPetscLib, rnd::PetscRandom, type::PetscRandomType )

    @chk ccall(
               (:PetscRandomSetType, $petsc_library),
               PetscErrorCode,
               (PetscRandom, PetscRandomType),
               rnd, type,
              )


	return nothing
end 

"""
	type::PetscRandomType = PetscRandomGetType(petsclib::PetscLibType,rnd::PetscRandom) 
Gets the type name (as a string) from the `PetscRandom`.

Not Collective

Input Parameter:
- `rnd` - The random number generator context

Output Parameter:
- `type` - The type name

Level: intermediate

-seealso: `PetscRandom`, `PetscRandomType`, `PetscRandomSetType()`, `PetscRandomCreate()`

# External Links
$(_doc_external("Sys/PetscRandomGetType"))
"""
function PetscRandomGetType(petsclib::PetscLibType, rnd::PetscRandom) end

@for_petsc function PetscRandomGetType(petsclib::$UnionPetscLib, rnd::PetscRandom )
	type_ = Ref{PetscRandomType}()

    @chk ccall(
               (:PetscRandomGetType, $petsc_library),
               PetscErrorCode,
               (PetscRandom, Ptr{PetscRandomType}),
               rnd, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	PetscRandomRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Adds a new `PetscRandom` implementation

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - The name of a new user-defined creation routine
- `function` - The creation routine

Level: advanced

-seealso: `PetscRandom`, `PetscRandomRegisterAll()`, `PetscRandomRegisterDestroy()`

# External Links
$(_doc_external("Sys/PetscRandomRegister"))
"""
function PetscRandomRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function PetscRandomRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:PetscRandomRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	val::PetscScalar = PetscRandomGetValue(petsclib::PetscLibType,r::PetscRandom) 
Generates a random number.  Call this after first calling
`PetscRandomCreate()`.

Not Collective

Input Parameter:
- `r` - the random number generator context

Output Parameter:
- `val` - the value

Level: intermediate

-seealso: `PetscRandom`, `PetscRandomCreate()`, `PetscRandomDestroy()`, `VecSetRandom()`, `PetscRandomGetValueReal()`, `PetscRandomSetInterval()`

# External Links
$(_doc_external("Sys/PetscRandomGetValue"))
"""
function PetscRandomGetValue(petsclib::PetscLibType, r::PetscRandom) end

@for_petsc function PetscRandomGetValue(petsclib::$UnionPetscLib, r::PetscRandom )
	val_ = Ref{$PetscScalar}()

    @chk ccall(
               (:PetscRandomGetValue, $petsc_library),
               PetscErrorCode,
               (PetscRandom, Ptr{$PetscScalar}),
               r, val_,
              )

	val = val_[]

	return val
end 

"""
	val::PetscReal = PetscRandomGetValueReal(petsclib::PetscLibType,r::PetscRandom) 
Generates a real random number.  Call this after first calling
`PetscRandomCreate()`.

Not Collective

Input Parameter:
- `r` - the random number generator context

Output Parameter:
- `val` - the value

Level: intermediate

-seealso: `PetscRandom`, `PetscRandomCreate()`, `PetscRandomDestroy()`, `VecSetRandom()`, `PetscRandomGetValue()`

# External Links
$(_doc_external("Sys/PetscRandomGetValueReal"))
"""
function PetscRandomGetValueReal(petsclib::PetscLibType, r::PetscRandom) end

@for_petsc function PetscRandomGetValueReal(petsclib::$UnionPetscLib, r::PetscRandom )
	val_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscRandomGetValueReal, $petsc_library),
               PetscErrorCode,
               (PetscRandom, Ptr{$PetscReal}),
               r, val_,
              )

	val = val_[]

	return val
end 

"""
	val::PetscScalar = PetscRandomGetValues(petsclib::PetscLibType,r::PetscRandom, n::PetscInt) 
Generates a sequence of random numbers.  Call this after first calling
`PetscRandomCreate()`.

Not Collective

Input Parameters:
- `r` - the random number generator context
- `n` - number of random numbers to generate

Output Parameter:
- `val` - the array to hold the values

Level: intermediate

-seealso: `PetscRandom`, `PetscRandomCreate()`, `PetscRandomDestroy()`, `VecSetRandom()`, `PetscRandomGetValue()`

# External Links
$(_doc_external("Sys/PetscRandomGetValues"))
"""
function PetscRandomGetValues(petsclib::PetscLibType, r::PetscRandom, n::PetscInt) end

@for_petsc function PetscRandomGetValues(petsclib::$UnionPetscLib, r::PetscRandom, n::$PetscInt )
	val_ = Ref{$PetscScalar}()

    @chk ccall(
               (:PetscRandomGetValues, $petsc_library),
               PetscErrorCode,
               (PetscRandom, $PetscInt, Ptr{$PetscScalar}),
               r, n, val_,
              )

	val = val_[]

	return val
end 

"""
	val::PetscReal = PetscRandomGetValuesReal(petsclib::PetscLibType,r::PetscRandom, n::PetscInt) 
Generates a sequence of real random numbers.  Call this after first calling
`PetscRandomCreate()`.

Not Collective

Input Parameters:
- `r` - the random number generator context
- `n` - number of random numbers to generate

Output Parameter:
- `val` - the array to hold the values

Level: intermediate

-seealso: `PetscRandom`, `PetscRandomCreate()`, `PetscRandomDestroy()`, `VecSetRandom()`, `PetscRandomGetValues()`

# External Links
$(_doc_external("Sys/PetscRandomGetValuesReal"))
"""
function PetscRandomGetValuesReal(petsclib::PetscLibType, r::PetscRandom, n::PetscInt) end

@for_petsc function PetscRandomGetValuesReal(petsclib::$UnionPetscLib, r::PetscRandom, n::$PetscInt )
	val_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscRandomGetValuesReal, $petsc_library),
               PetscErrorCode,
               (PetscRandom, $PetscInt, Ptr{$PetscReal}),
               r, n, val_,
              )

	val = val_[]

	return val
end 

"""
	low::PetscScalar,high::PetscScalar = PetscRandomGetInterval(petsclib::PetscLibType,r::PetscRandom) 
Gets the interval over which the random numbers
will be distributed.  By default, this interval is [0,1).

Not Collective

Input Parameter:
- `r` - the random number generator context

Output Parameters:
- `low`  - The lower bound of the interval
- `high` - The upper bound of the interval

Level: intermediate

-seealso: `PetscRandom`, `PetscRandomCreate()`, `PetscRandomSetInterval()`

# External Links
$(_doc_external("Sys/PetscRandomGetInterval"))
"""
function PetscRandomGetInterval(petsclib::PetscLibType, r::PetscRandom) end

@for_petsc function PetscRandomGetInterval(petsclib::$UnionPetscLib, r::PetscRandom )
	low_ = Ref{$PetscScalar}()
	high_ = Ref{$PetscScalar}()

    @chk ccall(
               (:PetscRandomGetInterval, $petsc_library),
               PetscErrorCode,
               (PetscRandom, Ptr{$PetscScalar}, Ptr{$PetscScalar}),
               r, low_, high_,
              )

	low = low_[]
	high = high_[]

	return low,high
end 

"""
	PetscRandomSetInterval(petsclib::PetscLibType,r::PetscRandom, low::PetscScalar, high::PetscScalar) 
Sets the interval over which the random numbers
will be distributed.  By default, this interval is [0,1).

Not Collective

Input Parameters:
- `r`    - the random number generator context
- `low`  - The lower bound of the interval
- `high` - The upper bound of the interval

Level: intermediate

-seealso: `PetscRandomCreate()`, `PetscRandomGetInterval()`

# External Links
$(_doc_external("Sys/PetscRandomSetInterval"))
"""
function PetscRandomSetInterval(petsclib::PetscLibType, r::PetscRandom, low::PetscScalar, high::PetscScalar) end

@for_petsc function PetscRandomSetInterval(petsclib::$UnionPetscLib, r::PetscRandom, low::$PetscScalar, high::$PetscScalar )

    @chk ccall(
               (:PetscRandomSetInterval, $petsc_library),
               PetscErrorCode,
               (PetscRandom, $PetscScalar, $PetscScalar),
               r, low, high,
              )


	return nothing
end 

"""
	PetscRandomFinalizePackage(petsclib::PetscLibType) 
This function frees everything in the `PetscRandom` package. It is
called from `PetscFinalize()`.

Level: developer

-seealso: `PetscFinalize()`

# External Links
$(_doc_external("Sys/PetscRandomFinalizePackage"))
"""
function PetscRandomFinalizePackage(petsclib::PetscLibType) end

@for_petsc function PetscRandomFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscRandomFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscRandomInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `PetscRandom` package. It is called
from PetscDLLibraryRegister_petsc() when using dynamic libraries, and on the first call to `PetscRandomCreate()`
when using shared or static libraries.

Level: developer

-seealso: `PetscInitialize()`

# External Links
$(_doc_external("Sys/PetscRandomInitializePackage"))
"""
function PetscRandomInitializePackage(petsclib::PetscLibType) end

@for_petsc function PetscRandomInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscRandomInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

