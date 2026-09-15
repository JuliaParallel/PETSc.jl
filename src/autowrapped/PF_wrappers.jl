"""
	y::PetscScalar = PFApply(petsclib::PetscLibType,pf::AbstractPF, n::PetscInt, x::Vector{PetscScalar}) 
Applies the mathematical function to an array of values.

Collective

Input Parameters:
- `pf` - the function context
- `n`  - number of pointwise function evaluations to perform, each pointwise function evaluation
is a function of dimin variables and computes dimout variables where dimin and dimout are defined
in the call to `PFCreate()`
- `x`  - input array

Output Parameter:
- `y` - output array

Level: beginner

-seealso: `PF`, `PFApplyVec()`, `PFCreate()`, `PFDestroy()`, `PFSetType()`, `PFSet()`

# External Links
$(_doc_external("PF/PFApply"))
"""
function PFApply(petsclib::PetscLibType, pf::AbstractPF, n::Integer, x::AbstractVector{<:Number})
    error("PFApply: no generated method for these argument types")
end

@for_petsc function PFApply(petsclib::$UnionPetscLib, pf::AbstractPF, n::$PetscInt, x::Vector{$PetscScalar} )
	y_ = Ref{$PetscScalar}()

    @chk ccall(
               (:PFApply, $petsc_library),
               PetscErrorCode,
               (CPF, $PetscInt, Ptr{$PetscScalar}, Ptr{$PetscScalar}),
               pf, n, x, y_,
              )

	y = y_[]

	return y
end 

"""
	PFApplyVec(petsclib::PetscLibType,pf::AbstractPF, x::AbstractPetscVec, y::AbstractPetscVec) 
Applies the mathematical function to a vector

Collective

Input Parameters:
- `pf` - the function context
- `x`  - input vector (or `NULL` for the vector (0,1, .... N-1)

Output Parameter:
- `y` - output vector

Level: beginner

-seealso: `PF`, `PFApply()`, `PFCreate()`, `PFDestroy()`, `PFSetType()`, `PFSet()`

# External Links
$(_doc_external("PF/PFApplyVec"))
"""
function PFApplyVec(petsclib::PetscLibType, pf::AbstractPF, x::AbstractPetscVec, y::AbstractPetscVec)
    error("PFApplyVec: no generated method for these argument types")
end

@for_petsc function PFApplyVec(petsclib::$UnionPetscLib, pf::AbstractPF, x::AbstractPetscVec, y::AbstractPetscVec )

    @chk ccall(
               (:PFApplyVec, $petsc_library),
               PetscErrorCode,
               (CPF, CVec, CVec),
               pf, x, y,
              )


	return nothing
end 

"""
	pf::PF = PFCreate(petsclib::PetscLibType,comm::MPI_Comm, dimin::PetscInt, dimout::PetscInt) 
Creates a mathematical function context.

Collective

Input Parameters:
- `comm`   - MPI communicator
- `dimin`  - dimension of the space you are mapping from
- `dimout` - dimension of the space you are mapping to

Output Parameter:
- `pf` - the function context

Level: developer

-seealso: `PF`, `PFSet()`, `PFApply()`, `PFDestroy()`, `PFApplyVec()`

# External Links
$(_doc_external("PF/PFCreate"))
"""
function PFCreate(petsclib::PetscLibType, comm::MPI_Comm, dimin::Integer, dimout::Integer)
    error("PFCreate: no generated method for these argument types")
end

@for_petsc function PFCreate(petsclib::$UnionPetscLib, comm::MPI_Comm, dimin::$PetscInt, dimout::$PetscInt )
	pf_ = Ref{CPF}()

    @chk ccall(
               (:PFCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, Ptr{CPF}),
               comm, dimin, dimout, pf_,
              )

	pf = PF(pf_[], petsclib)

	return pf
end 

"""
	PFDestroy(petsclib::PetscLibType,pf::AbstractPF) 
Destroys `PF` context that was created with `PFCreate()`.

Collective

Input Parameter:
- `pf` - the function context

Level: beginner

-seealso: `PF`, `PFCreate()`, `PFSet()`, `PFSetType()`

# External Links
$(_doc_external("PF/PFDestroy"))
"""
function PFDestroy(petsclib::PetscLibType, pf::AbstractPF)
    error("PFDestroy: no generated method for these argument types")
end

@for_petsc function PFDestroy(petsclib::$UnionPetscLib, pf::AbstractPF )
	pf_ = Ref(pf.ptr)

    @chk ccall(
               (:PFDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{CPF},),
               pf_,
              )

	pf.ptr = C_NULL

	return nothing
end 

"""
	PFFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the PETSc `PF` package. It is
called from `PetscFinalize()`.

Level: developer

-seealso: `PF`, `PetscFinalize()`

# External Links
$(_doc_external("PF/PFFinalizePackage"))
"""
function PFFinalizePackage(petsclib::PetscLibType)
    error("PFFinalizePackage: no generated method for these argument types")
end

@for_petsc function PFFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PFFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	type::PFType = PFGetType(petsclib::PetscLibType,pf::AbstractPF) 
Gets the `PFType` name (as a string) from the `PF`
context.

Not Collective

Input Parameter:
- `pf` - the function context

Output Parameter:
- `type` - name of function

Level: intermediate

-seealso: `PF`, `PFSetType()`

# External Links
$(_doc_external("PF/PFGetType"))
"""
function PFGetType(petsclib::PetscLibType, pf::AbstractPF)
    error("PFGetType: no generated method for these argument types")
end

@for_petsc function PFGetType(petsclib::$UnionPetscLib, pf::AbstractPF )
	type_ = Ref{PFType}()

    @chk ccall(
               (:PFGetType, $petsc_library),
               PetscErrorCode,
               (CPF, Ptr{PFType}),
               pf, type_,
              )

	type = type_[] == C_NULL ? "" : unsafe_string(type_[])

	return type
end 

"""
	PFInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `PF` package. It is called
from PetscDLLibraryRegister_petscvec() when using dynamic libraries, and on the first call to `PFCreate()`
when using shared or static libraries.

Level: developer

-seealso: `PF`, `PetscInitialize()`

# External Links
$(_doc_external("PF/PFInitializePackage"))
"""
function PFInitializePackage(petsclib::PetscLibType)
    error("PFInitializePackage: no generated method for these argument types")
end

@for_petsc function PFInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PFInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PFRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Adds a method to the mathematical function package.

Not Collective

Input Parameters:
- `sname`    - name of a new user-defined solver
- `function` - routine to create method context

-seealso: `PF`, `PFRegisterAll()`, `PFRegisterDestroy()`

# External Links
$(_doc_external("PF/PFRegister"))
"""
function PFRegister(petsclib::PetscLibType, sname::String, fnc::external)
    error("PFRegister: no generated method for these argument types")
end

@for_petsc function PFRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:PFRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	PFSet(petsclib::PetscLibType,pf::AbstractPF, apply::external, applyvec::external, view::external, destroy::external, ctx::Ptr{Cvoid}) 
Sets the C/C++/Fortran functions to be used by the PF function

Collective

Input Parameters:
- `pf`       - the function context
- `apply`    - function to apply to an array
- `applyvec` - function to apply to a Vec
- `view`     - function that prints information about the `PF`
- `destroy`  - function to free the private function context
- `ctx`      - private function context

Level: beginner

-seealso: `PF`, `PFCreate()`, `PFDestroy()`, `PFSetType()`, `PFApply()`, `PFApplyVec()`

# External Links
$(_doc_external("PF/PFSet"))
"""
function PFSet(petsclib::PetscLibType, pf::AbstractPF, apply::external, applyvec::external, view::external, destroy::external, ctx::Ptr{Cvoid})
    error("PFSet: no generated method for these argument types")
end

@for_petsc function PFSet(petsclib::$UnionPetscLib, pf::AbstractPF, apply::external, applyvec::external, view::external, destroy::external, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:PFSet, $petsc_library),
               PetscErrorCode,
               (CPF, external, external, external, external, Ptr{Cvoid}),
               pf, apply, applyvec, view, destroy, ctx,
              )


	return nothing
end 

"""
	PFSetFromOptions(petsclib::PetscLibType,pf::AbstractPF) 
Sets `PF` options from the options database.

Collective

Input Parameters:
- `pf` - the mathematical function context

Level: intermediate

-seealso: `PF`

# External Links
$(_doc_external("PF/PFSetFromOptions"))
"""
function PFSetFromOptions(petsclib::PetscLibType, pf::AbstractPF)
    error("PFSetFromOptions: no generated method for these argument types")
end

@for_petsc function PFSetFromOptions(petsclib::$UnionPetscLib, pf::AbstractPF )

    @chk ccall(
               (:PFSetFromOptions, $petsc_library),
               PetscErrorCode,
               (CPF,),
               pf,
              )


	return nothing
end 

"""
	PFSetType(petsclib::PetscLibType,pf::AbstractPF, type::PFType, ctx::Ptr{Cvoid}) 
Builds `PF` for a particular function

Collective

Input Parameters:
- `pf`   - the function context.
- `type` - a known method
- `ctx`  - optional type dependent context

Options Database Key:
- `-pf_type <type>` - Sets PF type

Level: intermediate

-seealso: `PF`, `PFSet()`, `PFRegister()`, `PFCreate()`, `DMDACreatePF()`

# External Links
$(_doc_external("PF/PFSetType"))
"""
function PFSetType(petsclib::PetscLibType, pf::AbstractPF, type::PFType, ctx::Ptr{Cvoid})
    error("PFSetType: no generated method for these argument types")
end

@for_petsc function PFSetType(petsclib::$UnionPetscLib, pf::AbstractPF, type::PFType, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:PFSetType, $petsc_library),
               PetscErrorCode,
               (CPF, PFType, Ptr{Cvoid}),
               pf, type, ctx,
              )


	return nothing
end 

"""
	PFStringSetFunction(petsclib::PetscLibType,pf::AbstractPF, string::String) 
Creates a function from a string

Collective

Input Parameters:
- `pf`     - the function object
- `string` - the string that defines the function

Level: intermediate

-seealso: `PFSetFromOptions()`

# External Links
$(_doc_external("PF/PFStringSetFunction"))
"""
function PFStringSetFunction(petsclib::PetscLibType, pf::AbstractPF, string::String)
    error("PFStringSetFunction: no generated method for these argument types")
end

@for_petsc function PFStringSetFunction(petsclib::$UnionPetscLib, pf::AbstractPF, string::String )

    @chk ccall(
               (:PFStringSetFunction, $petsc_library),
               PetscErrorCode,
               (CPF, Ptr{Cchar}),
               pf, string,
              )


	return nothing
end 

"""
	PFView(petsclib::PetscLibType,pf::AbstractPF, viewer::PetscViewer) 
Prints information about a mathematical function

Collective unless `viewer` is `PETSC_VIEWER_STDOUT_SELF`

Input Parameters:
- `pf`     - the `PF` context
- `viewer` - optional visualization context

Level: developer

-seealso: `PF`, `PetscViewerCreate()`, `PetscViewerASCIIOpen()`

# External Links
$(_doc_external("PF/PFView"))
"""
function PFView(petsclib::PetscLibType, pf::AbstractPF, viewer::PetscViewer)
    error("PFView: no generated method for these argument types")
end

@for_petsc function PFView(petsclib::$UnionPetscLib, pf::AbstractPF, viewer::PetscViewer )

    @chk ccall(
               (:PFView, $petsc_library),
               PetscErrorCode,
               (CPF, PetscViewer),
               pf, viewer,
              )


	return nothing
end 

"""
	PFViewFromOptions(petsclib::PetscLibType,A::AbstractPF, obj, name::String) 
View a `PF` based on options set in the options database

Collective

Input Parameters:
- `A`    - the `PF` context
- `obj`  - Optional object that provides the prefix used to search the options database
- `name` - command line option

Level: intermediate

-seealso: `PF`, `PFView`, `PetscObjectViewFromOptions()`, `PFCreate()`

# External Links
$(_doc_external("PF/PFViewFromOptions"))
"""
function PFViewFromOptions(petsclib::PetscLibType, A::AbstractPF, obj, name::String)
    error("PFViewFromOptions: no generated method for these argument types")
end

@for_petsc function PFViewFromOptions(petsclib::$UnionPetscLib, A::AbstractPF, obj, name::String )

    @chk ccall(
               (:PFViewFromOptions, $petsc_library),
               PetscErrorCode,
               (CPF, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

