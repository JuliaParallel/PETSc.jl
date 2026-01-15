# autodefined type arguments for class ------
mutable struct _n_PetscMatlabEngine end
const PetscMatlabEngine = Ptr{_n_PetscMatlabEngine}
# -------------------------------------------------------

"""
	mengine::PetscMatlabEngine = PetscMatlabEngineCreate(petsclib::PetscLibType,comm::MPI_Comm, host::String) 
Creates a MATLAB engine object

Not Collective

Input Parameters:
- `comm` - a separate MATLAB engine is started for each process in the communicator
- `host` - name of machine where MATLAB engine is to be run (usually NULL)

Output Parameter:
- `mengine` - the resulting object

Options Database Keys:
- `-matlab_engine_graphics` - allow the MATLAB engine to display graphics
- `-matlab_engine_host`     - hostname, machine to run the MATLAB engine on
- `-info`                   - print out all requests to MATLAB and all if its responses (for debugging)

Level: advanced

-seealso: `PetscMatlabEngineDestroy()`, `PetscMatlabEnginePut()`, `PetscMatlabEngineGet()`,
`PetscMatlabEngineEvaluate()`, `PetscMatlabEngineGetOutput()`, `PetscMatlabEnginePrintOutput()`,
`PETSC_MATLAB_ENGINE_()`, `PetscMatlabEnginePutArray()`, `PetscMatlabEngineGetArray()`, `PetscMatlabEngine`

# External Links
$(_doc_external("Sys/PetscMatlabEngineCreate"))
"""
function PetscMatlabEngineCreate(petsclib::PetscLibType, comm::MPI_Comm, host::String) end

@for_petsc function PetscMatlabEngineCreate(petsclib::$UnionPetscLib, comm::MPI_Comm, host::String )
	mengine_ = Ref{PetscMatlabEngine}()

    @chk ccall(
               (:PetscMatlabEngineCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{PetscMatlabEngine}),
               comm, host, mengine_,
              )

	mengine = mengine_[]

	return mengine
end 

"""
	PetscMatlabEngineDestroy(petsclib::PetscLibType,v::PetscMatlabEngine) 
Shuts down a MATLAB engine.

Collective

Input Parameter:
- `v` - the engine

Level: advanced

-seealso: `PetscMatlabEngineCreate()`, `PetscMatlabEnginePut()`, `PetscMatlabEngineGet()`,
`PetscMatlabEngineEvaluate()`, `PetscMatlabEngineGetOutput()`, `PetscMatlabEnginePrintOutput()`,
`PETSC_MATLAB_ENGINE_()`, `PetscMatlabEnginePutArray()`, `PetscMatlabEngineGetArray()`, `PetscMatlabEngine`

# External Links
$(_doc_external("Sys/PetscMatlabEngineDestroy"))
"""
function PetscMatlabEngineDestroy(petsclib::PetscLibType, v::PetscMatlabEngine) end

@for_petsc function PetscMatlabEngineDestroy(petsclib::$UnionPetscLib, v::PetscMatlabEngine )

    @chk ccall(
               (:PetscMatlabEngineDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscMatlabEngine},),
               v,
              )


	return nothing
end 

"""
	PetscMatlabEngineGetOutput(petsclib::PetscLibType,mengine::PetscMatlabEngine, string::String) 
Gets a string buffer where the MATLAB output is
printed

Not Collective

Input Parameter:
- `mengine` - the MATLAB engine

Output Parameter:
- `string` - buffer where MATLAB output is printed

Level: advanced

-seealso: `PetscMatlabEngineDestroy()`, `PetscMatlabEnginePut()`, `PetscMatlabEngineGet()`,
`PetscMatlabEngineEvaluate()`, `PetscMatlabEngineCreate()`, `PetscMatlabEnginePrintOutput()`,
`PETSC_MATLAB_ENGINE_()`, `PetscMatlabEnginePutArray()`, `PetscMatlabEngineGetArray()`, `PetscMatlabEngine`

# External Links
$(_doc_external("Sys/PetscMatlabEngineGetOutput"))
"""
function PetscMatlabEngineGetOutput(petsclib::PetscLibType, mengine::PetscMatlabEngine, string::String) end

@for_petsc function PetscMatlabEngineGetOutput(petsclib::$UnionPetscLib, mengine::PetscMatlabEngine, string::String )
	string_ = Ref(pointer(string))

    @chk ccall(
               (:PetscMatlabEngineGetOutput, $petsc_library),
               PetscErrorCode,
               (PetscMatlabEngine, Ptr{Ptr{Cchar}}),
               mengine, string_,
              )


	return nothing
end 

"""
	PetscMatlabEnginePrintOutput(petsclib::PetscLibType,mengine::PetscMatlabEngine, fd::Libc.FILE) 
prints the output from MATLAB to an ASCII file

Collective

Input Parameters:
- `mengine` - the MATLAB engine
- `fd`      - the file

Level: advanced

-seealso: `PetscMatlabEngineDestroy()`, `PetscMatlabEnginePut()`, `PetscMatlabEngineGet()`,
`PetscMatlabEngineEvaluate()`, `PetscMatlabEngineGetOutput()`, `PetscMatlabEngineCreate()`,
`PETSC_MATLAB_ENGINE_()`, `PetscMatlabEnginePutArray()`, `PetscMatlabEngineGetArray()`, `PetscMatlabEngine`

# External Links
$(_doc_external("Sys/PetscMatlabEnginePrintOutput"))
"""
function PetscMatlabEnginePrintOutput(petsclib::PetscLibType, mengine::PetscMatlabEngine, fd::Libc.FILE) end

@for_petsc function PetscMatlabEnginePrintOutput(petsclib::$UnionPetscLib, mengine::PetscMatlabEngine, fd::Libc.FILE )

    @chk ccall(
               (:PetscMatlabEnginePrintOutput, $petsc_library),
               PetscErrorCode,
               (PetscMatlabEngine, Ptr{Libc.FILE}),
               mengine, fd,
              )


	return nothing
end 

"""
	PetscMatlabEnginePut(petsclib::PetscLibType,mengine::PetscMatlabEngine, obj::PetscObject) 
Puts a PETSc object, such as a `Mat` or `Vec` into the MATLAB space. For parallel objects,
each processor's part is put in a separate  MATLAB process.

Collective

Input Parameters:
- `mengine` - the MATLAB engine
- `obj`     - the PETSc object, for example Vec

Level: advanced

-seealso: `PetscMatlabEngineDestroy()`, `PetscMatlabEngineCreate()`, `PetscMatlabEngineGet()`,
`PetscMatlabEngineEvaluate()`, `PetscMatlabEngineGetOutput()`, `PetscMatlabEnginePrintOutput()`,
`PETSC_MATLAB_ENGINE_()`, `PetscMatlabEnginePutArray()`, `PetscMatlabEngineGetArray()`, `PetscMatlabEngine`

# External Links
$(_doc_external("Sys/PetscMatlabEnginePut"))
"""
function PetscMatlabEnginePut(petsclib::PetscLibType, mengine::PetscMatlabEngine, obj::PetscObject) end

@for_petsc function PetscMatlabEnginePut(petsclib::$UnionPetscLib, mengine::PetscMatlabEngine, obj::PetscObject )

    @chk ccall(
               (:PetscMatlabEnginePut, $petsc_library),
               PetscErrorCode,
               (PetscMatlabEngine, PetscObject),
               mengine, obj,
              )


	return nothing
end 

"""
	PetscMatlabEngineGet(petsclib::PetscLibType,mengine::PetscMatlabEngine, obj::PetscObject) 
Gets a variable from MATLAB into a PETSc object.

Collective

Input Parameters:
- `mengine` - the MATLAB engine
- `obj`     - the PETSc object, for example a `Vec`

Level: advanced

-seealso: `PetscMatlabEngineDestroy()`, `PetscMatlabEnginePut()`, `PetscMatlabEngineCreate()`,
`PetscMatlabEngineEvaluate()`, `PetscMatlabEngineGetOutput()`, `PetscMatlabEnginePrintOutput()`,
`PETSC_MATLAB_ENGINE_()`, `PetscMatlabEnginePutArray()`, `PetscMatlabEngineGetArray()`, `PetscMatlabEngine`

# External Links
$(_doc_external("Sys/PetscMatlabEngineGet"))
"""
function PetscMatlabEngineGet(petsclib::PetscLibType, mengine::PetscMatlabEngine, obj::PetscObject) end

@for_petsc function PetscMatlabEngineGet(petsclib::$UnionPetscLib, mengine::PetscMatlabEngine, obj::PetscObject )

    @chk ccall(
               (:PetscMatlabEngineGet, $petsc_library),
               PetscErrorCode,
               (PetscMatlabEngine, PetscObject),
               mengine, obj,
              )


	return nothing
end 

"""
	PetscMatlabEnginePutArray(petsclib::PetscLibType,mengine::PetscMatlabEngine, m::Cint, n::Cint, array::Vector{PetscScalar}, name::String) 
Puts an array into the MATLAB space, treating it as a Fortran style (column major ordering) array. For parallel objects,
each processors part is put in a separate  MATLAB process.

Collective

Input Parameters:
- `mengine` - the MATLAB engine
- `m`       - the x dimension of the array
- `n`       - the y dimension of the array
- `array`   - the array (represented in one dimension)
- `name`    - the name of the array

Level: advanced

-seealso: `PetscMatlabEngineDestroy()`, `PetscMatlabEngineCreate()`, `PetscMatlabEngineGet()`,
`PetscMatlabEngineEvaluate()`, `PetscMatlabEngineGetOutput()`, `PetscMatlabEnginePrintOutput()`,
`PETSC_MATLAB_ENGINE_()`, `PetscMatlabEnginePut()`, `PetscMatlabEngineGetArray()`, `PetscMatlabEngine`

# External Links
$(_doc_external("Sys/PetscMatlabEnginePutArray"))
"""
function PetscMatlabEnginePutArray(petsclib::PetscLibType, mengine::PetscMatlabEngine, m::Cint, n::Cint, array::Vector{PetscScalar}, name::String) end

@for_petsc function PetscMatlabEnginePutArray(petsclib::$UnionPetscLib, mengine::PetscMatlabEngine, m::Cint, n::Cint, array::Vector{$PetscScalar}, name::String )

    @chk ccall(
               (:PetscMatlabEnginePutArray, $petsc_library),
               PetscErrorCode,
               (PetscMatlabEngine, Cint, Cint, Ptr{$PetscScalar}, Ptr{Cchar}),
               mengine, m, n, array, name,
              )


	return nothing
end 

"""
	PetscMatlabEngineGetArray(petsclib::PetscLibType,mengine::PetscMatlabEngine, m::Cint, n::Cint, array::Vector{PetscScalar}, name::String) 
Gets a variable from MATLAB into an array

Not Collective

Input Parameters:
- `mengine` - the MATLAB engine
- `m`       - the x dimension of the array
- `n`       - the y dimension of the array
- `array`   - the array (represented in one dimension), much be large enough to hold all the data
- `name`    - the name of the array

Level: advanced

-seealso: `PetscMatlabEngineDestroy()`, `PetscMatlabEnginePut()`, `PetscMatlabEngineCreate()`,
`PetscMatlabEngineEvaluate()`, `PetscMatlabEngineGetOutput()`, `PetscMatlabEnginePrintOutput()`,
`PETSC_MATLAB_ENGINE_()`, `PetscMatlabEnginePutArray()`, `PetscMatlabEngineGet()`, `PetscMatlabEngine`

# External Links
$(_doc_external("Sys/PetscMatlabEngineGetArray"))
"""
function PetscMatlabEngineGetArray(petsclib::PetscLibType, mengine::PetscMatlabEngine, m::Cint, n::Cint, array::Vector{PetscScalar}, name::String) end

@for_petsc function PetscMatlabEngineGetArray(petsclib::$UnionPetscLib, mengine::PetscMatlabEngine, m::Cint, n::Cint, array::Vector{$PetscScalar}, name::String )

    @chk ccall(
               (:PetscMatlabEngineGetArray, $petsc_library),
               PetscErrorCode,
               (PetscMatlabEngine, Cint, Cint, Ptr{$PetscScalar}, Ptr{Cchar}),
               mengine, m, n, array, name,
              )


	return nothing
end 

