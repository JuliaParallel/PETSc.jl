"""
	adm::PetscDM,ax::PetscVec = DMAdaptorAdapt(petsclib::PetscLibType, adaptor::DMAdaptor, x::AbstractPetscVec, strategy::DMAdaptationStrategy) 
Creates a new `DM` that is adapted to the problem

Not Collective

Input Parameters:
- `adaptor`  - The `DMAdaptor` object
- `x`        - The global approximate solution
- `strategy` - The adaptation strategy, see `DMAdaptationStrategy`

Output Parameters:
- `adm` - The adapted `DM`
- `ax`  - The adapted solution

Options Database Keys:
- `-snes_adapt (initial|sequential|multigrid)` - adaption strategy, see `DMAdaptationStrategy`
- `-adapt_gradient_view`                       - View the Clement interpolant of the solution gradient
- `-adapt_hessian_view`                        - View the Clement interpolant of the solution Hessian
- `-adapt_metric_view`                         - View the metric tensor for adaptive mesh refinement

Level: intermediate

See also: `DMAdaptor`, `DMAdaptationStrategy`, `DMAdaptorSetSolver()`, `DMAdaptorCreate()`

# External Links
$(_doc_external("DM/DMAdaptorAdapt"))
"""
function DMAdaptorAdapt(petsclib::PetscLibType, adaptor::DMAdaptor, x::AbstractPetscVec, strategy::DMAdaptationStrategy)
    error("DMAdaptorAdapt: no generated method for these argument types")
end

@for_petsc function DMAdaptorAdapt(petsclib::$UnionPetscLib, adaptor::DMAdaptor, x::AbstractPetscVec, strategy::DMAdaptationStrategy )
	adm_ = Ref{CDM}()
	ax_ = Ref{CVec}()

    @chk ccall(
               (:DMAdaptorAdapt, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, CVec, DMAdaptationStrategy, Ptr{CDM}, Ptr{CVec}),
               adaptor, x, strategy, adm_, ax_,
              )

	adm = PetscDM(adm_[], petsclib)
	ax = PetscVec(ax_[], petsclib)

	return adm,ax
end 

"""
	adaptor::DMAdaptor = DMAdaptorCreate(petsclib::PetscLibType, comm::MPI_Comm) 
Create a `DMAdaptor` object. Its purpose is to construct a adaptation `DMLabel` or metric `Vec` that can be used to modify the `DM`.

Collective

Input Parameter:
- `comm` - The communicator for the `DMAdaptor` object

Output Parameter:
- `adaptor` - The `DMAdaptor` object

Level: beginner

See also: `DM`, `DMAdaptor`, `DMAdaptorDestroy()`, `DMAdaptorAdapt()`, `PetscConvEst`, `PetscConvEstCreate()`

# External Links
$(_doc_external("DM/DMAdaptorCreate"))
"""
function DMAdaptorCreate(petsclib::PetscLibType, comm::MPI_Comm)
    error("DMAdaptorCreate: no generated method for these argument types")
end

@for_petsc function DMAdaptorCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	adaptor_ = Ref{DMAdaptor}()

    @chk ccall(
               (:DMAdaptorCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{DMAdaptor}),
               comm, adaptor_,
              )

	adaptor = adaptor_[]

	return adaptor
end 

"""
	DMAdaptorDestroy(petsclib::PetscLibType, adaptor::Union{DMAdaptor, Ref{DMAdaptor}}) 
Destroys a `DMAdaptor` object

Collective

Input Parameter:
- `adaptor` - The `DMAdaptor` object

Level: beginner

See also: `DM`, `DMAdaptor`, `DMAdaptorCreate()`, `DMAdaptorAdapt()`

# External Links
$(_doc_external("DM/DMAdaptorDestroy"))
"""
function DMAdaptorDestroy(petsclib::PetscLibType, adaptor::Union{DMAdaptor, Ref{DMAdaptor}})
    error("DMAdaptorDestroy: no generated method for these argument types")
end

@for_petsc function DMAdaptorDestroy(petsclib::$UnionPetscLib, adaptor::Union{DMAdaptor, Ref{DMAdaptor}} )
	adaptor_ = adaptor isa Base.RefValue ? adaptor : Ref{DMAdaptor}(adaptor)

    @chk ccall(
               (:DMAdaptorDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{DMAdaptor},),
               adaptor_,
              )


	return nothing
end 

"""
	criterion::DMAdaptationCriterion = DMAdaptorGetCriterion(petsclib::PetscLibType, adaptor::DMAdaptor) 
Get the adaptation criterion

Not Collective

Input Parameter:
- `adaptor` - the `DMAdaptor`

Output Parameter:
- `criterion` - the criterion for adaptation

Level: advanced

See also: `DMAdaptor`, `DMAdaptorSetCriterion()`, `DMAdaptationCriterion`

# External Links
$(_doc_external("DM/DMAdaptorGetCriterion"))
"""
function DMAdaptorGetCriterion(petsclib::PetscLibType, adaptor::DMAdaptor)
    error("DMAdaptorGetCriterion: no generated method for these argument types")
end

@for_petsc function DMAdaptorGetCriterion(petsclib::$UnionPetscLib, adaptor::DMAdaptor )
	criterion_ = Ref{DMAdaptationCriterion}()

    @chk ccall(
               (:DMAdaptorGetCriterion, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, Ptr{DMAdaptationCriterion}),
               adaptor, criterion_,
              )

	criterion = criterion_[]

	return criterion
end 

"""
	DMAdaptorGetMixedSetupFunction(petsclib::PetscLibType, adaptor::DMAdaptor, noname::Ptr{Cvoid}) 
Get the function setting up the mixed problem, if it exists

Not Collective

Input Parameter:
- `adaptor` - the `DMAdaptor`

Output Parameter:
- `setupFunc` - the function setting up the mixed problem, or `NULL`

Level: advanced

See also: `DMAdaptor`, `DMAdaptorSetMixedSetupFunction()`, `DMAdaptorAdapt()`

# External Links
$(_doc_external("DM/DMAdaptorGetMixedSetupFunction"))
"""
function DMAdaptorGetMixedSetupFunction(petsclib::PetscLibType, adaptor::DMAdaptor, noname::Ptr{Cvoid})
    error("DMAdaptorGetMixedSetupFunction: no generated method for these argument types")
end

@for_petsc function DMAdaptorGetMixedSetupFunction(petsclib::$UnionPetscLib, adaptor::DMAdaptor, noname::Ptr{Cvoid} )

    @chk ccall(
               (:DMAdaptorGetMixedSetupFunction, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, Ptr{Cvoid}),
               adaptor, noname,
              )


	return nothing
end 

"""
	num::PetscInt = DMAdaptorGetSequenceLength(petsclib::PetscLibType, adaptor::DMAdaptor) 
Gets the number of sequential adaptations used by an adapter

Not Collective

Input Parameter:
- `adaptor` - The `DMAdaptor` object

Output Parameter:
- `num` - The number of adaptations

Level: intermediate

See also: `DMAdaptor`, `DMAdaptorSetSequenceLength()`, `DMAdaptorCreate()`, `DMAdaptorAdapt()`

# External Links
$(_doc_external("DM/DMAdaptorGetSequenceLength"))
"""
function DMAdaptorGetSequenceLength(petsclib::PetscLibType, adaptor::DMAdaptor)
    error("DMAdaptorGetSequenceLength: no generated method for these argument types")
end

@for_petsc function DMAdaptorGetSequenceLength(petsclib::$UnionPetscLib, adaptor::DMAdaptor )
	num_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMAdaptorGetSequenceLength, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, Ptr{$PetscInt}),
               adaptor, num_,
              )

	num = num_[]

	return num
end 

"""
	snes::SNES = DMAdaptorGetSolver(petsclib::PetscLibType, adaptor::DMAdaptor) 
Gets the solver used to produce discrete solutions

Not Collective

Input Parameter:
- `adaptor` - The `DMAdaptor` object

Output Parameter:
- `snes` - The solver

Level: intermediate

See also: `DM`, `DMAdaptor`, `DMAdaptorSetSolver()`, `DMAdaptorCreate()`, `DMAdaptorAdapt()`

# External Links
$(_doc_external("DM/DMAdaptorGetSolver"))
"""
function DMAdaptorGetSolver(petsclib::PetscLibType, adaptor::DMAdaptor)
    error("DMAdaptorGetSolver: no generated method for these argument types")
end

@for_petsc function DMAdaptorGetSolver(petsclib::$UnionPetscLib, adaptor::DMAdaptor )
	snes_ = Ref{CSNES}()

    @chk ccall(
               (:DMAdaptorGetSolver, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, Ptr{CSNES}),
               adaptor, snes_,
              )

	snes = SNES(snes_[], petsclib)

	return snes
end 

"""
	DMAdaptorGetTransferFunction(petsclib::PetscLibType, adaptor::DMAdaptor, noname::Ptr{Cvoid}) 
Get the callback used by a `DMAdaptor` to transfer a solution vector from an old `DM` to the adapted `DM`

Not Collective

Input Parameter:
- `adaptor` - the `DMAdaptor` object

Output Parameter:
- `tfunc` - pointer to the transfer callback

Calling sequence of `tfunc`:
- `adaptor` - the `DMAdaptor` object
- `dm`      - the current `DM`
- `xin`     - the current solution
- `newdm`   - the adapted `DM`
- `xout`    - the transferred solution on `newdm`
- `ctx`     - application context, set with `DMSetApplicationContext()`

Level: developer

See also: `DMAdaptor`, `DMAdaptorSetTransferFunction()`, `DMAdaptorAdapt()`

# External Links
$(_doc_external("DM/DMAdaptorGetTransferFunction"))
"""
function DMAdaptorGetTransferFunction(petsclib::PetscLibType, adaptor::DMAdaptor, noname::Ptr{Cvoid})
    error("DMAdaptorGetTransferFunction: no generated method for these argument types")
end

@for_petsc function DMAdaptorGetTransferFunction(petsclib::$UnionPetscLib, adaptor::DMAdaptor, noname::Ptr{Cvoid} )

    @chk ccall(
               (:DMAdaptorGetTransferFunction, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, Ptr{Cvoid}),
               adaptor, noname,
              )


	return nothing
end 

"""
	type::DMAdaptorType = DMAdaptorGetType(petsclib::PetscLibType, adaptor::DMAdaptor) 
Gets the type name (as a string) from the adaptor.

Not Collective

Input Parameter:
- `adaptor` - The `DMAdaptor`

Output Parameter:
- `type` - The `DMAdaptorType` name

Level: intermediate

See also: `DM`, `DMPLEX`, `DMAdaptor`, `DMAdaptorType`, `DMAdaptorSetType()`, `DMAdaptorCreate()`

# External Links
$(_doc_external("DM/DMAdaptorGetType"))
"""
function DMAdaptorGetType(petsclib::PetscLibType, adaptor::DMAdaptor)
    error("DMAdaptorGetType: no generated method for these argument types")
end

@for_petsc function DMAdaptorGetType(petsclib::$UnionPetscLib, adaptor::DMAdaptor )
	type_ = Ref{DMAdaptorType}()

    @chk ccall(
               (:DMAdaptorGetType, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, Ptr{DMAdaptorType}),
               adaptor, type_,
              )

	type = type_[] == C_NULL ? "" : unsafe_string(type_[])

	return type
end 

"""
	DMAdaptorMonitor(petsclib::PetscLibType, adaptor::DMAdaptor, it::PetscInt, odm::AbstractPetscDM, adm::AbstractPetscDM, Nf::PetscInt, enorms::Vector{PetscReal}, error::AbstractPetscVec) 
runs the user provided monitor routines, if they exist

Collective

Input Parameters:
- `adaptor` - the `DMAdaptor`
- `it`      - iteration number
- `odm`     - the original `DM`
- `adm`     - the adapted `DM`
- `Nf`      - the number of fields
- `enorms`  - the 2-norm error values for each field
- `error`   - `Vec` of cellwise errors

Level: developer

See also: `DMAdaptorMonitorSet()`

# External Links
$(_doc_external("DM/DMAdaptorMonitor"))
"""
function DMAdaptorMonitor(petsclib::PetscLibType, adaptor::DMAdaptor, it::Integer, odm::AbstractPetscDM, adm::AbstractPetscDM, Nf::Integer, enorms::AbstractVector{<:Number}, error::AbstractPetscVec)
    error("DMAdaptorMonitor: no generated method for these argument types")
end

@for_petsc function DMAdaptorMonitor(petsclib::$UnionPetscLib, adaptor::DMAdaptor, it::$PetscInt, odm::AbstractPetscDM, adm::AbstractPetscDM, Nf::$PetscInt, enorms::Vector{$PetscReal}, error::AbstractPetscVec )

    @chk ccall(
               (:DMAdaptorMonitor, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, $PetscInt, CDM, CDM, $PetscInt, Ptr{$PetscReal}, CVec),
               adaptor, it, odm, adm, Nf, enorms, error,
              )


	return nothing
end 

"""
	DMAdaptorMonitorCancel(petsclib::PetscLibType, adaptor::DMAdaptor) 
Clears all monitors for a `DMAdaptor` object.

Logically Collective

Input Parameter:
- `adaptor` - the `DMAdaptor`

Options Database Key:
- `-dm_adaptor_monitor_cancel` - Cancels all monitors that have been hardwired into a code by calls to `DMAdaptorMonitorSet()`, but does not cancel those set via the options database.

Level: intermediate

See also: `DMAdaptorMonitorError()`, `DMAdaptorMonitorSet()`, `DMAdaptor`

# External Links
$(_doc_external("DM/DMAdaptorMonitorCancel"))
"""
function DMAdaptorMonitorCancel(petsclib::PetscLibType, adaptor::DMAdaptor)
    error("DMAdaptorMonitorCancel: no generated method for these argument types")
end

@for_petsc function DMAdaptorMonitorCancel(petsclib::$UnionPetscLib, adaptor::DMAdaptor )

    @chk ccall(
               (:DMAdaptorMonitorCancel, $petsc_library),
               PetscErrorCode,
               (DMAdaptor,),
               adaptor,
              )


	return nothing
end 

"""
	DMAdaptorMonitorError(petsclib::PetscLibType, adaptor::DMAdaptor, n::PetscInt, odm::AbstractPetscDM, adm::AbstractPetscDM, Nf::PetscInt, enorms::Vector{PetscReal}, error::AbstractPetscVec, vf::Vector{PetscViewerAndFormat}) 
Prints the error norm at each iteration of an adaptation loop.

Collective

Input Parameters:
- `adaptor` - the `DMAdaptor`
- `n`       - iteration number
- `odm`     - the original `DM`
- `adm`     - the adapted `DM`
- `Nf`      - number of fields
- `enorms`  - 2-norm error values for each field (may be estimated).
- `error`   - `Vec` of cellwise errors
- `vf`      - The viewer context

Options Database Key:
- `-adaptor_monitor_error` - Activates `DMAdaptorMonitorError()`

Level: intermediate

See also: `DMAdaptor`, `DMAdaptorMonitorSet()`, `DMAdaptorMonitorErrorDraw()`, `DMAdaptorMonitorErrorDrawLG()`

# External Links
$(_doc_external("DM/DMAdaptorMonitorError"))
"""
function DMAdaptorMonitorError(petsclib::PetscLibType, adaptor::DMAdaptor, n::Integer, odm::AbstractPetscDM, adm::AbstractPetscDM, Nf::Integer, enorms::AbstractVector{<:Number}, error::AbstractPetscVec, vf::Vector{PetscViewerAndFormat})
    error("DMAdaptorMonitorError: no generated method for these argument types")
end

@for_petsc function DMAdaptorMonitorError(petsclib::$UnionPetscLib, adaptor::DMAdaptor, n::$PetscInt, odm::AbstractPetscDM, adm::AbstractPetscDM, Nf::$PetscInt, enorms::Vector{$PetscReal}, error::AbstractPetscVec, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:DMAdaptorMonitorError, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, $PetscInt, CDM, CDM, $PetscInt, Ptr{$PetscReal}, CVec, Ptr{PetscViewerAndFormat}),
               adaptor, n, odm, adm, Nf, enorms, error, vf,
              )


	return nothing
end 

"""
	DMAdaptorMonitorErrorDraw(petsclib::PetscLibType, adaptor::DMAdaptor, n::PetscInt, odm::AbstractPetscDM, adm::AbstractPetscDM, Nf::PetscInt, enorms::Vector{PetscReal}, error::AbstractPetscVec, vf::Vector{PetscViewerAndFormat}) 
Plots the error at each iteration of an iterative solver.

Collective

Input Parameters:
- `adaptor` - the `DMAdaptor`
- `n`       - iteration number
- `odm`     - the original `DM`
- `adm`     - the adapted `DM`
- `Nf`      - number of fields
- `enorms`  - 2-norm error values for each field (may be estimated).
- `error`   - `Vec` of cellwise errors
- `vf`      - The viewer context

Options Database Key:
- `-adaptor_monitor_error draw` - Activates `DMAdaptorMonitorErrorDraw()`

Level: intermediate

See also: `PETSCVIEWERDRAW`, `DMAdaptor`, `DMAdaptorMonitorSet()`, `DMAdaptorMonitorErrorDrawLG()`

# External Links
$(_doc_external("DM/DMAdaptorMonitorErrorDraw"))
"""
function DMAdaptorMonitorErrorDraw(petsclib::PetscLibType, adaptor::DMAdaptor, n::Integer, odm::AbstractPetscDM, adm::AbstractPetscDM, Nf::Integer, enorms::AbstractVector{<:Number}, error::AbstractPetscVec, vf::Vector{PetscViewerAndFormat})
    error("DMAdaptorMonitorErrorDraw: no generated method for these argument types")
end

@for_petsc function DMAdaptorMonitorErrorDraw(petsclib::$UnionPetscLib, adaptor::DMAdaptor, n::$PetscInt, odm::AbstractPetscDM, adm::AbstractPetscDM, Nf::$PetscInt, enorms::Vector{$PetscReal}, error::AbstractPetscVec, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:DMAdaptorMonitorErrorDraw, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, $PetscInt, CDM, CDM, $PetscInt, Ptr{$PetscReal}, CVec, Ptr{PetscViewerAndFormat}),
               adaptor, n, odm, adm, Nf, enorms, error, vf,
              )


	return nothing
end 

"""
	DMAdaptorMonitorErrorDrawLG(petsclib::PetscLibType, adaptor::DMAdaptor, n::PetscInt, odm::AbstractPetscDM, adm::AbstractPetscDM, Nf::PetscInt, enorms::Vector{PetscReal}, error::AbstractPetscVec, vf::Vector{PetscViewerAndFormat}) 
Plots the error norm at each iteration of an adaptive loop.

Collective

Input Parameters:
- `adaptor` - the `DMAdaptor`
- `n`       - iteration number
- `odm`     - the original `DM`
- `adm`     - the adapted `DM`
- `Nf`      - number of fields
- `enorms`  - 2-norm error values for each field (may be estimated).
- `error`   - `Vec` of cellwise errors
- `vf`      - The viewer context, obtained via `DMAdaptorMonitorErrorDrawLGCreate()`

Options Database Key:
- `-adaptor_error draw::draw_lg` - Activates `DMAdaptorMonitorErrorDrawLG()`

Level: intermediate

See also: `PETSCVIEWERDRAW`, `DMAdaptor`, `DMAdaptorMonitorSet()`, `DMAdaptorMonitorErrorDraw()`, `DMAdaptorMonitorError()`,
`DMAdaptorMonitorTrueResidualDrawLGCreate()`

# External Links
$(_doc_external("DM/DMAdaptorMonitorErrorDrawLG"))
"""
function DMAdaptorMonitorErrorDrawLG(petsclib::PetscLibType, adaptor::DMAdaptor, n::Integer, odm::AbstractPetscDM, adm::AbstractPetscDM, Nf::Integer, enorms::AbstractVector{<:Number}, error::AbstractPetscVec, vf::Vector{PetscViewerAndFormat})
    error("DMAdaptorMonitorErrorDrawLG: no generated method for these argument types")
end

@for_petsc function DMAdaptorMonitorErrorDrawLG(petsclib::$UnionPetscLib, adaptor::DMAdaptor, n::$PetscInt, odm::AbstractPetscDM, adm::AbstractPetscDM, Nf::$PetscInt, enorms::Vector{$PetscReal}, error::AbstractPetscVec, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:DMAdaptorMonitorErrorDrawLG, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, $PetscInt, CDM, CDM, $PetscInt, Ptr{$PetscReal}, CVec, Ptr{PetscViewerAndFormat}),
               adaptor, n, odm, adm, Nf, enorms, error, vf,
              )


	return nothing
end 

"""
	vf::Ptr{PetscViewerAndFormat} = DMAdaptorMonitorErrorDrawLGCreate(petsclib::PetscLibType, viewer::PetscViewer, format::PetscViewerFormat, ctx::Ptr{Cvoid}) 
Creates the context for the error plotter `DMAdaptorMonitorErrorDrawLG()`

Collective

Input Parameters:
- `viewer` - The `PetscViewer`
- `format` - The viewer format
- `ctx`    - An optional application context

Output Parameter:
- `vf` - The viewer context

Level: intermediate

See also: `PETSCVIEWERDRAW`, `PetscViewerMonitorGLSetUp()`, `DMAdaptor`, `DMAdaptorMonitorSet()`, `DMAdaptorMonitorErrorDrawLG()`

# External Links
$(_doc_external("DM/DMAdaptorMonitorErrorDrawLGCreate"))
"""
function DMAdaptorMonitorErrorDrawLGCreate(petsclib::PetscLibType, viewer::PetscViewer, format::PetscViewerFormat, ctx::Ptr{Cvoid})
    error("DMAdaptorMonitorErrorDrawLGCreate: no generated method for these argument types")
end

@for_petsc function DMAdaptorMonitorErrorDrawLGCreate(petsclib::$UnionPetscLib, viewer::PetscViewer, format::PetscViewerFormat, ctx::Ptr{Cvoid} )
	vf_ = Ref{Ptr{PetscViewerAndFormat}}()

    @chk ccall(
               (:DMAdaptorMonitorErrorDrawLGCreate, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscViewerFormat, Ptr{Cvoid}, Ptr{Ptr{PetscViewerAndFormat}}),
               viewer, format, ctx, vf_,
              )

	vf = vf_[]

	return vf
end 

"""
	DMAdaptorMonitorRegister(petsclib::PetscLibType, name::String, vtype::PetscViewerType, format::PetscViewerFormat, monitor::external, create::external, destroy::external) 
Registers a mesh adaptation monitor routine that may be accessed with `DMAdaptorMonitorSetFromOptions()`

Not Collective

Input Parameters:
- `name`    - name of a new monitor routine
- `vtype`   - A `PetscViewerType` for the output
- `format`  - A `PetscViewerFormat` for the output
- `monitor` - Monitor routine
- `create`  - Creation routine, or `NULL`
- `destroy` - Destruction routine, or `NULL`

Level: advanced

See also: `DMAdaptor`, `DMAdaptorMonitorSet()`, `DMAdaptorMonitorRegisterAll()`, `DMAdaptorMonitorSetFromOptions()`

# External Links
$(_doc_external("DM/DMAdaptorMonitorRegister"))
"""
function DMAdaptorMonitorRegister(petsclib::PetscLibType, name::String, vtype::PetscViewerType, format::PetscViewerFormat, monitor::external, create::external, destroy::external)
    error("DMAdaptorMonitorRegister: no generated method for these argument types")
end

@for_petsc function DMAdaptorMonitorRegister(petsclib::$UnionPetscLib, name::String, vtype::PetscViewerType, format::PetscViewerFormat, monitor::external, create::external, destroy::external )

    @chk ccall(
               (:DMAdaptorMonitorRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, PetscViewerType, PetscViewerFormat, external, external, external),
               name, vtype, format, monitor, create, destroy,
              )


	return nothing
end 

"""
	DMAdaptorMonitorRegisterAll(petsclib::PetscLibType) 
Registers all of the mesh adaptation monitors in the `SNES` package.

Not Collective

Level: advanced

See also: `SNES`, `DM`, `DMAdaptorMonitorRegister()`, `DMAdaptorRegister()`

# External Links
$(_doc_external("DM/DMAdaptorMonitorRegisterAll"))
"""
function DMAdaptorMonitorRegisterAll(petsclib::PetscLibType)
    error("DMAdaptorMonitorRegisterAll: no generated method for these argument types")
end

@for_petsc function DMAdaptorMonitorRegisterAll(petsclib::$UnionPetscLib)

    @chk ccall(
               (:DMAdaptorMonitorRegisterAll, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	DMAdaptorMonitorRegisterDestroy(petsclib::PetscLibType) 
This function destroys the registered monitors for `DMAdaptor`. It is called from `PetscFinalize()`.

Not collective

Level: developer

See also: `DM`, `DMPLEX`, `DMAdaptorMonitorRegisterAll()`, `DMAdaptor`, `PetscFinalize()`

# External Links
$(_doc_external("DM/DMAdaptorMonitorRegisterDestroy"))
"""
function DMAdaptorMonitorRegisterDestroy(petsclib::PetscLibType)
    error("DMAdaptorMonitorRegisterDestroy: no generated method for these argument types")
end

@for_petsc function DMAdaptorMonitorRegisterDestroy(petsclib::$UnionPetscLib)

    @chk ccall(
               (:DMAdaptorMonitorRegisterDestroy, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	DMAdaptorMonitorSet(petsclib::PetscLibType, adaptor::DMAdaptor, monitor::external, ctx::Ptr{Cvoid}, monitordestroy::Ptr{Cvoid}) 
Sets an ADDITIONAL function to be called at every iteration to monitor
the error etc.

Logically Collective

Input Parameters:
- `adaptor`        - the `DMAdaptor`
- `monitor`        - pointer to function (if this is `NULL`, it turns off monitoring
- `ctx`            - [optional] context for private data for the monitor routine (use `NULL` if no context is needed)
- `monitordestroy` - [optional] routine that frees monitor context (may be `NULL`), see `PetscCtxDestroyFn` for its calling sequence

Calling sequence of `monitor`:
- `adaptor` - the `DMAdaptor`
- `it`      - iteration number
- `odm`     - the original `DM`
- `adm`     - the adapted `DM`
- `Nf`      - number of fields
- `enorms`  - (estimated) 2-norm of the error for each field
- `error`   - `Vec` of cellwise errors
- `ctx`     - optional monitoring context, as set by `DMAdaptorMonitorSet()`

Options Database Keys:
- `-adaptor_monitor_size`                - sets `DMAdaptorMonitorSize()`
- `-adaptor_monitor_error`               - sets `DMAdaptorMonitorError()`
- `-adaptor_monitor_error draw`          - sets `DMAdaptorMonitorErrorDraw()` and plots error
- `-adaptor_monitor_error draw::draw_lg` - sets `DMAdaptorMonitorErrorDrawLG()` and plots error
- `-dm_adaptor_monitor_cancel`           - Cancels all monitors that have been hardwired into a code by calls to `DMAdaptorMonitorSet()`, but does not cancel those set via the options database.

Level: beginner

See also: `DMAdaptorMonitorError()`, `DMAdaptor`, `PetscCtxDestroyFn`

# External Links
$(_doc_external("DM/DMAdaptorMonitorSet"))
"""
function DMAdaptorMonitorSet(petsclib::PetscLibType, adaptor::DMAdaptor, monitor::external, ctx::Ptr{Cvoid}, monitordestroy::Ptr{Cvoid})
    error("DMAdaptorMonitorSet: no generated method for these argument types")
end

@for_petsc function DMAdaptorMonitorSet(petsclib::$UnionPetscLib, adaptor::DMAdaptor, monitor::external, ctx::Ptr{Cvoid}, monitordestroy::Ptr{Cvoid} )

    @chk ccall(
               (:DMAdaptorMonitorSet, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, external, Ptr{Cvoid}, Ptr{Cvoid}),
               adaptor, monitor, ctx, monitordestroy,
              )


	return nothing
end 

"""
	DMAdaptorMonitorSetFromOptions(petsclib::PetscLibType, adaptor::DMAdaptor, opt::String, name::String, ctx::Ptr{Cvoid}) 
Sets a monitor function and viewer appropriate for the type indicated by the user in the options database

Collective

Input Parameters:
- `adaptor` - `DMadaptor` object you wish to monitor
- `opt`     - the command line option for this monitor
- `name`    - the monitor type one is seeking
- `ctx`     - An optional application context for the monitor, or `NULL`

Level: developer

See also: `DMAdaptorMonitorRegister()`, `DMAdaptorMonitorSet()`, `PetscOptionsGetViewer()`

# External Links
$(_doc_external("DM/DMAdaptorMonitorSetFromOptions"))
"""
function DMAdaptorMonitorSetFromOptions(petsclib::PetscLibType, adaptor::DMAdaptor, opt::String, name::String, ctx::Ptr{Cvoid})
    error("DMAdaptorMonitorSetFromOptions: no generated method for these argument types")
end

@for_petsc function DMAdaptorMonitorSetFromOptions(petsclib::$UnionPetscLib, adaptor::DMAdaptor, opt::String, name::String, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:DMAdaptorMonitorSetFromOptions, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cvoid}),
               adaptor, opt, name, ctx,
              )


	return nothing
end 

"""
	DMAdaptorMonitorSize(petsclib::PetscLibType, adaptor::DMAdaptor, n::PetscInt, odm::AbstractPetscDM, adm::AbstractPetscDM, Nf::PetscInt, enorms::Vector{PetscReal}, error::AbstractPetscVec, vf::Vector{PetscViewerAndFormat}) 
Prints the mesh sizes at each iteration of an adaptation loop.

Collective

Input Parameters:
- `adaptor` - the `DMAdaptor`
- `n`       - iteration number
- `odm`     - the original `DM`
- `adm`     - the adapted `DM`
- `Nf`      - number of fields
- `enorms`  - 2-norm error values for each field (may be estimated).
- `error`   - `Vec` of cellwise errors
- `vf`      - The viewer context

Options Database Key:
- `-adaptor_monitor_size` - Activates `DMAdaptorMonitorSize()`

Level: intermediate

See also: `DMAdaptor`, `DMAdaptorMonitorSet()`, `DMAdaptorMonitorError()`, `DMAdaptorMonitorErrorDraw()`, `DMAdaptorMonitorErrorDrawLG()`

# External Links
$(_doc_external("DM/DMAdaptorMonitorSize"))
"""
function DMAdaptorMonitorSize(petsclib::PetscLibType, adaptor::DMAdaptor, n::Integer, odm::AbstractPetscDM, adm::AbstractPetscDM, Nf::Integer, enorms::AbstractVector{<:Number}, error::AbstractPetscVec, vf::Vector{PetscViewerAndFormat})
    error("DMAdaptorMonitorSize: no generated method for these argument types")
end

@for_petsc function DMAdaptorMonitorSize(petsclib::$UnionPetscLib, adaptor::DMAdaptor, n::$PetscInt, odm::AbstractPetscDM, adm::AbstractPetscDM, Nf::$PetscInt, enorms::Vector{$PetscReal}, error::AbstractPetscVec, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:DMAdaptorMonitorSize, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, $PetscInt, CDM, CDM, $PetscInt, Ptr{$PetscReal}, CVec, Ptr{PetscViewerAndFormat}),
               adaptor, n, odm, adm, Nf, enorms, error, vf,
              )


	return nothing
end 

"""
	DMAdaptorRegister(petsclib::PetscLibType, name::String, noname::Ptr{Cvoid}) 
Adds a new adaptor component implementation

Not Collective

Input Parameters:
- `name`        - The name of a new user-defined creation routine
- `create_func` - The creation routine

See also: `DM`, `DMPLEX`, `DMAdaptor`, `DMAdaptorRegisterAll()`, `DMAdaptorRegisterDestroy()`

# External Links
$(_doc_external("DM/DMAdaptorRegister"))
"""
function DMAdaptorRegister(petsclib::PetscLibType, name::String, noname::Ptr{Cvoid})
    error("DMAdaptorRegister: no generated method for these argument types")
end

@for_petsc function DMAdaptorRegister(petsclib::$UnionPetscLib, name::String, noname::Ptr{Cvoid} )

    @chk ccall(
               (:DMAdaptorRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cvoid}),
               name, noname,
              )


	return nothing
end 

"""
	DMAdaptorRegisterAll(petsclib::PetscLibType) 
Registers all of the adaptor components in the `DM` package.

Not Collective

Level: advanced

See also: `DM`, `DMPLEX`, `DMAdaptorType`, `DMRegisterAll()`, `DMAdaptorRegisterDestroy()`

# External Links
$(_doc_external("DM/DMAdaptorRegisterAll"))
"""
function DMAdaptorRegisterAll(petsclib::PetscLibType)
    error("DMAdaptorRegisterAll: no generated method for these argument types")
end

@for_petsc function DMAdaptorRegisterAll(petsclib::$UnionPetscLib)

    @chk ccall(
               (:DMAdaptorRegisterAll, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	DMAdaptorRegisterDestroy(petsclib::PetscLibType) 
This function destroys the registered `DMAdaptorType`. It is called from `PetscFinalize()`.

Not collective

Level: developer

See also: `DM`, `DMPLEX`, `DMAdaptorRegisterAll()`, `DMAdaptorType`, `PetscFinalize()`

# External Links
$(_doc_external("DM/DMAdaptorRegisterDestroy"))
"""
function DMAdaptorRegisterDestroy(petsclib::PetscLibType)
    error("DMAdaptorRegisterDestroy: no generated method for these argument types")
end

@for_petsc function DMAdaptorRegisterDestroy(petsclib::$UnionPetscLib)

    @chk ccall(
               (:DMAdaptorRegisterDestroy, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	DMAdaptorSetCriterion(petsclib::PetscLibType, adaptor::DMAdaptor, criterion::DMAdaptationCriterion) 
Set the adaptation criterion

Not Collective

Input Parameters:
- `adaptor`   - the `DMAdaptor`
- `criterion` - the adaptation criterion

Level: advanced

See also: `DMAdaptor`, `DMAdaptorGetCriterion()`, `DMAdaptationCriterion`

# External Links
$(_doc_external("DM/DMAdaptorSetCriterion"))
"""
function DMAdaptorSetCriterion(petsclib::PetscLibType, adaptor::DMAdaptor, criterion::DMAdaptationCriterion)
    error("DMAdaptorSetCriterion: no generated method for these argument types")
end

@for_petsc function DMAdaptorSetCriterion(petsclib::$UnionPetscLib, adaptor::DMAdaptor, criterion::DMAdaptationCriterion )

    @chk ccall(
               (:DMAdaptorSetCriterion, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, DMAdaptationCriterion),
               adaptor, criterion,
              )


	return nothing
end 

"""
	DMAdaptorSetFromOptions(petsclib::PetscLibType, adaptor::DMAdaptor) 
Sets properties of a `DMAdaptor` object from values in the options database

Collective

Input Parameter:
- `adaptor` - The `DMAdaptor` object

Options Database Keys:
- `-adaptor_monitor_size`              - Monitor the mesh size
- `-adaptor_monitor_error`             - Monitor the solution error
- `-adaptor_sequence_num num`          - Number of adaptations to generate an optimal grid
- `-adaptor_target_num num`            - Set the target number of vertices N_adapt, -1 for automatic determination
- `-adaptor_refinement_factor r`       - Set r such that N_adapt = r^dim N_orig
- `-adaptor_mixed_setup_function func` - Set the function func that sets up the mixed problem

Level: beginner

See also: `DM`, `DMAdaptor`, `DMAdaptorCreate()`, `DMAdaptorAdapt()`

# External Links
$(_doc_external("DM/DMAdaptorSetFromOptions"))
"""
function DMAdaptorSetFromOptions(petsclib::PetscLibType, adaptor::DMAdaptor)
    error("DMAdaptorSetFromOptions: no generated method for these argument types")
end

@for_petsc function DMAdaptorSetFromOptions(petsclib::$UnionPetscLib, adaptor::DMAdaptor )

    @chk ccall(
               (:DMAdaptorSetFromOptions, $petsc_library),
               PetscErrorCode,
               (DMAdaptor,),
               adaptor,
              )


	return nothing
end 

"""
	DMAdaptorSetMixedSetupFunction(petsclib::PetscLibType, adaptor::DMAdaptor, setupFunc::external) 
Set the function setting up the mixed problem

Not Collective

Input Parameters:
- `adaptor`   - the `DMAdaptor`
- `setupFunc` - the function setting up the mixed problem

Calling sequence of setupFunc:
- `adaptor` - the `DMAdaptor`
- `dm`      - the `DM`

Level: advanced

See also: `DMAdaptor`, `DMAdaptorGetMixedSetupFunction()`, `DMAdaptorAdapt()`

# External Links
$(_doc_external("DM/DMAdaptorSetMixedSetupFunction"))
"""
function DMAdaptorSetMixedSetupFunction(petsclib::PetscLibType, adaptor::DMAdaptor, setupFunc::external)
    error("DMAdaptorSetMixedSetupFunction: no generated method for these argument types")
end

@for_petsc function DMAdaptorSetMixedSetupFunction(petsclib::$UnionPetscLib, adaptor::DMAdaptor, setupFunc::external )

    @chk ccall(
               (:DMAdaptorSetMixedSetupFunction, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, external),
               adaptor, setupFunc,
              )


	return nothing
end 

"""
	DMAdaptorSetOptionsPrefix(petsclib::PetscLibType, adaptor::DMAdaptor, prefix::String) 
Sets the prefix used for searching for all `DMAdaptor` options in the database.

Logically Collective

Input Parameters:
- `adaptor` - the `DMAdaptor`
- `prefix`  - the prefix to prepend to all option names

Level: advanced

See also: `DMAdaptor`, `SNESSetOptionsPrefix()`, `DMAdaptorSetFromOptions()`

# External Links
$(_doc_external("DM/DMAdaptorSetOptionsPrefix"))
"""
function DMAdaptorSetOptionsPrefix(petsclib::PetscLibType, adaptor::DMAdaptor, prefix::String)
    error("DMAdaptorSetOptionsPrefix: no generated method for these argument types")
end

@for_petsc function DMAdaptorSetOptionsPrefix(petsclib::$UnionPetscLib, adaptor::DMAdaptor, prefix::String )

    @chk ccall(
               (:DMAdaptorSetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, Ptr{Cchar}),
               adaptor, prefix,
              )


	return nothing
end 

"""
	DMAdaptorSetSequenceLength(petsclib::PetscLibType, adaptor::DMAdaptor, num::PetscInt) 
Sets the number of sequential adaptations

Not Collective

Input Parameters:
- `adaptor` - The `DMAdaptor` object
- `num`     - The number of adaptations

Level: intermediate

See also: `DMAdaptorGetSequenceLength()`, `DMAdaptorCreate()`, `DMAdaptorAdapt()`

# External Links
$(_doc_external("DM/DMAdaptorSetSequenceLength"))
"""
function DMAdaptorSetSequenceLength(petsclib::PetscLibType, adaptor::DMAdaptor, num::Integer)
    error("DMAdaptorSetSequenceLength: no generated method for these argument types")
end

@for_petsc function DMAdaptorSetSequenceLength(petsclib::$UnionPetscLib, adaptor::DMAdaptor, num::$PetscInt )

    @chk ccall(
               (:DMAdaptorSetSequenceLength, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, $PetscInt),
               adaptor, num,
              )


	return nothing
end 

"""
	DMAdaptorSetSolver(petsclib::PetscLibType, adaptor::DMAdaptor, snes::AbstractSNES) 
Sets the solver used to produce discrete solutions

Not Collective

Input Parameters:
- `adaptor` - The `DMAdaptor` object
- `snes`    - The solver, this MUST have an attached `DM`/`PetscDS`, so that the exact solution can be computed

Level: intermediate

See also: `DMAdaptor`, `DMAdaptorGetSolver()`, `DMAdaptorCreate()`, `DMAdaptorAdapt()`

# External Links
$(_doc_external("DM/DMAdaptorSetSolver"))
"""
function DMAdaptorSetSolver(petsclib::PetscLibType, adaptor::DMAdaptor, snes::AbstractSNES)
    error("DMAdaptorSetSolver: no generated method for these argument types")
end

@for_petsc function DMAdaptorSetSolver(petsclib::$UnionPetscLib, adaptor::DMAdaptor, snes::AbstractSNES )

    @chk ccall(
               (:DMAdaptorSetSolver, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, CSNES),
               adaptor, snes,
              )


	return nothing
end 

"""
	DMAdaptorSetTransferFunction(petsclib::PetscLibType, adaptor::DMAdaptor, tfunc::external) 
Set the callback used by a `DMAdaptor` to transfer a solution vector from an old `DM` to the adapted `DM`

Logically Collective

Input Parameters:
- `adaptor` - the `DMAdaptor` object
- `tfunc`   - the transfer callback

Calling sequence of `tfunc`:
- `adaptor` - the `DMAdaptor` object
- `dm`      - the current `DM`
- `xin`     - the current solution
- `newdm`   - the adapted `DM`
- `xout`    - the transferred solution on `newdm`
- `ctx`     - application context, set with `DMSetApplicationContext()`

Level: developer

See also: `DMAdaptor`, `DMAdaptorGetTransferFunction()`, `DMAdaptorAdapt()`

# External Links
$(_doc_external("DM/DMAdaptorSetTransferFunction"))
"""
function DMAdaptorSetTransferFunction(petsclib::PetscLibType, adaptor::DMAdaptor, tfunc::external)
    error("DMAdaptorSetTransferFunction: no generated method for these argument types")
end

@for_petsc function DMAdaptorSetTransferFunction(petsclib::$UnionPetscLib, adaptor::DMAdaptor, tfunc::external )

    @chk ccall(
               (:DMAdaptorSetTransferFunction, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, external),
               adaptor, tfunc,
              )


	return nothing
end 

"""
	DMAdaptorSetType(petsclib::PetscLibType, adaptor::DMAdaptor, method::DMAdaptorType) 
Sets the particular implementation for a adaptor.

Collective

Input Parameters:
- `adaptor` - The `DMAdaptor`
- `method`  - The name of the adaptor type

Options Database Key:
- `-adaptor_type type` - Sets the adaptor type; see `DMAdaptorType`

Level: intermediate

See also: `DM`, `DMPLEX`, `DMAdaptor`, `DMAdaptorType`, `DMAdaptorGetType()`, `DMAdaptorCreate()`

# External Links
$(_doc_external("DM/DMAdaptorSetType"))
"""
function DMAdaptorSetType(petsclib::PetscLibType, adaptor::DMAdaptor, method::DMAdaptorType)
    error("DMAdaptorSetType: no generated method for these argument types")
end

@for_petsc function DMAdaptorSetType(petsclib::$UnionPetscLib, adaptor::DMAdaptor, method::DMAdaptorType )

    @chk ccall(
               (:DMAdaptorSetType, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, DMAdaptorType),
               adaptor, method,
              )


	return nothing
end 

"""
	DMAdaptorSetUp(petsclib::PetscLibType, adaptor::DMAdaptor) 
After the solver is specified, creates data structures for controlling adaptivity

Collective

Input Parameter:
- `adaptor` - The `DMAdaptor` object

Level: beginner

See also: `DMAdaptor`, `DMAdaptorCreate()`, `DMAdaptorAdapt()`

# External Links
$(_doc_external("DM/DMAdaptorSetUp"))
"""
function DMAdaptorSetUp(petsclib::PetscLibType, adaptor::DMAdaptor)
    error("DMAdaptorSetUp: no generated method for these argument types")
end

@for_petsc function DMAdaptorSetUp(petsclib::$UnionPetscLib, adaptor::DMAdaptor )

    @chk ccall(
               (:DMAdaptorSetUp, $petsc_library),
               PetscErrorCode,
               (DMAdaptor,),
               adaptor,
              )


	return nothing
end 

"""
	DMAdaptorView(petsclib::PetscLibType, adaptor::DMAdaptor, viewer::PetscViewer) 
Views a `DMAdaptor` object

Collective

Input Parameters:
- `adaptor` - The `DMAdaptor` object
- `viewer`  - The `PetscViewer` object

Level: beginner

See also: `DM`, `DMAdaptor`, `DMAdaptorCreate()`, `DMAdaptorAdapt()`

# External Links
$(_doc_external("DM/DMAdaptorView"))
"""
function DMAdaptorView(petsclib::PetscLibType, adaptor::DMAdaptor, viewer::PetscViewer)
    error("DMAdaptorView: no generated method for these argument types")
end

@for_petsc function DMAdaptorView(petsclib::$UnionPetscLib, adaptor::DMAdaptor, viewer::PetscViewer )

    @chk ccall(
               (:DMAdaptorView, $petsc_library),
               PetscErrorCode,
               (DMAdaptor, PetscViewer),
               adaptor, viewer,
              )


	return nothing
end 

"""
	field::DMField = DMFieldCreateDA(petsclib::PetscLibType, dm::AbstractPetscDM, nc::PetscInt, cornerValues::Vector{PetscScalar}) 
Create a `DMField` of type `DMFIELDDA` that represents a multilinear field on a `DMDA` given by its values at the corners of the reference element.

Collective

Input Parameters:
- `dm`           - the `DMDA` on which the field lives
- `nc`           - the number of components of the field
- `cornerValues` - array of length `nc * (1 << dim)` holding the field values at each corner of the reference element, ordered by lexicographic corner index

Output Parameter:
- `field` - the newly created `DMField`

Level: intermediate

See also: `DMField`, `DMFIELDDA`, `DMDA`, `DMFieldCreate()`, `DMFieldCreateDS()`, `DMFieldCreateShell()`

# External Links
$(_doc_external("DM/DMFieldCreateDA"))
"""
function DMFieldCreateDA(petsclib::PetscLibType, dm::AbstractPetscDM, nc::Integer, cornerValues::AbstractVector{<:Number})
    error("DMFieldCreateDA: no generated method for these argument types")
end

@for_petsc function DMFieldCreateDA(petsclib::$UnionPetscLib, dm::AbstractPetscDM, nc::$PetscInt, cornerValues::Vector{$PetscScalar} )
	field_ = Ref{DMField}()

    @chk ccall(
               (:DMFieldCreateDA, $petsc_library),
               PetscErrorCode,
               (CDM, $PetscInt, Ptr{$PetscScalar}, Ptr{DMField}),
               dm, nc, cornerValues, field_,
              )

	field = field_[]

	return field
end 

"""
	field::DMField = DMFieldCreateDS(petsclib::PetscLibType, dm::AbstractPetscDM, fieldNum::PetscInt, vec::AbstractPetscVec) 
Create a `DMField` of type `DMFIELDDS` for a `PetscDS` field on a `DM`.

Collective

Input Parameters:
- `dm`       - the `DM` carrying the discretization
- `fieldNum` - the field number within the `DM`'s `PetscDS`
- `vec`      - local vector holding the coefficients

Output Parameter:
- `field` - the newly created `DMField`

Level: intermediate

See also: `DMField`, `DMFIELDDS`, `DMFieldCreateDSWithDG()`, `DMFieldCreate()`, `PetscDS`

# External Links
$(_doc_external("DM/DMFieldCreateDS"))
"""
function DMFieldCreateDS(petsclib::PetscLibType, dm::AbstractPetscDM, fieldNum::Integer, vec::AbstractPetscVec)
    error("DMFieldCreateDS: no generated method for these argument types")
end

@for_petsc function DMFieldCreateDS(petsclib::$UnionPetscLib, dm::AbstractPetscDM, fieldNum::$PetscInt, vec::AbstractPetscVec )
	field_ = Ref{DMField}()

    @chk ccall(
               (:DMFieldCreateDS, $petsc_library),
               PetscErrorCode,
               (CDM, $PetscInt, CVec, Ptr{DMField}),
               dm, fieldNum, vec, field_,
              )

	field = field_[]

	return field
end 

"""
	field::DMField = DMFieldCreateDSWithDG(petsclib::PetscLibType, dm::AbstractPetscDM, dmDG::AbstractPetscDM, fieldNum::PetscInt, vec::AbstractPetscVec, vecDG::AbstractPetscVec) 
Create a `DMField` of type `DMFIELDDS` for a `PetscDS` field, optionally paired with a matching discontinuous-Galerkin representation on a companion `DM`.

Collective

Input Parameters:
- `dm`       - the `DM` carrying the primary (continuous) discretization
- `dmDG`     - optional `DM` carrying a matching discontinuous-Galerkin discretization, or `NULL`
- `fieldNum` - the field number within the `DM`'s `PetscDS`
- `vec`      - local vector holding the coefficients on `dm`
- `vecDG`    - local vector holding the coefficients on `dmDG`, or `NULL` if `dmDG` is `NULL`

Output Parameter:
- `field` - the newly created `DMField`

Level: intermediate

See also: `DMField`, `DMFIELDDS`, `DMFieldCreateDS()`, `DMFieldCreate()`, `PetscDS`, `PetscFE`

# External Links
$(_doc_external("DM/DMFieldCreateDSWithDG"))
"""
function DMFieldCreateDSWithDG(petsclib::PetscLibType, dm::AbstractPetscDM, dmDG::AbstractPetscDM, fieldNum::Integer, vec::AbstractPetscVec, vecDG::AbstractPetscVec)
    error("DMFieldCreateDSWithDG: no generated method for these argument types")
end

@for_petsc function DMFieldCreateDSWithDG(petsclib::$UnionPetscLib, dm::AbstractPetscDM, dmDG::AbstractPetscDM, fieldNum::$PetscInt, vec::AbstractPetscVec, vecDG::AbstractPetscVec )
	field_ = Ref{DMField}()

    @chk ccall(
               (:DMFieldCreateDSWithDG, $petsc_library),
               PetscErrorCode,
               (CDM, CDM, $PetscInt, CVec, CVec, Ptr{DMField}),
               dm, dmDG, fieldNum, vec, vecDG, field_,
              )

	field = field_[]

	return field
end 

"""
	quad::PetscQuadrature = DMFieldCreateDefaultFaceQuadrature(petsclib::PetscLibType, field::DMField, pointIS::AbstractIS) 
Creates a quadrature sufficient to integrate the field on all faces of the selected cells via pullback onto the reference element

Not Collective

Input Parameters:
- `field`   - the `DMField` object
- `pointIS` - the index set of points over which we wish to integrate the field over faces

Output Parameter:
- `quad` - a `PetscQuadrature` object

Level: developer

See also: `DMFieldCreateDefaultQuadrature()`, `DMField`, `PetscQuadrature`, `IS`, `DMFieldEvaluteFE()`, `DMFieldGetDegree()`

# External Links
$(_doc_external("DM/DMFieldCreateDefaultFaceQuadrature"))
"""
function DMFieldCreateDefaultFaceQuadrature(petsclib::PetscLibType, field::DMField, pointIS::AbstractIS)
    error("DMFieldCreateDefaultFaceQuadrature: no generated method for these argument types")
end

@for_petsc function DMFieldCreateDefaultFaceQuadrature(petsclib::$UnionPetscLib, field::DMField, pointIS::AbstractIS )
	quad_ = Ref{PetscQuadrature}()

    @chk ccall(
               (:DMFieldCreateDefaultFaceQuadrature, $petsc_library),
               PetscErrorCode,
               (DMField, CIS, Ptr{PetscQuadrature}),
               field, pointIS, quad_,
              )

	quad = quad_[]

	return quad
end 

"""
	quad::PetscQuadrature = DMFieldCreateDefaultQuadrature(petsclib::PetscLibType, field::DMField, pointIS::AbstractIS) 
Creates a quadrature sufficient to integrate the field on the selected
points via pullback onto the reference element

Not Collective

Input Parameters:
- `field`   - the `DMField` object
- `pointIS` - the index set of points over which we wish to integrate the field

Output Parameter:
- `quad` - a `PetscQuadrature` object

Level: developer

See also: `DMFieldCreateDefaultFaceQuadrature()`, `DMField`, `PetscQuadrature`, `IS`, `DMFieldEvaluteFE()`, `DMFieldGetDegree()`

# External Links
$(_doc_external("DM/DMFieldCreateDefaultQuadrature"))
"""
function DMFieldCreateDefaultQuadrature(petsclib::PetscLibType, field::DMField, pointIS::AbstractIS)
    error("DMFieldCreateDefaultQuadrature: no generated method for these argument types")
end

@for_petsc function DMFieldCreateDefaultQuadrature(petsclib::$UnionPetscLib, field::DMField, pointIS::AbstractIS )
	quad_ = Ref{PetscQuadrature}()

    @chk ccall(
               (:DMFieldCreateDefaultQuadrature, $petsc_library),
               PetscErrorCode,
               (DMField, CIS, Ptr{PetscQuadrature}),
               field, pointIS, quad_,
              )

	quad = quad_[]

	return quad
end 

"""
	geom::Ptr{PetscFEGeom} = DMFieldCreateFEGeom(petsclib::PetscLibType, field::DMField, pointIS::AbstractIS, quad::PetscQuadrature, mode::PetscFEGeomMode) 
Compute and create the geometric factors of a coordinate field

Not Collective

Input Parameters:
- `field`   - the `DMField` object
- `pointIS` - the index set of points over which we wish to integrate the field
- `quad`    - the quadrature points at which to evaluate the geometric factors
- `mode`    - Type of geometry data to store

Output Parameter:
- `geom` - the geometric factors

Level: developer

See also: `DMField`, `PetscQuadrature`, `IS`, `PetscFEGeom`, `DMFieldEvaluateFE()`, `DMFieldCreateDefaulteQuadrature()`, `DMFieldGetDegree()`

# External Links
$(_doc_external("DM/DMFieldCreateFEGeom"))
"""
function DMFieldCreateFEGeom(petsclib::PetscLibType, field::DMField, pointIS::AbstractIS, quad::PetscQuadrature, mode::PetscFEGeomMode)
    error("DMFieldCreateFEGeom: no generated method for these argument types")
end

@for_petsc function DMFieldCreateFEGeom(petsclib::$UnionPetscLib, field::DMField, pointIS::AbstractIS, quad::PetscQuadrature, mode::PetscFEGeomMode )
	geom_ = Ref{Ptr{PetscFEGeom}}()

    @chk ccall(
               (:DMFieldCreateFEGeom, $petsc_library),
               PetscErrorCode,
               (DMField, CIS, PetscQuadrature, PetscFEGeomMode, Ptr{Ptr{PetscFEGeom}}),
               field, pointIS, quad, mode, geom_,
              )

	geom = geom_[]

	return geom
end 

"""
	field::DMField = DMFieldCreateShell(petsclib::PetscLibType, dm::AbstractPetscDM, numComponents::PetscInt, continuity::DMFieldContinuity, ctx::Ptr{Cvoid}) 
Create a `DMFIELDSHELL`, a `DMField` whose evaluation is implemented entirely by user-supplied callbacks.

Collective

Input Parameters:
- `dm`            - the `DM` on which the field lives
- `numComponents` - the number of components of the field
- `continuity`    - the continuity of the field (e.g. `DMFIELD_VERTEX`)
- `ctx`           - optional application context returned by `DMFieldShellGetContext()`

Output Parameter:
- `field` - the newly created `DMField` of type `DMFIELDSHELL`

Level: intermediate

See also: `DMField`, `DMFIELDSHELL`, `DMFieldShellGetContext()`, `DMFieldShellSetEvaluate()`, `DMFieldShellSetEvaluateFE()`, `DMFieldShellSetEvaluateFV()`, `DMFieldShellSetDestroy()`

# External Links
$(_doc_external("DM/DMFieldCreateShell"))
"""
function DMFieldCreateShell(petsclib::PetscLibType, dm::AbstractPetscDM, numComponents::Integer, continuity::DMFieldContinuity, ctx::Ptr{Cvoid})
    error("DMFieldCreateShell: no generated method for these argument types")
end

@for_petsc function DMFieldCreateShell(petsclib::$UnionPetscLib, dm::AbstractPetscDM, numComponents::$PetscInt, continuity::DMFieldContinuity, ctx::Ptr{Cvoid} )
	field_ = Ref{DMField}()

    @chk ccall(
               (:DMFieldCreateShell, $petsc_library),
               PetscErrorCode,
               (CDM, $PetscInt, DMFieldContinuity, Ptr{Cvoid}, Ptr{DMField}),
               dm, numComponents, continuity, ctx, field_,
              )

	field = field_[]

	return field
end 

"""
	DMFieldDestroy(petsclib::PetscLibType, field::Union{DMField, Ref{DMField}}) 
destroy a `DMField`

Collective

Input Parameter:
- `field` - address of `DMField`

Level: advanced

See also: `DMField`, `DMFieldCreate()`

# External Links
$(_doc_external("DM/DMFieldDestroy"))
"""
function DMFieldDestroy(petsclib::PetscLibType, field::Union{DMField, Ref{DMField}})
    error("DMFieldDestroy: no generated method for these argument types")
end

@for_petsc function DMFieldDestroy(petsclib::$UnionPetscLib, field::Union{DMField, Ref{DMField}} )
	field_ = field isa Base.RefValue ? field : Ref{DMField}(field)

    @chk ccall(
               (:DMFieldDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{DMField},),
               field_,
              )


	return nothing
end 

"""
	B::Ptr{Cvoid},D::Ptr{Cvoid},H::Ptr{Cvoid} = DMFieldEvaluate(petsclib::PetscLibType, field::DMField, points::AbstractPetscVec, datatype::PetscDataType) 
Evaluate the field and its derivatives on a set of points

Collective

Input Parameters:
- `field`    - The `DMField` object
- `points`   - The points at which to evaluate the field.  Should have size d x n,
where d is the coordinate dimension of the manifold and n is the number
of points
- `datatype` - The PetscDataType of the output arrays: either `PETSC_REAL` or `PETSC_SCALAR`.
If the field is complex and datatype is `PETSC_REAL`, the real part of the
field is returned.

Output Parameters:
- `B` - pointer to data of size c * n * sizeof(datatype), where c is the number of components in the field.
If B is not `NULL`, the values of the field are written in this array, varying first by component,
then by point.
- `D` - pointer to data of size d * c * n * sizeof(datatype).
If `D` is not `NULL`, the values of the field's spatial derivatives are written in this array,
varying first by the partial derivative component, then by field component, then by point.
- `H` - pointer to data of size d * d * c * n * sizeof(datatype).
If `H` is not `NULL`, the values of the field's second spatial derivatives are written in this array,
varying first by the second partial derivative component, then by field component, then by point.

Level: intermediate

See also: `DMField`, `DMFieldGetDM()`, `DMFieldGetNumComponents()`, `DMFieldEvaluateFE()`, `DMFieldEvaluateFV()`, `PetscDataType`

# External Links
$(_doc_external("DM/DMFieldEvaluate"))
"""
function DMFieldEvaluate(petsclib::PetscLibType, field::DMField, points::AbstractPetscVec, datatype::PetscDataType)
    error("DMFieldEvaluate: no generated method for these argument types")
end

@for_petsc function DMFieldEvaluate(petsclib::$UnionPetscLib, field::DMField, points::AbstractPetscVec, datatype::PetscDataType )
	B_ = Ref{Ptr{Cvoid}}()
	D_ = Ref{Ptr{Cvoid}}()
	H_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:DMFieldEvaluate, $petsc_library),
               PetscErrorCode,
               (DMField, CVec, PetscDataType, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
               field, points, datatype, B_, D_, H_,
              )

	B = B_[]
	D = D_[]
	H = H_[]

	return B,D,H
end 

"""
	B::Ptr{Cvoid},D::Ptr{Cvoid},H::Ptr{Cvoid} = DMFieldEvaluateFE(petsclib::PetscLibType, field::DMField, cellIS::AbstractIS, points::PetscQuadrature, datatype::PetscDataType) 
Evaluate the field and its derivatives on a set of points mapped from
quadrature points on a reference point.  The derivatives are taken with respect to the
reference coordinates.

Not Collective

Input Parameters:
- `field`    - The `DMField` object
- `cellIS`   - Index set for cells on which to evaluate the field
- `points`   - The quadature containing the points in the reference cell at which to evaluate the field.
- `datatype` - The PetscDataType of the output arrays: either `PETSC_REAL` or `PETSC_SCALAR`.
If the field is complex and datatype is `PETSC_REAL`, the real part of the
field is returned.

Output Parameters:
- `B` - pointer to data of size c * n * sizeof(datatype), where c is the number of components in the field.
If B is not `NULL`, the values of the field are written in this array, varying first by component,
then by point.
- `D` - pointer to data of size d * c * n * sizeof(datatype).
If D is not `NULL`, the values of the field's spatial derivatives are written in this array,
varying first by the partial derivative component, then by field component, then by point.
- `H` - pointer to data of size d * d * c * n * sizeof(datatype).
If H is not `NULL`, the values of the field's second spatial derivatives are written in this array,
varying first by the second partial derivative component, then by field component, then by point.

Level: intermediate

See also: `DMField`, `DM`, `DMFieldGetNumComponents()`, `DMFieldEvaluate()`, `DMFieldEvaluateFV()`

# External Links
$(_doc_external("DM/DMFieldEvaluateFE"))
"""
function DMFieldEvaluateFE(petsclib::PetscLibType, field::DMField, cellIS::AbstractIS, points::PetscQuadrature, datatype::PetscDataType)
    error("DMFieldEvaluateFE: no generated method for these argument types")
end

@for_petsc function DMFieldEvaluateFE(petsclib::$UnionPetscLib, field::DMField, cellIS::AbstractIS, points::PetscQuadrature, datatype::PetscDataType )
	B_ = Ref{Ptr{Cvoid}}()
	D_ = Ref{Ptr{Cvoid}}()
	H_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:DMFieldEvaluateFE, $petsc_library),
               PetscErrorCode,
               (DMField, CIS, PetscQuadrature, PetscDataType, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
               field, cellIS, points, datatype, B_, D_, H_,
              )

	B = B_[]
	D = D_[]
	H = H_[]

	return B,D,H
end 

"""
	B::Ptr{Cvoid},D::Ptr{Cvoid},H::Ptr{Cvoid} = DMFieldEvaluateFV(petsclib::PetscLibType, field::DMField, cellIS::AbstractIS, datatype::PetscDataType) 
Evaluate the mean of a field and its finite volume derivatives on a set of points.

Not Collective

Input Parameters:
- `field`    - The `DMField` object
- `cellIS`   - Index set for cells on which to evaluate the field
- `datatype` - The PetscDataType of the output arrays: either `PETSC_REAL` or `PETSC_SCALAR`.
If the field is complex and datatype is `PETSC_REAL`, the real part of the
field is returned.

Output Parameters:
- `B` - pointer to data of size c * n * sizeof(datatype), where c is the number of components in the field.
If B is not `NULL`, the values of the field are written in this array, varying first by component,
then by point.
- `D` - pointer to data of size d * c * n * sizeof(datatype).
If D is not `NULL`, the values of the field's spatial derivatives are written in this array,
varying first by the partial derivative component, then by field component, then by point.
- `H` - pointer to data of size d * d * c * n * sizeof(datatype).
If H is not `NULL`, the values of the field's second spatial derivatives are written in this array,
varying first by the second partial derivative component, then by field component, then by point.

Level: intermediate

See also: `DMField`, `IS`, `DMFieldGetNumComponents()`, `DMFieldEvaluate()`, `DMFieldEvaluateFE()`, `PetscDataType`

# External Links
$(_doc_external("DM/DMFieldEvaluateFV"))
"""
function DMFieldEvaluateFV(petsclib::PetscLibType, field::DMField, cellIS::AbstractIS, datatype::PetscDataType)
    error("DMFieldEvaluateFV: no generated method for these argument types")
end

@for_petsc function DMFieldEvaluateFV(petsclib::$UnionPetscLib, field::DMField, cellIS::AbstractIS, datatype::PetscDataType )
	B_ = Ref{Ptr{Cvoid}}()
	D_ = Ref{Ptr{Cvoid}}()
	H_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:DMFieldEvaluateFV, $petsc_library),
               PetscErrorCode,
               (DMField, CIS, PetscDataType, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
               field, cellIS, datatype, B_, D_, H_,
              )

	B = B_[]
	D = D_[]
	H = H_[]

	return B,D,H
end 

"""
	DMFieldFinalizePackage(petsclib::PetscLibType) 
Finalize `DMField` package, it is called from `PetscFinalize()`

Logically Collective

Level: developer

See also: `DMFieldInitializePackage()`

# External Links
$(_doc_external("DM/DMFieldFinalizePackage"))
"""
function DMFieldFinalizePackage(petsclib::PetscLibType)
    error("DMFieldFinalizePackage: no generated method for these argument types")
end

@for_petsc function DMFieldFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:DMFieldFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	dm::PetscDM = DMFieldGetDM(petsclib::PetscLibType, field::DMField) 
Returns the `DM` for the manifold over which the field is defined.

Not Collective

Input Parameter:
- `field` - The `DMField` object

Output Parameter:
- `dm` - The `DM` object

Level: intermediate

See also: `DMField`, `DM`, `DMFieldEvaluate()`

# External Links
$(_doc_external("DM/DMFieldGetDM"))
"""
function DMFieldGetDM(petsclib::PetscLibType, field::DMField)
    error("DMFieldGetDM: no generated method for these argument types")
end

@for_petsc function DMFieldGetDM(petsclib::$UnionPetscLib, field::DMField )
	dm_ = Ref{CDM}()

    @chk ccall(
               (:DMFieldGetDM, $petsc_library),
               PetscErrorCode,
               (DMField, Ptr{CDM}),
               field, dm_,
              )

	dm = PetscDM(dm_[], petsclib)

	return dm
end 

"""
	minDegree::PetscInt,maxDegree::PetscInt = DMFieldGetDegree(petsclib::PetscLibType, field::DMField, cellIS::AbstractIS) 
Get the polynomial degree of a field when pulled back onto the
reference element

Not Collective

Input Parameters:
- `field`  - the `DMField` object
- `cellIS` - the index set of points over which we want know the invariance

Output Parameters:
- `minDegree` - the degree of the largest polynomial space contained in the field on each element
- `maxDegree` - the largest degree of the smallest polynomial space containing the field on any element

Level: intermediate

See also: `DMField`, `IS`, `DMFieldEvaluateFE()`

# External Links
$(_doc_external("DM/DMFieldGetDegree"))
"""
function DMFieldGetDegree(petsclib::PetscLibType, field::DMField, cellIS::AbstractIS)
    error("DMFieldGetDegree: no generated method for these argument types")
end

@for_petsc function DMFieldGetDegree(petsclib::$UnionPetscLib, field::DMField, cellIS::AbstractIS )
	minDegree_ = Ref{$PetscInt}()
	maxDegree_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMFieldGetDegree, $petsc_library),
               PetscErrorCode,
               (DMField, CIS, Ptr{$PetscInt}, Ptr{$PetscInt}),
               field, cellIS, minDegree_, maxDegree_,
              )

	minDegree = minDegree_[]
	maxDegree = maxDegree_[]

	return minDegree,maxDegree
end 

"""
	nc::PetscInt = DMFieldGetNumComponents(petsclib::PetscLibType, field::DMField) 
Returns the number of components in the field

Not Collective

Input Parameter:
- `field` - The `DMField` object

Output Parameter:
- `nc` - The number of field components

Level: intermediate

See also: `DMField`, `DMFieldEvaluate()`

# External Links
$(_doc_external("DM/DMFieldGetNumComponents"))
"""
function DMFieldGetNumComponents(petsclib::PetscLibType, field::DMField)
    error("DMFieldGetNumComponents: no generated method for these argument types")
end

@for_petsc function DMFieldGetNumComponents(petsclib::$UnionPetscLib, field::DMField )
	nc_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMFieldGetNumComponents, $petsc_library),
               PetscErrorCode,
               (DMField, Ptr{$PetscInt}),
               field, nc_,
              )

	nc = nc_[]

	return nc
end 

"""
	type::DMFieldType = DMFieldGetType(petsclib::PetscLibType, field::DMField) 
Gets the `DMFieldType` name (as a string) from the `DMField`.

Not Collective

Input Parameter:
- `field` - The `DMField` context

Output Parameter:
- `type` - The `DMFieldType` name

Level: advanced

See also: `DMField`, `DMFieldSetType()`, `DMFieldType`, `PetscObjectTypeCompare()`, `PetscObjectTypeCompareAny()`

# External Links
$(_doc_external("DM/DMFieldGetType"))
"""
function DMFieldGetType(petsclib::PetscLibType, field::DMField)
    error("DMFieldGetType: no generated method for these argument types")
end

@for_petsc function DMFieldGetType(petsclib::$UnionPetscLib, field::DMField )
	type_ = Ref{DMFieldType}()

    @chk ccall(
               (:DMFieldGetType, $petsc_library),
               PetscErrorCode,
               (DMField, Ptr{DMFieldType}),
               field, type_,
              )

	type = type_[] == C_NULL ? "" : unsafe_string(type_[])

	return type
end 

"""
	DMFieldInitializePackage(petsclib::PetscLibType) 
Initialize `DMField` package

Logically Collective

Level: developer

See also: `DMFieldFinalizePackage()`

# External Links
$(_doc_external("DM/DMFieldInitializePackage"))
"""
function DMFieldInitializePackage(petsclib::PetscLibType)
    error("DMFieldInitializePackage: no generated method for these argument types")
end

@for_petsc function DMFieldInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:DMFieldInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	DMFieldRegister(petsclib::PetscLibType, sname::String, fnc::external) 
Adds an implementation of the `DMField` object.

Not collective, No Fortran Support

Input Parameters:
- `sname`    - name of a new user-defined implementation
- `function` - routine to create method context

See also: `DMField`, `DMFieldRegisterAll()`, `DMFieldRegisterDestroy()`

# External Links
$(_doc_external("DM/DMFieldRegister"))
"""
function DMFieldRegister(petsclib::PetscLibType, sname::String, fnc::external)
    error("DMFieldRegister: no generated method for these argument types")
end

@for_petsc function DMFieldRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:DMFieldRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	DMFieldSetType(petsclib::PetscLibType, field::DMField, type::DMFieldType) 
set the `DMField` implementation

Collective

Input Parameters:
- `field` - the `DMField` context
- `type`  - a known method, see `DMFieldType`

Level: advanced

See also: `DMField`, `DMFieldGetType()`, `DMFieldType`

# External Links
$(_doc_external("DM/DMFieldSetType"))
"""
function DMFieldSetType(petsclib::PetscLibType, field::DMField, type::DMFieldType)
    error("DMFieldSetType: no generated method for these argument types")
end

@for_petsc function DMFieldSetType(petsclib::$UnionPetscLib, field::DMField, type::DMFieldType )

    @chk ccall(
               (:DMFieldSetType, $petsc_library),
               PetscErrorCode,
               (DMField, DMFieldType),
               field, type,
              )


	return nothing
end 

"""
	B::Ptr{Cvoid},D::Ptr{Cvoid},H::Ptr{Cvoid} = DMFieldShellEvaluateFEDefault(petsclib::PetscLibType, field::DMField, pointIS::AbstractIS, quad::PetscQuadrature, type::PetscDataType) 
Default finite-element evaluation for a `DMFIELDSHELL` that maps the quadrature points to real space using the coordinate `DMField` and then calls `DMFieldEvaluate()`.

Not Collective

Input Parameters:
- `field`   - the `DMField` of type `DMFIELDSHELL`
- `pointIS` - the `IS` of mesh points at which to evaluate
- `quad`    - the reference-element quadrature
- `type`    - `PETSC_SCALAR` or `PETSC_REAL`

Output Parameters:
- `B` - values at quadrature points, or `NULL`
- `D` - derivatives at quadrature points, or `NULL`
- `H` - Hessians at quadrature points, or `NULL`

Level: developer

See also: `DMField`, `DMFIELDSHELL`, `DMFieldShellSetEvaluateFE()`, `DMFieldShellEvaluateFVDefault()`, `DMFieldEvaluate()`

# External Links
$(_doc_external("DM/DMFieldShellEvaluateFEDefault"))
"""
function DMFieldShellEvaluateFEDefault(petsclib::PetscLibType, field::DMField, pointIS::AbstractIS, quad::PetscQuadrature, type::PetscDataType)
    error("DMFieldShellEvaluateFEDefault: no generated method for these argument types")
end

@for_petsc function DMFieldShellEvaluateFEDefault(petsclib::$UnionPetscLib, field::DMField, pointIS::AbstractIS, quad::PetscQuadrature, type::PetscDataType )
	B_ = Ref{Ptr{Cvoid}}()
	D_ = Ref{Ptr{Cvoid}}()
	H_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:DMFieldShellEvaluateFEDefault, $petsc_library),
               PetscErrorCode,
               (DMField, CIS, PetscQuadrature, PetscDataType, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
               field, pointIS, quad, type, B_, D_, H_,
              )

	B = B_[]
	D = D_[]
	H = H_[]

	return B,D,H
end 

"""
	B::Ptr{Cvoid},D::Ptr{Cvoid},H::Ptr{Cvoid} = DMFieldShellEvaluateFVDefault(petsclib::PetscLibType, field::DMField, pointIS::AbstractIS, type::PetscDataType) 
Default finite-volume evaluation for a `DMFIELDSHELL` that samples at cell centroids using the coordinate `DMField`'s default quadrature and calls `DMFieldEvaluate()`.

Not Collective

Input Parameters:
- `field`   - the `DMField` of type `DMFIELDSHELL`
- `pointIS` - the `IS` of mesh cells at which to evaluate
- `type`    - `PETSC_SCALAR` or `PETSC_REAL`

Output Parameters:
- `B` - cell-averaged values, or `NULL`
- `D` - cell-averaged derivatives, or `NULL`
- `H` - cell-averaged Hessians, or `NULL`

Level: developer

See also: `DMField`, `DMFIELDSHELL`, `DMFieldShellSetEvaluateFV()`, `DMFieldShellEvaluateFEDefault()`, `DMFieldEvaluate()`

# External Links
$(_doc_external("DM/DMFieldShellEvaluateFVDefault"))
"""
function DMFieldShellEvaluateFVDefault(petsclib::PetscLibType, field::DMField, pointIS::AbstractIS, type::PetscDataType)
    error("DMFieldShellEvaluateFVDefault: no generated method for these argument types")
end

@for_petsc function DMFieldShellEvaluateFVDefault(petsclib::$UnionPetscLib, field::DMField, pointIS::AbstractIS, type::PetscDataType )
	B_ = Ref{Ptr{Cvoid}}()
	D_ = Ref{Ptr{Cvoid}}()
	H_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:DMFieldShellEvaluateFVDefault, $petsc_library),
               PetscErrorCode,
               (DMField, CIS, PetscDataType, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
               field, pointIS, type, B_, D_, H_,
              )

	B = B_[]
	D = D_[]
	H = H_[]

	return B,D,H
end 

"""
	ctx::Ptr{Cvoid} = DMFieldShellGetContext(petsclib::PetscLibType, field::DMField) 
Retrieve the user-supplied context associated with a `DMFIELDSHELL`.

Not Collective

Input Parameter:
- `field` - the `DMField` of type `DMFIELDSHELL`

Output Parameter:
- `ctx` - the context pointer that was passed to `DMFieldCreateShell()`

Level: intermediate

See also: `DMField`, `DMFIELDSHELL`, `DMFieldCreateShell()`

# External Links
$(_doc_external("DM/DMFieldShellGetContext"))
"""
function DMFieldShellGetContext(petsclib::PetscLibType, field::DMField)
    error("DMFieldShellGetContext: no generated method for these argument types")
end

@for_petsc function DMFieldShellGetContext(petsclib::$UnionPetscLib, field::DMField )
	ctx_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:DMFieldShellGetContext, $petsc_library),
               PetscErrorCode,
               (DMField, Ptr{Cvoid}),
               field, ctx_,
              )

	ctx = ctx_[]

	return ctx
end 

"""
	DMFieldShellSetCreateDefaultQuadrature(petsclib::PetscLibType, field::DMField, create::external) 
Register the routine that supplies a default `PetscQuadrature` sufficient to integrate a `DMFIELDSHELL` exactly over a set of mesh points.

Logically Collective

Input Parameters:
- `field`  - the `DMField` of type `DMFIELDSHELL`
- `create` - callback that returns a newly created `PetscQuadrature` for the given point `IS`

Calling sequence of `create`:
- `f`    - the `DMField` of type `DMFIELDSHELL`
- `is`   - the `IS` of mesh points over which the field will be integrated
- `quad` - the newly created `PetscQuadrature`

Level: intermediate

See also: `DMField`, `DMFIELDSHELL`, `DMFieldCreateShell()`, `DMFieldCreateDefaultQuadrature()`

# External Links
$(_doc_external("DM/DMFieldShellSetCreateDefaultQuadrature"))
"""
function DMFieldShellSetCreateDefaultQuadrature(petsclib::PetscLibType, field::DMField, create::external)
    error("DMFieldShellSetCreateDefaultQuadrature: no generated method for these argument types")
end

@for_petsc function DMFieldShellSetCreateDefaultQuadrature(petsclib::$UnionPetscLib, field::DMField, create::external )

    @chk ccall(
               (:DMFieldShellSetCreateDefaultQuadrature, $petsc_library),
               PetscErrorCode,
               (DMField, external),
               field, create,
              )


	return nothing
end 

"""
	DMFieldShellSetDestroy(petsclib::PetscLibType, field::DMField, destroy::external) 
Register a destroy callback that will be invoked when a `DMFIELDSHELL` is destroyed.

Logically Collective

Input Parameters:
- `field`   - the `DMField` of type `DMFIELDSHELL`
- `destroy` - the destroy routine, called before the shell's own data is freed

Calling sequence of `destroy`:
- `field` - the `DMField` of type `DMFIELDSHELL` being destroyed

Level: intermediate

See also: `DMField`, `DMFIELDSHELL`, `DMFieldCreateShell()`, `DMFieldDestroy()`

# External Links
$(_doc_external("DM/DMFieldShellSetDestroy"))
"""
function DMFieldShellSetDestroy(petsclib::PetscLibType, field::DMField, destroy::external)
    error("DMFieldShellSetDestroy: no generated method for these argument types")
end

@for_petsc function DMFieldShellSetDestroy(petsclib::$UnionPetscLib, field::DMField, destroy::external )

    @chk ccall(
               (:DMFieldShellSetDestroy, $petsc_library),
               PetscErrorCode,
               (DMField, external),
               field, destroy,
              )


	return nothing
end 

"""
	DMFieldShellSetEvaluate(petsclib::PetscLibType, field::DMField, evaluate::external) 
Register the routine that evaluates a `DMFIELDSHELL` at an arbitrary set of real-space points supplied as a `Vec` of coordinates.

Logically Collective

Input Parameters:
- `field`    - the `DMField` of type `DMFIELDSHELL`
- `evaluate` - the evaluation callback

Calling sequence of `evaluate`:
- `field` - the `DMField` of type `DMFIELDSHELL`
- `u`     - the points at which to evaluate the field, as a `Vec` of coordinates of size d x n
- `dtype` - `PETSC_SCALAR` or `PETSC_REAL`
- `B`     - array of field values at each point, or `NULL`
- `D`     - array of field spatial derivatives at each point, or `NULL`
- `H`     - array of field spatial Hessians at each point, or `NULL`

Level: intermediate

See also: `DMField`, `DMFIELDSHELL`, `DMFieldCreateShell()`, `DMFieldEvaluate()`, `DMFieldShellSetEvaluateFE()`, `DMFieldShellSetEvaluateFV()`

# External Links
$(_doc_external("DM/DMFieldShellSetEvaluate"))
"""
function DMFieldShellSetEvaluate(petsclib::PetscLibType, field::DMField, evaluate::external)
    error("DMFieldShellSetEvaluate: no generated method for these argument types")
end

@for_petsc function DMFieldShellSetEvaluate(petsclib::$UnionPetscLib, field::DMField, evaluate::external )

    @chk ccall(
               (:DMFieldShellSetEvaluate, $petsc_library),
               PetscErrorCode,
               (DMField, external),
               field, evaluate,
              )


	return nothing
end 

"""
	DMFieldShellSetEvaluateFE(petsclib::PetscLibType, field::DMField, evaluateFE::external) 
Register the routine that evaluates a `DMFIELDSHELL` at finite-element quadrature points over a set of mesh points.

Logically Collective

Input Parameters:
- `field`      - the `DMField` of type `DMFIELDSHELL`
- `evaluateFE` - the FE evaluation callback

Calling sequence of `evaluateFE`:
- `field` - the `DMField` of type `DMFIELDSHELL`
- `is`    - the `IS` of mesh cells on which to evaluate the field
- `quad`  - the reference-cell `PetscQuadrature` supplying the evaluation points
- `dtype` - `PETSC_SCALAR` or `PETSC_REAL`
- `B`     - array of field values at each quadrature point, or `NULL`
- `D`     - array of field reference derivatives at each quadrature point, or `NULL`
- `H`     - array of field reference Hessians at each quadrature point, or `NULL`

Level: intermediate

See also: `DMField`, `DMFIELDSHELL`, `DMFieldCreateShell()`, `DMFieldEvaluateFE()`, `DMFieldShellEvaluateFEDefault()`, `DMFieldShellSetEvaluateFV()`

# External Links
$(_doc_external("DM/DMFieldShellSetEvaluateFE"))
"""
function DMFieldShellSetEvaluateFE(petsclib::PetscLibType, field::DMField, evaluateFE::external)
    error("DMFieldShellSetEvaluateFE: no generated method for these argument types")
end

@for_petsc function DMFieldShellSetEvaluateFE(petsclib::$UnionPetscLib, field::DMField, evaluateFE::external )

    @chk ccall(
               (:DMFieldShellSetEvaluateFE, $petsc_library),
               PetscErrorCode,
               (DMField, external),
               field, evaluateFE,
              )


	return nothing
end 

"""
	DMFieldShellSetEvaluateFV(petsclib::PetscLibType, field::DMField, evaluateFV::external) 
Register the routine that evaluates a `DMFIELDSHELL` as cell averages over a set of mesh cells.

Logically Collective

Input Parameters:
- `field`      - the `DMField` of type `DMFIELDSHELL`
- `evaluateFV` - the FV evaluation callback

Calling sequence of `evaluateFV`:
- `field` - the `DMField` of type `DMFIELDSHELL`
- `is`    - the `IS` of mesh cells on which to evaluate the field
- `dtype` - `PETSC_SCALAR` or `PETSC_REAL`
- `B`     - array of cell-averaged field values, or `NULL`
- `D`     - array of cell-averaged field derivatives, or `NULL`
- `H`     - array of cell-averaged field Hessians, or `NULL`

Level: intermediate

See also: `DMField`, `DMFIELDSHELL`, `DMFieldCreateShell()`, `DMFieldEvaluateFV()`, `DMFieldShellEvaluateFVDefault()`, `DMFieldShellSetEvaluateFE()`

# External Links
$(_doc_external("DM/DMFieldShellSetEvaluateFV"))
"""
function DMFieldShellSetEvaluateFV(petsclib::PetscLibType, field::DMField, evaluateFV::external)
    error("DMFieldShellSetEvaluateFV: no generated method for these argument types")
end

@for_petsc function DMFieldShellSetEvaluateFV(petsclib::$UnionPetscLib, field::DMField, evaluateFV::external )

    @chk ccall(
               (:DMFieldShellSetEvaluateFV, $petsc_library),
               PetscErrorCode,
               (DMField, external),
               field, evaluateFV,
              )


	return nothing
end 

"""
	DMFieldShellSetGetDegree(petsclib::PetscLibType, field::DMField, getDegree::external) 
Register the routine that reports the polynomial degree bounds of a `DMFIELDSHELL` over a set of mesh points.

Logically Collective

Input Parameters:
- `field`     - the `DMField` of type `DMFIELDSHELL`
- `getDegree` - callback that returns the minimum and maximum polynomial degrees of the field over the given point `IS`

Calling sequence of `getDegree`:
- `field`     - the `DMField` of type `DMFIELDSHELL`
- `is`        - the `IS` of mesh points over which the degree bounds are requested
- `minDegree` - the degree of the largest polynomial space contained in the field on each element
- `maxDegree` - the largest degree of the smallest polynomial space containing the field on any element

Level: intermediate

See also: `DMField`, `DMFIELDSHELL`, `DMFieldCreateShell()`, `DMFieldGetDegree()`

# External Links
$(_doc_external("DM/DMFieldShellSetGetDegree"))
"""
function DMFieldShellSetGetDegree(petsclib::PetscLibType, field::DMField, getDegree::external)
    error("DMFieldShellSetGetDegree: no generated method for these argument types")
end

@for_petsc function DMFieldShellSetGetDegree(petsclib::$UnionPetscLib, field::DMField, getDegree::external )

    @chk ccall(
               (:DMFieldShellSetGetDegree, $petsc_library),
               PetscErrorCode,
               (DMField, external),
               field, getDegree,
              )


	return nothing
end 

"""
	DMFieldView(petsclib::PetscLibType, field::DMField, viewer::PetscViewer) 
view a `DMField`

Collective

Input Parameters:
- `field`  - `DMField`
- `viewer` - viewer to display field, for example `PETSC_VIEWER_STDOUT_WORLD`

Level: advanced

See also: `DMField`, `DMFieldCreate()`

# External Links
$(_doc_external("DM/DMFieldView"))
"""
function DMFieldView(petsclib::PetscLibType, field::DMField, viewer::PetscViewer)
    error("DMFieldView: no generated method for these argument types")
end

@for_petsc function DMFieldView(petsclib::$UnionPetscLib, field::DMField, viewer::PetscViewer )

    @chk ccall(
               (:DMFieldView, $petsc_library),
               PetscErrorCode,
               (DMField, PetscViewer),
               field, viewer,
              )


	return nothing
end 

"""
	DMLabelAddStrata(petsclib::PetscLibType, label::DMLabel, numStrata::PetscInt, stratumValues::Vector{PetscInt}) 
Adds new stratum values in a `DMLabel`

Not Collective

Input Parameters:
- `label`         - The `DMLabel`
- `numStrata`     - The number of stratum values
- `stratumValues` - The stratum values

Level: beginner

See also: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelDestroy()`

# External Links
$(_doc_external("DMLabel/DMLabelAddStrata"))
"""
function DMLabelAddStrata(petsclib::PetscLibType, label::DMLabel, numStrata::Integer, stratumValues::AbstractVector{<:Number})
    error("DMLabelAddStrata: no generated method for these argument types")
end

@for_petsc function DMLabelAddStrata(petsclib::$UnionPetscLib, label::DMLabel, numStrata::$PetscInt, stratumValues::Vector{$PetscInt} )

    @chk ccall(
               (:DMLabelAddStrata, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, Ptr{$PetscInt}),
               label, numStrata, stratumValues,
              )


	return nothing
end 

"""
	DMLabelAddStrataIS(petsclib::PetscLibType, label::DMLabel, valueIS::AbstractIS) 
Adds new stratum values in a `DMLabel`

Not Collective

Input Parameters:
- `label`   - The `DMLabel`
- `valueIS` - Index set with stratum values

Level: beginner

See also: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelDestroy()`

# External Links
$(_doc_external("DMLabel/DMLabelAddStrataIS"))
"""
function DMLabelAddStrataIS(petsclib::PetscLibType, label::DMLabel, valueIS::AbstractIS)
    error("DMLabelAddStrataIS: no generated method for these argument types")
end

@for_petsc function DMLabelAddStrataIS(petsclib::$UnionPetscLib, label::DMLabel, valueIS::AbstractIS )

    @chk ccall(
               (:DMLabelAddStrataIS, $petsc_library),
               PetscErrorCode,
               (DMLabel, CIS),
               label, valueIS,
              )


	return nothing
end 

"""
	DMLabelAddStratum(petsclib::PetscLibType, label::DMLabel, value::PetscInt) 
Adds a new stratum value in a `DMLabel`

Input Parameters:
- `label` - The `DMLabel`
- `value` - The stratum value

Level: beginner

See also: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelDestroy()`

# External Links
$(_doc_external("DMLabel/DMLabelAddStratum"))
"""
function DMLabelAddStratum(petsclib::PetscLibType, label::DMLabel, value::Integer)
    error("DMLabelAddStratum: no generated method for these argument types")
end

@for_petsc function DMLabelAddStratum(petsclib::$UnionPetscLib, label::DMLabel, value::$PetscInt )

    @chk ccall(
               (:DMLabelAddStratum, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt),
               label, value,
              )


	return nothing
end 

"""
	DMLabelClearStratum(petsclib::PetscLibType, label::DMLabel, value::PetscInt) 
Remove a stratum

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `value` - the stratum value

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("DMLabel/DMLabelClearStratum"))
"""
function DMLabelClearStratum(petsclib::PetscLibType, label::DMLabel, value::Integer)
    error("DMLabelClearStratum: no generated method for these argument types")
end

@for_petsc function DMLabelClearStratum(petsclib::$UnionPetscLib, label::DMLabel, value::$PetscInt )

    @chk ccall(
               (:DMLabelClearStratum, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt),
               label, value,
              )


	return nothing
end 

"""
	DMLabelClearValue(petsclib::PetscLibType, label::DMLabel, point::PetscInt, value::PetscInt) 
Clear the value a label assigns to a point

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `point` - the point
- `value` - The point value

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`

# External Links
$(_doc_external("DMLabel/DMLabelClearValue"))
"""
function DMLabelClearValue(petsclib::PetscLibType, label::DMLabel, point::Integer, value::Integer)
    error("DMLabelClearValue: no generated method for these argument types")
end

@for_petsc function DMLabelClearValue(petsclib::$UnionPetscLib, label::DMLabel, point::$PetscInt, value::$PetscInt )

    @chk ccall(
               (:DMLabelClearValue, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, $PetscInt),
               label, point, value,
              )


	return nothing
end 

"""
	equal::PetscBool,message::String = DMLabelCompare(petsclib::PetscLibType, comm::MPI_Comm, l0::DMLabel, l1::DMLabel) 
Compare two `DMLabel` objects

Collective; No Fortran Support

Input Parameters:
- `comm` - Comm over which to compare labels
- `l0`   - First `DMLabel`
- `l1`   - Second `DMLabel`

Output Parameters:
- `equal`   - (Optional) Flag whether the two labels are equal
- `message` - (Optional) Message describing the difference

Level: intermediate

See also: `DMLabel`, `DM`, `DMCompareLabels()`, `DMLabelGetNumValues()`, `DMLabelGetDefaultValue()`, `DMLabelGetNonEmptyStratumValuesIS()`, `DMLabelGetStratumIS()`

# External Links
$(_doc_external("DMLabel/DMLabelCompare"))
"""
function DMLabelCompare(petsclib::PetscLibType, comm::MPI_Comm, l0::DMLabel, l1::DMLabel)
    error("DMLabelCompare: no generated method for these argument types")
end

@for_petsc function DMLabelCompare(petsclib::$UnionPetscLib, comm::MPI_Comm, l0::DMLabel, l1::DMLabel )
	equal_ = Ref{PetscBool}()
	message_ = Ref{Ptr{Cchar}}()

    @chk ccall(
               (:DMLabelCompare, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, DMLabel, DMLabel, Ptr{PetscBool}, Ptr{Ptr{Cchar}}),
               comm, l0, l1, equal_, message_,
              )

	equal = equal_[]
	message = unsafe_string(message_[])

	return equal,message
end 

"""
	DMLabelComputeIndex(petsclib::PetscLibType, label::DMLabel) 
Create an index structure for membership determination, automatically determining the bounds

Not Collective

Input Parameter:
- `label` - The `DMLabel`

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelHasPoint()`, `DMLabelCreateIndex()`, `DMLabelDestroyIndex()`, `DMLabelGetValue()`, `DMLabelSetValue()`

# External Links
$(_doc_external("DMLabel/DMLabelComputeIndex"))
"""
function DMLabelComputeIndex(petsclib::PetscLibType, label::DMLabel)
    error("DMLabelComputeIndex: no generated method for these argument types")
end

@for_petsc function DMLabelComputeIndex(petsclib::$UnionPetscLib, label::DMLabel )

    @chk ccall(
               (:DMLabelComputeIndex, $petsc_library),
               PetscErrorCode,
               (DMLabel,),
               label,
              )


	return nothing
end 

"""
	section::PetscSection,is::IS = DMLabelConvertToSection(petsclib::PetscLibType, label::DMLabel) 
Make a `PetscSection`/`IS` pair that encodes the label

Not Collective

Input Parameter:
- `label` - the `DMLabel`

Output Parameters:
- `section` - the section giving offsets for each stratum
- `is`      - An `IS` containing all the label points

Level: developer

See also: `DMLabel`, `DM`, `DMLabelDistribute()`

# External Links
$(_doc_external("DMLabel/DMLabelConvertToSection"))
"""
function DMLabelConvertToSection(petsclib::PetscLibType, label::DMLabel)
    error("DMLabelConvertToSection: no generated method for these argument types")
end

@for_petsc function DMLabelConvertToSection(petsclib::$UnionPetscLib, label::DMLabel )
	section_ = Ref{PetscSection}()
	is_ = Ref{CIS}()

    @chk ccall(
               (:DMLabelConvertToSection, $petsc_library),
               PetscErrorCode,
               (DMLabel, Ptr{PetscSection}, Ptr{CIS}),
               label, section_, is_,
              )

	section = section_[]
	is = IS(is_[], petsclib)

	return section,is
end 

"""
	label::DMLabel = DMLabelCreate(petsclib::PetscLibType, comm::MPI_Comm, name::String) 
Create a `DMLabel` object, which is a multimap

Collective

Input Parameters:
- `comm` - The communicator, usually `PETSC_COMM_SELF`
- `name` - The label name

Output Parameter:
- `label` - The `DMLabel`

Level: beginner

See also: `DMLabel`, `DM`, `DMLabelDestroy()`

# External Links
$(_doc_external("DMLabel/DMLabelCreate"))
"""
function DMLabelCreate(petsclib::PetscLibType, comm::MPI_Comm, name::String)
    error("DMLabelCreate: no generated method for these argument types")
end

@for_petsc function DMLabelCreate(petsclib::$UnionPetscLib, comm::MPI_Comm, name::String )
	label_ = Ref{DMLabel}()

    @chk ccall(
               (:DMLabelCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{DMLabel}),
               comm, name, label_,
              )

	label = label_[]

	return label
end 

"""
	DMLabelCreateIndex(petsclib::PetscLibType, label::DMLabel, pStart::PetscInt, pEnd::PetscInt) 
Create an index structure for membership determination

Not Collective

Input Parameters:
- `label`  - The `DMLabel`
- `pStart` - The smallest point
- `pEnd`   - The largest point + 1

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelHasPoint()`, `DMLabelComputeIndex()`, `DMLabelDestroyIndex()`, `DMLabelGetValue()`, `DMLabelSetValue()`

# External Links
$(_doc_external("DMLabel/DMLabelCreateIndex"))
"""
function DMLabelCreateIndex(petsclib::PetscLibType, label::DMLabel, pStart::Integer, pEnd::Integer)
    error("DMLabelCreateIndex: no generated method for these argument types")
end

@for_petsc function DMLabelCreateIndex(petsclib::$UnionPetscLib, label::DMLabel, pStart::$PetscInt, pEnd::$PetscInt )

    @chk ccall(
               (:DMLabelCreateIndex, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, $PetscInt),
               label, pStart, pEnd,
              )


	return nothing
end 

"""
	DMLabelDestroy(petsclib::PetscLibType, label::Union{DMLabel, Ref{DMLabel}}) 
Destroys a `DMLabel`

Collective

Input Parameter:
- `label` - The `DMLabel`

Level: beginner

See also: `DMLabel`, `DM`, `DMLabelReset()`, `DMLabelCreate()`

# External Links
$(_doc_external("DMLabel/DMLabelDestroy"))
"""
function DMLabelDestroy(petsclib::PetscLibType, label::Union{DMLabel, Ref{DMLabel}})
    error("DMLabelDestroy: no generated method for these argument types")
end

@for_petsc function DMLabelDestroy(petsclib::$UnionPetscLib, label::Union{DMLabel, Ref{DMLabel}} )
	label_ = label isa Base.RefValue ? label : Ref{DMLabel}(label)

    @chk ccall(
               (:DMLabelDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{DMLabel},),
               label_,
              )


	return nothing
end 

"""
	DMLabelDestroyIndex(petsclib::PetscLibType, label::DMLabel) 
Destroy the index structure

Not Collective

Input Parameter:
- `label` - the `DMLabel`

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelHasPoint()`, `DMLabelCreateIndex()`, `DMLabelGetValue()`, `DMLabelSetValue()`

# External Links
$(_doc_external("DMLabel/DMLabelDestroyIndex"))
"""
function DMLabelDestroyIndex(petsclib::PetscLibType, label::DMLabel)
    error("DMLabelDestroyIndex: no generated method for these argument types")
end

@for_petsc function DMLabelDestroyIndex(petsclib::$UnionPetscLib, label::DMLabel )

    @chk ccall(
               (:DMLabelDestroyIndex, $petsc_library),
               PetscErrorCode,
               (DMLabel,),
               label,
              )


	return nothing
end 

"""
	labelNew::DMLabel = DMLabelDistribute(petsclib::PetscLibType, label::DMLabel, sf::PetscSF) 
Create a new label pushed forward over the `PetscSF`

Collective

Input Parameters:
- `label` - the `DMLabel`
- `sf`    - the map from old to new distribution

Output Parameter:
- `labelNew` - the new redistributed label

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("DMLabel/DMLabelDistribute"))
"""
function DMLabelDistribute(petsclib::PetscLibType, label::DMLabel, sf::PetscSF)
    error("DMLabelDistribute: no generated method for these argument types")
end

@for_petsc function DMLabelDistribute(petsclib::$UnionPetscLib, label::DMLabel, sf::PetscSF )
	labelNew_ = Ref{DMLabel}()

    @chk ccall(
               (:DMLabelDistribute, $petsc_library),
               PetscErrorCode,
               (DMLabel, PetscSF, Ptr{DMLabel}),
               label, sf, labelNew_,
              )

	labelNew = labelNew_[]

	return labelNew
end 

"""
	labelnew::DMLabel = DMLabelDuplicate(petsclib::PetscLibType, label::DMLabel) 
Duplicates a `DMLabel`

Collective

Input Parameter:
- `label` - The `DMLabel`

Output Parameter:
- `labelnew` - new label

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelDestroy()`

# External Links
$(_doc_external("DMLabel/DMLabelDuplicate"))
"""
function DMLabelDuplicate(petsclib::PetscLibType, label::DMLabel)
    error("DMLabelDuplicate: no generated method for these argument types")
end

@for_petsc function DMLabelDuplicate(petsclib::$UnionPetscLib, label::DMLabel )
	labelnew_ = Ref{DMLabel}()

    @chk ccall(
               (:DMLabelDuplicate, $petsc_library),
               PetscErrorCode,
               (DMLabel, Ptr{DMLabel}),
               label, labelnew_,
              )

	labelnew = labelnew_[]

	return labelnew
end 

"""
	olabel::DMLabel = DMLabelEphemeralGetLabel(petsclib::PetscLibType, label::DMLabel) 
Get the base label for this ephemeral label

Not Collective

Input Parameter:
- `label` - the `DMLabel`

Output Parameter:
- `olabel` - the base label for this ephemeral label

Level: intermediate

See also: `DMLabelEphemeralSetLabel()`, `DMLabelEphemeralGetTransform()`, `DMLabelSetType()`

# External Links
$(_doc_external("DMLabel/DMLabelEphemeralGetLabel"))
"""
function DMLabelEphemeralGetLabel(petsclib::PetscLibType, label::DMLabel)
    error("DMLabelEphemeralGetLabel: no generated method for these argument types")
end

@for_petsc function DMLabelEphemeralGetLabel(petsclib::$UnionPetscLib, label::DMLabel )
	olabel_ = Ref{DMLabel}()

    @chk ccall(
               (:DMLabelEphemeralGetLabel, $petsc_library),
               PetscErrorCode,
               (DMLabel, Ptr{DMLabel}),
               label, olabel_,
              )

	olabel = olabel_[]

	return olabel
end 

"""
	tr::DMPlexTransform = DMLabelEphemeralGetTransform(petsclib::PetscLibType, label::DMLabel) 
Get the transform for this ephemeral label

Not Collective

Input Parameter:
- `label` - the `DMLabel`

Output Parameter:
- `tr` - the transform for this ephemeral label

Level: intermediate

See also: `DMLabelEphemeralSetTransform()`, `DMLabelEphemeralGetLabel()`, `DMLabelSetType()`

# External Links
$(_doc_external("DMLabel/DMLabelEphemeralGetTransform"))
"""
function DMLabelEphemeralGetTransform(petsclib::PetscLibType, label::DMLabel)
    error("DMLabelEphemeralGetTransform: no generated method for these argument types")
end

@for_petsc function DMLabelEphemeralGetTransform(petsclib::$UnionPetscLib, label::DMLabel )
	tr_ = Ref{DMPlexTransform}()

    @chk ccall(
               (:DMLabelEphemeralGetTransform, $petsc_library),
               PetscErrorCode,
               (DMLabel, Ptr{DMPlexTransform}),
               label, tr_,
              )

	tr = tr_[]

	return tr
end 

"""
	DMLabelEphemeralSetLabel(petsclib::PetscLibType, label::DMLabel, olabel::DMLabel) 
Set the base label for this ephemeral label

Not Collective

Input Parameters:
- `label`  - the `DMLabel`
- `olabel` - the base label for this ephemeral label

Level: intermediate

See also: `DMLabelEphemeralGetLabel()`, `DMLabelEphemeralSetTransform()`, `DMLabelSetType()`

# External Links
$(_doc_external("DMLabel/DMLabelEphemeralSetLabel"))
"""
function DMLabelEphemeralSetLabel(petsclib::PetscLibType, label::DMLabel, olabel::DMLabel)
    error("DMLabelEphemeralSetLabel: no generated method for these argument types")
end

@for_petsc function DMLabelEphemeralSetLabel(petsclib::$UnionPetscLib, label::DMLabel, olabel::DMLabel )

    @chk ccall(
               (:DMLabelEphemeralSetLabel, $petsc_library),
               PetscErrorCode,
               (DMLabel, DMLabel),
               label, olabel,
              )


	return nothing
end 

"""
	DMLabelEphemeralSetTransform(petsclib::PetscLibType, label::DMLabel, tr::DMPlexTransform) 
Set the transform for this ephemeral label

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `tr`    - the transform for this ephemeral label

Level: intermediate

See also: `DMLabelEphemeralGetTransform()`, `DMLabelEphemeralSetLabel()`, `DMLabelSetType()`

# External Links
$(_doc_external("DMLabel/DMLabelEphemeralSetTransform"))
"""
function DMLabelEphemeralSetTransform(petsclib::PetscLibType, label::DMLabel, tr::DMPlexTransform)
    error("DMLabelEphemeralSetTransform: no generated method for these argument types")
end

@for_petsc function DMLabelEphemeralSetTransform(petsclib::$UnionPetscLib, label::DMLabel, tr::DMPlexTransform )

    @chk ccall(
               (:DMLabelEphemeralSetTransform, $petsc_library),
               PetscErrorCode,
               (DMLabel, DMPlexTransform),
               label, tr,
              )


	return nothing
end 

"""
	DMLabelFilter(petsclib::PetscLibType, label::DMLabel, start::PetscInt, end_::PetscInt) 
Remove all points outside of [`start`, `end`)

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `start` - the first point kept
- `end`   - one more than the last point kept

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("DMLabel/DMLabelFilter"))
"""
function DMLabelFilter(petsclib::PetscLibType, label::DMLabel, start::Integer, end_::Integer)
    error("DMLabelFilter: no generated method for these argument types")
end

@for_petsc function DMLabelFilter(petsclib::$UnionPetscLib, label::DMLabel, start::$PetscInt, end_::$PetscInt )

    @chk ccall(
               (:DMLabelFilter, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, $PetscInt),
               label, start, end_,
              )


	return nothing
end 

"""
	labelNew::DMLabel = DMLabelGather(petsclib::PetscLibType, label::DMLabel, sf::PetscSF) 
Gather all label values from leafs into roots

Collective

Input Parameters:
- `label` - the `DMLabel`
- `sf`    - the `PetscSF` communication map

Output Parameter:
- `labelNew` - the new `DMLabel` with localised leaf values

Level: developer

See also: `DMLabel`, `DM`, `DMLabelDistribute()`

# External Links
$(_doc_external("DMLabel/DMLabelGather"))
"""
function DMLabelGather(petsclib::PetscLibType, label::DMLabel, sf::PetscSF)
    error("DMLabelGather: no generated method for these argument types")
end

@for_petsc function DMLabelGather(petsclib::$UnionPetscLib, label::DMLabel, sf::PetscSF )
	labelNew_ = Ref{DMLabel}()

    @chk ccall(
               (:DMLabelGather, $petsc_library),
               PetscErrorCode,
               (DMLabel, PetscSF, Ptr{DMLabel}),
               label, sf, labelNew_,
              )

	labelNew = labelNew_[]

	return labelNew
end 

"""
	pStart::PetscInt,pEnd::PetscInt = DMLabelGetBounds(petsclib::PetscLibType, label::DMLabel) 
Return the smallest and largest point in the label

Not Collective

Input Parameter:
- `label` - the `DMLabel`

Output Parameters:
- `pStart` - The smallest point
- `pEnd`   - The largest point + 1

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelHasPoint()`, `DMLabelCreateIndex()`, `DMLabelGetValue()`, `DMLabelSetValue()`

# External Links
$(_doc_external("DMLabel/DMLabelGetBounds"))
"""
function DMLabelGetBounds(petsclib::PetscLibType, label::DMLabel)
    error("DMLabelGetBounds: no generated method for these argument types")
end

@for_petsc function DMLabelGetBounds(petsclib::$UnionPetscLib, label::DMLabel )
	pStart_ = Ref{$PetscInt}()
	pEnd_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMLabelGetBounds, $petsc_library),
               PetscErrorCode,
               (DMLabel, Ptr{$PetscInt}, Ptr{$PetscInt}),
               label, pStart_, pEnd_,
              )

	pStart = pStart_[]
	pEnd = pEnd_[]

	return pStart,pEnd
end 

"""
	defaultValue::PetscInt = DMLabelGetDefaultValue(petsclib::PetscLibType, label::DMLabel) 
Get the default value returned by `DMLabelGetValue()` if a point has not been explicitly given a value.
When a label is created, it is initialized to -1.

Not Collective

Input Parameter:
- `label` - a `DMLabel` object

Output Parameter:
- `defaultValue` - the default value

Level: beginner

See also: `DMLabel`, `DM`, `DMLabelSetDefaultValue()`, `DMLabelGetValue()`, `DMLabelSetValue()`

# External Links
$(_doc_external("DMLabel/DMLabelGetDefaultValue"))
"""
function DMLabelGetDefaultValue(petsclib::PetscLibType, label::DMLabel)
    error("DMLabelGetDefaultValue: no generated method for these argument types")
end

@for_petsc function DMLabelGetDefaultValue(petsclib::$UnionPetscLib, label::DMLabel )
	defaultValue_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMLabelGetDefaultValue, $petsc_library),
               PetscErrorCode,
               (DMLabel, Ptr{$PetscInt}),
               label, defaultValue_,
              )

	defaultValue = defaultValue_[]

	return defaultValue
end 

"""
	values::IS = DMLabelGetNonEmptyStratumValuesIS(petsclib::PetscLibType, label::DMLabel) 
Get an `IS` of all values that the `DMlabel` takes

Not Collective

Input Parameter:
- `label` - the `DMLabel`

Output Parameter:
- `values` - the value `IS`

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelGetValueIS()`, `DMLabelGetValueISGlobal()`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("DMLabel/DMLabelGetNonEmptyStratumValuesIS"))
"""
function DMLabelGetNonEmptyStratumValuesIS(petsclib::PetscLibType, label::DMLabel)
    error("DMLabelGetNonEmptyStratumValuesIS: no generated method for these argument types")
end

@for_petsc function DMLabelGetNonEmptyStratumValuesIS(petsclib::$UnionPetscLib, label::DMLabel )
	values_ = Ref{CIS}()

    @chk ccall(
               (:DMLabelGetNonEmptyStratumValuesIS, $petsc_library),
               PetscErrorCode,
               (DMLabel, Ptr{CIS}),
               label, values_,
              )

	values = IS(values_[], petsclib)

	return values
end 

"""
	numValues::PetscInt = DMLabelGetNumValues(petsclib::PetscLibType, label::DMLabel) 
Get the number of values that the `DMLabel` takes

Not Collective

Input Parameter:
- `label` - the `DMLabel`

Output Parameter:
- `numValues` - the number of values

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("DMLabel/DMLabelGetNumValues"))
"""
function DMLabelGetNumValues(petsclib::PetscLibType, label::DMLabel)
    error("DMLabelGetNumValues: no generated method for these argument types")
end

@for_petsc function DMLabelGetNumValues(petsclib::$UnionPetscLib, label::DMLabel )
	numValues_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMLabelGetNumValues, $petsc_library),
               PetscErrorCode,
               (DMLabel, Ptr{$PetscInt}),
               label, numValues_,
              )

	numValues = numValues_[]

	return numValues
end 

"""
	start::PetscInt,end_::PetscInt = DMLabelGetStratumBounds(petsclib::PetscLibType, label::DMLabel, value::PetscInt) 
Get the largest and smallest point of a stratum

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `value` - the stratum value

Output Parameters:
- `start` - the smallest point in the stratum
- `end`   - the largest point in the stratum

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("DMLabel/DMLabelGetStratumBounds"))
"""
function DMLabelGetStratumBounds(petsclib::PetscLibType, label::DMLabel, value::Integer)
    error("DMLabelGetStratumBounds: no generated method for these argument types")
end

@for_petsc function DMLabelGetStratumBounds(petsclib::$UnionPetscLib, label::DMLabel, value::$PetscInt )
	start_ = Ref{$PetscInt}()
	end__ = Ref{$PetscInt}()

    @chk ccall(
               (:DMLabelGetStratumBounds, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
               label, value, start_, end__,
              )

	start = start_[]
	end_ = end__[]

	return start,end_
end 

"""
	points::IS = DMLabelGetStratumIS(petsclib::PetscLibType, label::DMLabel, value::PetscInt) 
Get an `IS` with the stratum points

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `value` - the stratum value

Output Parameter:
- `points` - The stratum points

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("DMLabel/DMLabelGetStratumIS"))
"""
function DMLabelGetStratumIS(petsclib::PetscLibType, label::DMLabel, value::Integer)
    error("DMLabelGetStratumIS: no generated method for these argument types")
end

@for_petsc function DMLabelGetStratumIS(petsclib::$UnionPetscLib, label::DMLabel, value::$PetscInt )
	points_ = Ref{CIS}()

    @chk ccall(
               (:DMLabelGetStratumIS, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, Ptr{CIS}),
               label, value, points_,
              )

	points = IS(points_[], petsclib)

	return points
end 

"""
	index::PetscInt = DMLabelGetStratumPointIndex(petsclib::PetscLibType, label::DMLabel, value::PetscInt, p::PetscInt) 
Get the index of a point in a given stratum

Not Collective

Input Parameters:
- `label` - The `DMLabel`
- `value` - The label value
- `p`     - A point with this value

Output Parameter:
- `index` - The index of this point in the stratum, or -1 if the point is not in the stratum or the stratum does not exist

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelGetValueIndex()`, `DMLabelGetStratumIS()`, `DMLabelCreate()`

# External Links
$(_doc_external("DMLabel/DMLabelGetStratumPointIndex"))
"""
function DMLabelGetStratumPointIndex(petsclib::PetscLibType, label::DMLabel, value::Integer, p::Integer)
    error("DMLabelGetStratumPointIndex: no generated method for these argument types")
end

@for_petsc function DMLabelGetStratumPointIndex(petsclib::$UnionPetscLib, label::DMLabel, value::$PetscInt, p::$PetscInt )
	index_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMLabelGetStratumPointIndex, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, $PetscInt, Ptr{$PetscInt}),
               label, value, p, index_,
              )

	index = index_[]

	return index
end 

"""
	size::PetscInt = DMLabelGetStratumSize(petsclib::PetscLibType, label::DMLabel, value::PetscInt) 
Get the size of a stratum

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `value` - the stratum value

Output Parameter:
- `size` - The number of points in the stratum

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("DMLabel/DMLabelGetStratumSize"))
"""
function DMLabelGetStratumSize(petsclib::PetscLibType, label::DMLabel, value::Integer)
    error("DMLabelGetStratumSize: no generated method for these argument types")
end

@for_petsc function DMLabelGetStratumSize(petsclib::$UnionPetscLib, label::DMLabel, value::$PetscInt )
	size_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMLabelGetStratumSize, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, Ptr{$PetscInt}),
               label, value, size_,
              )

	size = size_[]

	return size
end 

"""
	type::DMLabelType = DMLabelGetType(petsclib::PetscLibType, label::DMLabel) 
Gets the type name (as a string) from the label.

Not Collective

Input Parameter:
- `label` - The `DMLabel`

Output Parameter:
- `type` - The `DMLabel` type name

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelSetType()`, `DMLabelCreate()`

# External Links
$(_doc_external("DMLabel/DMLabelGetType"))
"""
function DMLabelGetType(petsclib::PetscLibType, label::DMLabel)
    error("DMLabelGetType: no generated method for these argument types")
end

@for_petsc function DMLabelGetType(petsclib::$UnionPetscLib, label::DMLabel )
	type_ = Ref{DMLabelType}()

    @chk ccall(
               (:DMLabelGetType, $petsc_library),
               PetscErrorCode,
               (DMLabel, Ptr{DMLabelType}),
               label, type_,
              )

	type = type_[] == C_NULL ? "" : unsafe_string(type_[])

	return type
end 

"""
	value::PetscInt = DMLabelGetValue(petsclib::PetscLibType, label::DMLabel, point::PetscInt) 
Return the value a label assigns to a point, or the label's default value (which is initially -1, and can be changed with
`DMLabelSetDefaultValue()`)

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `point` - the point

Output Parameter:
- `value` - The point value, or the default value (-1 by default)

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelSetValue()`, `DMLabelClearValue()`, `DMLabelGetDefaultValue()`, `DMLabelSetDefaultValue()`

# External Links
$(_doc_external("DMLabel/DMLabelGetValue"))
"""
function DMLabelGetValue(petsclib::PetscLibType, label::DMLabel, point::Integer)
    error("DMLabelGetValue: no generated method for these argument types")
end

@for_petsc function DMLabelGetValue(petsclib::$UnionPetscLib, label::DMLabel, point::$PetscInt )
	value_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMLabelGetValue, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, Ptr{$PetscInt}),
               label, point, value_,
              )

	value = value_[]

	return value
end 

"""
	minValue::PetscInt,maxValue::PetscInt = DMLabelGetValueBounds(petsclib::PetscLibType, label::DMLabel) 
Return the smallest and largest value in the label

Not Collective

Input Parameter:
- `label` - the `DMLabel`

Output Parameters:
- `minValue` - The smallest value
- `maxValue` - The largest value

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelGetBounds()`, `DMLabelGetValue()`, `DMLabelSetValue()`

# External Links
$(_doc_external("DMLabel/DMLabelGetValueBounds"))
"""
function DMLabelGetValueBounds(petsclib::PetscLibType, label::DMLabel)
    error("DMLabelGetValueBounds: no generated method for these argument types")
end

@for_petsc function DMLabelGetValueBounds(petsclib::$UnionPetscLib, label::DMLabel )
	minValue_ = Ref{$PetscInt}()
	maxValue_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMLabelGetValueBounds, $petsc_library),
               PetscErrorCode,
               (DMLabel, Ptr{$PetscInt}, Ptr{$PetscInt}),
               label, minValue_, maxValue_,
              )

	minValue = minValue_[]
	maxValue = maxValue_[]

	return minValue,maxValue
end 

"""
	values::IS = DMLabelGetValueIS(petsclib::PetscLibType, label::DMLabel) 
Get an `IS` of all values that the `DMlabel` takes

Not Collective

Input Parameter:
- `label` - the `DMLabel`

Output Parameter:
- `values` - the value `IS`

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelGetNonEmptyStratumValuesIS()`, `DMLabelGetValueISGlobal()`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("DMLabel/DMLabelGetValueIS"))
"""
function DMLabelGetValueIS(petsclib::PetscLibType, label::DMLabel)
    error("DMLabelGetValueIS: no generated method for these argument types")
end

@for_petsc function DMLabelGetValueIS(petsclib::$UnionPetscLib, label::DMLabel )
	values_ = Ref{CIS}()

    @chk ccall(
               (:DMLabelGetValueIS, $petsc_library),
               PetscErrorCode,
               (DMLabel, Ptr{CIS}),
               label, values_,
              )

	values = IS(values_[], petsclib)

	return values
end 

"""
	values::IS = DMLabelGetValueISGlobal(petsclib::PetscLibType, comm::MPI_Comm, label::DMLabel, get_nonempty::PetscBool) 
Get an `IS` of all values that the `DMlabel` takes across all ranks

Collective

Input Parameter:
- `comm`         - MPI communicator to collect values
- `label`        - the `DMLabel`, may be `NULL` for ranks in `comm` which do not have the corresponding `DMLabel`
- `get_nonempty` - whether to get nonempty stratum values (akin to `DMLabelGetNonEmptyStratumValuesIS()`)

Output Parameter:
- `values` - the value `IS`

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelGetValueIS()`, `DMLabelGetNonEmptyStratumValuesIS()`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("DMLabel/DMLabelGetValueISGlobal"))
"""
function DMLabelGetValueISGlobal(petsclib::PetscLibType, comm::MPI_Comm, label::DMLabel, get_nonempty::PetscBool)
    error("DMLabelGetValueISGlobal: no generated method for these argument types")
end

@for_petsc function DMLabelGetValueISGlobal(petsclib::$UnionPetscLib, comm::MPI_Comm, label::DMLabel, get_nonempty::PetscBool )
	values_ = Ref{CIS}()

    @chk ccall(
               (:DMLabelGetValueISGlobal, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, DMLabel, PetscBool, Ptr{CIS}),
               comm, label, get_nonempty, values_,
              )

	values = IS(values_[], petsclib)

	return values
end 

"""
	index::PetscInt = DMLabelGetValueIndex(petsclib::PetscLibType, label::DMLabel, value::PetscInt) 
Get the index of a given value in the list of values for the `DMlabel`, or -1 if it is not present

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `value` - the value

Output Parameter:
- `index` - the index of value in the list of values

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelGetValueIS()`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("DMLabel/DMLabelGetValueIndex"))
"""
function DMLabelGetValueIndex(petsclib::PetscLibType, label::DMLabel, value::Integer)
    error("DMLabelGetValueIndex: no generated method for these argument types")
end

@for_petsc function DMLabelGetValueIndex(petsclib::$UnionPetscLib, label::DMLabel, value::$PetscInt )
	index_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMLabelGetValueIndex, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, Ptr{$PetscInt}),
               label, value, index_,
              )

	index = index_[]

	return index
end 

"""
	contains::PetscBool = DMLabelHasPoint(petsclib::PetscLibType, label::DMLabel, point::PetscInt) 
Determine whether a label assigns a value to a point

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `point` - the point

Output Parameter:
- `contains` - Flag indicating whether the label maps this point to a value

Level: developer

See also: `DMLabel`, `DM`, `DMLabelCreateIndex()`, `DMLabelGetValue()`, `DMLabelSetValue()`

# External Links
$(_doc_external("DMLabel/DMLabelHasPoint"))
"""
function DMLabelHasPoint(petsclib::PetscLibType, label::DMLabel, point::Integer)
    error("DMLabelHasPoint: no generated method for these argument types")
end

@for_petsc function DMLabelHasPoint(petsclib::$UnionPetscLib, label::DMLabel, point::$PetscInt )
	contains_ = Ref{PetscBool}()

    @chk ccall(
               (:DMLabelHasPoint, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, Ptr{PetscBool}),
               label, point, contains_,
              )

	contains = contains_[]

	return contains
end 

"""
	exists::PetscBool = DMLabelHasStratum(petsclib::PetscLibType, label::DMLabel, value::PetscInt) 
Determine whether points exist with the given value

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `value` - the stratum value

Output Parameter:
- `exists` - Flag saying whether points exist

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("DMLabel/DMLabelHasStratum"))
"""
function DMLabelHasStratum(petsclib::PetscLibType, label::DMLabel, value::Integer)
    error("DMLabelHasStratum: no generated method for these argument types")
end

@for_petsc function DMLabelHasStratum(petsclib::$UnionPetscLib, label::DMLabel, value::$PetscInt )
	exists_ = Ref{PetscBool}()

    @chk ccall(
               (:DMLabelHasStratum, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, Ptr{PetscBool}),
               label, value, exists_,
              )

	exists = exists_[]

	return exists
end 

"""
	contains::PetscBool = DMLabelHasValue(petsclib::PetscLibType, label::DMLabel, value::PetscInt) 
Determine whether a label assigns the value to any point

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `value` - the value

Output Parameter:
- `contains` - Flag indicating whether the label maps this value to any point

Level: developer

See also: `DMLabel`, `DM`, `DMLabelHasPoint()`, `DMLabelGetValue()`, `DMLabelSetValue()`

# External Links
$(_doc_external("DMLabel/DMLabelHasValue"))
"""
function DMLabelHasValue(petsclib::PetscLibType, label::DMLabel, value::Integer)
    error("DMLabelHasValue: no generated method for these argument types")
end

@for_petsc function DMLabelHasValue(petsclib::$UnionPetscLib, label::DMLabel, value::$PetscInt )
	contains_ = Ref{PetscBool}()

    @chk ccall(
               (:DMLabelHasValue, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, Ptr{PetscBool}),
               label, value, contains_,
              )

	contains = contains_[]

	return contains
end 

"""
	DMLabelInsertIS(petsclib::PetscLibType, label::DMLabel, is::AbstractIS, value::PetscInt) 
Set all points in the `IS` to a value

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `is`    - the point `IS`
- `value` - The point value

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("DMLabel/DMLabelInsertIS"))
"""
function DMLabelInsertIS(petsclib::PetscLibType, label::DMLabel, is::AbstractIS, value::Integer)
    error("DMLabelInsertIS: no generated method for these argument types")
end

@for_petsc function DMLabelInsertIS(petsclib::$UnionPetscLib, label::DMLabel, is::AbstractIS, value::$PetscInt )

    @chk ccall(
               (:DMLabelInsertIS, $petsc_library),
               PetscErrorCode,
               (DMLabel, CIS, $PetscInt),
               label, is, value,
              )


	return nothing
end 

"""
	labelNew::DMLabel = DMLabelPermute(petsclib::PetscLibType, label::DMLabel, permutation::AbstractIS) 
Create a new label with permuted points

Not Collective

Input Parameters:
- `label`       - the `DMLabel`
- `permutation` - the point permutation

Output Parameter:
- `labelNew` - the new label containing the permuted points

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("DMLabel/DMLabelPermute"))
"""
function DMLabelPermute(petsclib::PetscLibType, label::DMLabel, permutation::AbstractIS)
    error("DMLabelPermute: no generated method for these argument types")
end

@for_petsc function DMLabelPermute(petsclib::$UnionPetscLib, label::DMLabel, permutation::AbstractIS )
	labelNew_ = Ref{DMLabel}()

    @chk ccall(
               (:DMLabelPermute, $petsc_library),
               PetscErrorCode,
               (DMLabel, CIS, Ptr{DMLabel}),
               label, permutation, labelNew_,
              )

	labelNew = labelNew_[]

	return labelNew
end 

"""
	DMLabelPermuteValues(petsclib::PetscLibType, label::DMLabel, permutation::AbstractIS) 
Permute the values in a label

Not collective

Input Parameters:
- `label`       - the `DMLabel`
- `permutation` - the value permutation, permutation[old value] = new value

Output Parameter:
- `label` - the `DMLabel` now with permuted values

See also: `DMLabelRewriteValues()`, `DMLabel`, `DM`, `DMLabelPermute()`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("DMLabel/DMLabelPermuteValues"))
"""
function DMLabelPermuteValues(petsclib::PetscLibType, label::DMLabel, permutation::AbstractIS)
    error("DMLabelPermuteValues: no generated method for these argument types")
end

@for_petsc function DMLabelPermuteValues(petsclib::$UnionPetscLib, label::DMLabel, permutation::AbstractIS )

    @chk ccall(
               (:DMLabelPermuteValues, $petsc_library),
               PetscErrorCode,
               (DMLabel, CIS),
               label, permutation,
              )


	return nothing
end 

"""
	DMLabelPropagateBegin(petsclib::PetscLibType, label::DMLabel, sf::PetscSF) 
Setup a cycle of label propagation

Collective

Input Parameters:
- `label` - The `DMLabel` to propagate across processes
- `sf`    - The `PetscSF` describing parallel layout of the label points

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelPropagateEnd()`, `DMLabelPropagatePush()`

# External Links
$(_doc_external("DMLabel/DMLabelPropagateBegin"))
"""
function DMLabelPropagateBegin(petsclib::PetscLibType, label::DMLabel, sf::PetscSF)
    error("DMLabelPropagateBegin: no generated method for these argument types")
end

@for_petsc function DMLabelPropagateBegin(petsclib::$UnionPetscLib, label::DMLabel, sf::PetscSF )

    @chk ccall(
               (:DMLabelPropagateBegin, $petsc_library),
               PetscErrorCode,
               (DMLabel, PetscSF),
               label, sf,
              )


	return nothing
end 

"""
	DMLabelPropagateEnd(petsclib::PetscLibType, label::DMLabel, pointSF::PetscSF) 
Tear down a cycle of label propagation

Collective

Input Parameters:
- `label`   - The `DMLabel` to propagate across processes
- `pointSF` - The `PetscSF` describing parallel layout of the label points

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelPropagateBegin()`, `DMLabelPropagatePush()`

# External Links
$(_doc_external("DMLabel/DMLabelPropagateEnd"))
"""
function DMLabelPropagateEnd(petsclib::PetscLibType, label::DMLabel, pointSF::PetscSF)
    error("DMLabelPropagateEnd: no generated method for these argument types")
end

@for_petsc function DMLabelPropagateEnd(petsclib::$UnionPetscLib, label::DMLabel, pointSF::PetscSF )

    @chk ccall(
               (:DMLabelPropagateEnd, $petsc_library),
               PetscErrorCode,
               (DMLabel, PetscSF),
               label, pointSF,
              )


	return nothing
end 

"""
	DMLabelPropagatePush(petsclib::PetscLibType, label::DMLabel, pointSF::PetscSF, markPoint::external, ctx::Ptr{Cvoid}) 
Tear down a cycle of label propagation

Collective

Input Parameters:
- `label`     - The `DMLabel` to propagate across processes
- `pointSF`   - The `PetscSF` describing parallel layout of the label points
- `markPoint` - An optional callback that is called when a point is marked, or `NULL`
- `ctx`       - An optional application context for the callback, or `NULL`

Calling sequence of `markPoint`:
- `label` - The `DMLabel`
- `p`     - The point being marked
- `val`   - The label value for `p`
- `ctx`   - An optional application context

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelPropagateBegin()`, `DMLabelPropagateEnd()`

# External Links
$(_doc_external("DMLabel/DMLabelPropagatePush"))
"""
function DMLabelPropagatePush(petsclib::PetscLibType, label::DMLabel, pointSF::PetscSF, markPoint::external, ctx::Ptr{Cvoid})
    error("DMLabelPropagatePush: no generated method for these argument types")
end

@for_petsc function DMLabelPropagatePush(petsclib::$UnionPetscLib, label::DMLabel, pointSF::PetscSF, markPoint::external, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:DMLabelPropagatePush, $petsc_library),
               PetscErrorCode,
               (DMLabel, PetscSF, external, Ptr{Cvoid}),
               label, pointSF, markPoint, ctx,
              )


	return nothing
end 

"""
	DMLabelRegister(petsclib::PetscLibType, name::String, noname::Ptr{Cvoid}) 
Adds a new label component implementation

Not Collective

Input Parameters:
- `name`        - The name of a new user-defined creation routine
- `create_func` - The creation routine itself

See also: `DMLabel`, `DM`, `DMLabelType`, `DMLabelRegisterAll()`, `DMLabelRegisterDestroy()`

# External Links
$(_doc_external("DMLabel/DMLabelRegister"))
"""
function DMLabelRegister(petsclib::PetscLibType, name::String, noname::Ptr{Cvoid})
    error("DMLabelRegister: no generated method for these argument types")
end

@for_petsc function DMLabelRegister(petsclib::$UnionPetscLib, name::String, noname::Ptr{Cvoid} )

    @chk ccall(
               (:DMLabelRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cvoid}),
               name, noname,
              )


	return nothing
end 

"""
	DMLabelRegisterAll(petsclib::PetscLibType) 
Registers all of the `DMLabel` implementations in the `DM` package.

Not Collective

Level: advanced

See also: `DMLabel`, `DM`, `DMRegisterAll()`, `DMLabelRegisterDestroy()`

# External Links
$(_doc_external("DMLabel/DMLabelRegisterAll"))
"""
function DMLabelRegisterAll(petsclib::PetscLibType)
    error("DMLabelRegisterAll: no generated method for these argument types")
end

@for_petsc function DMLabelRegisterAll(petsclib::$UnionPetscLib)

    @chk ccall(
               (:DMLabelRegisterAll, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	DMLabelRegisterDestroy(petsclib::PetscLibType) 
This function destroys the `DMLabel` registry. It is called from `PetscFinalize()`.

Level: developer

See also: `DMLabel`, `DM`, `PetscInitialize()`

# External Links
$(_doc_external("DMLabel/DMLabelRegisterDestroy"))
"""
function DMLabelRegisterDestroy(petsclib::PetscLibType)
    error("DMLabelRegisterDestroy: no generated method for these argument types")
end

@for_petsc function DMLabelRegisterDestroy(petsclib::$UnionPetscLib)

    @chk ccall(
               (:DMLabelRegisterDestroy, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	DMLabelReset(petsclib::PetscLibType, label::DMLabel) 
Destroys internal data structures in a `DMLabel`

Not Collective

Input Parameter:
- `label` - The `DMLabel`

Level: beginner

See also: `DMLabel`, `DM`, `DMLabelDestroy()`, `DMLabelCreate()`

# External Links
$(_doc_external("DMLabel/DMLabelReset"))
"""
function DMLabelReset(petsclib::PetscLibType, label::DMLabel)
    error("DMLabelReset: no generated method for these argument types")
end

@for_petsc function DMLabelReset(petsclib::$UnionPetscLib, label::DMLabel )

    @chk ccall(
               (:DMLabelReset, $petsc_library),
               PetscErrorCode,
               (DMLabel,),
               label,
              )


	return nothing
end 

"""
	DMLabelRewriteValues(petsclib::PetscLibType, label::DMLabel, permutation::AbstractIS) 
Permute the values in a label, but some may be omitted

Not collective

Input Parameters:
- `label`       - the `DMLabel`
- `permutation` - the value permutation, permutation[old value] = new value, but some maybe omitted

Output Parameter:
- `label` - the `DMLabel` now with permuted values

See also: `DMLabelPermuteValues()`, `DMLabel`, `DM`, `DMLabelPermute()`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("DMLabel/DMLabelRewriteValues"))
"""
function DMLabelRewriteValues(petsclib::PetscLibType, label::DMLabel, permutation::AbstractIS)
    error("DMLabelRewriteValues: no generated method for these argument types")
end

@for_petsc function DMLabelRewriteValues(petsclib::$UnionPetscLib, label::DMLabel, permutation::AbstractIS )

    @chk ccall(
               (:DMLabelRewriteValues, $petsc_library),
               PetscErrorCode,
               (DMLabel, CIS),
               label, permutation,
              )


	return nothing
end 

"""
	DMLabelSetDefaultValue(petsclib::PetscLibType, label::DMLabel, defaultValue::PetscInt) 
Set the default value returned by `DMLabelGetValue()` if a point has not been explicitly given a value.
When a label is created, it is initialized to -1.

Not Collective

Input Parameter:
- `label` - a `DMLabel` object

Output Parameter:
- `defaultValue` - the default value

Level: beginner

See also: `DMLabel`, `DM`, `DMLabelGetDefaultValue()`, `DMLabelGetValue()`, `DMLabelSetValue()`

# External Links
$(_doc_external("DMLabel/DMLabelSetDefaultValue"))
"""
function DMLabelSetDefaultValue(petsclib::PetscLibType, label::DMLabel, defaultValue::Integer)
    error("DMLabelSetDefaultValue: no generated method for these argument types")
end

@for_petsc function DMLabelSetDefaultValue(petsclib::$UnionPetscLib, label::DMLabel, defaultValue::$PetscInt )

    @chk ccall(
               (:DMLabelSetDefaultValue, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt),
               label, defaultValue,
              )


	return nothing
end 

"""
	DMLabelSetStratumBounds(petsclib::PetscLibType, label::DMLabel, value::PetscInt, pStart::PetscInt, pEnd::PetscInt) 
Efficiently give a contiguous set of points a given label value

Not Collective

Input Parameters:
- `label`  - The `DMLabel`
- `value`  - The label value for all points
- `pStart` - The first point
- `pEnd`   - A point beyond all marked points

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelSetStratumIS()`, `DMLabelGetStratumIS()`

# External Links
$(_doc_external("DMLabel/DMLabelSetStratumBounds"))
"""
function DMLabelSetStratumBounds(petsclib::PetscLibType, label::DMLabel, value::Integer, pStart::Integer, pEnd::Integer)
    error("DMLabelSetStratumBounds: no generated method for these argument types")
end

@for_petsc function DMLabelSetStratumBounds(petsclib::$UnionPetscLib, label::DMLabel, value::$PetscInt, pStart::$PetscInt, pEnd::$PetscInt )

    @chk ccall(
               (:DMLabelSetStratumBounds, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, $PetscInt, $PetscInt),
               label, value, pStart, pEnd,
              )


	return nothing
end 

"""
	DMLabelSetStratumIS(petsclib::PetscLibType, label::DMLabel, value::PetscInt, is::AbstractIS) 
Set the stratum points using an `IS`

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `value` - the stratum value
- `is`    - The stratum points

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("DMLabel/DMLabelSetStratumIS"))
"""
function DMLabelSetStratumIS(petsclib::PetscLibType, label::DMLabel, value::Integer, is::AbstractIS)
    error("DMLabelSetStratumIS: no generated method for these argument types")
end

@for_petsc function DMLabelSetStratumIS(petsclib::$UnionPetscLib, label::DMLabel, value::$PetscInt, is::AbstractIS )

    @chk ccall(
               (:DMLabelSetStratumIS, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, CIS),
               label, value, is,
              )


	return nothing
end 

"""
	DMLabelSetType(petsclib::PetscLibType, label::DMLabel, method::DMLabelType) 
Sets the particular implementation for a label.

Collective

Input Parameters:
- `label`  - The label
- `method` - The name of the label type

Options Database Key:
- `-dm_label_type type` - Sets the label type; see `DMLabelType`

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelGetType()`, `DMLabelCreate()`, `DMLabelType`

# External Links
$(_doc_external("DMLabel/DMLabelSetType"))
"""
function DMLabelSetType(petsclib::PetscLibType, label::DMLabel, method::DMLabelType)
    error("DMLabelSetType: no generated method for these argument types")
end

@for_petsc function DMLabelSetType(petsclib::$UnionPetscLib, label::DMLabel, method::DMLabelType )

    @chk ccall(
               (:DMLabelSetType, $petsc_library),
               PetscErrorCode,
               (DMLabel, DMLabelType),
               label, method,
              )


	return nothing
end 

"""
	DMLabelSetUp(petsclib::PetscLibType, label::DMLabel) 
SetUp a `DMLabel` object

Collective

Input Parameters:
- `label` - The `DMLabel`

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelDestroy()`

# External Links
$(_doc_external("DMLabel/DMLabelSetUp"))
"""
function DMLabelSetUp(petsclib::PetscLibType, label::DMLabel)
    error("DMLabelSetUp: no generated method for these argument types")
end

@for_petsc function DMLabelSetUp(petsclib::$UnionPetscLib, label::DMLabel )

    @chk ccall(
               (:DMLabelSetUp, $petsc_library),
               PetscErrorCode,
               (DMLabel,),
               label,
              )


	return nothing
end 

"""
	DMLabelSetValue(petsclib::PetscLibType, label::DMLabel, point::PetscInt, value::PetscInt) 
Set the value a label assigns to a point.  If the value is the same as the label's default value (which is initially -1, and can
be changed with `DMLabelSetDefaultValue()` to something different), then this function will do nothing.

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `point` - the point
- `value` - The point value

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelGetValue()`, `DMLabelClearValue()`, `DMLabelGetDefaultValue()`, `DMLabelSetDefaultValue()`

# External Links
$(_doc_external("DMLabel/DMLabelSetValue"))
"""
function DMLabelSetValue(petsclib::PetscLibType, label::DMLabel, point::Integer, value::Integer)
    error("DMLabelSetValue: no generated method for these argument types")
end

@for_petsc function DMLabelSetValue(petsclib::$UnionPetscLib, label::DMLabel, point::$PetscInt, value::$PetscInt )

    @chk ccall(
               (:DMLabelSetValue, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, $PetscInt),
               label, point, value,
              )


	return nothing
end 

"""
	contains::PetscBool = DMLabelStratumHasPoint(petsclib::PetscLibType, label::DMLabel, value::PetscInt, point::PetscInt) 
Return true if the stratum contains a point

Not Collective

Input Parameters:
- `label` - the `DMLabel`
- `value` - the stratum value
- `point` - the point

Output Parameter:
- `contains` - true if the stratum contains the point

Level: intermediate

See also: `DMLabel`, `DM`, `DMLabelCreate()`, `DMLabelSetValue()`, `DMLabelClearValue()`

# External Links
$(_doc_external("DMLabel/DMLabelStratumHasPoint"))
"""
function DMLabelStratumHasPoint(petsclib::PetscLibType, label::DMLabel, value::Integer, point::Integer)
    error("DMLabelStratumHasPoint: no generated method for these argument types")
end

@for_petsc function DMLabelStratumHasPoint(petsclib::$UnionPetscLib, label::DMLabel, value::$PetscInt, point::$PetscInt )
	contains_ = Ref{PetscBool}()

    @chk ccall(
               (:DMLabelStratumHasPoint, $petsc_library),
               PetscErrorCode,
               (DMLabel, $PetscInt, $PetscInt, Ptr{PetscBool}),
               label, value, point, contains_,
              )

	contains = contains_[]

	return contains
end 

"""
	DMLabelView(petsclib::PetscLibType, label::DMLabel, viewer::PetscViewer) 
View the label

Collective

Input Parameters:
- `label`  - The `DMLabel`
- `viewer` - The `PetscViewer`

Level: intermediate

See also: `DMLabel`, `PetscViewer`, `DM`, `DMLabelCreate()`, `DMLabelDestroy()`

# External Links
$(_doc_external("DMLabel/DMLabelView"))
"""
function DMLabelView(petsclib::PetscLibType, label::DMLabel, viewer::PetscViewer)
    error("DMLabelView: no generated method for these argument types")
end

@for_petsc function DMLabelView(petsclib::$UnionPetscLib, label::DMLabel, viewer::PetscViewer )

    @chk ccall(
               (:DMLabelView, $petsc_library),
               PetscErrorCode,
               (DMLabel, PetscViewer),
               label, viewer,
              )


	return nothing
end 

"""
	DMLabelViewFromOptions(petsclib::PetscLibType, label::DMLabel, obj, name::String) 
View a `DMLabel` in a particular way based on a request in the options database

Collective

Input Parameters:
- `label` - the `DMLabel` object
- `obj`   - optional object that provides the prefix for the options database (if `NULL` then the prefix in `obj` is used)
- `name`  - option string that is used to activate viewing

Options Database Key:
- `-name [viewertype][:...]` - option name and values. See `PetscObjectViewFromOptions()` for the possible arguments

Level: intermediate

See also: `DMLabel`, `DMLabelView()`, `PetscObjectViewFromOptions()`, `DMLabelCreate()`

# External Links
$(_doc_external("DMLabel/DMLabelViewFromOptions"))
"""
function DMLabelViewFromOptions(petsclib::PetscLibType, label::DMLabel, obj, name::String)
    error("DMLabelViewFromOptions: no generated method for these argument types")
end

@for_petsc function DMLabelViewFromOptions(petsclib::$UnionPetscLib, label::DMLabel, obj, name::String )

    @chk ccall(
               (:DMLabelViewFromOptions, $petsc_library),
               PetscErrorCode,
               (DMLabel, PetscObject, Ptr{Cchar}),
               label, obj, name,
              )


	return nothing
end 

"""
	DMNetworkMonitorAdd(petsclib::PetscLibType, monitor::DMNetworkMonitor, name::String, element::PetscInt, nodes::PetscInt, start::PetscInt, blocksize::PetscInt, xmin::PetscReal, xmax::PetscReal, ymin::PetscReal, ymax::PetscReal, hold::PetscBool) 
Adds a new viewer to a `DMNetworkMonitor`

Collective

Input Parameters:
- `monitor`   - the monitor
- `name`      - name of viewer
- `element`   - vertex / edge number
- `nodes`     - number of nodes
- `start`     - variable starting offset
- `blocksize` - variable blocksize
- `xmin`      - xmin (or `PETSC_DECIDE`) for viewer
- `xmax`      - xmax (or `PETSC_DECIDE`) for viewer
- `ymin`      - ymin for viewer
- `ymax`      - ymax for viewer
- `hold`      - determines if plot limits should be held

Level: intermediate

See also: `DM`, `DMNETWORK`, `DMNetworkMonitor`, `DMNetworkMonitorCreate()`, `DMNetworkMonitorDestroy()`

# External Links
$(_doc_external("DMNetwork/DMNetworkMonitorAdd"))
"""
function DMNetworkMonitorAdd(petsclib::PetscLibType, monitor::DMNetworkMonitor, name::String, element::Integer, nodes::Integer, start::Integer, blocksize::Integer, xmin::Real, xmax::Real, ymin::Real, ymax::Real, hold::PetscBool)
    error("DMNetworkMonitorAdd: no generated method for these argument types")
end

@for_petsc function DMNetworkMonitorAdd(petsclib::$UnionPetscLib, monitor::DMNetworkMonitor, name::String, element::$PetscInt, nodes::$PetscInt, start::$PetscInt, blocksize::$PetscInt, xmin::$PetscReal, xmax::$PetscReal, ymin::$PetscReal, ymax::$PetscReal, hold::PetscBool )

    @chk ccall(
               (:DMNetworkMonitorAdd, $petsc_library),
               PetscErrorCode,
               (DMNetworkMonitor, Ptr{Cchar}, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscReal, $PetscReal, $PetscReal, $PetscReal, PetscBool),
               monitor, name, element, nodes, start, blocksize, xmin, xmax, ymin, ymax, hold,
              )


	return nothing
end 

"""
	monitorptr::DMNetworkMonitor = DMNetworkMonitorCreate(petsclib::PetscLibType, network::AbstractPetscDM) 
Creates a network monitor context

Collective

Input Parameter:
- `network` - network to monitor

Output Parameter:
- `monitorptr` - the `DMNetworkMonitor` object

Level: intermediate

See also: `DM`, `DMNETWORK`, `DMNetworkMonitor`, `DMNetworkMonitorDestroy()`, `DMNetworkMonitorAdd()`

# External Links
$(_doc_external("DMNetwork/DMNetworkMonitorCreate"))
"""
function DMNetworkMonitorCreate(petsclib::PetscLibType, network::AbstractPetscDM)
    error("DMNetworkMonitorCreate: no generated method for these argument types")
end

@for_petsc function DMNetworkMonitorCreate(petsclib::$UnionPetscLib, network::AbstractPetscDM )
	monitorptr_ = Ref{DMNetworkMonitor}()

    @chk ccall(
               (:DMNetworkMonitorCreate, $petsc_library),
               PetscErrorCode,
               (CDM, Ptr{DMNetworkMonitor}),
               network, monitorptr_,
              )

	monitorptr = monitorptr_[]

	return monitorptr
end 

"""
	DMNetworkMonitorDestroy(petsclib::PetscLibType, monitor::Union{DMNetworkMonitor, Ref{DMNetworkMonitor}}) 
Destroys a network monitor and all associated viewers

Collective

Input Parameter:
- `monitor` - monitor to destroy

Level: intermediate

See also: `DM`, `DMNETWORK`, `DMNetworkMonitor`, `DMNetworkMonitorCreate()`, `DMNetworkMonitorAdd()`

# External Links
$(_doc_external("DMNetwork/DMNetworkMonitorDestroy"))
"""
function DMNetworkMonitorDestroy(petsclib::PetscLibType, monitor::Union{DMNetworkMonitor, Ref{DMNetworkMonitor}})
    error("DMNetworkMonitorDestroy: no generated method for these argument types")
end

@for_petsc function DMNetworkMonitorDestroy(petsclib::$UnionPetscLib, monitor::Union{DMNetworkMonitor, Ref{DMNetworkMonitor}} )
	monitor_ = monitor isa Base.RefValue ? monitor : Ref{DMNetworkMonitor}(monitor)

    @chk ccall(
               (:DMNetworkMonitorDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{DMNetworkMonitor},),
               monitor_,
              )


	return nothing
end 

"""
	DMNetworkMonitorPop(petsclib::PetscLibType, monitor::DMNetworkMonitor) 
Removes the most recently added viewer to a `DMNetworkMonitor`

Collective

Input Parameter:
- `monitor` - the monitor

Level: intermediate

See also: `DM`, `DMNETWORK`, `DMNetworkMonitor`, `DMNetworkMonitorCreate()`, `DMNetworkMonitorDestroy()`

# External Links
$(_doc_external("DMNetwork/DMNetworkMonitorPop"))
"""
function DMNetworkMonitorPop(petsclib::PetscLibType, monitor::DMNetworkMonitor)
    error("DMNetworkMonitorPop: no generated method for these argument types")
end

@for_petsc function DMNetworkMonitorPop(petsclib::$UnionPetscLib, monitor::DMNetworkMonitor )

    @chk ccall(
               (:DMNetworkMonitorPop, $petsc_library),
               PetscErrorCode,
               (DMNetworkMonitor,),
               monitor,
              )


	return nothing
end 

"""
	DMNetworkMonitorView(petsclib::PetscLibType, monitor::DMNetworkMonitor, x::AbstractPetscVec) 
A `DMNETWORK` specific monitor function for `TSMonitorSet()`

Collective, No Fortran support

Input Parameters:
- `monitor` - `DMNetworkMonitor` object
- `x`       - `TS` solution vector

Level: intermediate

See also: `DM`, `DMNETWORK`, `DMNetworkMonitor`, `DMNetworkMonitorCreate()`, `DMNetworkMonitorDestroy()`, `DMNetworkMonitorAdd()`

# External Links
$(_doc_external("DMNetwork/DMNetworkMonitorView"))
"""
function DMNetworkMonitorView(petsclib::PetscLibType, monitor::DMNetworkMonitor, x::AbstractPetscVec)
    error("DMNetworkMonitorView: no generated method for these argument types")
end

@for_petsc function DMNetworkMonitorView(petsclib::$UnionPetscLib, monitor::DMNetworkMonitor, x::AbstractPetscVec )

    @chk ccall(
               (:DMNetworkMonitorView, $petsc_library),
               PetscErrorCode,
               (DMNetworkMonitor, CVec),
               monitor, x,
              )


	return nothing
end 

"""
	p::PetscInt = DMPlexPointQueueBack(petsclib::PetscLibType, queue::DMPlexPoCintQueue) 
Return, without removing, the mesh point at the back of a `DMPlexPointQueue`.

Not Collective

Input Parameter:
- `queue` - the queue

Output Parameter:
- `p` - the mesh point at the back of the queue

Level: developer

See also: `DMPLEX`, `DMPlexPointQueue`, `DMPlexPointQueueFront()`, `DMPlexPointQueueEnqueue()`, `DMPlexPointQueueEmpty()`

# External Links
$(_doc_external("DMPlex/DMPlexPointQueueBack"))
"""
function DMPlexPointQueueBack(petsclib::PetscLibType, queue::DMPlexPoCintQueue)
    error("DMPlexPointQueueBack: no generated method for these argument types")
end

@for_petsc function DMPlexPointQueueBack(petsclib::$UnionPetscLib, queue::DMPlexPoCintQueue )
	p_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMPlexPointQueueBack, $petsc_library),
               PetscErrorCode,
               (DMPlexPoCintQueue, Ptr{$PetscInt}),
               queue, p_,
              )

	p = p_[]

	return p
end 

"""
	queue::DMPlexPoCintQueue = DMPlexPointQueueCreate(petsclib::PetscLibType, size::PetscInt) 
Create a `DMPlexPointQueue`, a simple FIFO queue of `PetscInt` mesh points used by `DMPLEX` traversal routines.

Not Collective

Input Parameter:
- `size` - the initial capacity of the queue

Output Parameter:
- `queue` - the newly created `DMPlexPointQueue`

Level: developer

See also: `DMPLEX`, `DMPlexPointQueue`, `DMPlexPointQueueDestroy()`, `DMPlexPointQueueEnqueue()`, `DMPlexPointQueueDequeue()`

# External Links
$(_doc_external("DMPlex/DMPlexPointQueueCreate"))
"""
function DMPlexPointQueueCreate(petsclib::PetscLibType, size::Integer)
    error("DMPlexPointQueueCreate: no generated method for these argument types")
end

@for_petsc function DMPlexPointQueueCreate(petsclib::$UnionPetscLib, size::$PetscInt )
	queue_ = Ref{DMPlexPoCintQueue}()

    @chk ccall(
               (:DMPlexPointQueueCreate, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{DMPlexPoCintQueue}),
               size, queue_,
              )

	queue = queue_[]

	return queue
end 

"""
	p::PetscInt = DMPlexPointQueueDequeue(petsclib::PetscLibType, queue::DMPlexPoCintQueue) 
Remove and return the mesh point at the front of a `DMPlexPointQueue`.

Not Collective

Input Parameter:
- `queue` - the queue

Output Parameter:
- `p` - the mesh point that was at the front of the queue

Level: developer

See also: `DMPLEX`, `DMPlexPointQueue`, `DMPlexPointQueueEnqueue()`, `DMPlexPointQueueFront()`, `DMPlexPointQueueEmpty()`

# External Links
$(_doc_external("DMPlex/DMPlexPointQueueDequeue"))
"""
function DMPlexPointQueueDequeue(petsclib::PetscLibType, queue::DMPlexPoCintQueue)
    error("DMPlexPointQueueDequeue: no generated method for these argument types")
end

@for_petsc function DMPlexPointQueueDequeue(petsclib::$UnionPetscLib, queue::DMPlexPoCintQueue )
	p_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMPlexPointQueueDequeue, $petsc_library),
               PetscErrorCode,
               (DMPlexPoCintQueue, Ptr{$PetscInt}),
               queue, p_,
              )

	p = p_[]

	return p
end 

"""
	DMPlexPointQueueDestroy(petsclib::PetscLibType, queue::Union{DMPlexPoCintQueue, Ref{DMPlexPoCintQueue}}) 
Destroy a `DMPlexPointQueue` previously created with `DMPlexPointQueueCreate()`.

Not Collective

Input Parameter:
- `queue` - the queue to destroy; set to `NULL` on return

Level: developer

See also: `DMPLEX`, `DMPlexPointQueue`, `DMPlexPointQueueCreate()`

# External Links
$(_doc_external("DMPlex/DMPlexPointQueueDestroy"))
"""
function DMPlexPointQueueDestroy(petsclib::PetscLibType, queue::Union{DMPlexPoCintQueue, Ref{DMPlexPoCintQueue}})
    error("DMPlexPointQueueDestroy: no generated method for these argument types")
end

@for_petsc function DMPlexPointQueueDestroy(petsclib::$UnionPetscLib, queue::Union{DMPlexPoCintQueue, Ref{DMPlexPoCintQueue}} )
	queue_ = queue isa Base.RefValue ? queue : Ref{DMPlexPoCintQueue}(queue)

    @chk ccall(
               (:DMPlexPointQueueDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{DMPlexPoCintQueue},),
               queue_,
              )


	return nothing
end 

"""
	empty::PetscBool = DMPlexPointQueueEmptyCollective(petsclib::PetscLibType, obj, queue::DMPlexPoCintQueue) 
Collectively determine whether a `DMPlexPointQueue` is empty on every rank of a communicator.

Collective

Input Parameters:
- `obj`   - a `PetscObject` whose communicator is used for the reduction
- `queue` - the queue

Output Parameter:
- `empty` - `PETSC_TRUE` if the queue is empty on every rank, `PETSC_FALSE` otherwise

Level: developer

See also: `DMPLEX`, `DMPlexPointQueue`, `DMPlexPointQueueEmpty()`

# External Links
$(_doc_external("DMPlex/DMPlexPointQueueEmptyCollective"))
"""
function DMPlexPointQueueEmptyCollective(petsclib::PetscLibType, obj, queue::DMPlexPoCintQueue)
    error("DMPlexPointQueueEmptyCollective: no generated method for these argument types")
end

@for_petsc function DMPlexPointQueueEmptyCollective(petsclib::$UnionPetscLib, obj, queue::DMPlexPoCintQueue )
	empty_ = Ref{PetscBool}()

    @chk ccall(
               (:DMPlexPointQueueEmptyCollective, $petsc_library),
               PetscErrorCode,
               (PetscObject, DMPlexPoCintQueue, Ptr{PetscBool}),
               obj, queue, empty_,
              )

	empty = empty_[]

	return empty
end 

"""
	DMPlexPointQueueEnqueue(petsclib::PetscLibType, queue::DMPlexPoCintQueue, p::PetscInt) 
Add a mesh point to the back of a `DMPlexPointQueue`.

Not Collective

Input Parameters:
- `queue` - the queue
- `p`     - the mesh point to enqueue

Level: developer

See also: `DMPLEX`, `DMPlexPointQueue`, `DMPlexPointQueueDequeue()`, `DMPlexPointQueueBack()`

# External Links
$(_doc_external("DMPlex/DMPlexPointQueueEnqueue"))
"""
function DMPlexPointQueueEnqueue(petsclib::PetscLibType, queue::DMPlexPoCintQueue, p::Integer)
    error("DMPlexPointQueueEnqueue: no generated method for these argument types")
end

@for_petsc function DMPlexPointQueueEnqueue(petsclib::$UnionPetscLib, queue::DMPlexPoCintQueue, p::$PetscInt )

    @chk ccall(
               (:DMPlexPointQueueEnqueue, $petsc_library),
               PetscErrorCode,
               (DMPlexPoCintQueue, $PetscInt),
               queue, p,
              )


	return nothing
end 

"""
	DMPlexPointQueueEnsureSize(petsclib::PetscLibType, queue::DMPlexPoCintQueue) 
Ensure that a `DMPlexPointQueue` has room for at least one more entry, doubling its capacity if it is full.

Not Collective

Input Parameter:
- `queue` - the queue

Level: developer

See also: `DMPLEX`, `DMPlexPointQueue`, `DMPlexPointQueueCreate()`, `DMPlexPointQueueEnqueue()`

# External Links
$(_doc_external("DMPlex/DMPlexPointQueueEnsureSize"))
"""
function DMPlexPointQueueEnsureSize(petsclib::PetscLibType, queue::DMPlexPoCintQueue)
    error("DMPlexPointQueueEnsureSize: no generated method for these argument types")
end

@for_petsc function DMPlexPointQueueEnsureSize(petsclib::$UnionPetscLib, queue::DMPlexPoCintQueue )

    @chk ccall(
               (:DMPlexPointQueueEnsureSize, $petsc_library),
               PetscErrorCode,
               (DMPlexPoCintQueue,),
               queue,
              )


	return nothing
end 

"""
	p::PetscInt = DMPlexPointQueueFront(petsclib::PetscLibType, queue::DMPlexPoCintQueue) 
Return, without removing, the mesh point at the front of a `DMPlexPointQueue`.

Not Collective

Input Parameter:
- `queue` - the queue

Output Parameter:
- `p` - the mesh point at the front of the queue

Level: developer

See also: `DMPLEX`, `DMPlexPointQueue`, `DMPlexPointQueueBack()`, `DMPlexPointQueueDequeue()`, `DMPlexPointQueueEmpty()`

# External Links
$(_doc_external("DMPlex/DMPlexPointQueueFront"))
"""
function DMPlexPointQueueFront(petsclib::PetscLibType, queue::DMPlexPoCintQueue)
    error("DMPlexPointQueueFront: no generated method for these argument types")
end

@for_petsc function DMPlexPointQueueFront(petsclib::$UnionPetscLib, queue::DMPlexPoCintQueue )
	p_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMPlexPointQueueFront, $petsc_library),
               PetscErrorCode,
               (DMPlexPoCintQueue, Ptr{$PetscInt}),
               queue, p_,
              )

	p = p_[]

	return p
end 

"""
	rdm::PetscDM = DMPlexTransformAdaptLabel(petsclib::PetscLibType, dm::AbstractPetscDM, metric::AbstractPetscVec, adaptLabel::DMLabel, rgLabel::DMLabel) 
Adapt a `DMPLEX` using a `DMPlexTransform` driven by a `DMLabel` marking cells to be refined or coarsened.

Collective

Input Parameters:
- `dm`         - the input `DMPLEX`
- `metric`     - unused; present to conform to the `DMAdaptor` label-based interface
- `adaptLabel` - a `DMLabel` marking cells with `DM_ADAPT_REFINE`, `DM_ADAPT_COARSEN`, etc.
- `rgLabel`    - unused region-tag label; present to conform to the `DMAdaptor` interface

Output Parameter:
- `rdm` - the adapted `DMPLEX`

Level: developer

See also: `DMPLEX`, `DMPlexTransform`, `DMAdaptLabel()`, `DMPlexTransformApply()`, `DMPlexTransformCreate()`, `DMLabel`

# External Links
$(_doc_external("DMPlex/DMPlexTransformAdaptLabel"))
"""
function DMPlexTransformAdaptLabel(petsclib::PetscLibType, dm::AbstractPetscDM, metric::AbstractPetscVec, adaptLabel::DMLabel, rgLabel::DMLabel)
    error("DMPlexTransformAdaptLabel: no generated method for these argument types")
end

@for_petsc function DMPlexTransformAdaptLabel(petsclib::$UnionPetscLib, dm::AbstractPetscDM, metric::AbstractPetscVec, adaptLabel::DMLabel, rgLabel::DMLabel )
	rdm_ = Ref{CDM}()

    @chk ccall(
               (:DMPlexTransformAdaptLabel, $petsc_library),
               PetscErrorCode,
               (CDM, CVec, DMLabel, DMLabel, Ptr{CDM}),
               dm, metric, adaptLabel, rgLabel, rdm_,
              )

	rdm = PetscDM(rdm_[], petsclib)

	return rdm
end 

"""
	trdm::PetscDM = DMPlexTransformApply(petsclib::PetscLibType, tr::DMPlexTransform, dm::AbstractPetscDM) 
Execute the transformation, producing another `DM`

Collective

Input Parameters:
- `tr` - The `DMPlexTransform` object
- `dm` - The original `DM`

Output Parameter:
- `trdm` - The transformed `DM`

Level: intermediate

Options Database Keys:
- `-dm_plex_transform_label_match_strata`    - Only label points of the same stratum as the producing point
- `-dm_plex_transform_label_replica_inc num` - Increment for the label value to be multiplied by the replica number
- `-dm_plex_transform_active name`           - Name for active mesh label

See also: [](plex_transform_table), `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformCreate()`, `DMPlexTransformSetDM()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformApply"))
"""
function DMPlexTransformApply(petsclib::PetscLibType, tr::DMPlexTransform, dm::AbstractPetscDM)
    error("DMPlexTransformApply: no generated method for these argument types")
end

@for_petsc function DMPlexTransformApply(petsclib::$UnionPetscLib, tr::DMPlexTransform, dm::AbstractPetscDM )
	trdm_ = Ref{CDM}()

    @chk ccall(
               (:DMPlexTransformApply, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, CDM, Ptr{CDM}),
               tr, dm, trdm_,
              )

	trdm = PetscDM(trdm_[], petsclib)

	return trdm
end 

"""
	rt::PetscInt,Nt::PetscInt,target::Vector{DMPolytopeType},size::Vector{PetscInt},cone::Ptr{PetscInt},ornt::Ptr{PetscInt} = DMPlexTransformCellTransform(petsclib::PetscLibType, tr::DMPlexTransform, source::DMPolytopeType, p::PetscInt) 
Describes the transform of a given source cell into a set of other target cells. These produced cells become the new mesh.

Input Parameters:
- `tr`     - The `DMPlexTransform` object
- `source` - The source cell type
- `p`      - The source point, which can also determine the refine type

Output Parameters:
- `rt`     - The refine type for this point
- `Nt`     - The number of types produced by this point
- `target` - An array of length `Nt` giving the types produced
- `size`   - An array of length `Nt` giving the number of cells of each type produced
- `cone`   - An array of length `Nt`*size[t]*coneSize[t] giving the cell type for each point in the cone of each produced point
- `ornt`   - An array of length `Nt`*size[t]*coneSize[t] giving the orientation for each point in the cone of each produced point

Level: advanced

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPolytopeType`, `DMPlexTransformApply()`, `DMPlexTransformCreate()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformCellTransform"))
"""
function DMPlexTransformCellTransform(petsclib::PetscLibType, tr::DMPlexTransform, source::DMPolytopeType, p::Integer)
    error("DMPlexTransformCellTransform: no generated method for these argument types")
end

@for_petsc function DMPlexTransformCellTransform(petsclib::$UnionPetscLib, tr::DMPlexTransform, source::DMPolytopeType, p::$PetscInt )
	rt_ = Ref{$PetscInt}()
	Nt_ = Ref{$PetscInt}()
	target_ = Ref{Ptr{DMPolytopeType}}()
	size_ = Ref{Ptr{$PetscInt}}()
	cone_ = Ref{Ptr{$PetscInt}}()
	ornt_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:DMPlexTransformCellTransform, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, DMPolytopeType, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{DMPolytopeType}}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}),
               tr, source, p, rt_, Nt_, target_, size_, cone_, ornt_,
              )

	rt = rt_[]
	Nt = Nt_[]
	cone = cone_[]
	ornt = ornt_[]
	target = target_[] == C_NULL ? DMPolytopeType[] : unsafe_wrap(Array, target_[], Nt; own = false)
	size = size_[] == C_NULL ? $PetscInt[] : unsafe_wrap(Array, size_[], Nt; own = false)

	return rt,Nt,target,size,cone,ornt
end 

"""
	rt::PetscInt,Nt::PetscInt,target::Vector{DMPolytopeType},size::Vector{PetscInt},cone::Ptr{PetscInt},ornt::Ptr{PetscInt} = DMPlexTransformCellTransformIdentity(petsclib::PetscLibType, tr::DMPlexTransform, source::DMPolytopeType, p::PetscInt) 
Default `celltransform` implementation for transforms that reproduce the input mesh

Not Collective

Input Parameters:
- `tr`     - The `DMPlexTransform`
- `source` - The cell type of the source point
- `p`      - The source point

Output Parameters:
- `rt`     - Refinement type of the source point (set to 0), or `NULL`
- `Nt`     - Number of target cell types produced (always 1)
- `target` - Array of produced cell types (a single-element array containing `source`)
- `size`   - Array of replica counts for each produced type (a single-element array containing 1)
- `cone`   - Cone description used by `DMPlexTransformGetCone()`; encodes that the replica takes the entire parent cone
- `ornt`   - Orientation array associated with `cone`; all zero for identity

Level: developer

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPolytopeType`, `DMPlexTransformCellTransform()`, `DMPlexTransformGetSubcellOrientationIdentity()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformCellTransformIdentity"))
"""
function DMPlexTransformCellTransformIdentity(petsclib::PetscLibType, tr::DMPlexTransform, source::DMPolytopeType, p::Integer)
    error("DMPlexTransformCellTransformIdentity: no generated method for these argument types")
end

@for_petsc function DMPlexTransformCellTransformIdentity(petsclib::$UnionPetscLib, tr::DMPlexTransform, source::DMPolytopeType, p::$PetscInt )
	rt_ = Ref{$PetscInt}()
	Nt_ = Ref{$PetscInt}()
	target_ = Ref{Ptr{DMPolytopeType}}()
	size_ = Ref{Ptr{$PetscInt}}()
	cone_ = Ref{Ptr{$PetscInt}}()
	ornt_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:DMPlexTransformCellTransformIdentity, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, DMPolytopeType, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{DMPolytopeType}}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}),
               tr, source, p, rt_, Nt_, target_, size_, cone_, ornt_,
              )

	rt = rt_[]
	Nt = Nt_[]
	cone = cone_[]
	ornt = ornt_[]
	target = target_[] == C_NULL ? DMPolytopeType[] : unsafe_wrap(Array, target_[], Nt; own = false)
	size = size_[] == C_NULL ? $PetscInt[] : unsafe_wrap(Array, size_[], Nt; own = false)

	return rt,Nt,target,size,cone,ornt
end 

"""
	useTensor::PetscBool = DMPlexTransformCohesiveExtrudeGetTensor(petsclib::PetscLibType, tr::DMPlexTransform) 
Get the flag to use tensor cells

Not Collective

Input Parameter:
- `tr` - The `DMPlexTransform`

Output Parameter:
- `useTensor` - The flag to use tensor cells

See also: `DMPlexTransform`, `DMPlexTransformCohesiveExtrudeSetTensor()`, `DMPlexTransformExtrudeGetTensor()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformCohesiveExtrudeGetTensor"))
"""
function DMPlexTransformCohesiveExtrudeGetTensor(petsclib::PetscLibType, tr::DMPlexTransform)
    error("DMPlexTransformCohesiveExtrudeGetTensor: no generated method for these argument types")
end

@for_petsc function DMPlexTransformCohesiveExtrudeGetTensor(petsclib::$UnionPetscLib, tr::DMPlexTransform )
	useTensor_ = Ref{PetscBool}()

    @chk ccall(
               (:DMPlexTransformCohesiveExtrudeGetTensor, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{PetscBool}),
               tr, useTensor_,
              )

	useTensor = useTensor_[]

	return useTensor
end 

"""
	unsplit::DMLabel = DMPlexTransformCohesiveExtrudeGetUnsplit(petsclib::PetscLibType, tr::DMPlexTransform) 
Get a new label marking the unsplit points in the transformed mesh

Not Collective

Input Parameter:
- `tr` - The `DMPlexTransform`

Output Parameter:
- `unsplit` - A new `DMLabel` marking the unsplit points in the transformed mesh

Level: intermediate

See also: `DMPlexTransform`, `DMPlexTransformGetTransformTypes()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformCohesiveExtrudeGetUnsplit"))
"""
function DMPlexTransformCohesiveExtrudeGetUnsplit(petsclib::PetscLibType, tr::DMPlexTransform)
    error("DMPlexTransformCohesiveExtrudeGetUnsplit: no generated method for these argument types")
end

@for_petsc function DMPlexTransformCohesiveExtrudeGetUnsplit(petsclib::$UnionPetscLib, tr::DMPlexTransform )
	unsplit_ = Ref{DMLabel}()

    @chk ccall(
               (:DMPlexTransformCohesiveExtrudeGetUnsplit, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{DMLabel}),
               tr, unsplit_,
              )

	unsplit = unsplit_[]

	return unsplit
end 

"""
	width::PetscReal = DMPlexTransformCohesiveExtrudeGetWidth(petsclib::PetscLibType, tr::DMPlexTransform) 
Get the width of extruded cells

Not Collective

Input Parameter:
- `tr` - The `DMPlexTransform`

Output Parameter:
- `width` - The width of extruded cells, or 0.

Level: intermediate

See also: `DMPlexTransform`, `DMPlexTransformCohesiveExtrudeSetWidth()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformCohesiveExtrudeGetWidth"))
"""
function DMPlexTransformCohesiveExtrudeGetWidth(petsclib::PetscLibType, tr::DMPlexTransform)
    error("DMPlexTransformCohesiveExtrudeGetWidth: no generated method for these argument types")
end

@for_petsc function DMPlexTransformCohesiveExtrudeGetWidth(petsclib::$UnionPetscLib, tr::DMPlexTransform )
	width_ = Ref{$PetscReal}()

    @chk ccall(
               (:DMPlexTransformCohesiveExtrudeGetWidth, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{$PetscReal}),
               tr, width_,
              )

	width = width_[]

	return width
end 

"""
	DMPlexTransformCohesiveExtrudeSetTensor(petsclib::PetscLibType, tr::DMPlexTransform, useTensor::PetscBool) 
Set the flag to use tensor cells

Not Collective

Input Parameters:
- `tr`        - The `DMPlexTransform`
- `useTensor` - The flag for tensor cells

See also: `DMPlexTransform`, `DMPlexTransformCohesiveExtrudeGetTensor()`, `DMPlexTransformExtrudeSetTensor()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformCohesiveExtrudeSetTensor"))
"""
function DMPlexTransformCohesiveExtrudeSetTensor(petsclib::PetscLibType, tr::DMPlexTransform, useTensor::PetscBool)
    error("DMPlexTransformCohesiveExtrudeSetTensor: no generated method for these argument types")
end

@for_petsc function DMPlexTransformCohesiveExtrudeSetTensor(petsclib::$UnionPetscLib, tr::DMPlexTransform, useTensor::PetscBool )

    @chk ccall(
               (:DMPlexTransformCohesiveExtrudeSetTensor, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, PetscBool),
               tr, useTensor,
              )


	return nothing
end 

"""
	DMPlexTransformCohesiveExtrudeSetWidth(petsclib::PetscLibType, tr::DMPlexTransform, width::PetscReal) 
Set the width of extruded cells

Not Collective

Input Parameters:
- `tr`    - The `DMPlexTransform`
- `width` - The width of the extruded cells, or 0.

Level: intermediate

See also: `DMPlexTransform`, `DMPlexTransformCohesiveExtrudeGetWidth()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformCohesiveExtrudeSetWidth"))
"""
function DMPlexTransformCohesiveExtrudeSetWidth(petsclib::PetscLibType, tr::DMPlexTransform, width::Real)
    error("DMPlexTransformCohesiveExtrudeSetWidth: no generated method for these argument types")
end

@for_petsc function DMPlexTransformCohesiveExtrudeSetWidth(petsclib::$UnionPetscLib, tr::DMPlexTransform, width::$PetscReal )

    @chk ccall(
               (:DMPlexTransformCohesiveExtrudeSetWidth, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, $PetscReal),
               tr, width,
              )


	return nothing
end 

"""
	tr::DMPlexTransform = DMPlexTransformCreate(petsclib::PetscLibType, comm::MPI_Comm) 
Creates an empty transform object. The type can then be set with `DMPlexTransformSetType()`.

Collective

Input Parameter:
- `comm` - The communicator for the transform object

Output Parameter:
- `tr` - The transform object

Level: beginner

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformType`, `DMPlexTransformSetType()`, `DMPLEXREFINEREGULAR`, `DMPLEXTRANSFORMFILTER`

# External Links
$(_doc_external("DMPlex/DMPlexTransformCreate"))
"""
function DMPlexTransformCreate(petsclib::PetscLibType, comm::MPI_Comm)
    error("DMPlexTransformCreate: no generated method for these argument types")
end

@for_petsc function DMPlexTransformCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	tr_ = Ref{DMPlexTransform}()

    @chk ccall(
               (:DMPlexTransformCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{DMPlexTransform}),
               comm, tr_,
              )

	tr = tr_[]

	return tr
end 

"""
	DMPlexTransformCreateDiscLabels(petsclib::PetscLibType, tr::DMPlexTransform, rdm::AbstractPetscDM) 
Refine the labels which define field and discrete system regions on the transformed `DM`

Not Collective

Input Parameters:
- `tr`  - The `DMPlexTransform`
- `rdm` - The refined `DM` produced by the transform

Level: developer

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformApply()`, `DMSetField()`, `DMSetRegionNumDS()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformCreateDiscLabels"))
"""
function DMPlexTransformCreateDiscLabels(petsclib::PetscLibType, tr::DMPlexTransform, rdm::AbstractPetscDM)
    error("DMPlexTransformCreateDiscLabels: no generated method for these argument types")
end

@for_petsc function DMPlexTransformCreateDiscLabels(petsclib::$UnionPetscLib, tr::DMPlexTransform, rdm::AbstractPetscDM )

    @chk ccall(
               (:DMPlexTransformCreateDiscLabels, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, CDM),
               tr, rdm,
              )


	return nothing
end 

"""
	DMPlexTransformDestroy(petsclib::PetscLibType, tr::Union{DMPlexTransform, Ref{DMPlexTransform}}) 
Destroys a `DMPlexTransform`

Collective

Input Parameter:
- `tr` - the transform object to destroy

Level: beginner

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformView()`, `DMPlexTransformCreate()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformDestroy"))
"""
function DMPlexTransformDestroy(petsclib::PetscLibType, tr::Union{DMPlexTransform, Ref{DMPlexTransform}})
    error("DMPlexTransformDestroy: no generated method for these argument types")
end

@for_petsc function DMPlexTransformDestroy(petsclib::$UnionPetscLib, tr::Union{DMPlexTransform, Ref{DMPlexTransform}} )
	tr_ = tr isa Base.RefValue ? tr : Ref{DMPlexTransform}(tr)

    @chk ccall(
               (:DMPlexTransformDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{DMPlexTransform},),
               tr_,
              )


	return nothing
end 

"""
	layers::PetscInt = DMPlexTransformExtrudeGetLayers(petsclib::PetscLibType, tr::DMPlexTransform) 
Get the number of extruded layers.

Not Collective

Input Parameter:
- `tr` - The `DMPlexTransform`

Output Parameter:
- `layers` - The number of layers

Level: intermediate

See also: `DMPlexTransform`, `DMPlexTransformExtrudeSetLayers()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformExtrudeGetLayers"))
"""
function DMPlexTransformExtrudeGetLayers(petsclib::PetscLibType, tr::DMPlexTransform)
    error("DMPlexTransformExtrudeGetLayers: no generated method for these argument types")
end

@for_petsc function DMPlexTransformExtrudeGetLayers(petsclib::$UnionPetscLib, tr::DMPlexTransform )
	layers_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMPlexTransformExtrudeGetLayers, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{$PetscInt}),
               tr, layers_,
              )

	layers = layers_[]

	return layers
end 

"""
	DMPlexTransformExtrudeGetNormal(petsclib::PetscLibType, tr::DMPlexTransform, normal::Vector{PetscReal}) 
Get the extrusion normal vector

Not Collective

Input Parameter:
- `tr` - The `DMPlexTransform`

Output Parameter:
- `normal` - The extrusion direction

See also: `DMPlexTransform`, `DMPlexTransformExtrudeSetNormal()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformExtrudeGetNormal"))
"""
function DMPlexTransformExtrudeGetNormal(petsclib::PetscLibType, tr::DMPlexTransform, normal::AbstractVector{<:Number})
    error("DMPlexTransformExtrudeGetNormal: no generated method for these argument types")
end

@for_petsc function DMPlexTransformExtrudeGetNormal(petsclib::$UnionPetscLib, tr::DMPlexTransform, normal::Vector{$PetscReal} )

    @chk ccall(
               (:DMPlexTransformExtrudeGetNormal, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{$PetscReal}),
               tr, normal,
              )


	return nothing
end 

"""
	periodic::PetscBool = DMPlexTransformExtrudeGetPeriodic(petsclib::PetscLibType, tr::DMPlexTransform) 
Get the flag to extrude periodically from the initial surface

Not Collective

Input Parameter:
- `tr` - The `DMPlexTransform`

Output Parameter:
- `periodic` - The flag to extrude periodically

Level: intermediate

See also: `DMPlexTransform`, `DMPlexTransformExtrudeSetPeriodic()`, `DMPlexTransformExtrudeGetSymmetric()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformExtrudeGetPeriodic"))
"""
function DMPlexTransformExtrudeGetPeriodic(petsclib::PetscLibType, tr::DMPlexTransform)
    error("DMPlexTransformExtrudeGetPeriodic: no generated method for these argument types")
end

@for_petsc function DMPlexTransformExtrudeGetPeriodic(petsclib::$UnionPetscLib, tr::DMPlexTransform )
	periodic_ = Ref{PetscBool}()

    @chk ccall(
               (:DMPlexTransformExtrudeGetPeriodic, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{PetscBool}),
               tr, periodic_,
              )

	periodic = periodic_[]

	return periodic
end 

"""
	symmetric::PetscBool = DMPlexTransformExtrudeGetSymmetric(petsclib::PetscLibType, tr::DMPlexTransform) 
Get the flag to extrude symmetrically from the initial surface

Not Collective

Input Parameter:
- `tr` - The `DMPlexTransform`

Output Parameter:
- `symmetric` - The flag to extrude symmetrically

Level: intermediate

See also: `DMPlexTransform`, `DMPlexTransformExtrudeSetSymmetric()`, `DMPlexTransformExtrudeGetPeriodic()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformExtrudeGetSymmetric"))
"""
function DMPlexTransformExtrudeGetSymmetric(petsclib::PetscLibType, tr::DMPlexTransform)
    error("DMPlexTransformExtrudeGetSymmetric: no generated method for these argument types")
end

@for_petsc function DMPlexTransformExtrudeGetSymmetric(petsclib::$UnionPetscLib, tr::DMPlexTransform )
	symmetric_ = Ref{PetscBool}()

    @chk ccall(
               (:DMPlexTransformExtrudeGetSymmetric, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{PetscBool}),
               tr, symmetric_,
              )

	symmetric = symmetric_[]

	return symmetric
end 

"""
	useTensor::PetscBool = DMPlexTransformExtrudeGetTensor(petsclib::PetscLibType, tr::DMPlexTransform) 
Get the flag to use tensor cells

Not Collective

Input Parameter:
- `tr` - The `DMPlexTransform`

Output Parameter:
- `useTensor` - The flag to use tensor cells

See also: `DMPlexTransform`, `DMPlexTransformExtrudeSetTensor()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformExtrudeGetTensor"))
"""
function DMPlexTransformExtrudeGetTensor(petsclib::PetscLibType, tr::DMPlexTransform)
    error("DMPlexTransformExtrudeGetTensor: no generated method for these argument types")
end

@for_petsc function DMPlexTransformExtrudeGetTensor(petsclib::$UnionPetscLib, tr::DMPlexTransform )
	useTensor_ = Ref{PetscBool}()

    @chk ccall(
               (:DMPlexTransformExtrudeGetTensor, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{PetscBool}),
               tr, useTensor_,
              )

	useTensor = useTensor_[]

	return useTensor
end 

"""
	thickness::PetscReal = DMPlexTransformExtrudeGetThickness(petsclib::PetscLibType, tr::DMPlexTransform) 
Get the total thickness of the layers

Not Collective

Input Parameter:
- `tr` - The `DMPlexTransform`

Output Parameter:
- `thickness` - The total thickness of the layers

Level: intermediate

See also: `DMPlexTransform`, `DMPlexTransformExtrudeSetThickness()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformExtrudeGetThickness"))
"""
function DMPlexTransformExtrudeGetThickness(petsclib::PetscLibType, tr::DMPlexTransform)
    error("DMPlexTransformExtrudeGetThickness: no generated method for these argument types")
end

@for_petsc function DMPlexTransformExtrudeGetThickness(petsclib::$UnionPetscLib, tr::DMPlexTransform )
	thickness_ = Ref{$PetscReal}()

    @chk ccall(
               (:DMPlexTransformExtrudeGetThickness, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{$PetscReal}),
               tr, thickness_,
              )

	thickness = thickness_[]

	return thickness
end 

"""
	DMPlexTransformExtrudeSetLayers(petsclib::PetscLibType, tr::DMPlexTransform, layers::PetscInt) 
Set the number of extruded layers.

Not Collective

Input Parameters:
- `tr`     - The `DMPlexTransform`
- `layers` - The number of layers

Level: intermediate

See also: `DMPlexTransform`, `DMPlexTransformExtrudeGetLayers()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformExtrudeSetLayers"))
"""
function DMPlexTransformExtrudeSetLayers(petsclib::PetscLibType, tr::DMPlexTransform, layers::Integer)
    error("DMPlexTransformExtrudeSetLayers: no generated method for these argument types")
end

@for_petsc function DMPlexTransformExtrudeSetLayers(petsclib::$UnionPetscLib, tr::DMPlexTransform, layers::$PetscInt )

    @chk ccall(
               (:DMPlexTransformExtrudeSetLayers, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, $PetscInt),
               tr, layers,
              )


	return nothing
end 

"""
	DMPlexTransformExtrudeSetNormal(petsclib::PetscLibType, tr::DMPlexTransform, normal::Vector{PetscReal}) 
Set the extrusion normal

Not Collective

Input Parameters:
- `tr`     - The `DMPlexTransform`
- `normal` - The extrusion direction

Level: intermediate

See also: `DMPlexTransform`, `DMPlexTransformExtrudeGetNormal()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformExtrudeSetNormal"))
"""
function DMPlexTransformExtrudeSetNormal(petsclib::PetscLibType, tr::DMPlexTransform, normal::AbstractVector{<:Number})
    error("DMPlexTransformExtrudeSetNormal: no generated method for these argument types")
end

@for_petsc function DMPlexTransformExtrudeSetNormal(petsclib::$UnionPetscLib, tr::DMPlexTransform, normal::Vector{$PetscReal} )

    @chk ccall(
               (:DMPlexTransformExtrudeSetNormal, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{$PetscReal}),
               tr, normal,
              )


	return nothing
end 

"""
	DMPlexTransformExtrudeSetNormalFunction(petsclib::PetscLibType, tr::DMPlexTransform, normalFunc::Ptr{Cvoid}) 
Set a function to determine the extrusion normal

Not Collective

Input Parameters:
- `tr`         - The `DMPlexTransform`
- `normalFunc` - A function determining the extrusion direction, see `PetscSimplePointFn` for the calling sequence

Level: intermediate

See also: `DMPlexTransform`, `DMPlexTransformExtrudeGetNormal()`, `PetscSimplePointFn`

# External Links
$(_doc_external("DMPlex/DMPlexTransformExtrudeSetNormalFunction"))
"""
function DMPlexTransformExtrudeSetNormalFunction(petsclib::PetscLibType, tr::DMPlexTransform, normalFunc::Ptr{Cvoid})
    error("DMPlexTransformExtrudeSetNormalFunction: no generated method for these argument types")
end

@for_petsc function DMPlexTransformExtrudeSetNormalFunction(petsclib::$UnionPetscLib, tr::DMPlexTransform, normalFunc::Ptr{Cvoid} )

    @chk ccall(
               (:DMPlexTransformExtrudeSetNormalFunction, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{Cvoid}),
               tr, normalFunc,
              )


	return nothing
end 

"""
	DMPlexTransformExtrudeSetPeriodic(petsclib::PetscLibType, tr::DMPlexTransform, periodic::PetscBool) 
Set the flag to extrude periodically from the initial surface

Not Collective

Input Parameters:
- `tr`       - The `DMPlexTransform`
- `periodic` - The flag to extrude periodically

Level: intermediate

See also: `DMPlexTransform`, `DMPlexTransformExtrudeGetPeriodic()`, `DMPlexTransformExtrudeSetSymmetric()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformExtrudeSetPeriodic"))
"""
function DMPlexTransformExtrudeSetPeriodic(petsclib::PetscLibType, tr::DMPlexTransform, periodic::PetscBool)
    error("DMPlexTransformExtrudeSetPeriodic: no generated method for these argument types")
end

@for_petsc function DMPlexTransformExtrudeSetPeriodic(petsclib::$UnionPetscLib, tr::DMPlexTransform, periodic::PetscBool )

    @chk ccall(
               (:DMPlexTransformExtrudeSetPeriodic, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, PetscBool),
               tr, periodic,
              )


	return nothing
end 

"""
	DMPlexTransformExtrudeSetSymmetric(petsclib::PetscLibType, tr::DMPlexTransform, symmetric::PetscBool) 
Set the flag to extrude symmetrically from the initial surface

Not Collective

Input Parameters:
- `tr`        - The `DMPlexTransform`
- `symmetric` - The flag to extrude symmetrically

Level: intermediate

See also: `DMPlexTransform`, `DMPlexTransformExtrudeGetSymmetric()`, `DMPlexTransformExtrudeSetPeriodic()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformExtrudeSetSymmetric"))
"""
function DMPlexTransformExtrudeSetSymmetric(petsclib::PetscLibType, tr::DMPlexTransform, symmetric::PetscBool)
    error("DMPlexTransformExtrudeSetSymmetric: no generated method for these argument types")
end

@for_petsc function DMPlexTransformExtrudeSetSymmetric(petsclib::$UnionPetscLib, tr::DMPlexTransform, symmetric::PetscBool )

    @chk ccall(
               (:DMPlexTransformExtrudeSetSymmetric, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, PetscBool),
               tr, symmetric,
              )


	return nothing
end 

"""
	DMPlexTransformExtrudeSetTensor(petsclib::PetscLibType, tr::DMPlexTransform, useTensor::PetscBool) 
Set the flag to use tensor cells

Not Collective

Input Parameters:
- `tr`        - The `DMPlexTransform`
- `useTensor` - The flag for tensor cells

See also: `DMPlexTransform`, `DMPlexTransformExtrudeGetTensor()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformExtrudeSetTensor"))
"""
function DMPlexTransformExtrudeSetTensor(petsclib::PetscLibType, tr::DMPlexTransform, useTensor::PetscBool)
    error("DMPlexTransformExtrudeSetTensor: no generated method for these argument types")
end

@for_petsc function DMPlexTransformExtrudeSetTensor(petsclib::$UnionPetscLib, tr::DMPlexTransform, useTensor::PetscBool )

    @chk ccall(
               (:DMPlexTransformExtrudeSetTensor, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, PetscBool),
               tr, useTensor,
              )


	return nothing
end 

"""
	DMPlexTransformExtrudeSetThickness(petsclib::PetscLibType, tr::DMPlexTransform, thickness::PetscReal) 
Set the total thickness of the layers

Not Collective

Input Parameters:
- `tr`        - The `DMPlexTransform`
- `thickness` - The total thickness of the layers

Level: intermediate

See also: `DMPlexTransform`, `DMPlexTransformExtrudeGetThickness()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformExtrudeSetThickness"))
"""
function DMPlexTransformExtrudeSetThickness(petsclib::PetscLibType, tr::DMPlexTransform, thickness::Real)
    error("DMPlexTransformExtrudeSetThickness: no generated method for these argument types")
end

@for_petsc function DMPlexTransformExtrudeSetThickness(petsclib::$UnionPetscLib, tr::DMPlexTransform, thickness::$PetscReal )

    @chk ccall(
               (:DMPlexTransformExtrudeSetThickness, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, $PetscReal),
               tr, thickness,
              )


	return nothing
end 

"""
	DMPlexTransformExtrudeSetThicknesses(petsclib::PetscLibType, tr::DMPlexTransform, Nth::PetscInt, thicknesses::Vector{PetscReal}) 
Set the thickness of each layer

Not Collective

Input Parameters:
- `tr`          - The `DMPlexTransform`
- `Nth`         - The number of thicknesses
- `thicknesses` - The array of thicknesses

Level: intermediate

See also: `DMPlexTransform`, `DMPlexTransformExtrudeSetThickness()`, `DMPlexTransformExtrudeGetThickness()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformExtrudeSetThicknesses"))
"""
function DMPlexTransformExtrudeSetThicknesses(petsclib::PetscLibType, tr::DMPlexTransform, Nth::Integer, thicknesses::AbstractVector{<:Number})
    error("DMPlexTransformExtrudeSetThicknesses: no generated method for these argument types")
end

@for_petsc function DMPlexTransformExtrudeSetThicknesses(petsclib::$UnionPetscLib, tr::DMPlexTransform, Nth::$PetscInt, thicknesses::Vector{$PetscReal} )

    @chk ccall(
               (:DMPlexTransformExtrudeSetThicknesses, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, $PetscInt, Ptr{$PetscReal}),
               tr, Nth, thicknesses,
              )


	return nothing
end 

"""
	active::DMLabel = DMPlexTransformGetActive(petsclib::PetscLibType, tr::DMPlexTransform) 
Get the `DMLabel` marking the active points for the transform

Input Parameter:
- `tr` - The `DMPlexTransform` object

Output Parameter:
- `active` - The `DMLabel` indicating which points will be transformed

Level: intermediate

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformSetActive()`, `DMPlexTransformApply()`, `DMPlexTransformCreate()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformGetActive"))
"""
function DMPlexTransformGetActive(petsclib::PetscLibType, tr::DMPlexTransform)
    error("DMPlexTransformGetActive: no generated method for these argument types")
end

@for_petsc function DMPlexTransformGetActive(petsclib::$UnionPetscLib, tr::DMPlexTransform )
	active_ = Ref{DMLabel}()

    @chk ccall(
               (:DMPlexTransformGetActive, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{DMLabel}),
               tr, active_,
              )

	active = active_[]

	return active
end 

"""
	celltype::DMPolytopeType = DMPlexTransformGetCellType(petsclib::PetscLibType, tr::DMPlexTransform, cell::PetscInt) 
Return the cell type for a point in the transformed mesh

Not Collective

Input Parameters:
- `tr`   - The `DMPlexTransform`
- `cell` - The point number in the transformed mesh

Output Parameter:
- `celltype` - The `DMPolytopeType` of the point

Level: developer

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPolytopeType`, `DMPlexTransformGetChart()`, `DMPlexTransformGetCellTypeStratum()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformGetCellType"))
"""
function DMPlexTransformGetCellType(petsclib::PetscLibType, tr::DMPlexTransform, cell::Integer)
    error("DMPlexTransformGetCellType: no generated method for these argument types")
end

@for_petsc function DMPlexTransformGetCellType(petsclib::$UnionPetscLib, tr::DMPlexTransform, cell::$PetscInt )
	celltype_ = Ref{DMPolytopeType}()

    @chk ccall(
               (:DMPlexTransformGetCellType, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, $PetscInt, Ptr{DMPolytopeType}),
               tr, cell, celltype_,
              )

	celltype = celltype_[]

	return celltype
end 

"""
	start::PetscInt,end_::PetscInt = DMPlexTransformGetCellTypeStratum(petsclib::PetscLibType, tr::DMPlexTransform, celltype::DMPolytopeType) 
Return the point range for a given cell type in the transformed mesh

Not Collective

Input Parameters:
- `tr`       - The `DMPlexTransform`
- `celltype` - The `DMPolytopeType` of the requested stratum

Output Parameters:
- `start` - The first point of the stratum, or `NULL` if not needed
- `end`   - One past the last point of the stratum, or `NULL` if not needed

Level: developer

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPolytopeType`, `DMPlexTransformGetCellType()`, `DMPlexTransformGetChart()`, `DMPlexGetDepthStratum()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformGetCellTypeStratum"))
"""
function DMPlexTransformGetCellTypeStratum(petsclib::PetscLibType, tr::DMPlexTransform, celltype::DMPolytopeType)
    error("DMPlexTransformGetCellTypeStratum: no generated method for these argument types")
end

@for_petsc function DMPlexTransformGetCellTypeStratum(petsclib::$UnionPetscLib, tr::DMPlexTransform, celltype::DMPolytopeType )
	start_ = Ref{$PetscInt}()
	end__ = Ref{$PetscInt}()

    @chk ccall(
               (:DMPlexTransformGetCellTypeStratum, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, DMPolytopeType, Ptr{$PetscInt}, Ptr{$PetscInt}),
               tr, celltype, start_, end__,
              )

	start = start_[]
	end_ = end__[]

	return start,end_
end 

"""
	Nv::PetscInt,trVerts::Ptr{PetscScalar} = DMPlexTransformGetCellVertices(petsclib::PetscLibType, tr::DMPlexTransform, ct::DMPolytopeType) 
Get the set of transformed vertices lying in the closure of a reference cell of given type

Input Parameters:
- `tr` - The `DMPlexTransform` object
- `ct` - The cell type

Output Parameters:
- `Nv`      - The number of transformed vertices in the closure of the reference cell of given type
- `trVerts` - The coordinates of these vertices in the reference cell

Level: developer

See also: `DMPLEX`, `DMPlexTransform`, `DMPolytopeType`, `DMPlexTransformGetSubcellVertices()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformGetCellVertices"))
"""
function DMPlexTransformGetCellVertices(petsclib::PetscLibType, tr::DMPlexTransform, ct::DMPolytopeType)
    error("DMPlexTransformGetCellVertices: no generated method for these argument types")
end

@for_petsc function DMPlexTransformGetCellVertices(petsclib::$UnionPetscLib, tr::DMPlexTransform, ct::DMPolytopeType )
	Nv_ = Ref{$PetscInt}()
	trVerts_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:DMPlexTransformGetCellVertices, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, DMPolytopeType, Ptr{$PetscInt}, Ptr{Ptr{$PetscScalar}}),
               tr, ct, Nv_, trVerts_,
              )

	Nv = Nv_[]
	trVerts = trVerts_[]

	return Nv,trVerts
end 

"""
	pStart::PetscInt,pEnd::PetscInt = DMPlexTransformGetChart(petsclib::PetscLibType, tr::DMPlexTransform) 
Get the chart `[pStart, pEnd)` for the points produced by the transform

Not Collective

Input Parameter:
- `tr` - The `DMPlexTransform`

Output Parameters:
- `pStart` - The first point in the transformed mesh, or `NULL` if not needed
- `pEnd`   - One past the last point in the transformed mesh, or `NULL` if not needed

Level: developer

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformApply()`, `DMPlexTransformGetCellType()`, `DMPlexTransformGetCellTypeStratum()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformGetChart"))
"""
function DMPlexTransformGetChart(petsclib::PetscLibType, tr::DMPlexTransform)
    error("DMPlexTransformGetChart: no generated method for these argument types")
end

@for_petsc function DMPlexTransformGetChart(petsclib::$UnionPetscLib, tr::DMPlexTransform )
	pStart_ = Ref{$PetscInt}()
	pEnd_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMPlexTransformGetChart, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{$PetscInt}, Ptr{$PetscInt}),
               tr, pStart_, pEnd_,
              )

	pStart = pStart_[]
	pEnd = pEnd_[]

	return pStart,pEnd
end 

"""
	cone::Ptr{PetscInt},ornt::Ptr{PetscInt} = DMPlexTransformGetCone(petsclib::PetscLibType, tr::DMPlexTransform, q::PetscInt) 
Return the cone of a point in the transformed mesh

Not Collective

Input Parameters:
- `tr` - The `DMPlexTransform`
- `q`  - The point number in the transformed mesh

Output Parameters:
- `cone` - The cone points, obtained from an internal work array, or `NULL` if not requested
- `ornt` - The orientations of the cone points, obtained from an internal work array, or `NULL` if not requested

Level: developer

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformRestoreCone()`, `DMPlexTransformGetConeOriented()`, `DMPlexTransformGetConeSize()`, `DMPlexGetCone()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformGetCone"))
"""
function DMPlexTransformGetCone(petsclib::PetscLibType, tr::DMPlexTransform, q::Integer)
    error("DMPlexTransformGetCone: no generated method for these argument types")
end

@for_petsc function DMPlexTransformGetCone(petsclib::$UnionPetscLib, tr::DMPlexTransform, q::$PetscInt )
	cone_ = Ref{Ptr{$PetscInt}}()
	ornt_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:DMPlexTransformGetCone, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, $PetscInt, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}),
               tr, q, cone_, ornt_,
              )

	cone = cone_[]
	ornt = ornt_[]

	return cone,ornt
end 

"""
	cone::Ptr{PetscInt},ornt::Ptr{PetscInt} = DMPlexTransformGetConeOriented(petsclib::PetscLibType, tr::DMPlexTransform, q::PetscInt, po::PetscInt) 
Return the cone of a point in the transformed mesh, computed using a specified parent orientation

Not Collective

Input Parameters:
- `tr` - The `DMPlexTransform`
- `q`  - The point number in the transformed mesh
- `po` - The orientation of the parent cell in the original mesh to use when producing the cone

Output Parameters:
- `cone` - The cone points, obtained from an internal work array
- `ornt` - The orientations of the cone points, obtained from an internal work array

Level: developer

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformGetCone()`, `DMPlexTransformRestoreCone()`, `DMPlexTransformGetConeSize()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformGetConeOriented"))
"""
function DMPlexTransformGetConeOriented(petsclib::PetscLibType, tr::DMPlexTransform, q::Integer, po::Integer)
    error("DMPlexTransformGetConeOriented: no generated method for these argument types")
end

@for_petsc function DMPlexTransformGetConeOriented(petsclib::$UnionPetscLib, tr::DMPlexTransform, q::$PetscInt, po::$PetscInt )
	cone_ = Ref{Ptr{$PetscInt}}()
	ornt_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:DMPlexTransformGetConeOriented, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, $PetscInt, $PetscInt, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}),
               tr, q, po, cone_, ornt_,
              )

	cone = cone_[]
	ornt = ornt_[]

	return cone,ornt
end 

"""
	coneSize::PetscInt = DMPlexTransformGetConeSize(petsclib::PetscLibType, tr::DMPlexTransform, q::PetscInt) 
Return the cone size of a point in the transformed mesh

Not Collective

Input Parameters:
- `tr` - The `DMPlexTransform`
- `q`  - The point number in the transformed mesh

Output Parameter:
- `coneSize` - The number of points in the cone of `q`

Level: developer

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformGetCone()`, `DMPlexTransformGetCellType()`, `DMPlexGetConeSize()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformGetConeSize"))
"""
function DMPlexTransformGetConeSize(petsclib::PetscLibType, tr::DMPlexTransform, q::Integer)
    error("DMPlexTransformGetConeSize: no generated method for these argument types")
end

@for_petsc function DMPlexTransformGetConeSize(petsclib::$UnionPetscLib, tr::DMPlexTransform, q::$PetscInt )
	coneSize_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMPlexTransformGetConeSize, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, $PetscInt, Ptr{$PetscInt}),
               tr, q, coneSize_,
              )

	coneSize = coneSize_[]

	return coneSize
end 

"""
	dm::PetscDM = DMPlexTransformGetDM(petsclib::PetscLibType, tr::DMPlexTransform) 
Get the base `DM` for the transform

Input Parameter:
- `tr` - The `DMPlexTransform` object

Output Parameter:
- `dm` - The original `DM` which will be transformed

Level: intermediate

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformSetDM()`, `DMPlexTransformApply()`, `DMPlexTransformCreate()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformGetDM"))
"""
function DMPlexTransformGetDM(petsclib::PetscLibType, tr::DMPlexTransform)
    error("DMPlexTransformGetDM: no generated method for these argument types")
end

@for_petsc function DMPlexTransformGetDM(petsclib::$UnionPetscLib, tr::DMPlexTransform )
	dm_ = Ref{CDM}()

    @chk ccall(
               (:DMPlexTransformGetDM, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{CDM}),
               tr, dm_,
              )

	dm = PetscDM(dm_[], petsclib)

	return dm
end 

"""
	depth::PetscInt = DMPlexTransformGetDepth(petsclib::PetscLibType, tr::DMPlexTransform) 
Return the topological depth of the transformed mesh

Not Collective

Input Parameter:
- `tr` - The `DMPlexTransform`

Output Parameter:
- `depth` - The depth of the transformed mesh

Level: developer

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformGetDepthStratum()`, `DMPlexGetDepth()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformGetDepth"))
"""
function DMPlexTransformGetDepth(petsclib::PetscLibType, tr::DMPlexTransform)
    error("DMPlexTransformGetDepth: no generated method for these argument types")
end

@for_petsc function DMPlexTransformGetDepth(petsclib::$UnionPetscLib, tr::DMPlexTransform )
	depth_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMPlexTransformGetDepth, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{$PetscInt}),
               tr, depth_,
              )

	depth = depth_[]

	return depth
end 

"""
	start::PetscInt,end_::PetscInt = DMPlexTransformGetDepthStratum(petsclib::PetscLibType, tr::DMPlexTransform, depth::PetscInt) 
Return the point range for a given depth in the transformed mesh

Not Collective

Input Parameters:
- `tr`    - The `DMPlexTransform`
- `depth` - The requested depth in the transformed mesh

Output Parameters:
- `start` - The first point at the given depth, or `NULL` if not needed
- `end`   - One past the last point at the given depth, or `NULL` if not needed

Level: developer

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformGetDepth()`, `DMPlexGetDepthStratum()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformGetDepthStratum"))
"""
function DMPlexTransformGetDepthStratum(petsclib::PetscLibType, tr::DMPlexTransform, depth::Integer)
    error("DMPlexTransformGetDepthStratum: no generated method for these argument types")
end

@for_petsc function DMPlexTransformGetDepthStratum(petsclib::$UnionPetscLib, tr::DMPlexTransform, depth::$PetscInt )
	start_ = Ref{$PetscInt}()
	end__ = Ref{$PetscInt}()

    @chk ccall(
               (:DMPlexTransformGetDepthStratum, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
               tr, depth, start_, end__,
              )

	start = start_[]
	end_ = end__[]

	return start,end_
end 

"""
	match::PetscBool = DMPlexTransformGetMatchStrata(petsclib::PetscLibType, tr::DMPlexTransform) 
Get the flag which determines what points get added to the transformed labels

Not Collective

Input Parameter:
- `tr` - The `DMPlexTransform`

Output Parameter:
- `match` - If `PETSC_TRUE`, only add produced points at the same stratum as the original point to new labels

Level: intermediate

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformSetMatchStrata()`, `DMPlexGetPointDepth()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformGetMatchStrata"))
"""
function DMPlexTransformGetMatchStrata(petsclib::PetscLibType, tr::DMPlexTransform)
    error("DMPlexTransformGetMatchStrata: no generated method for these argument types")
end

@for_petsc function DMPlexTransformGetMatchStrata(petsclib::$UnionPetscLib, tr::DMPlexTransform )
	match_ = Ref{PetscBool}()

    @chk ccall(
               (:DMPlexTransformGetMatchStrata, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{PetscBool}),
               tr, match_,
              )

	match = match_[]

	return match
end 

"""
	ct::DMPolytopeType,ctNew::DMPolytopeType,p::PetscInt,r::PetscInt = DMPlexTransformGetSourcePoint(petsclib::PetscLibType, tr::DMPlexTransform, pNew::PetscInt) 
Get the number of a point in the original mesh based on information from the transformed mesh.

Not Collective

Input Parameters:
- `tr`   - The `DMPlexTransform`
- `pNew` - The new point number

Output Parameters:
- `ct`    - The type of the original point which produces the new point
- `ctNew` - The type of the new point
- `p`     - The original point which produces the new point
- `r`     - The replica number of the new point, meaning it is the rth point of type ctNew produced from p

Level: developer

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPolytopeType`, `DMPlexTransformGetTargetPoint()`, `DMPlexTransformCellTransform()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformGetSourcePoint"))
"""
function DMPlexTransformGetSourcePoint(petsclib::PetscLibType, tr::DMPlexTransform, pNew::Integer)
    error("DMPlexTransformGetSourcePoint: no generated method for these argument types")
end

@for_petsc function DMPlexTransformGetSourcePoint(petsclib::$UnionPetscLib, tr::DMPlexTransform, pNew::$PetscInt )
	ct_ = Ref{DMPolytopeType}()
	ctNew_ = Ref{DMPolytopeType}()
	p_ = Ref{$PetscInt}()
	r_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMPlexTransformGetSourcePoint, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, $PetscInt, Ptr{DMPolytopeType}, Ptr{DMPolytopeType}, Ptr{$PetscInt}, Ptr{$PetscInt}),
               tr, pNew, ct_, ctNew_, p_, r_,
              )

	ct = ct_[]
	ctNew = ctNew_[]
	p = p_[]
	r = r_[]

	return ct,ctNew,p,r
end 

"""
	rnew::PetscInt,onew::PetscInt = DMPlexTransformGetSubcellOrientation(petsclib::PetscLibType, tr::DMPlexTransform, sct::DMPolytopeType, sp::PetscInt, so::PetscInt, tct::DMPolytopeType, r::PetscInt, o::PetscInt) 
Transform the replica number and orientation for a target point according to the group action for the source point

Not Collective

Input Parameters:
- `tr`  - The `DMPlexTransform`
- `sct` - The source point cell type, from whom the new cell is being produced
- `sp`  - The source point
- `so`  - The orientation of the source point in its enclosing parent
- `tct` - The target point cell type
- `r`   - The replica number requested for the produced cell type
- `o`   - The orientation of the replica

Output Parameters:
- `rnew` - The replica number, given the orientation of the parent
- `onew` - The replica orientation, given the orientation of the parent

Level: advanced

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPolytopeType`, `DMPlexTransformCellTransform()`, `DMPlexTransformApply()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformGetSubcellOrientation"))
"""
function DMPlexTransformGetSubcellOrientation(petsclib::PetscLibType, tr::DMPlexTransform, sct::DMPolytopeType, sp::Integer, so::Integer, tct::DMPolytopeType, r::Integer, o::Integer)
    error("DMPlexTransformGetSubcellOrientation: no generated method for these argument types")
end

@for_petsc function DMPlexTransformGetSubcellOrientation(petsclib::$UnionPetscLib, tr::DMPlexTransform, sct::DMPolytopeType, sp::$PetscInt, so::$PetscInt, tct::DMPolytopeType, r::$PetscInt, o::$PetscInt )
	rnew_ = Ref{$PetscInt}()
	onew_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMPlexTransformGetSubcellOrientation, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, DMPolytopeType, $PetscInt, $PetscInt, DMPolytopeType, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
               tr, sct, sp, so, tct, r, o, rnew_, onew_,
              )

	rnew = rnew_[]
	onew = onew_[]

	return rnew,onew
end 

"""
	rnew::PetscInt,onew::PetscInt = DMPlexTransformGetSubcellOrientationIdentity(petsclib::PetscLibType, tr::DMPlexTransform, sct::DMPolytopeType, sp::PetscInt, so::PetscInt, tct::DMPolytopeType, r::PetscInt, o::PetscInt) 
Default `getsubcellorientation` implementation for transforms that reproduce the input mesh

Not Collective

Input Parameters:
- `tr`  - The `DMPlexTransform`
- `sct` - The source point cell type
- `sp`  - The source point
- `so`  - The orientation of the source point in its enclosing parent
- `tct` - The target point cell type
- `r`   - The replica number requested for the produced cell type
- `o`   - The orientation of the replica

Output Parameters:
- `rnew` - The replica number, given the orientation of the parent (returns `r`)
- `onew` - The replica orientation composed with the source orientation

Level: developer

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformGetSubcellOrientation()`, `DMPlexTransformCellTransformIdentity()`, `DMPolytopeTypeComposeOrientation()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformGetSubcellOrientationIdentity"))
"""
function DMPlexTransformGetSubcellOrientationIdentity(petsclib::PetscLibType, tr::DMPlexTransform, sct::DMPolytopeType, sp::Integer, so::Integer, tct::DMPolytopeType, r::Integer, o::Integer)
    error("DMPlexTransformGetSubcellOrientationIdentity: no generated method for these argument types")
end

@for_petsc function DMPlexTransformGetSubcellOrientationIdentity(petsclib::$UnionPetscLib, tr::DMPlexTransform, sct::DMPolytopeType, sp::$PetscInt, so::$PetscInt, tct::DMPolytopeType, r::$PetscInt, o::$PetscInt )
	rnew_ = Ref{$PetscInt}()
	onew_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMPlexTransformGetSubcellOrientationIdentity, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, DMPolytopeType, $PetscInt, $PetscInt, DMPolytopeType, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
               tr, sct, sp, so, tct, r, o, rnew_, onew_,
              )

	rnew = rnew_[]
	onew = onew_[]

	return rnew,onew
end 

"""
	subVerts::Ptr{PetscInt} = DMPlexTransformGetSubcellVertices(petsclib::PetscLibType, tr::DMPlexTransform, ct::DMPolytopeType, rct::DMPolytopeType, r::PetscInt) 
Get the set of transformed vertices defining a subcell in the reference cell of given type

Input Parameters:
- `tr`  - The `DMPlexTransform` object
- `ct`  - The cell type
- `rct` - The subcell type
- `r`   - The subcell index

Output Parameter:
- `subVerts` - The indices of these vertices in the set of vertices returned by `DMPlexTransformGetCellVertices()`

Level: developer

See also: `DMPLEX`, `DMPlexTransform`, `DMPolytopeType`, `DMPlexTransformGetCellVertices()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformGetSubcellVertices"))
"""
function DMPlexTransformGetSubcellVertices(petsclib::PetscLibType, tr::DMPlexTransform, ct::DMPolytopeType, rct::DMPolytopeType, r::Integer)
    error("DMPlexTransformGetSubcellVertices: no generated method for these argument types")
end

@for_petsc function DMPlexTransformGetSubcellVertices(petsclib::$UnionPetscLib, tr::DMPlexTransform, ct::DMPolytopeType, rct::DMPolytopeType, r::$PetscInt )
	subVerts_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:DMPlexTransformGetSubcellVertices, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, DMPolytopeType, DMPolytopeType, $PetscInt, Ptr{Ptr{$PetscInt}}),
               tr, ct, rct, r, subVerts_,
              )

	subVerts = subVerts_[]

	return subVerts
end 

"""
	pNew::PetscInt = DMPlexTransformGetTargetPoint(petsclib::PetscLibType, tr::DMPlexTransform, ct::DMPolytopeType, ctNew::DMPolytopeType, p::PetscInt, r::PetscInt) 
Get the number of a point in the transformed mesh based on information from the original mesh.

Not Collective

Input Parameters:
- `tr`    - The `DMPlexTransform`
- `ct`    - The type of the original point which produces the new point
- `ctNew` - The type of the new point
- `p`     - The original point which produces the new point
- `r`     - The replica number of the new point, meaning it is the rth point of type `ctNew` produced from `p`

Output Parameter:
- `pNew` - The new point number

Level: developer

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPolytopeType`, `DMPlexTransformGetSourcePoint()`, `DMPlexTransformCellTransform()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformGetTargetPoint"))
"""
function DMPlexTransformGetTargetPoint(petsclib::PetscLibType, tr::DMPlexTransform, ct::DMPolytopeType, ctNew::DMPolytopeType, p::Integer, r::Integer)
    error("DMPlexTransformGetTargetPoint: no generated method for these argument types")
end

@for_petsc function DMPlexTransformGetTargetPoint(petsclib::$UnionPetscLib, tr::DMPlexTransform, ct::DMPolytopeType, ctNew::DMPolytopeType, p::$PetscInt, r::$PetscInt )
	pNew_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMPlexTransformGetTargetPoint, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, DMPolytopeType, DMPolytopeType, $PetscInt, $PetscInt, Ptr{$PetscInt}),
               tr, ct, ctNew, p, r, pNew_,
              )

	pNew = pNew_[]

	return pNew
end 

"""
	trType::DMLabel = DMPlexTransformGetTransformTypes(petsclib::PetscLibType, tr::DMPlexTransform) 
Get the `DMLabel` marking the transform type of each point for the transform

Input Parameter:
- `tr` - The `DMPlexTransform` object

Output Parameter:
- `trType` - The `DMLabel` indicating the transform type for each point

Level: intermediate

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexSetTransformType()`, `DMPlexTransformGetActive()`, `DMPlexTransformApply()`, `DMPlexTransformCreate()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformGetTransformTypes"))
"""
function DMPlexTransformGetTransformTypes(petsclib::PetscLibType, tr::DMPlexTransform)
    error("DMPlexTransformGetTransformTypes: no generated method for these argument types")
end

@for_petsc function DMPlexTransformGetTransformTypes(petsclib::$UnionPetscLib, tr::DMPlexTransform )
	trType_ = Ref{DMLabel}()

    @chk ccall(
               (:DMPlexTransformGetTransformTypes, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{DMLabel}),
               tr, trType_,
              )

	trType = trType_[]

	return trType
end 

"""
	type::DMPlexTransformType = DMPlexTransformGetType(petsclib::PetscLibType, tr::DMPlexTransform) 
Gets the type name (as a string) from the transform.

Not Collective

Input Parameter:
- `tr` - The `DMPlexTransform`

Output Parameter:
- `type` - The `DMPlexTransformType` name

Level: intermediate

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformType`, `DMPlexTransformSetType()`, `DMPlexTransformCreate()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformGetType"))
"""
function DMPlexTransformGetType(petsclib::PetscLibType, tr::DMPlexTransform)
    error("DMPlexTransformGetType: no generated method for these argument types")
end

@for_petsc function DMPlexTransformGetType(petsclib::$UnionPetscLib, tr::DMPlexTransform )
	type_ = Ref{DMPlexTransformType}()

    @chk ccall(
               (:DMPlexTransformGetType, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, Ptr{DMPlexTransformType}),
               tr, type_,
              )

	type = type_[] == C_NULL ? "" : unsafe_string(type_[])

	return type
end 

"""
	DMPlexTransformMapCoordinates(petsclib::PetscLibType, tr::DMPlexTransform, pct::DMPolytopeType, ct::DMPolytopeType, p::PetscInt, r::PetscInt, Nv::PetscInt, dE::PetscInt, in::Vector{PetscScalar}, out::Vector{PetscScalar}) 
Calculate new coordinates for produced points

Not collective

Input Parameters:
- `tr`  - The `DMPlexTransform`
- `pct` - The cell type of the parent, from whom the new cell is being produced
- `ct`  - The type being produced
- `p`   - The original point
- `r`   - The replica number requested for the produced cell type
- `Nv`  - Number of vertices in the closure of the parent cell
- `dE`  - Spatial dimension
- `in`  - array of size Nv*dE, holding coordinates of the vertices in the closure of the parent cell

Output Parameter:
- `out` - The coordinates of the new vertices

Level: intermediate

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPolytopeType`, `DMPlexTransformApply()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformMapCoordinates"))
"""
function DMPlexTransformMapCoordinates(petsclib::PetscLibType, tr::DMPlexTransform, pct::DMPolytopeType, ct::DMPolytopeType, p::Integer, r::Integer, Nv::Integer, dE::Integer, in::AbstractVector{<:Number}, out::AbstractVector{<:Number})
    error("DMPlexTransformMapCoordinates: no generated method for these argument types")
end

@for_petsc function DMPlexTransformMapCoordinates(petsclib::$UnionPetscLib, tr::DMPlexTransform, pct::DMPolytopeType, ct::DMPolytopeType, p::$PetscInt, r::$PetscInt, Nv::$PetscInt, dE::$PetscInt, in::Vector{$PetscScalar}, out::Vector{$PetscScalar} )

    @chk ccall(
               (:DMPlexTransformMapCoordinates, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, DMPolytopeType, DMPolytopeType, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscScalar}, Ptr{$PetscScalar}),
               tr, pct, ct, p, r, Nv, dE, in, out,
              )


	return nothing
end 

"""
	DMPlexTransformRegister(petsclib::PetscLibType, name::String, noname::Ptr{Cvoid}) 
Adds a new transform component implementation

Not Collective

Input Parameters:
- `name`        - The name of a new user-defined creation routine
- `create_func` - The creation routine

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformRegisterAll()`, `DMPlexTransformRegisterDestroy()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformRegister"))
"""
function DMPlexTransformRegister(petsclib::PetscLibType, name::String, noname::Ptr{Cvoid})
    error("DMPlexTransformRegister: no generated method for these argument types")
end

@for_petsc function DMPlexTransformRegister(petsclib::$UnionPetscLib, name::String, noname::Ptr{Cvoid} )

    @chk ccall(
               (:DMPlexTransformRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{Cvoid}),
               name, noname,
              )


	return nothing
end 

"""
	DMPlexTransformRegisterAll(petsclib::PetscLibType) 
Registers all of the transform components in the `DM` package.

Not Collective

Level: advanced

See also: `DM`, `DMPLEX`, `DMPlexTransformType`, `DMRegisterAll()`, `DMPlexTransformRegisterDestroy()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformRegisterAll"))
"""
function DMPlexTransformRegisterAll(petsclib::PetscLibType)
    error("DMPlexTransformRegisterAll: no generated method for these argument types")
end

@for_petsc function DMPlexTransformRegisterAll(petsclib::$UnionPetscLib)

    @chk ccall(
               (:DMPlexTransformRegisterAll, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	DMPlexTransformRegisterDestroy(petsclib::PetscLibType) 
This function destroys the registered `DMPlexTransformType`. It is called from `PetscFinalize()`.

Not collective

Level: developer

See also: `DM`, `DMPLEX`, `DMRegisterAll()`, `DMPlexTransformType`, `PetscInitialize()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformRegisterDestroy"))
"""
function DMPlexTransformRegisterDestroy(petsclib::PetscLibType)
    error("DMPlexTransformRegisterDestroy: no generated method for these argument types")
end

@for_petsc function DMPlexTransformRegisterDestroy(petsclib::$UnionPetscLib)

    @chk ccall(
               (:DMPlexTransformRegisterDestroy, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	DMPlexTransformRestoreCone(petsclib::PetscLibType, tr::DMPlexTransform, q::PetscInt, cone::Union{Ptr, AbstractArray{PetscInt}}, ornt::Union{Ptr, AbstractArray{PetscInt}}) 
Return the work arrays produced by `DMPlexTransformGetCone()` or `DMPlexTransformGetConeOriented()`

Not Collective

Input Parameters:
- `tr`   - The `DMPlexTransform`
- `q`    - The point number in the transformed mesh
- `cone` - The cone points to release, or `NULL`
- `ornt` - The orientations to release, or `NULL`

Level: developer

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformGetCone()`, `DMPlexTransformGetConeOriented()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformRestoreCone"))
"""
function DMPlexTransformRestoreCone(petsclib::PetscLibType, tr::DMPlexTransform, q::Integer, cone::Union{Ptr, AbstractArray{<:Number}}, ornt::Union{Ptr, AbstractArray{<:Number}})
    error("DMPlexTransformRestoreCone: no generated method for these argument types")
end

@for_petsc function DMPlexTransformRestoreCone(petsclib::$UnionPetscLib, tr::DMPlexTransform, q::$PetscInt, cone::Union{Ptr, AbstractArray{$PetscInt}}, ornt::Union{Ptr, AbstractArray{$PetscInt}} )
	cone_ = Ref{Ptr{$PetscInt}}(cone isa Ptr ? cone : pointer(cone))
	ornt_ = Ref{Ptr{$PetscInt}}(ornt isa Ptr ? ornt : pointer(ornt))

    @chk ccall(
               (:DMPlexTransformRestoreCone, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, $PetscInt, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{$PetscInt}}),
               tr, q, cone_, ornt_,
              )


	return nothing
end 

"""
	DMPlexTransformSetActive(petsclib::PetscLibType, tr::DMPlexTransform, active::DMLabel) 
Set the `DMLabel` marking the active points for the transform

Input Parameters:
- `tr`     - The `DMPlexTransform` object
- `active` - The `DMLabel` indicating which points will be transformed

Level: intermediate

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformGetActive()`, `DMPlexTransformApply()`, `DMPlexTransformCreate()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformSetActive"))
"""
function DMPlexTransformSetActive(petsclib::PetscLibType, tr::DMPlexTransform, active::DMLabel)
    error("DMPlexTransformSetActive: no generated method for these argument types")
end

@for_petsc function DMPlexTransformSetActive(petsclib::$UnionPetscLib, tr::DMPlexTransform, active::DMLabel )

    @chk ccall(
               (:DMPlexTransformSetActive, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, DMLabel),
               tr, active,
              )


	return nothing
end 

"""
	DMPlexTransformSetDM(petsclib::PetscLibType, tr::DMPlexTransform, dm::AbstractPetscDM) 
Set the base `DM` for the transform

Input Parameters:
- `tr` - The `DMPlexTransform` object
- `dm` - The original `DM` which will be transformed

Level: intermediate

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformGetDM()`, `DMPlexTransformApply()`, `DMPlexTransformCreate()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformSetDM"))
"""
function DMPlexTransformSetDM(petsclib::PetscLibType, tr::DMPlexTransform, dm::AbstractPetscDM)
    error("DMPlexTransformSetDM: no generated method for these argument types")
end

@for_petsc function DMPlexTransformSetDM(petsclib::$UnionPetscLib, tr::DMPlexTransform, dm::AbstractPetscDM )

    @chk ccall(
               (:DMPlexTransformSetDM, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, CDM),
               tr, dm,
              )


	return nothing
end 

"""
	DMPlexTransformSetDimensions(petsclib::PetscLibType, tr::DMPlexTransform, dm::AbstractPetscDM, trdm::AbstractPetscDM) 
Set the dimensions for the transformed `DM`

Input Parameters:
- `tr` - The `DMPlexTransform` object
- `dm` - The original `DM`

Output Parameter:
- `trdm` - The transformed `DM`

Level: advanced

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformApply()`, `DMPlexTransformCreate()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformSetDimensions"))
"""
function DMPlexTransformSetDimensions(petsclib::PetscLibType, tr::DMPlexTransform, dm::AbstractPetscDM, trdm::AbstractPetscDM)
    error("DMPlexTransformSetDimensions: no generated method for these argument types")
end

@for_petsc function DMPlexTransformSetDimensions(petsclib::$UnionPetscLib, tr::DMPlexTransform, dm::AbstractPetscDM, trdm::AbstractPetscDM )

    @chk ccall(
               (:DMPlexTransformSetDimensions, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, CDM, CDM),
               tr, dm, trdm,
              )


	return nothing
end 

"""
	DMPlexTransformSetFromOptions(petsclib::PetscLibType, tr::DMPlexTransform) 
Sets parameters in a transform from values in the options database

Collective

Input Parameter:
- `tr` - the `DMPlexTransform` object to set options for

Options Database Keys:
- `-dm_plex_transform_type type`               - Set the transform type, e.g. refine_regular
- `-dm_plex_transform_label_match_strata`      - Only label points of the same stratum as the producing point
- `-dm_plex_transform_label_replica_inc inc`   - Increment for the label value to be multiplied by the replica number, so that the new label value is oldValue + r * inc
- `-dm_plex_transform_active name`             - Name for active mesh label
- `-dm_plex_transform_active_values v0,v1,...` - Values in the active label

Level: intermediate

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformView()`, `DMPlexTransformCreate()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformSetFromOptions"))
"""
function DMPlexTransformSetFromOptions(petsclib::PetscLibType, tr::DMPlexTransform)
    error("DMPlexTransformSetFromOptions: no generated method for these argument types")
end

@for_petsc function DMPlexTransformSetFromOptions(petsclib::$UnionPetscLib, tr::DMPlexTransform )

    @chk ccall(
               (:DMPlexTransformSetFromOptions, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform,),
               tr,
              )


	return nothing
end 

"""
	DMPlexTransformSetMatchStrata(petsclib::PetscLibType, tr::DMPlexTransform, match::PetscBool) 
Set the flag which determines what points get added to the transformed labels

Not Collective

Input Parameters:
- `tr`    - The `DMPlexTransform`
- `match` - If `PETSC_TRUE`, only add produced points at the same stratum as the original point to new labels

Level: intermediate

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformGetMatchStrata()`, `DMPlexGetPointDepth()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformSetMatchStrata"))
"""
function DMPlexTransformSetMatchStrata(petsclib::PetscLibType, tr::DMPlexTransform, match::PetscBool)
    error("DMPlexTransformSetMatchStrata: no generated method for these argument types")
end

@for_petsc function DMPlexTransformSetMatchStrata(petsclib::$UnionPetscLib, tr::DMPlexTransform, match::PetscBool )

    @chk ccall(
               (:DMPlexTransformSetMatchStrata, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, PetscBool),
               tr, match,
              )


	return nothing
end 

"""
	DMPlexTransformSetTransformTypes(petsclib::PetscLibType, tr::DMPlexTransform, trType::DMLabel) 
Set the `DMLabel` marking the transform type of each point for the transform

Input Parameters:
- `tr`     - The `DMPlexTransform` object
- `trType` - The original `DM` which will be transformed

Level: intermediate

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformGetTransformTypes()`, `DMPlexTransformGetActive())`, `DMPlexTransformApply()`, `DMPlexTransformCreate()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformSetTransformTypes"))
"""
function DMPlexTransformSetTransformTypes(petsclib::PetscLibType, tr::DMPlexTransform, trType::DMLabel)
    error("DMPlexTransformSetTransformTypes: no generated method for these argument types")
end

@for_petsc function DMPlexTransformSetTransformTypes(petsclib::$UnionPetscLib, tr::DMPlexTransform, trType::DMLabel )

    @chk ccall(
               (:DMPlexTransformSetTransformTypes, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, DMLabel),
               tr, trType,
              )


	return nothing
end 

"""
	DMPlexTransformSetType(petsclib::PetscLibType, tr::DMPlexTransform, method::DMPlexTransformType) 
Sets the particular implementation for a transform.

Collective

Input Parameters:
- `tr`     - The transform
- `method` - The name of the transform type

Options Database Key:
- `-dm_plex_transform_type type` - Sets the transform type; see `DMPlexTransformType`

Level: intermediate

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformType`, `DMPlexTransformGetType()`, `DMPlexTransformCreate()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformSetType"))
"""
function DMPlexTransformSetType(petsclib::PetscLibType, tr::DMPlexTransform, method::DMPlexTransformType)
    error("DMPlexTransformSetType: no generated method for these argument types")
end

@for_petsc function DMPlexTransformSetType(petsclib::$UnionPetscLib, tr::DMPlexTransform, method::DMPlexTransformType )

    @chk ccall(
               (:DMPlexTransformSetType, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, DMPlexTransformType),
               tr, method,
              )


	return nothing
end 

"""
	DMPlexTransformSetUp(petsclib::PetscLibType, tr::DMPlexTransform) 
Create the tables that drive the transform

Input Parameter:
- `tr` - The `DMPlexTransform` object

Level: intermediate

See also: [](plex_transform_table), `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformApply()`, `DMPlexTransformCreate()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformSetUp"))
"""
function DMPlexTransformSetUp(petsclib::PetscLibType, tr::DMPlexTransform)
    error("DMPlexTransformSetUp: no generated method for these argument types")
end

@for_petsc function DMPlexTransformSetUp(petsclib::$UnionPetscLib, tr::DMPlexTransform )

    @chk ccall(
               (:DMPlexTransformSetUp, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform,),
               tr,
              )


	return nothing
end 

"""
	DMPlexTransformView(petsclib::PetscLibType, tr::DMPlexTransform, v::PetscViewer) 
Views a `DMPlexTransform`

Collective

Input Parameters:
- `tr` - the `DMPlexTransform` object to view
- `v`  - the viewer

Level: beginner

See also: `DM`, `DMPLEX`, `DMPlexTransform`, `DMPlexTransformType`, `PetscViewer`, `DMPlexTransformDestroy()`, `DMPlexTransformCreate()`

# External Links
$(_doc_external("DMPlex/DMPlexTransformView"))
"""
function DMPlexTransformView(petsclib::PetscLibType, tr::DMPlexTransform, v::PetscViewer)
    error("DMPlexTransformView: no generated method for these argument types")
end

@for_petsc function DMPlexTransformView(petsclib::$UnionPetscLib, tr::DMPlexTransform, v::PetscViewer )

    @chk ccall(
               (:DMPlexTransformView, $petsc_library),
               PetscErrorCode,
               (DMPlexTransform, PetscViewer),
               tr, v,
              )


	return nothing
end 

"""
	celldm::DMSwarmCellDM = DMSwarmCellDMCreate(petsclib::PetscLibType, dm::AbstractPetscDM, Nf::PetscInt, dmFields::String, Nfc::PetscInt, coordFields::String) 
create a `DMSwarmCellDM`

Collective

Input Parameters:
- `dm`          - The background `DM` for the `DMSwarm`
- `Nf`          - The number of swarm fields defined over `dm`
- `dmFields`    - The swarm field names for the `dm` fields
- `Nfc`         - The number of swarm fields to use for coordinates over `dm`
- `coordFields` - The swarm field names for the `dm` coordinate fields

Output Parameter:
- `celldm` - The new `DMSwarmCellDM`

Level: advanced

See also: `DMSwarmCellDM`, `DMSWARM`, `DMSetType()`

# External Links
$(_doc_external("DMSwarm/DMSwarmCellDMCreate"))
"""
function DMSwarmCellDMCreate(petsclib::PetscLibType, dm::AbstractPetscDM, Nf::Integer, dmFields::String, Nfc::Integer, coordFields::String)
    error("DMSwarmCellDMCreate: no generated method for these argument types")
end

@for_petsc function DMSwarmCellDMCreate(petsclib::$UnionPetscLib, dm::AbstractPetscDM, Nf::$PetscInt, dmFields::String, Nfc::$PetscInt, coordFields::String )
	dmFields_ = Ref{Ptr{Cchar}}(dmFields isa Ptr ? dmFields : pointer(dmFields))
	coordFields_ = Ref{Ptr{Cchar}}(coordFields isa Ptr ? coordFields : pointer(coordFields))
	celldm_ = Ref{DMSwarmCellDM}()

    @chk ccall(
               (:DMSwarmCellDMCreate, $petsc_library),
               PetscErrorCode,
               (CDM, $PetscInt, Ptr{Ptr{Cchar}}, $PetscInt, Ptr{Ptr{Cchar}}, Ptr{DMSwarmCellDM}),
               dm, Nf, dmFields_, Nfc, coordFields_, celldm_,
              )

	celldm = celldm_[]

	return celldm
end 

"""
	DMSwarmCellDMDestroy(petsclib::PetscLibType, celldm::Union{DMSwarmCellDM, Ref{DMSwarmCellDM}}) 
destroy a `DMSwarmCellDM`

Collective

Input Parameter:
- `celldm` - address of `DMSwarmCellDM`

Level: advanced

See also: `DMSwarmCellDM`, `DMSwarmCellDMCreate()`

# External Links
$(_doc_external("DMSwarm/DMSwarmCellDMDestroy"))
"""
function DMSwarmCellDMDestroy(petsclib::PetscLibType, celldm::Union{DMSwarmCellDM, Ref{DMSwarmCellDM}})
    error("DMSwarmCellDMDestroy: no generated method for these argument types")
end

@for_petsc function DMSwarmCellDMDestroy(petsclib::$UnionPetscLib, celldm::Union{DMSwarmCellDM, Ref{DMSwarmCellDM}} )
	celldm_ = celldm isa Base.RefValue ? celldm : Ref{DMSwarmCellDM}(celldm)

    @chk ccall(
               (:DMSwarmCellDMDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{DMSwarmCellDM},),
               celldm_,
              )


	return nothing
end 

"""
	bs::PetscInt = DMSwarmCellDMGetBlockSize(petsclib::PetscLibType, celldm::DMSwarmCellDM, sw::AbstractPetscDM) 
Returns the total blocksize for the `DM` fields

Not Collective

Input Parameters:
- `celldm` - The `DMSwarmCellDM` object
- `sw`     - The `DMSwarm` object

Output Parameter:
- `bs` - The total block size

Level: intermediate

See also: `DMSwarmCellDM`, `DM`, `DMSwarmSetCellDM()`

# External Links
$(_doc_external("DMSwarm/DMSwarmCellDMGetBlockSize"))
"""
function DMSwarmCellDMGetBlockSize(petsclib::PetscLibType, celldm::DMSwarmCellDM, sw::AbstractPetscDM)
    error("DMSwarmCellDMGetBlockSize: no generated method for these argument types")
end

@for_petsc function DMSwarmCellDMGetBlockSize(petsclib::$UnionPetscLib, celldm::DMSwarmCellDM, sw::AbstractPetscDM )
	bs_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMSwarmCellDMGetBlockSize, $petsc_library),
               PetscErrorCode,
               (DMSwarmCellDM, CDM, Ptr{$PetscInt}),
               celldm, sw, bs_,
              )

	bs = bs_[]

	return bs
end 

"""
	cellid::Ptr{Cchar} = DMSwarmCellDMGetCellID(petsclib::PetscLibType, celldm::DMSwarmCellDM) 
Returns the cell id field name for the `DMSwarm`

Not Collective

Input Parameter:
- `celldm` - The `DMSwarmCellDM` object

Output Parameters:
- `cellid` - The cell id field name in the `DMSWARM`

Level: intermediate

See also: `DMSwarmCellDM`, `DM`, `DMSwarmSetCellDM()`

# External Links
$(_doc_external("DMSwarm/DMSwarmCellDMGetCellID"))
"""
function DMSwarmCellDMGetCellID(petsclib::PetscLibType, celldm::DMSwarmCellDM)
    error("DMSwarmCellDMGetCellID: no generated method for these argument types")
end

@for_petsc function DMSwarmCellDMGetCellID(petsclib::$UnionPetscLib, celldm::DMSwarmCellDM )
	cellid_ = Ref{Ptr{Cchar}}()

    @chk ccall(
               (:DMSwarmCellDMGetCellID, $petsc_library),
               PetscErrorCode,
               (DMSwarmCellDM, Ptr{Ptr{Cchar}}),
               celldm, cellid_,
              )

	cellid = cellid_[]

	return cellid
end 

"""
	Nfc::PetscInt,names::Vector{String} = DMSwarmCellDMGetCoordinateFields(petsclib::PetscLibType, celldm::DMSwarmCellDM) 
Returns the `DM` coordinate fields for the `DMSwarm`

Not Collective

Input Parameter:
- `celldm` - The `DMSwarmCellDM` object

Output Parameters:
- `Nfc`   - The number of coordinate fields
- `names` - The array of coordinate field names in the `DMSWARM`

Level: intermediate

See also: `DMSwarmCellDM`, `DM`, `DMSwarmSetCellDM()`

# External Links
$(_doc_external("DMSwarm/DMSwarmCellDMGetCoordinateFields"))
"""
function DMSwarmCellDMGetCoordinateFields(petsclib::PetscLibType, celldm::DMSwarmCellDM)
    error("DMSwarmCellDMGetCoordinateFields: no generated method for these argument types")
end

@for_petsc function DMSwarmCellDMGetCoordinateFields(petsclib::$UnionPetscLib, celldm::DMSwarmCellDM )
	Nfc_ = Ref{$PetscInt}()
	names_ = Ref{Ptr{Ptr{Cchar}}}()

    @chk ccall(
               (:DMSwarmCellDMGetCoordinateFields, $petsc_library),
               PetscErrorCode,
               (DMSwarmCellDM, Ptr{$PetscInt}, Ptr{Ptr{Ptr{Cchar}}}),
               celldm, Nfc_, names_,
              )

	Nfc = Nfc_[]
	names = names_[] == C_NULL ? String[] : [unsafe_string(p) for p in unsafe_wrap(Array, names_[], Nfc; own = false)]

	return Nfc,names
end 

"""
	dm::PetscDM = DMSwarmCellDMGetDM(petsclib::PetscLibType, celldm::DMSwarmCellDM) 
Returns the background `DM` for the `DMSwarm`

Not Collective

Input Parameter:
- `celldm` - The `DMSwarmCellDM` object

Output Parameter:
- `dm` - The `DM` object

Level: intermediate

See also: `DMSwarmCellDM`, `DM`, `DMSwarmSetCellDM()`

# External Links
$(_doc_external("DMSwarm/DMSwarmCellDMGetDM"))
"""
function DMSwarmCellDMGetDM(petsclib::PetscLibType, celldm::DMSwarmCellDM)
    error("DMSwarmCellDMGetDM: no generated method for these argument types")
end

@for_petsc function DMSwarmCellDMGetDM(petsclib::$UnionPetscLib, celldm::DMSwarmCellDM )
	dm_ = Ref{CDM}()

    @chk ccall(
               (:DMSwarmCellDMGetDM, $petsc_library),
               PetscErrorCode,
               (DMSwarmCellDM, Ptr{CDM}),
               celldm, dm_,
              )

	dm = PetscDM(dm_[], petsclib)

	return dm
end 

"""
	Nf::PetscInt,names::Vector{String} = DMSwarmCellDMGetFields(petsclib::PetscLibType, celldm::DMSwarmCellDM) 
Returns the `DM` fields for the `DMSwarm`

Not Collective

Input Parameter:
- `celldm` - The `DMSwarmCellDM` object

Output Parameters:
- `Nf`    - The number of fields
- `names` - The array of field names in the `DMSWARM`

Level: intermediate

See also: `DMSwarmCellDM`, `DM`, `DMSwarmSetCellDM()`

# External Links
$(_doc_external("DMSwarm/DMSwarmCellDMGetFields"))
"""
function DMSwarmCellDMGetFields(petsclib::PetscLibType, celldm::DMSwarmCellDM)
    error("DMSwarmCellDMGetFields: no generated method for these argument types")
end

@for_petsc function DMSwarmCellDMGetFields(petsclib::$UnionPetscLib, celldm::DMSwarmCellDM )
	Nf_ = Ref{$PetscInt}()
	names_ = Ref{Ptr{Ptr{Cchar}}}()

    @chk ccall(
               (:DMSwarmCellDMGetFields, $petsc_library),
               PetscErrorCode,
               (DMSwarmCellDM, Ptr{$PetscInt}, Ptr{Ptr{Ptr{Cchar}}}),
               celldm, Nf_, names_,
              )

	Nf = Nf_[]
	names = names_[] == C_NULL ? String[] : [unsafe_string(p) for p in unsafe_wrap(Array, names_[], Nf; own = false)]

	return Nf,names
end 

"""
	sort::DMSwarmSort = DMSwarmCellDMGetSort(petsclib::PetscLibType, celldm::DMSwarmCellDM) 
Returns the sort context over the active `DMSwarmCellDM` for the `DMSwarm`

Not Collective

Input Parameter:
- `celldm` - The `DMSwarmCellDM` object

Output Parameter:
- `sort` - The `DMSwarmSort` object

Level: intermediate

See also: `DMSwarmCellDM`, `DM`, `DMSwarmSetCellDM()`

# External Links
$(_doc_external("DMSwarm/DMSwarmCellDMGetSort"))
"""
function DMSwarmCellDMGetSort(petsclib::PetscLibType, celldm::DMSwarmCellDM)
    error("DMSwarmCellDMGetSort: no generated method for these argument types")
end

@for_petsc function DMSwarmCellDMGetSort(petsclib::$UnionPetscLib, celldm::DMSwarmCellDM )
	sort_ = Ref{DMSwarmSort}()

    @chk ccall(
               (:DMSwarmCellDMGetSort, $petsc_library),
               PetscErrorCode,
               (DMSwarmCellDM, Ptr{DMSwarmSort}),
               celldm, sort_,
              )

	sort = sort_[]

	return sort
end 

"""
	DMSwarmCellDMSetSort(petsclib::PetscLibType, celldm::DMSwarmCellDM, sort::DMSwarmSort) 
Sets the sort context over the active `DMSwarmCellDM` for the `DMSwarm`

Not Collective

Input Parameters:
- `celldm` - The `DMSwarmCellDM` object
- `sort`   - The `DMSwarmSort` object

Level: intermediate

See also: `DMSwarmCellDM`, `DM`, `DMSwarmSetCellDM()`

# External Links
$(_doc_external("DMSwarm/DMSwarmCellDMSetSort"))
"""
function DMSwarmCellDMSetSort(petsclib::PetscLibType, celldm::DMSwarmCellDM, sort::DMSwarmSort)
    error("DMSwarmCellDMSetSort: no generated method for these argument types")
end

@for_petsc function DMSwarmCellDMSetSort(petsclib::$UnionPetscLib, celldm::DMSwarmCellDM, sort::DMSwarmSort )

    @chk ccall(
               (:DMSwarmCellDMSetSort, $petsc_library),
               PetscErrorCode,
               (DMSwarmCellDM, DMSwarmSort),
               celldm, sort,
              )


	return nothing
end 

"""
	DMSwarmCellDMView(petsclib::PetscLibType, celldm::DMSwarmCellDM, viewer::PetscViewer) 
view a `DMSwarmCellDM`

Collective

Input Parameters:
- `celldm` - `DMSwarmCellDM`
- `viewer` - viewer to display field, for example `PETSC_VIEWER_STDOUT_WORLD`

Level: advanced

See also: `DMSwarmCellDM`, `DMSwarmCellDMCreate()`

# External Links
$(_doc_external("DMSwarm/DMSwarmCellDMView"))
"""
function DMSwarmCellDMView(petsclib::PetscLibType, celldm::DMSwarmCellDM, viewer::PetscViewer)
    error("DMSwarmCellDMView: no generated method for these argument types")
end

@for_petsc function DMSwarmCellDMView(petsclib::$UnionPetscLib, celldm::DMSwarmCellDM, viewer::PetscViewer )

    @chk ccall(
               (:DMSwarmCellDMView, $petsc_library),
               PetscErrorCode,
               (DMSwarmCellDM, PetscViewer),
               celldm, viewer,
              )


	return nothing
end 

"""
	gfield::DMSwarmDataField = DMSwarmDataBucketGetDMSwarmDataFieldByName(petsclib::PetscLibType, db::DMSwarmDataBucket, name::String) 
Return the `DMSwarmDataField` handle registered in a `DMSwarmDataBucket` under a given name.

Not Collective

Input Parameters:
- `db`   - the `DMSwarmDataBucket`
- `name` - the field name

Output Parameter:
- `gfield` - the `DMSwarmDataField` handle

Level: developer

See also: `DMSwarmDataBucket`, `DMSwarmDataField`, `DMSwarmDataBucketGetDMSwarmDataFieldIdByName()`, `DMSwarmDataBucketQueryDMSwarmDataFieldByName()`

# External Links
$(_doc_external("DMSwarm/DMSwarmDataBucketGetDMSwarmDataFieldByName"))
"""
function DMSwarmDataBucketGetDMSwarmDataFieldByName(petsclib::PetscLibType, db::DMSwarmDataBucket, name::String)
    error("DMSwarmDataBucketGetDMSwarmDataFieldByName: no generated method for these argument types")
end

@for_petsc function DMSwarmDataBucketGetDMSwarmDataFieldByName(petsclib::$UnionPetscLib, db::DMSwarmDataBucket, name::String )
	gfield_ = Ref{DMSwarmDataField}()

    @chk ccall(
               (:DMSwarmDataBucketGetDMSwarmDataFieldByName, $petsc_library),
               PetscErrorCode,
               (DMSwarmDataBucket, Ptr{Cchar}, Ptr{DMSwarmDataField}),
               db, name, gfield_,
              )

	gfield = gfield_[]

	return gfield
end 

"""
	idx::PetscInt = DMSwarmDataBucketGetDMSwarmDataFieldIdByName(petsclib::PetscLibType, db::DMSwarmDataBucket, name::String) 
Return the index of a `DMSwarmDataField` within a `DMSwarmDataBucket` given its name.

Not Collective

Input Parameters:
- `db`   - the `DMSwarmDataBucket`
- `name` - the field name

Output Parameter:
- `idx` - the index of the field within the bucket

Level: developer

See also: `DMSwarmDataBucket`, `DMSwarmDataField`, `DMSwarmDataBucketGetDMSwarmDataFieldByName()`, `DMSwarmDataBucketQueryDMSwarmDataFieldByName()`

# External Links
$(_doc_external("DMSwarm/DMSwarmDataBucketGetDMSwarmDataFieldIdByName"))
"""
function DMSwarmDataBucketGetDMSwarmDataFieldIdByName(petsclib::PetscLibType, db::DMSwarmDataBucket, name::String)
    error("DMSwarmDataBucketGetDMSwarmDataFieldIdByName: no generated method for these argument types")
end

@for_petsc function DMSwarmDataBucketGetDMSwarmDataFieldIdByName(petsclib::$UnionPetscLib, db::DMSwarmDataBucket, name::String )
	idx_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMSwarmDataBucketGetDMSwarmDataFieldIdByName, $petsc_library),
               PetscErrorCode,
               (DMSwarmDataBucket, Ptr{Cchar}, Ptr{$PetscInt}),
               db, name, idx_,
              )

	idx = idx_[]

	return idx
end 

"""
	found::PetscBool = DMSwarmDataBucketQueryDMSwarmDataFieldByName(petsclib::PetscLibType, db::DMSwarmDataBucket, name::String) 
Test whether a `DMSwarmDataBucket` contains a `DMSwarmDataField` with the given name.

Not Collective

Input Parameters:
- `db`   - the `DMSwarmDataBucket`
- `name` - the field name to look up

Output Parameter:
- `found` - `PETSC_TRUE` if a field with the given name is registered, otherwise `PETSC_FALSE`

Level: developer

See also: `DMSwarmDataBucket`, `DMSwarmDataField`, `DMSwarmDataBucketGetDMSwarmDataFieldByName()`, `DMSwarmDataBucketGetDMSwarmDataFieldIdByName()`

# External Links
$(_doc_external("DMSwarm/DMSwarmDataBucketQueryDMSwarmDataFieldByName"))
"""
function DMSwarmDataBucketQueryDMSwarmDataFieldByName(petsclib::PetscLibType, db::DMSwarmDataBucket, name::String)
    error("DMSwarmDataBucketQueryDMSwarmDataFieldByName: no generated method for these argument types")
end

@for_petsc function DMSwarmDataBucketQueryDMSwarmDataFieldByName(petsclib::$UnionPetscLib, db::DMSwarmDataBucket, name::String )
	found_ = Ref{PetscBool}()

    @chk ccall(
               (:DMSwarmDataBucketQueryDMSwarmDataFieldByName, $petsc_library),
               PetscErrorCode,
               (DMSwarmDataBucket, Ptr{Cchar}, Ptr{PetscBool}),
               db, name, found_,
              )

	found = found_[]

	return found
end 

"""
	data::Ptr{Cvoid} = DMSwarmDataFieldGetEntries(petsclib::PetscLibType, gfield::DMSwarmDataField) 
Return a pointer to the raw contiguous storage backing a `DMSwarmDataField`.

Not Collective

Input Parameter:
- `gfield` - the `DMSwarmDataField`

Output Parameter:
- `data` - pointer to the raw entries; must be released with `DMSwarmDataFieldRestoreEntries()`

Level: developer

See also: `DMSwarmDataField`, `DMSwarmDataFieldRestoreEntries()`, `DMSwarmDataFieldGetNumEntries()`, `DMSwarmDataFieldGetAtomicSize()`

# External Links
$(_doc_external("DMSwarm/DMSwarmDataFieldGetEntries"))
"""
function DMSwarmDataFieldGetEntries(petsclib::PetscLibType, gfield::DMSwarmDataField)
    error("DMSwarmDataFieldGetEntries: no generated method for these argument types")
end

@for_petsc function DMSwarmDataFieldGetEntries(petsclib::$UnionPetscLib, gfield::DMSwarmDataField )
	data_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:DMSwarmDataFieldGetEntries, $petsc_library),
               PetscErrorCode,
               (DMSwarmDataField, Ptr{Ptr{Cvoid}}),
               gfield, data_,
              )

	data = data_[]

	return data
end 

"""
	data::Ptr{Cvoid} = DMSwarmDataFieldRestoreEntries(petsclib::PetscLibType, gfield::DMSwarmDataField) 
Release a pointer obtained from `DMSwarmDataFieldGetEntries()`, clearing the caller's handle to `NULL`.

Not Collective

Input Parameter:
- `gfield` - the `DMSwarmDataField`

Output Parameter:
- `data` - pointer that will be set to `NULL`

Level: developer

See also: `DMSwarmDataField`, `DMSwarmDataFieldGetEntries()`

# External Links
$(_doc_external("DMSwarm/DMSwarmDataFieldRestoreEntries"))
"""
function DMSwarmDataFieldRestoreEntries(petsclib::PetscLibType, gfield::DMSwarmDataField)
    error("DMSwarmDataFieldRestoreEntries: no generated method for these argument types")
end

@for_petsc function DMSwarmDataFieldRestoreEntries(petsclib::$UnionPetscLib, gfield::DMSwarmDataField )
	data_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:DMSwarmDataFieldRestoreEntries, $petsc_library),
               PetscErrorCode,
               (DMSwarmDataField, Ptr{Ptr{Cvoid}}),
               gfield, data_,
              )

	data = data_[]

	return data
end 

"""
	DMSwarmSortDestroy(petsclib::PetscLibType, ctx::Union{DMSwarmSort, Ref{DMSwarmSort}}) 

# External Links
$(_doc_external("DMSwarm/DMSwarmSortDestroy"))
"""
function DMSwarmSortDestroy(petsclib::PetscLibType, ctx::Union{DMSwarmSort, Ref{DMSwarmSort}})
    error("DMSwarmSortDestroy: no generated method for these argument types")
end

@for_petsc function DMSwarmSortDestroy(petsclib::$UnionPetscLib, ctx::Union{DMSwarmSort, Ref{DMSwarmSort}} )
	ctx_ = ctx isa Base.RefValue ? ctx : Ref{DMSwarmSort}(ctx)

    @chk ccall(
               (:DMSwarmSortDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{DMSwarmSort},),
               ctx_,
              )


	return nothing
end 

"""
	DMSwarmSortGetAccess(petsclib::PetscLibType, sw::AbstractPetscDM) 
Setups up a `DMSWARM` point sort context for efficient traversal of points within a cell

Not Collective

Input Parameter:
- `sw` - a `DMSWARM` object

Level: advanced

See also: `DMSWARM`, `DMSwarmSetType()`, `DMSwarmSortRestoreAccess()`

# External Links
$(_doc_external("DMSwarm/DMSwarmSortGetAccess"))
"""
function DMSwarmSortGetAccess(petsclib::PetscLibType, sw::AbstractPetscDM)
    error("DMSwarmSortGetAccess: no generated method for these argument types")
end

@for_petsc function DMSwarmSortGetAccess(petsclib::$UnionPetscLib, sw::AbstractPetscDM )

    @chk ccall(
               (:DMSwarmSortGetAccess, $petsc_library),
               PetscErrorCode,
               (CDM,),
               sw,
              )


	return nothing
end 

"""
	isvalid::PetscBool = DMSwarmSortGetIsValid(petsclib::PetscLibType, sw::AbstractPetscDM) 
Gets the isvalid flag associated with a `DMSWARM` point sorting context

Not Collective

Input Parameter:
- `sw` - a `DMSWARM` object

Output Parameter:
- `isvalid` - flag indicating whether the sort context is up-to-date

Level: advanced

See also: `DMSWARM`, `DMSwarmSetType()`, `DMSwarmSortGetAccess()`

# External Links
$(_doc_external("DMSwarm/DMSwarmSortGetIsValid"))
"""
function DMSwarmSortGetIsValid(petsclib::PetscLibType, sw::AbstractPetscDM)
    error("DMSwarmSortGetIsValid: no generated method for these argument types")
end

@for_petsc function DMSwarmSortGetIsValid(petsclib::$UnionPetscLib, sw::AbstractPetscDM )
	isvalid_ = Ref{PetscBool}()

    @chk ccall(
               (:DMSwarmSortGetIsValid, $petsc_library),
               PetscErrorCode,
               (CDM, Ptr{PetscBool}),
               sw, isvalid_,
              )

	isvalid = isvalid_[]

	return isvalid
end 

"""
	npoints::PetscInt = DMSwarmSortGetNumberOfPointsPerCell(petsclib::PetscLibType, sw::AbstractPetscDM, cell::PetscInt) 
Returns the number of points in a cell

Not Collective

Input Parameters:
- `sw`   - a `DMSWARM` objects
- `cell` - the cell number in the cell `DM`

Output Parameter:
- `npoints` - the number of points in the cell

Level: advanced

See also: `DMSWARM`, `DMSwarmSetType()`, `DMSwarmSortGetAccess()`, `DMSwarmSortGetPointsPerCell()`

# External Links
$(_doc_external("DMSwarm/DMSwarmSortGetNumberOfPointsPerCell"))
"""
function DMSwarmSortGetNumberOfPointsPerCell(petsclib::PetscLibType, sw::AbstractPetscDM, cell::Integer)
    error("DMSwarmSortGetNumberOfPointsPerCell: no generated method for these argument types")
end

@for_petsc function DMSwarmSortGetNumberOfPointsPerCell(petsclib::$UnionPetscLib, sw::AbstractPetscDM, cell::$PetscInt )
	npoints_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMSwarmSortGetNumberOfPointsPerCell, $petsc_library),
               PetscErrorCode,
               (CDM, $PetscInt, Ptr{$PetscInt}),
               sw, cell, npoints_,
              )

	npoints = npoints_[]

	return npoints
end 

"""
	npoints::PetscInt = DMSwarmSortGetPointsPerCell(petsclib::PetscLibType, sw::AbstractPetscDM, cell::PetscInt, pidlist::PetscInt) 
Creates an array of point indices for all points in a cell

Not Collective

Input Parameters:
- `sw`      - a `DMSWARM` object
- `cell`    - the cell number in the cell `DM`
- `npoints` - the number of points in the cell
- `pidlist` - array of the indices identifying all points in cell e

Level: advanced

See also: `DMSWARM`, `DMSwarmSetType()`, `DMSwarmRestorePointsPerCell()`, `DMSwarmSortGetAccess()`, `DMSwarmSortGetNumberOfPointsPerCell()`

# External Links
$(_doc_external("DMSwarm/DMSwarmSortGetPointsPerCell"))
"""
function DMSwarmSortGetPointsPerCell(petsclib::PetscLibType, sw::AbstractPetscDM, cell::Integer, pidlist::Integer)
    error("DMSwarmSortGetPointsPerCell: no generated method for these argument types")
end

@for_petsc function DMSwarmSortGetPointsPerCell(petsclib::$UnionPetscLib, sw::AbstractPetscDM, cell::$PetscInt, pidlist::$PetscInt )
	npoints_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMSwarmSortGetPointsPerCell, $petsc_library),
               PetscErrorCode,
               (CDM, $PetscInt, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
               sw, cell, npoints_, pidlist,
              )

	npoints = npoints_[]

	return npoints
end 

"""
	ncells::PetscInt,npoints::PetscInt = DMSwarmSortGetSizes(petsclib::PetscLibType, sw::AbstractPetscDM) 
Gets the sizes associated with a `DMSWARM` point sorting context

Not Collective

Input Parameter:
- `sw` - a `DMSWARM` object

Output Parameters:
- `ncells`  - number of cells within the sort context (pass `NULL` to ignore)
- `npoints` - number of points used to create the sort context (pass `NULL` to ignore)

Level: advanced

See also: `DMSWARM`, `DMSwarmSetType()`, `DMSwarmSortGetAccess()`

# External Links
$(_doc_external("DMSwarm/DMSwarmSortGetSizes"))
"""
function DMSwarmSortGetSizes(petsclib::PetscLibType, sw::AbstractPetscDM)
    error("DMSwarmSortGetSizes: no generated method for these argument types")
end

@for_petsc function DMSwarmSortGetSizes(petsclib::$UnionPetscLib, sw::AbstractPetscDM )
	ncells_ = Ref{$PetscInt}()
	npoints_ = Ref{$PetscInt}()

    @chk ccall(
               (:DMSwarmSortGetSizes, $petsc_library),
               PetscErrorCode,
               (CDM, Ptr{$PetscInt}, Ptr{$PetscInt}),
               sw, ncells_, npoints_,
              )

	ncells = ncells_[]
	npoints = npoints_[]

	return ncells,npoints
end 

"""
	DMSwarmSortRestoreAccess(petsclib::PetscLibType, sw::AbstractPetscDM) 
Invalidates the `DMSWARM` point sorting context previously computed with `DMSwarmSortGetAccess()`

Not Collective

Input Parameter:
- `sw` - a `DMSWARM` object

Level: advanced

See also: `DMSWARM`, `DMSwarmSetType()`, `DMSwarmSortGetAccess()`

# External Links
$(_doc_external("DMSwarm/DMSwarmSortRestoreAccess"))
"""
function DMSwarmSortRestoreAccess(petsclib::PetscLibType, sw::AbstractPetscDM)
    error("DMSwarmSortRestoreAccess: no generated method for these argument types")
end

@for_petsc function DMSwarmSortRestoreAccess(petsclib::$UnionPetscLib, sw::AbstractPetscDM )

    @chk ccall(
               (:DMSwarmSortRestoreAccess, $petsc_library),
               PetscErrorCode,
               (CDM,),
               sw,
              )


	return nothing
end 

"""
	DMSwarmSortRestorePointsPerCell(petsclib::PetscLibType, dm::AbstractPetscDM, e::PetscInt, npoints::PetscInt, pidlist::PetscInt) 
Restores an array of point indices for all points in a cell

Not Collective

Input Parameters:
- `dm`      - a `DMSWARM` object
- `e`       - the index of the cell
- `npoints` - the number of points in the cell
- `pidlist` - array of the indices identifying all points in cell e

Level: advanced

See also: `DMSWARM`, `DMSwarmSetType()`, `DMSwarmSortGetPointsPerCell()`, `DMSwarmSortGetAccess()`, `DMSwarmSortGetNumberOfPointsPerCell()`

# External Links
$(_doc_external("DMSwarm/DMSwarmSortRestorePointsPerCell"))
"""
function DMSwarmSortRestorePointsPerCell(petsclib::PetscLibType, dm::AbstractPetscDM, e::Integer, npoints::Integer, pidlist::Integer)
    error("DMSwarmSortRestorePointsPerCell: no generated method for these argument types")
end

@for_petsc function DMSwarmSortRestorePointsPerCell(petsclib::$UnionPetscLib, dm::AbstractPetscDM, e::$PetscInt, npoints::$PetscInt, pidlist::$PetscInt )
	npoints_ = Ref{$PetscInt}(npoints)

    @chk ccall(
               (:DMSwarmSortRestorePointsPerCell, $petsc_library),
               PetscErrorCode,
               (CDM, $PetscInt, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
               dm, e, npoints_, pidlist,
              )


	return nothing
end 

