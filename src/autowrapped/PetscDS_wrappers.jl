"""
	bd::PetscInt = PetscDSAddBoundary(petsclib::PetscLibType,ds::PetscDS, type::DMBoundaryConditionType, name::String, label::DMLabel, Nv::PetscInt, values::Vector{PetscInt}, field::PetscInt, Nc::PetscInt, comps::Vector{PetscInt}, bcFunc::Ptr{Cvoid}, bcFunc_t::Ptr{Cvoid}, ctx::Ptr{Cvoid}) 
Add a boundary condition to the model.

Collective

Input Parameters:
- `ds`       - The `PetscDS` object
- `type`     - The type of condition, e.g. `DM_BC_ESSENTIAL`/`DM_BC_ESSENTIAL_FIELD` (Dirichlet), or `DM_BC_NATURAL` (Neumann)
- `name`     - The name for the boundary condition
- `label`    - The label defining constrained points
- `Nv`       - The number of `DMLabel` values for constrained points
- `values`   - An array of label values for constrained points
- `field`    - The field to constrain
- `Nc`       - The number of constrained field components (0 will constrain all fields)
- `comps`    - An array of constrained component numbers
- `bcFunc`   - A pointwise function giving boundary values
- `bcFunc_t` - A pointwise function giving the time derivative of the boundary values, or `NULL`
- `ctx`      - An optional application context for `bcFunc`

Output Parameter:
- `bd` - The boundary number

Options Database Keys:
- `-bc_NAME values`     - comma separated list of values for the boundary condition NAME
- `-bc_NAME_comp comps` - comma separated list of components for the boundary condition NAME

Level: developer

-seealso: `PetscDS`, `PetscWeakForm`, `DMLabel`, `DMBoundaryConditionType`, `PetscDSAddBoundaryByName()`, `PetscDSGetBoundary()`, `PetscDSSetResidual()`, `PetscDSSetBdResidual()`

# External Links
$(_doc_external("DT/PetscDSAddBoundary"))
"""
function PetscDSAddBoundary(petsclib::PetscLibType, ds::PetscDS, type::DMBoundaryConditionType, name::String, label::DMLabel, Nv::Integer, values::AbstractVector{<:Number}, field::Integer, Nc::Integer, comps::AbstractVector{<:Number}, bcFunc::Ptr{Cvoid}, bcFunc_t::Ptr{Cvoid}, ctx::Ptr{Cvoid})
    error("PetscDSAddBoundary: no generated method for these argument types")
end

@for_petsc function PetscDSAddBoundary(petsclib::$UnionPetscLib, ds::PetscDS, type::DMBoundaryConditionType, name::String, label::DMLabel, Nv::$PetscInt, values::Vector{$PetscInt}, field::$PetscInt, Nc::$PetscInt, comps::Vector{$PetscInt}, bcFunc::Ptr{Cvoid}, bcFunc_t::Ptr{Cvoid}, ctx::Ptr{Cvoid} )
	bd_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDSAddBoundary, $petsc_library),
               PetscErrorCode,
               (PetscDS, DMBoundaryConditionType, Ptr{Cchar}, DMLabel, $PetscInt, Ptr{$PetscInt}, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{$PetscInt}),
               ds, type, name, label, Nv, values, field, Nc, comps, bcFunc, bcFunc_t, ctx, bd_,
              )

	bd = bd_[]

	return bd
end 

"""
	bd::PetscInt = PetscDSAddBoundaryByName(petsclib::PetscLibType,ds::PetscDS, type::DMBoundaryConditionType, name::String, lname::String, Nv::PetscInt, values::Vector{PetscInt}, field::PetscInt, Nc::PetscInt, comps::Vector{PetscInt}, bcFunc::Ptr{Cvoid}, bcFunc_t::Ptr{Cvoid}, ctx::Ptr{Cvoid}) 
Add a boundary condition to the model.

Collective

Input Parameters:
- `ds`       - The `PetscDS` object
- `type`     - The type of condition, e.g. `DM_BC_ESSENTIAL`/`DM_BC_ESSENTIAL_FIELD` (Dirichlet), or `DM_BC_NATURAL` (Neumann)
- `name`     - The boundary condition name
- `lname`    - The name of the label defining constrained points
- `Nv`       - The number of `DMLabel` values for constrained points
- `values`   - An array of label values for constrained points
- `field`    - The field to constrain
- `Nc`       - The number of constrained field components (0 will constrain all fields)
- `comps`    - An array of constrained component numbers
- `bcFunc`   - A pointwise function giving boundary values
- `bcFunc_t` - A pointwise function giving the time derivative of the boundary values, or `NULL`
- `ctx`      - An optional application context for `bcFunc`

Output Parameter:
- `bd` - The boundary number

Options Database Keys:
- `-bc_NAME values`     - comma separated list of values for the boundary condition NAME
- `-bc_NAME_comp comps` - comma separated list of components for the boundary condition NAME

Calling Sequence of `bcFunc` and `bcFunc_t`:
If the type is `DM_BC_ESSENTIAL`
-seealso: `PetscDS`, `PetscWeakForm`, `DMLabel`, `DMBoundaryConditionType`, `PetscDSAddBoundary()`, `PetscDSGetBoundary()`, `PetscDSSetResidual()`, `PetscDSSetBdResidual()`

# External Links
$(_doc_external("DT/PetscDSAddBoundaryByName"))
"""
function PetscDSAddBoundaryByName(petsclib::PetscLibType, ds::PetscDS, type::DMBoundaryConditionType, name::String, lname::String, Nv::Integer, values::AbstractVector{<:Number}, field::Integer, Nc::Integer, comps::AbstractVector{<:Number}, bcFunc::Ptr{Cvoid}, bcFunc_t::Ptr{Cvoid}, ctx::Ptr{Cvoid})
    error("PetscDSAddBoundaryByName: no generated method for these argument types")
end

@for_petsc function PetscDSAddBoundaryByName(petsclib::$UnionPetscLib, ds::PetscDS, type::DMBoundaryConditionType, name::String, lname::String, Nv::$PetscInt, values::Vector{$PetscInt}, field::$PetscInt, Nc::$PetscInt, comps::Vector{$PetscInt}, bcFunc::Ptr{Cvoid}, bcFunc_t::Ptr{Cvoid}, ctx::Ptr{Cvoid} )
	bd_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDSAddBoundaryByName, $petsc_library),
               PetscErrorCode,
               (PetscDS, DMBoundaryConditionType, Ptr{Cchar}, Ptr{Cchar}, $PetscInt, Ptr{$PetscInt}, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{$PetscInt}),
               ds, type, name, lname, Nv, values, field, Nc, comps, bcFunc, bcFunc_t, ctx, bd_,
              )

	bd = bd_[]

	return bd
end 

"""
	PetscDSAddDiscretization(petsclib::PetscLibType,prob::PetscDS, disc) 
Adds a discretization object

Not Collective

Input Parameters:
- `prob` - The `PetscDS` object
- `disc` - The discretization object, this can be a `PetscFE` or `PetscFV`

Level: beginner

-seealso: `PetscWeakForm`, `PetscFE`, `PetscFV`, `PetscDSGetDiscretization()`, `PetscDSSetDiscretization()`, `PetscDSGetNumFields()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSAddDiscretization"))
"""
function PetscDSAddDiscretization(petsclib::PetscLibType, prob::PetscDS, disc)
    error("PetscDSAddDiscretization: no generated method for these argument types")
end

@for_petsc function PetscDSAddDiscretization(petsclib::$UnionPetscLib, prob::PetscDS, disc )

    @chk ccall(
               (:PetscDSAddDiscretization, $petsc_library),
               PetscErrorCode,
               (PetscDS, PetscObject),
               prob, disc,
              )


	return nothing
end 

"""
	PetscDSCopy(petsclib::PetscLibType,ds::PetscDS, minDegree::PetscInt, maxDegree::PetscInt, dmNew::AbstractPetscDM, dsNew::PetscDS) 
Copy the contents of a `PetscDS` into another `PetscDS` on a new `DM`.

Collective

Input Parameters:
- `ds`        - the source `PetscDS`
- `minDegree` - the minimum polynomial degree to consider when selecting discretizations, or `PETSC_DETERMINE`
- `maxDegree` - the maximum polynomial degree to consider when selecting discretizations, or `PETSC_DETERMINE`
- `dmNew`     - the target `DM` used to resolve boundary condition labels for the copied boundaries

Output Parameter:
- `dsNew` - the destination `PetscDS`

Level: developer

-seealso: `PetscDS`, `PetscDSCopyEquations()`, `PetscDSCopyConstants()`, `PetscDSCopyExactSolutions()`, `PetscDSCopyBounds()`, `PetscDSCopyBoundary()`

# External Links
$(_doc_external("DT/PetscDSCopy"))
"""
function PetscDSCopy(petsclib::PetscLibType, ds::PetscDS, minDegree::Integer, maxDegree::Integer, dmNew::AbstractPetscDM, dsNew::PetscDS)
    error("PetscDSCopy: no generated method for these argument types")
end

@for_petsc function PetscDSCopy(petsclib::$UnionPetscLib, ds::PetscDS, minDegree::$PetscInt, maxDegree::$PetscInt, dmNew::AbstractPetscDM, dsNew::PetscDS )

    @chk ccall(
               (:PetscDSCopy, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, $PetscInt, CDM, PetscDS),
               ds, minDegree, maxDegree, dmNew, dsNew,
              )


	return nothing
end 

"""
	PetscDSCopyBoundary(petsclib::PetscLibType,ds::PetscDS, numFields::PetscInt, fields::Vector{PetscInt}, newds::PetscDS) 
Copy all boundary condition objects to the new `PetscDS`

Not Collective

Input Parameters:
- `ds`        - The source `PetscDS` object
- `numFields` - The number of selected fields, or `PETSC_DEFAULT` for all fields
- `fields`    - The selected fields, or `NULL` for all fields

Output Parameter:
- `newds` - The target `PetscDS`, now with a copy of the boundary conditions

Level: intermediate

-seealso: `PetscDS`, `DMBoundary`, `PetscDSCopyEquations()`, `PetscDSSetResidual()`, `PetscDSSetJacobian()`, `PetscDSSetRiemannSolver()`, `PetscDSSetBdResidual()`, `PetscDSSetBdJacobian()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSCopyBoundary"))
"""
function PetscDSCopyBoundary(petsclib::PetscLibType, ds::PetscDS, numFields::Integer, fields::AbstractVector{<:Number}, newds::PetscDS)
    error("PetscDSCopyBoundary: no generated method for these argument types")
end

@for_petsc function PetscDSCopyBoundary(petsclib::$UnionPetscLib, ds::PetscDS, numFields::$PetscInt, fields::Vector{$PetscInt}, newds::PetscDS )

    @chk ccall(
               (:PetscDSCopyBoundary, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{$PetscInt}, PetscDS),
               ds, numFields, fields, newds,
              )


	return nothing
end 

"""
	PetscDSCopyBounds(petsclib::PetscLibType,ds::PetscDS, newds::PetscDS) 
Copy lower and upper solution bounds set with `PetscDSSetLowerBound()` and `PetscDSSetLowerBound()` to another `PetscDS`

Not Collective

Input Parameter:
- `ds` - The `PetscDS` object

Output Parameter:
- `newds` - The `PetscDS` copy

Level: intermediate

-seealso: `PetscDS`, `PetscDSCopyBoundary()`, `PetscDSCopyEquations()`, `PetscDSCopyExactSolutions()`, `PetscDSSetResidual()`, `PetscDSSetJacobian()`, `PetscDSSetRiemannSolver()`, `PetscDSSetBdResidual()`, `PetscDSSetBdJacobian()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSCopyBounds"))
"""
function PetscDSCopyBounds(petsclib::PetscLibType, ds::PetscDS, newds::PetscDS)
    error("PetscDSCopyBounds: no generated method for these argument types")
end

@for_petsc function PetscDSCopyBounds(petsclib::$UnionPetscLib, ds::PetscDS, newds::PetscDS )

    @chk ccall(
               (:PetscDSCopyBounds, $petsc_library),
               PetscErrorCode,
               (PetscDS, PetscDS),
               ds, newds,
              )


	return nothing
end 

"""
	PetscDSCopyConstants(petsclib::PetscLibType,prob::PetscDS, newprob::PetscDS) 
Copy all constants set with `PetscDSSetConstants()` to another `PetscDS`

Not Collective

Input Parameter:
- `prob` - The `PetscDS` object

Output Parameter:
- `newprob` - The `PetscDS` copy

Level: intermediate

-seealso: `PetscDS`, `PetscDSCopyBoundary()`, `PetscDSCopyEquations()`, `PetscDSSetResidual()`, `PetscDSSetJacobian()`, `PetscDSSetRiemannSolver()`, `PetscDSSetBdResidual()`, `PetscDSSetBdJacobian()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSCopyConstants"))
"""
function PetscDSCopyConstants(petsclib::PetscLibType, prob::PetscDS, newprob::PetscDS)
    error("PetscDSCopyConstants: no generated method for these argument types")
end

@for_petsc function PetscDSCopyConstants(petsclib::$UnionPetscLib, prob::PetscDS, newprob::PetscDS )

    @chk ccall(
               (:PetscDSCopyConstants, $petsc_library),
               PetscErrorCode,
               (PetscDS, PetscDS),
               prob, newprob,
              )


	return nothing
end 

"""
	PetscDSCopyEquations(petsclib::PetscLibType,prob::PetscDS, newprob::PetscDS) 
Copy all pointwise function pointers to another `PetscDS`

Not Collective

Input Parameter:
- `prob` - The `PetscDS` object

Output Parameter:
- `newprob` - The `PetscDS` copy

Level: intermediate

-seealso: `PetscDS`, `PetscDSCopyBoundary()`, `PetscDSSetResidual()`, `PetscDSSetJacobian()`, `PetscDSSetRiemannSolver()`, `PetscDSSetBdResidual()`, `PetscDSSetBdJacobian()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSCopyEquations"))
"""
function PetscDSCopyEquations(petsclib::PetscLibType, prob::PetscDS, newprob::PetscDS)
    error("PetscDSCopyEquations: no generated method for these argument types")
end

@for_petsc function PetscDSCopyEquations(petsclib::$UnionPetscLib, prob::PetscDS, newprob::PetscDS )

    @chk ccall(
               (:PetscDSCopyEquations, $petsc_library),
               PetscErrorCode,
               (PetscDS, PetscDS),
               prob, newprob,
              )


	return nothing
end 

"""
	PetscDSCopyExactSolutions(petsclib::PetscLibType,ds::PetscDS, newds::PetscDS) 
Copy all exact solutions set with `PetscDSSetExactSolution()` and `PetscDSSetExactSolutionTimeDerivative()` to another `PetscDS`

Not Collective

Input Parameter:
- `ds` - The `PetscDS` object

Output Parameter:
- `newds` - The `PetscDS` copy

Level: intermediate

-seealso: `PetscDS`, `PetscDSCopyBoundary()`, `PetscDSCopyEquations()`, `PetscDSCopyBounds()`, `PetscDSSetResidual()`, `PetscDSSetJacobian()`, `PetscDSSetRiemannSolver()`, `PetscDSSetBdResidual()`, `PetscDSSetBdJacobian()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSCopyExactSolutions"))
"""
function PetscDSCopyExactSolutions(petsclib::PetscLibType, ds::PetscDS, newds::PetscDS)
    error("PetscDSCopyExactSolutions: no generated method for these argument types")
end

@for_petsc function PetscDSCopyExactSolutions(petsclib::$UnionPetscLib, ds::PetscDS, newds::PetscDS )

    @chk ccall(
               (:PetscDSCopyExactSolutions, $petsc_library),
               PetscErrorCode,
               (PetscDS, PetscDS),
               ds, newds,
              )


	return nothing
end 

"""
	ds::PetscDS = PetscDSCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Creates an empty `PetscDS` object. The type can then be set with `PetscDSSetType()`.

Collective

Input Parameter:
- `comm` - The communicator for the `PetscDS` object

Output Parameter:
- `ds` - The `PetscDS` object

Level: beginner

-seealso: `PetscDS`, `PetscDSSetType()`, `PETSCDSBASIC`, `PetscDSType`, `PetscDSDestroy()`

# External Links
$(_doc_external("DT/PetscDSCreate"))
"""
function PetscDSCreate(petsclib::PetscLibType, comm::MPI_Comm)
    error("PetscDSCreate: no generated method for these argument types")
end

@for_petsc function PetscDSCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	ds_ = Ref{PetscDS}()

    @chk ccall(
               (:PetscDSCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscDS}),
               comm, ds_,
              )

	ds = ds_[]

	return ds
end 

"""
	PetscDSDestroy(petsclib::PetscLibType,ds::Union{PetscDS, Ref{PetscDS}}) 
Destroys a `PetscDS` object

Collective

Input Parameter:
- `ds` - the `PetscDS` object to destroy

Level: developer

-seealso: `PetscDSView()`

# External Links
$(_doc_external("DT/PetscDSDestroy"))
"""
function PetscDSDestroy(petsclib::PetscLibType, ds::Union{PetscDS, Ref{PetscDS}})
    error("PetscDSDestroy: no generated method for these argument types")
end

@for_petsc function PetscDSDestroy(petsclib::$UnionPetscLib, ds::Union{PetscDS, Ref{PetscDS}} )
	ds_ = ds isa Base.RefValue ? ds : Ref{PetscDS}(ds)

    @chk ccall(
               (:PetscDSDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscDS},),
               ds_,
              )


	return nothing
end 

"""
	PetscDSDestroyBoundary(petsclib::PetscLibType,ds::PetscDS) 
Remove all `DMBoundary` objects from the `PetscDS`

Not Collective

Input Parameter:
- `ds` - The `PetscDS` object

Level: intermediate

-seealso: `PetscDS`, `DMBoundary`, `PetscDSCopyBoundary()`, `PetscDSCopyEquations()`

# External Links
$(_doc_external("DT/PetscDSDestroyBoundary"))
"""
function PetscDSDestroyBoundary(petsclib::PetscLibType, ds::PetscDS)
    error("PetscDSDestroyBoundary: no generated method for these argument types")
end

@for_petsc function PetscDSDestroyBoundary(petsclib::$UnionPetscLib, ds::PetscDS )

    @chk ccall(
               (:PetscDSDestroyBoundary, $petsc_library),
               PetscErrorCode,
               (PetscDS,),
               ds,
              )


	return nothing
end 

"""
	PetscDSFinalizePackage(petsclib::PetscLibType) 
This function finalizes everything in the `PetscDS` package. It is called
from `PetscFinalize()`.

Level: developer

-seealso: `PetscInitialize()`

# External Links
$(_doc_external("DM/PetscDSFinalizePackage"))
"""
function PetscDSFinalizePackage(petsclib::PetscLibType)
    error("PetscDSFinalizePackage: no generated method for these argument types")
end

@for_petsc function PetscDSFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscDSFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	g0::Ptr{Cvoid},g1::Ptr{Cvoid},g2::Ptr{Cvoid},g3::Ptr{Cvoid} = PetscDSGetBdJacobian(petsclib::PetscLibType,ds::PetscDS, f::PetscInt, g::PetscInt) 
Get the pointwise boundary Jacobian function for given test and basis field

Not Collective

Input Parameters:
- `ds` - The `PetscDS`
- `f`  - The test field number
- `g`  - The field number

Output Parameters:
- `g0` - integrand for the test and basis function term, see `PetscBdPointJacFn`
- `g1` - integrand for the test function and basis function gradient term, see `PetscBdPointJacFn`
- `g2` - integrand for the test function gradient and basis function term, see `PetscBdPointJacFn`
- `g3` - integrand for the test function gradient and basis function gradient term, see `PetscBdPointJacFn`

Level: intermediate

-seealso: `PetscDS`, `PetscBdPointJacFn`, `PetscDSSetBdJacobian()`

# External Links
$(_doc_external("DT/PetscDSGetBdJacobian"))
"""
function PetscDSGetBdJacobian(petsclib::PetscLibType, ds::PetscDS, f::Integer, g::Integer)
    error("PetscDSGetBdJacobian: no generated method for these argument types")
end

@for_petsc function PetscDSGetBdJacobian(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt, g::$PetscInt )
	g0_ = Ref{Ptr{Cvoid}}()
	g1_ = Ref{Ptr{Cvoid}}()
	g2_ = Ref{Ptr{Cvoid}}()
	g3_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:PetscDSGetBdJacobian, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, $PetscInt, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}),
               ds, f, g, g0_, g1_, g2_, g3_,
              )

	g0 = g0_[]
	g1 = g1_[]
	g2 = g2_[]
	g3 = g3_[]

	return g0,g1,g2,g3
end 

"""
	g0::Ptr{Cvoid},g1::Ptr{Cvoid},g2::Ptr{Cvoid},g3::Ptr{Cvoid} = PetscDSGetBdJacobianPreconditioner(petsclib::PetscLibType,ds::PetscDS, f::PetscInt, g::PetscInt) 
Get the pointwise boundary Jacobian function for given test and basis field that constructs the
matrix used to construct the preconditioner

Not Collective; No Fortran Support

Input Parameters:
- `ds` - The `PetscDS`
- `f`  - The test field number
- `g`  - The field number

Output Parameters:
- `g0` - integrand for the test and basis function term, see `PetscBdPointJacFn`
- `g1` - integrand for the test function and basis function gradient term, see `PetscBdPointJacFn`
- `g2` - integrand for the test function gradient and basis function term, see `PetscBdPointJacFn`
- `g3` - integrand for the test function gradient and basis function gradient term, see `PetscBdPointJacFn`

Level: intermediate

-seealso: `PetscDS`, `PetscBdPointJacFn`, `PetscDSSetBdJacobianPreconditioner()`

# External Links
$(_doc_external("DT/PetscDSGetBdJacobianPreconditioner"))
"""
function PetscDSGetBdJacobianPreconditioner(petsclib::PetscLibType, ds::PetscDS, f::Integer, g::Integer)
    error("PetscDSGetBdJacobianPreconditioner: no generated method for these argument types")
end

@for_petsc function PetscDSGetBdJacobianPreconditioner(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt, g::$PetscInt )
	g0_ = Ref{Ptr{Cvoid}}()
	g1_ = Ref{Ptr{Cvoid}}()
	g2_ = Ref{Ptr{Cvoid}}()
	g3_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:PetscDSGetBdJacobianPreconditioner, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, $PetscInt, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}),
               ds, f, g, g0_, g1_, g2_, g3_,
              )

	g0 = g0_[]
	g1 = g1_[]
	g2 = g2_[]
	g3 = g3_[]

	return g0,g1,g2,g3
end 

"""
	f0::Ptr{Cvoid},f1::Ptr{Cvoid} = PetscDSGetBdResidual(petsclib::PetscLibType,ds::PetscDS, f::PetscInt) 
Get the pointwise boundary residual function for a given test field

Not Collective

Input Parameters:
- `ds` - The PetscDS
- `f`  - The test field number

Output Parameters:
- `f0` - boundary integrand for the test function term, see `PetscBdPointFn`
- `f1` - boundary integrand for the test function gradient term, see `PetscBdPointFn`

Level: intermediate

-seealso: `PetscDS`, `PetscBdPointFn`, `PetscDSSetBdResidual()`

# External Links
$(_doc_external("DT/PetscDSGetBdResidual"))
"""
function PetscDSGetBdResidual(petsclib::PetscLibType, ds::PetscDS, f::Integer)
    error("PetscDSGetBdResidual: no generated method for these argument types")
end

@for_petsc function PetscDSGetBdResidual(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt )
	f0_ = Ref{Ptr{Cvoid}}()
	f1_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:PetscDSGetBdResidual, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}),
               ds, f, f0_, f1_,
              )

	f0 = f0_[]
	f1 = f1_[]

	return f0,f1
end 

# override for PetscDSGetBoundary; C signature: PetscDSGetBoundary(PetscDS ds, PetscInt bd, PetscWeakForm* wf, DMBoundaryConditionType* type, char* name[], DMLabel* label, PetscInt* Nv, PetscInt* values[], PetscInt* field, PetscInt* Nc, PetscInt* comps[], PetscVoidFn** func, PetscVoidFn** func_t, void** ctx)
"""
	Nv::PetscInt,values::Vector{PetscInt},field::PetscInt,Nc::PetscInt,comps::Vector{PetscInt} = PetscDSGetBoundary(petsclib::PetscLibType,ds::PetscDS, bd::PetscInt, wf::PetscWeakForm, type::DMBoundaryConditionType, name::String, label::DMLabel, func::PetscVoidFn, func_t::PetscVoidFn, ctx::Cvoid) 
Gets a boundary condition from the model

Input Parameters:
- `ds` - The `PetscDS` object
- `bd` - The boundary condition number

Output Parameters:
- `wf`     - The `PetscWeakForm` holding the pointwise functions
- `type`   - The type of condition, e.g. `DM_BC_ESSENTIAL`/`DM_BC_ESSENTIAL_FIELD` (Dirichlet), or `DM_BC_NATURAL` (Neumann)
- `name`   - The boundary condition name
- `label`  - The label defining constrained points
- `Nv`     - The number of `DMLabel` ids for constrained points
- `values` - An array of ids for constrained points
- `field`  - The field to constrain
- `Nc`     - The number of constrained field components
- `comps`  - An array of constrained component numbers
- `func`   - A pointwise function giving boundary values
- `func_t` - A pointwise function giving the time derivative of the boundary values
- `ctx`    - An optional user context for `bcFunc`

Options Database Keys:
- `-bc_<boundary name> <num>`      - Overrides the boundary ids
- `-bc_<boundary name>_comp <num>` - Overrides the boundary components

Level: developer

-seealso: `PetscDS`, `PetscWeakForm`, `DMBoundaryConditionType`, `PetscDSAddBoundary()`, `DMLabel`

# External Links
$(_doc_external("Dm/PetscDSGetBoundary"))
"""
function PetscDSGetBoundary(petsclib::PetscLibType, ds::PetscDS, bd::PetscInt) end

@for_petsc function PetscDSGetBoundary(petsclib::$UnionPetscLib, ds::PetscDS, bd::$PetscInt)
	wf_ref    = Ref{PetscWeakForm}(C_NULL)
	type_ref  = Ref{DMBoundaryConditionType}()
	name_ref  = Ref{Ptr{Cchar}}(C_NULL)
	label_ref = Ref{DMLabel}(C_NULL)
	Nv_ref    = Ref{$PetscInt}(0)
	values_ref = Ref{Ptr{$PetscInt}}(C_NULL)
	field_ref = Ref{$PetscInt}(0)
	Nc_ref    = Ref{$PetscInt}(0)
	comps_ref = Ref{Ptr{$PetscInt}}(C_NULL)
	func_ref  = Ref{Ptr{Cvoid}}(C_NULL)
	func_t_ref = Ref{Ptr{Cvoid}}(C_NULL)
	ctx_ref   = Ref{Ptr{Cvoid}}(C_NULL)

    @chk ccall(
               (:PetscDSGetBoundary, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt,
                Ptr{PetscWeakForm}, Ptr{DMBoundaryConditionType}, Ptr{Ptr{Cchar}}, Ptr{DMLabel},
                Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}},
                Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}),
               ds, bd,
               wf_ref, type_ref, name_ref, label_ref,
               Nv_ref, values_ref, field_ref, Nc_ref, comps_ref,
               func_ref, func_t_ref, ctx_ref,
              )

	Nv = Nv_ref[]
	values = Nv > 0 ? unsafe_wrap(Array, values_ref[], Nv; own = false) : $PetscInt[]
	Nc = Nc_ref[]
	comps = Nc > 0 ? unsafe_wrap(Array, comps_ref[], Nc; own = false) : $PetscInt[]

	return (wf = wf_ref[], type = type_ref[], name = name_ref[] == C_NULL ? "" : unsafe_string(name_ref[]),
	        label = label_ref[], Nv = Nv, values = values, field = field_ref[],
	        Nc = Nc, comps = comps, func = func_ref[], func_t = func_t_ref[], ctx = ctx_ref[])
end

"""
	isCohesive::PetscBool = PetscDSGetCohesive(petsclib::PetscLibType,ds::PetscDS, f::PetscInt) 
Returns the flag indicating that a field is cohesive, meaning it is defined on the interior of a cohesive cell

Not Collective

Input Parameters:
- `ds` - The `PetscDS` object
- `f`  - The field index

Output Parameter:
- `isCohesive` - The flag

Level: developer

-seealso: `PetscDS`, `PetscDSSetCohesive()`, `PetscDSIsCohesive()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSGetCohesive"))
"""
function PetscDSGetCohesive(petsclib::PetscLibType, ds::PetscDS, f::Integer)
    error("PetscDSGetCohesive: no generated method for these argument types")
end

@for_petsc function PetscDSGetCohesive(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt )
	isCohesive_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscDSGetCohesive, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{PetscBool}),
               ds, f, isCohesive_,
              )

	isCohesive = isCohesive_[]

	return isCohesive
end 

"""
	offsets::Ptr{PetscInt} = PetscDSGetComponentDerivativeOffsets(petsclib::PetscLibType,prob::PetscDS) 
Returns the offset of each field derivative on an evaluation point

Not Collective

Input Parameter:
- `prob` - The `PetscDS` object

Output Parameter:
- `offsets` - The offsets

Level: beginner

-seealso: `PetscDS`, `PetscDSGetNumFields()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSGetComponentDerivativeOffsets"))
"""
function PetscDSGetComponentDerivativeOffsets(petsclib::PetscLibType, prob::PetscDS)
    error("PetscDSGetComponentDerivativeOffsets: no generated method for these argument types")
end

@for_petsc function PetscDSGetComponentDerivativeOffsets(petsclib::$UnionPetscLib, prob::PetscDS )
	offsets_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:PetscDSGetComponentDerivativeOffsets, $petsc_library),
               PetscErrorCode,
               (PetscDS, Ptr{Ptr{$PetscInt}}),
               prob, offsets_,
              )

	offsets = offsets_[]

	return offsets
end 

"""
	offsets::Ptr{PetscInt} = PetscDSGetComponentDerivativeOffsetsCohesive(petsclib::PetscLibType,ds::PetscDS, s::PetscInt) 
Returns the offset of each field derivative on an evaluation point

Not Collective

Input Parameters:
- `ds` - The `PetscDS` object
- `s`  - The cohesive side, 0 for negative, 1 for positive, 2 for cohesive

Output Parameter:
- `offsets` - The offsets

Level: beginner

-seealso: `PetscDS`, `PetscDSGetNumFields()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSGetComponentDerivativeOffsetsCohesive"))
"""
function PetscDSGetComponentDerivativeOffsetsCohesive(petsclib::PetscLibType, ds::PetscDS, s::Integer)
    error("PetscDSGetComponentDerivativeOffsetsCohesive: no generated method for these argument types")
end

@for_petsc function PetscDSGetComponentDerivativeOffsetsCohesive(petsclib::$UnionPetscLib, ds::PetscDS, s::$PetscInt )
	offsets_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:PetscDSGetComponentDerivativeOffsetsCohesive, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{Ptr{$PetscInt}}),
               ds, s, offsets_,
              )

	offsets = offsets_[]

	return offsets
end 

"""
	off::PetscInt = PetscDSGetComponentOffset(petsclib::PetscLibType,prob::PetscDS, f::PetscInt) 
Returns the offset of the given field on an evaluation point

Not Collective

Input Parameters:
- `prob` - The `PetscDS` object
- `f`    - The field number

Output Parameter:
- `off` - The offset

Level: beginner

-seealso: `PetscDS`, `PetscDSGetNumFields()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSGetComponentOffset"))
"""
function PetscDSGetComponentOffset(petsclib::PetscLibType, prob::PetscDS, f::Integer)
    error("PetscDSGetComponentOffset: no generated method for these argument types")
end

@for_petsc function PetscDSGetComponentOffset(petsclib::$UnionPetscLib, prob::PetscDS, f::$PetscInt )
	off_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDSGetComponentOffset, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{$PetscInt}),
               prob, f, off_,
              )

	off = off_[]

	return off
end 

"""
	offsets::Ptr{PetscInt} = PetscDSGetComponentOffsets(petsclib::PetscLibType,prob::PetscDS) 
Returns the offset of each field on an evaluation point

Not Collective

Input Parameter:
- `prob` - The `PetscDS` object

Output Parameter:
- `offsets` - The offsets

Level: beginner

-seealso: `PetscDS`, `PetscDSGetNumFields()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSGetComponentOffsets"))
"""
function PetscDSGetComponentOffsets(petsclib::PetscLibType, prob::PetscDS)
    error("PetscDSGetComponentOffsets: no generated method for these argument types")
end

@for_petsc function PetscDSGetComponentOffsets(petsclib::$UnionPetscLib, prob::PetscDS )
	offsets_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:PetscDSGetComponentOffsets, $petsc_library),
               PetscErrorCode,
               (PetscDS, Ptr{Ptr{$PetscInt}}),
               prob, offsets_,
              )

	offsets = offsets_[]

	return offsets
end 

"""
	offsets::Ptr{PetscInt} = PetscDSGetComponentOffsetsCohesive(petsclib::PetscLibType,ds::PetscDS, s::PetscInt) 
Returns the offset of each field on an evaluation point

Not Collective

Input Parameters:
- `ds` - The `PetscDS` object
- `s`  - The cohesive side, 0 for negative, 1 for positive, 2 for cohesive

Output Parameter:
- `offsets` - The offsets

Level: beginner

-seealso: `PetscDS`, `PetscDSGetNumFields()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSGetComponentOffsetsCohesive"))
"""
function PetscDSGetComponentOffsetsCohesive(petsclib::PetscLibType, ds::PetscDS, s::Integer)
    error("PetscDSGetComponentOffsetsCohesive: no generated method for these argument types")
end

@for_petsc function PetscDSGetComponentOffsetsCohesive(petsclib::$UnionPetscLib, ds::PetscDS, s::$PetscInt )
	offsets_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:PetscDSGetComponentOffsetsCohesive, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{Ptr{$PetscInt}}),
               ds, s, offsets_,
              )

	offsets = offsets_[]

	return offsets
end 

"""
	components::Ptr{PetscInt} = PetscDSGetComponents(petsclib::PetscLibType,prob::PetscDS) 
Returns the number of components for each field on an evaluation point

Not Collective

Input Parameter:
- `prob` - The `PetscDS` object

Output Parameter:
- `components` - The number of components

Level: beginner

-seealso: `PetscDS`, `PetscDSGetComponentOffsets()`, `PetscDSGetNumFields()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSGetComponents"))
"""
function PetscDSGetComponents(petsclib::PetscLibType, prob::PetscDS)
    error("PetscDSGetComponents: no generated method for these argument types")
end

@for_petsc function PetscDSGetComponents(petsclib::$UnionPetscLib, prob::PetscDS )
	components_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:PetscDSGetComponents, $petsc_library),
               PetscErrorCode,
               (PetscDS, Ptr{Ptr{$PetscInt}}),
               prob, components_,
              )

	components = components_[]

	return components
end 

"""
	numConstants::PetscInt,constants::Ptr{PetscScalar} = PetscDSGetConstants(petsclib::PetscLibType,ds::PetscDS) 
Returns the array of constants passed to point functions from a `PetscDS` object

Not Collective

Input Parameter:
- `ds` - The `PetscDS` object

Output Parameters:
- `numConstants` - The number of constants, or pass in `NULL` if not required
- `constants`    - The array of constants, `NULL` if there are none

Level: intermediate

-seealso: `PetscDS`, `PetscDSSetConstants()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSGetConstants"))
"""
function PetscDSGetConstants(petsclib::PetscLibType, ds::PetscDS)
    error("PetscDSGetConstants: no generated method for these argument types")
end

@for_petsc function PetscDSGetConstants(petsclib::$UnionPetscLib, ds::PetscDS )
	numConstants_ = Ref{$PetscInt}()
	constants_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:PetscDSGetConstants, $petsc_library),
               PetscErrorCode,
               (PetscDS, Ptr{$PetscInt}, Ptr{Ptr{$PetscScalar}}),
               ds, numConstants_, constants_,
              )

	numConstants = numConstants_[]
	constants = constants_[]

	return numConstants,constants
end 

"""
	PetscDSGetContext(petsclib::PetscLibType,ds::PetscDS, f::PetscInt, ctx::Ptr{Cvoid}) 
Returns the context that was passed by `PetscDSSetContext()`

Not Collective

Input Parameters:
- `ds`  - The `PetscDS`
- `f`   - The field number
- `ctx` - the context

Level: intermediate

-seealso: `PetscDS`, `PetscPointFn`, `PetscDSSetContext()`

# External Links
$(_doc_external("DT/PetscDSGetContext"))
"""
function PetscDSGetContext(petsclib::PetscLibType, ds::PetscDS, f::Integer, ctx::Ptr{Cvoid})
    error("PetscDSGetContext: no generated method for these argument types")
end

@for_petsc function PetscDSGetContext(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:PetscDSGetContext, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{Cvoid}),
               ds, f, ctx,
              )


	return nothing
end 

"""
	dimEmbed::PetscInt = PetscDSGetCoordinateDimension(petsclib::PetscLibType,prob::PetscDS) 
Returns the coordinate dimension of the `PetscDS`, meaning the dimension of the space into which the discretiaztions are embedded

Not Collective

Input Parameter:
- `prob` - The `PetscDS` object

Output Parameter:
- `dimEmbed` - The coordinate dimension

Level: beginner

-seealso: `PetscDS`, `PetscDSSetCoordinateDimension()`, `PetscDSGetSpatialDimension()`, `PetscDSGetNumFields()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSGetCoordinateDimension"))
"""
function PetscDSGetCoordinateDimension(petsclib::PetscLibType, prob::PetscDS)
    error("PetscDSGetCoordinateDimension: no generated method for these argument types")
end

@for_petsc function PetscDSGetCoordinateDimension(petsclib::$UnionPetscLib, prob::PetscDS )
	dimEmbed_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDSGetCoordinateDimension, $petsc_library),
               PetscErrorCode,
               (PetscDS, Ptr{$PetscInt}),
               prob, dimEmbed_,
              )

	dimEmbed = dimEmbed_[]

	return dimEmbed
end 

"""
	dimensions::Ptr{PetscInt} = PetscDSGetDimensions(petsclib::PetscLibType,prob::PetscDS) 
Returns the size of the approximation space for each field on an evaluation point

Not Collective

Input Parameter:
- `prob` - The `PetscDS` object

Output Parameter:
- `dimensions` - The number of dimensions

Level: beginner

-seealso: `PetscDS`, `PetscDSGetComponentOffsets()`, `PetscDSGetNumFields()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSGetDimensions"))
"""
function PetscDSGetDimensions(petsclib::PetscLibType, prob::PetscDS)
    error("PetscDSGetDimensions: no generated method for these argument types")
end

@for_petsc function PetscDSGetDimensions(petsclib::$UnionPetscLib, prob::PetscDS )
	dimensions_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:PetscDSGetDimensions, $petsc_library),
               PetscErrorCode,
               (PetscDS, Ptr{Ptr{$PetscInt}}),
               prob, dimensions_,
              )

	dimensions = dimensions_[]

	return dimensions
end 

"""
	disc::PetscObject = PetscDSGetDiscretization(petsclib::PetscLibType,prob::PetscDS, f::PetscInt) 
Returns the discretization object for the given field

Not Collective

Input Parameters:
- `prob` - The `PetscDS` object
- `f`    - The field number

Output Parameter:
- `disc` - The discretization object, this can be a `PetscFE` or a `PetscFV`

Level: beginner

-seealso: `PetscDS`, `PetscFE`, `PetscFV`, `PetscDSSetDiscretization()`, `PetscDSAddDiscretization()`, `PetscDSGetNumFields()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSGetDiscretization"))
"""
function PetscDSGetDiscretization(petsclib::PetscLibType, prob::PetscDS, f::Integer)
    error("PetscDSGetDiscretization: no generated method for these argument types")
end

@for_petsc function PetscDSGetDiscretization(petsclib::$UnionPetscLib, prob::PetscDS, f::$PetscInt )
	disc_ = Ref{PetscObject}()

    @chk ccall(
               (:PetscDSGetDiscretization, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{PetscObject}),
               prob, f, disc_,
              )

	disc = disc_[]

	return disc
end 

"""
	g0::Ptr{Cvoid},g1::Ptr{Cvoid},g2::Ptr{Cvoid},g3::Ptr{Cvoid} = PetscDSGetDynamicJacobian(petsclib::PetscLibType,ds::PetscDS, f::PetscInt, g::PetscInt) 
Get the pointwise dynamic Jacobian, dF/du_t, function for given test and basis field

Not Collective

Input Parameters:
- `ds` - The `PetscDS`
- `f`  - The test field number
- `g`  - The field number

Output Parameters:
- `g0` - integrand for the test and basis function term, see `PetscPointJacFn`
- `g1` - integrand for the test function and basis function gradient term, see `PetscPointJacFn`
- `g2` - integrand for the test function gradient and basis function term, see `PetscPointJacFn`
- `g3` - integrand for the test function gradient and basis function gradient term, see `PetscPointJacFn`

Level: intermediate

-seealso: `PetscDS`, `PetscDSSetJacobian()`, `PetscDSSetDynamicJacobian()`, `PetscPointJacFn`

# External Links
$(_doc_external("DT/PetscDSGetDynamicJacobian"))
"""
function PetscDSGetDynamicJacobian(petsclib::PetscLibType, ds::PetscDS, f::Integer, g::Integer)
    error("PetscDSGetDynamicJacobian: no generated method for these argument types")
end

@for_petsc function PetscDSGetDynamicJacobian(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt, g::$PetscInt )
	g0_ = Ref{Ptr{Cvoid}}()
	g1_ = Ref{Ptr{Cvoid}}()
	g2_ = Ref{Ptr{Cvoid}}()
	g3_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:PetscDSGetDynamicJacobian, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, $PetscInt, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}),
               ds, f, g, g0_, g1_, g2_, g3_,
              )

	g0 = g0_[]
	g1 = g1_[]
	g2 = g2_[]
	g3 = g3_[]

	return g0,g1,g2,g3
end 

"""
	u::Ptr{PetscScalar},u_t::Ptr{PetscScalar},u_x::Ptr{PetscScalar} = PetscDSGetEvaluationArrays(petsclib::PetscLibType,prob::PetscDS) 
Get scratch arrays used to evaluate fields, time derivatives, and field gradients at quadrature points.

Not Collective

Input Parameter:
- `prob` - the `PetscDS`

Output Parameters:
- `u`   - array for the field values, or `NULL` if not needed
- `u_t` - array for the field time derivatives, or `NULL` if not needed
- `u_x` - array for the field gradients, or `NULL` if not needed

Level: developer

-seealso: `PetscDS`, `PetscDSGetWeakFormArrays()`, `PetscDSGetWorkspace()`

# External Links
$(_doc_external("DT/PetscDSGetEvaluationArrays"))
"""
function PetscDSGetEvaluationArrays(petsclib::PetscLibType, prob::PetscDS)
    error("PetscDSGetEvaluationArrays: no generated method for these argument types")
end

@for_petsc function PetscDSGetEvaluationArrays(petsclib::$UnionPetscLib, prob::PetscDS )
	u_ = Ref{Ptr{$PetscScalar}}()
	u_t_ = Ref{Ptr{$PetscScalar}}()
	u_x_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:PetscDSGetEvaluationArrays, $petsc_library),
               PetscErrorCode,
               (PetscDS, Ptr{Ptr{$PetscScalar}}, Ptr{Ptr{$PetscScalar}}, Ptr{Ptr{$PetscScalar}}),
               prob, u_, u_t_, u_x_,
              )

	u = u_[]
	u_t = u_t_[]
	u_x = u_x_[]

	return u,u_t,u_x
end 

"""
	sol::Ptr{Cvoid},ctx::Ptr{Cvoid} = PetscDSGetExactSolution(petsclib::PetscLibType,prob::PetscDS, f::PetscInt) 
Get the pointwise exact solution function for a given test field

Not Collective

Input Parameters:
- `prob` - The `PetscDS`
- `f`    - The test field number

Output Parameters:
- `sol` - exact solution function for the test field, see `PetscPointExactSolutionFn`
- `ctx` - exact solution context

Level: intermediate

-seealso: `PetscDS`, `PetscPointExactSolutionFn`, `PetscDSSetExactSolution()`, `PetscDSGetExactSolutionTimeDerivative()`

# External Links
$(_doc_external("DT/PetscDSGetExactSolution"))
"""
function PetscDSGetExactSolution(petsclib::PetscLibType, prob::PetscDS, f::Integer)
    error("PetscDSGetExactSolution: no generated method for these argument types")
end

@for_petsc function PetscDSGetExactSolution(petsclib::$UnionPetscLib, prob::PetscDS, f::$PetscInt )
	sol_ = Ref{Ptr{Cvoid}}()
	ctx_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:PetscDSGetExactSolution, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}),
               prob, f, sol_, ctx_,
              )

	sol = sol_[]
	ctx = ctx_[]

	return sol,ctx
end 

"""
	sol::Ptr{Cvoid},ctx::Ptr{Cvoid} = PetscDSGetExactSolutionTimeDerivative(petsclib::PetscLibType,prob::PetscDS, f::PetscInt) 
Get the pointwise time derivative of the exact solution function for a given test field

Not Collective

Input Parameters:
- `prob` - The `PetscDS`
- `f`    - The test field number

Output Parameters:
- `sol` - time derivative of the exact solution for the test field, see `PetscPointExactSolutionFn`
- `ctx` - the exact solution context

Level: intermediate

-seealso: `PetscDS`, `PetscPointExactSolutionFn`, `PetscDSSetExactSolutionTimeDerivative()`, `PetscDSGetExactSolution()`

# External Links
$(_doc_external("DT/PetscDSGetExactSolutionTimeDerivative"))
"""
function PetscDSGetExactSolutionTimeDerivative(petsclib::PetscLibType, prob::PetscDS, f::Integer)
    error("PetscDSGetExactSolutionTimeDerivative: no generated method for these argument types")
end

@for_petsc function PetscDSGetExactSolutionTimeDerivative(petsclib::$UnionPetscLib, prob::PetscDS, f::$PetscInt )
	sol_ = Ref{Ptr{Cvoid}}()
	ctx_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:PetscDSGetExactSolutionTimeDerivative, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}),
               prob, f, sol_, ctx_,
              )

	sol = sol_[]
	ctx = ctx_[]

	return sol,ctx
end 

"""
	Tf::Ptr{PetscTabulation} = PetscDSGetFaceTabulation(petsclib::PetscLibType,prob::PetscDS) 
Return the basis tabulation at quadrature points on the faces

Not Collective

Input Parameter:
- `prob` - The `PetscDS` object

Output Parameter:
- `Tf` - The basis function and derivative tabulation on each local face at quadrature points for each field

Level: intermediate

-seealso: `PetscTabulation`, `PetscDS`, `PetscDSGetTabulation()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSGetFaceTabulation"))
"""
function PetscDSGetFaceTabulation(petsclib::PetscLibType, prob::PetscDS)
    error("PetscDSGetFaceTabulation: no generated method for these argument types")
end

@for_petsc function PetscDSGetFaceTabulation(petsclib::$UnionPetscLib, prob::PetscDS )
	Tf_ = Ref{Ptr{PetscTabulation}}()

    @chk ccall(
               (:PetscDSGetFaceTabulation, $petsc_library),
               PetscErrorCode,
               (PetscDS, Ptr{Ptr{PetscTabulation}}),
               prob, Tf_,
              )

	Tf = Tf_[]

	return Tf
end 

"""
	f::PetscInt = PetscDSGetFieldIndex(petsclib::PetscLibType,prob::PetscDS, disc) 
Returns the index of the given field

Not Collective

Input Parameters:
- `prob` - The `PetscDS` object
- `disc` - The discretization object

Output Parameter:
- `f` - The field number

Level: beginner

-seealso: `PetscDS`, `PetscGetDiscretization()`, `PetscDSGetNumFields()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSGetFieldIndex"))
"""
function PetscDSGetFieldIndex(petsclib::PetscLibType, prob::PetscDS, disc)
    error("PetscDSGetFieldIndex: no generated method for these argument types")
end

@for_petsc function PetscDSGetFieldIndex(petsclib::$UnionPetscLib, prob::PetscDS, disc )
	f_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDSGetFieldIndex, $petsc_library),
               PetscErrorCode,
               (PetscDS, PetscObject, Ptr{$PetscInt}),
               prob, disc, f_,
              )

	f = f_[]

	return f
end 

"""
	off::PetscInt = PetscDSGetFieldOffset(petsclib::PetscLibType,prob::PetscDS, f::PetscInt) 
Returns the offset of the given field in the full space basis

Not Collective

Input Parameters:
- `prob` - The `PetscDS` object
- `f`    - The field number

Output Parameter:
- `off` - The offset

Level: beginner

-seealso: `PetscDS`, `PetscDSGetFieldSize()`, `PetscDSGetNumFields()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSGetFieldOffset"))
"""
function PetscDSGetFieldOffset(petsclib::PetscLibType, prob::PetscDS, f::Integer)
    error("PetscDSGetFieldOffset: no generated method for these argument types")
end

@for_petsc function PetscDSGetFieldOffset(petsclib::$UnionPetscLib, prob::PetscDS, f::$PetscInt )
	off_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDSGetFieldOffset, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{$PetscInt}),
               prob, f, off_,
              )

	off = off_[]

	return off
end 

"""
	off::PetscInt = PetscDSGetFieldOffsetCohesive(petsclib::PetscLibType,ds::PetscDS, f::PetscInt) 
Returns the offset of the given field in the full space basis on a cohesive cell

Not Collective

Input Parameters:
- `ds` - The `PetscDS` object
- `f`  - The field number

Output Parameter:
- `off` - The offset

Level: beginner

-seealso: `PetscDS`, `PetscDSGetFieldSize()`, `PetscDSGetNumFields()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSGetFieldOffsetCohesive"))
"""
function PetscDSGetFieldOffsetCohesive(petsclib::PetscLibType, ds::PetscDS, f::Integer)
    error("PetscDSGetFieldOffsetCohesive: no generated method for these argument types")
end

@for_petsc function PetscDSGetFieldOffsetCohesive(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt )
	off_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDSGetFieldOffsetCohesive, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{$PetscInt}),
               ds, f, off_,
              )

	off = off_[]

	return off
end 

"""
	size::PetscInt = PetscDSGetFieldSize(petsclib::PetscLibType,prob::PetscDS, f::PetscInt) 
Returns the size of the given field in the full space basis

Not Collective

Input Parameters:
- `prob` - The `PetscDS` object
- `f`    - The field number

Output Parameter:
- `size` - The size

Level: beginner

-seealso: `PetscDS`, `PetscDSGetFieldOffset()`, `PetscDSGetNumFields()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSGetFieldSize"))
"""
function PetscDSGetFieldSize(petsclib::PetscLibType, prob::PetscDS, f::Integer)
    error("PetscDSGetFieldSize: no generated method for these argument types")
end

@for_petsc function PetscDSGetFieldSize(petsclib::$UnionPetscLib, prob::PetscDS, f::$PetscInt )
	size_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDSGetFieldSize, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{$PetscInt}),
               prob, f, size_,
              )

	size = size_[]

	return size
end 

"""
	forceQuad::PetscBool = PetscDSGetForceQuad(petsclib::PetscLibType,ds::PetscDS) 
Returns the flag to force matching quadratures among the field discretizations

Not collective

Input Parameter:
- `ds` - The `PetscDS` object

Output Parameter:
- `forceQuad` - The flag

Level: intermediate

-seealso: `PetscDS`, `PetscDSSetForceQuad()`, `PetscDSGetDiscretization()`, `PetscDSGetNumFields()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSGetForceQuad"))
"""
function PetscDSGetForceQuad(petsclib::PetscLibType, ds::PetscDS)
    error("PetscDSGetForceQuad: no generated method for these argument types")
end

@for_petsc function PetscDSGetForceQuad(petsclib::$UnionPetscLib, ds::PetscDS )
	forceQuad_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscDSGetForceQuad, $petsc_library),
               PetscErrorCode,
               (PetscDS, Ptr{PetscBool}),
               ds, forceQuad_,
              )

	forceQuad = forceQuad_[]

	return forceQuad
end 

"""
	subprob::PetscDS = PetscDSGetHeightSubspace(petsclib::PetscLibType,prob::PetscDS, height::PetscInt) 
Get the `PetscDS` for the trace subspace at a given height in the mesh.

Not Collective

Input Parameters:
- `prob`   - the `PetscDS`
- `height` - the height (0 for the ambient cell, 1 for faces, etc.)

Output Parameter:
- `subprob` - the `PetscDS` for the trace subspace; `prob` itself is returned when `height` is 0

Level: developer

-seealso: `PetscDS`, `PetscFE`, `PetscFEGetHeightSubspace()`, `PetscDSGetSpatialDimension()`

# External Links
$(_doc_external("DT/PetscDSGetHeightSubspace"))
"""
function PetscDSGetHeightSubspace(petsclib::PetscLibType, prob::PetscDS, height::Integer)
    error("PetscDSGetHeightSubspace: no generated method for these argument types")
end

@for_petsc function PetscDSGetHeightSubspace(petsclib::$UnionPetscLib, prob::PetscDS, height::$PetscInt )
	subprob_ = Ref{PetscDS}()

    @chk ccall(
               (:PetscDSGetHeightSubspace, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{PetscDS}),
               prob, height, subprob_,
              )

	subprob = subprob_[]

	return subprob
end 

"""
	implicit::PetscBool = PetscDSGetImplicit(petsclib::PetscLibType,prob::PetscDS, f::PetscInt) 
Returns the flag for implicit solve for this field. This is just a guide for `TSARKIMEX`

Not Collective

Input Parameters:
- `prob` - The `PetscDS` object
- `f`    - The field number

Output Parameter:
- `implicit` - The flag indicating what kind of solve to use for this field

Level: developer

-seealso: `TSARKIMEX`, `PetscDS`, `PetscDSSetImplicit()`, `PetscDSSetDiscretization()`, `PetscDSAddDiscretization()`, `PetscDSGetNumFields()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSGetImplicit"))
"""
function PetscDSGetImplicit(petsclib::PetscLibType, prob::PetscDS, f::Integer)
    error("PetscDSGetImplicit: no generated method for these argument types")
end

@for_petsc function PetscDSGetImplicit(petsclib::$UnionPetscLib, prob::PetscDS, f::$PetscInt )
	implicit_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscDSGetImplicit, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{PetscBool}),
               prob, f, implicit_,
              )

	implicit = implicit_[]

	return implicit
end 

"""
	g0::Ptr{Cvoid},g1::Ptr{Cvoid},g2::Ptr{Cvoid},g3::Ptr{Cvoid} = PetscDSGetJacobian(petsclib::PetscLibType,ds::PetscDS, f::PetscInt, g::PetscInt) 
Get the pointwise Jacobian function for given test and basis field

Not Collective

Input Parameters:
- `ds` - The `PetscDS`
- `f`  - The test field number
- `g`  - The field number

Output Parameters:
- `g0` - integrand for the test and basis function term, see `PetscPointJacFn`
- `g1` - integrand for the test function and basis function gradient term, see `PetscPointJacFn`
- `g2` - integrand for the test function gradient and basis function term, see `PetscPointJacFn`
- `g3` - integrand for the test function gradient and basis function gradient term, see `PetscPointJacFn`

Level: intermediate

-seealso: `PetscDS`, `PetscDSSetJacobian()`, `PetscPointJacFn`

# External Links
$(_doc_external("DT/PetscDSGetJacobian"))
"""
function PetscDSGetJacobian(petsclib::PetscLibType, ds::PetscDS, f::Integer, g::Integer)
    error("PetscDSGetJacobian: no generated method for these argument types")
end

@for_petsc function PetscDSGetJacobian(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt, g::$PetscInt )
	g0_ = Ref{Ptr{Cvoid}}()
	g1_ = Ref{Ptr{Cvoid}}()
	g2_ = Ref{Ptr{Cvoid}}()
	g3_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:PetscDSGetJacobian, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, $PetscInt, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}),
               ds, f, g, g0_, g1_, g2_, g3_,
              )

	g0 = g0_[]
	g1 = g1_[]
	g2 = g2_[]
	g3 = g3_[]

	return g0,g1,g2,g3
end 

"""
	g0::Ptr{Cvoid},g1::Ptr{Cvoid},g2::Ptr{Cvoid},g3::Ptr{Cvoid} = PetscDSGetJacobianPreconditioner(petsclib::PetscLibType,ds::PetscDS, f::PetscInt, g::PetscInt) 
Get the pointwise Jacobian function for given test and basis field that constructs the matrix used
to compute the preconditioner. If this is missing, the system matrix is used to build the preconditioner.

Not Collective

Input Parameters:
- `ds` - The `PetscDS`
- `f`  - The test field number
- `g`  - The field number

Output Parameters:
- `g0` - integrand for the test and basis function term, see `PetscPointJacFn`
- `g1` - integrand for the test function and basis function gradient term, see `PetscPointJacFn`
- `g2` - integrand for the test function gradient and basis function term, see `PetscPointJacFn`
- `g3` - integrand for the test function gradient and basis function gradient term, see `PetscPointJacFn`

Level: intermediate

-seealso: `PetscDS`, `PetscDSSetJacobianPreconditioner()`, `PetscDSGetJacobian()`, `PetscPointJacFn`

# External Links
$(_doc_external("DT/PetscDSGetJacobianPreconditioner"))
"""
function PetscDSGetJacobianPreconditioner(petsclib::PetscLibType, ds::PetscDS, f::Integer, g::Integer)
    error("PetscDSGetJacobianPreconditioner: no generated method for these argument types")
end

@for_petsc function PetscDSGetJacobianPreconditioner(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt, g::$PetscInt )
	g0_ = Ref{Ptr{Cvoid}}()
	g1_ = Ref{Ptr{Cvoid}}()
	g2_ = Ref{Ptr{Cvoid}}()
	g3_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:PetscDSGetJacobianPreconditioner, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, $PetscInt, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}),
               ds, f, g, g0_, g1_, g2_, g3_,
              )

	g0 = g0_[]
	g1 = g1_[]
	g2 = g2_[]
	g3 = g3_[]

	return g0,g1,g2,g3
end 

"""
	k::PetscInt = PetscDSGetJetDegree(petsclib::PetscLibType,ds::PetscDS, f::PetscInt) 
Returns the highest derivative for this field equation, or the k

Not Collective

Input Parameters:
- `ds` - The `PetscDS` object
- `f`  - The field number

Output Parameter:
- `k` - The highest derivative we need to tabulate

Level: developer

-seealso: `PetscDS`, `PetscDSSetJetDegree()`, `PetscDSSetDiscretization()`, `PetscDSAddDiscretization()`, `PetscDSGetNumFields()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSGetJetDegree"))
"""
function PetscDSGetJetDegree(petsclib::PetscLibType, ds::PetscDS, f::Integer)
    error("PetscDSGetJetDegree: no generated method for these argument types")
end

@for_petsc function PetscDSGetJetDegree(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt )
	k_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDSGetJetDegree, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{$PetscInt}),
               ds, f, k_,
              )

	k = k_[]

	return k
end 

"""
	lb::Ptr{Cvoid},ctx::Ptr{Cvoid} = PetscDSGetLowerBound(petsclib::PetscLibType,ds::PetscDS, f::PetscInt) 
Get the pointwise lower bound function for a given field

Not Collective

Input Parameters:
- `ds` - The PetscDS
- `f`  - The field number

Output Parameters:
- `lb`  - lower bound function for the field, see `PetscPointBoundFn`
- `ctx` - lower bound context that was set with `PetscDSSetLowerBound()`

Level: intermediate

-seealso: `PetscDS`, `PetscPointBoundFn`, `PetscDSSetLowerBound()`, `PetscDSGetUpperBound()`, `PetscDSGetExactSolution()`

# External Links
$(_doc_external("DT/PetscDSGetLowerBound"))
"""
function PetscDSGetLowerBound(petsclib::PetscLibType, ds::PetscDS, f::Integer)
    error("PetscDSGetLowerBound: no generated method for these argument types")
end

@for_petsc function PetscDSGetLowerBound(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt )
	lb_ = Ref{Ptr{Cvoid}}()
	ctx_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:PetscDSGetLowerBound, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}),
               ds, f, lb_, ctx_,
              )

	lb = lb_[]
	ctx = ctx_[]

	return lb,ctx
end 

"""
	numBd::PetscInt = PetscDSGetNumBoundary(petsclib::PetscLibType,ds::PetscDS) 
Get the number of registered boundary conditions

Input Parameter:
- `ds` - The `PetscDS` object

Output Parameter:
- `numBd` - The number of boundary conditions

Level: intermediate

-seealso: `PetscDS`, `PetscDSAddBoundary()`, `PetscDSGetBoundary()`

# External Links
$(_doc_external("DT/PetscDSGetNumBoundary"))
"""
function PetscDSGetNumBoundary(petsclib::PetscLibType, ds::PetscDS)
    error("PetscDSGetNumBoundary: no generated method for these argument types")
end

@for_petsc function PetscDSGetNumBoundary(petsclib::$UnionPetscLib, ds::PetscDS )
	numBd_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDSGetNumBoundary, $petsc_library),
               PetscErrorCode,
               (PetscDS, Ptr{$PetscInt}),
               ds, numBd_,
              )

	numBd = numBd_[]

	return numBd
end 

"""
	numCohesive::PetscInt = PetscDSGetNumCohesive(petsclib::PetscLibType,ds::PetscDS) 
Returns the number of cohesive fields, meaning those defined on the interior of a cohesive cell

Not Collective

Input Parameter:
- `ds` - The `PetscDS` object

Output Parameter:
- `numCohesive` - The number of cohesive fields

Level: developer

-seealso: `PetscDS`, `PetscDSSetCohesive()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSGetNumCohesive"))
"""
function PetscDSGetNumCohesive(petsclib::PetscLibType, ds::PetscDS)
    error("PetscDSGetNumCohesive: no generated method for these argument types")
end

@for_petsc function PetscDSGetNumCohesive(petsclib::$UnionPetscLib, ds::PetscDS )
	numCohesive_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDSGetNumCohesive, $petsc_library),
               PetscErrorCode,
               (PetscDS, Ptr{$PetscInt}),
               ds, numCohesive_,
              )

	numCohesive = numCohesive_[]

	return numCohesive
end 

"""
	Nf::PetscInt = PetscDSGetNumFields(petsclib::PetscLibType,prob::PetscDS) 
Returns the number of fields in the `PetscDS`

Not Collective

Input Parameter:
- `prob` - The `PetscDS` object

Output Parameter:
- `Nf` - The number of fields

Level: beginner

-seealso: `PetscDS`, `PetscDSGetSpatialDimension()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSGetNumFields"))
"""
function PetscDSGetNumFields(petsclib::PetscLibType, prob::PetscDS)
    error("PetscDSGetNumFields: no generated method for these argument types")
end

@for_petsc function PetscDSGetNumFields(petsclib::$UnionPetscLib, prob::PetscDS )
	Nf_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDSGetNumFields, $petsc_library),
               PetscErrorCode,
               (PetscDS, Ptr{$PetscInt}),
               prob, Nf_,
              )

	Nf = Nf_[]

	return Nf
end 

"""
	obj::Ptr{Cvoid} = PetscDSGetObjective(petsclib::PetscLibType,ds::PetscDS, f::PetscInt) 
Get the pointwise objective function for a given test field that was provided with `PetscDSSetObjective()`

Not Collective

Input Parameters:
- `ds` - The `PetscDS`
- `f`  - The test field number

Output Parameter:
- `obj` - integrand for the test function term, see `PetscPointFn`

Level: intermediate

-seealso: `PetscPointFn`, `PetscDS`, `PetscDSSetObjective()`, `PetscDSGetResidual()`

# External Links
$(_doc_external("DT/PetscDSGetObjective"))
"""
function PetscDSGetObjective(petsclib::PetscLibType, ds::PetscDS, f::Integer)
    error("PetscDSGetObjective: no generated method for these argument types")
end

@for_petsc function PetscDSGetObjective(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt )
	obj_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:PetscDSGetObjective, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{Ptr{Cvoid}}),
               ds, f, obj_,
              )

	obj = obj_[]

	return obj
end 

"""
	q::PetscQuadrature = PetscDSGetQuadrature(petsclib::PetscLibType,prob::PetscDS) 
Returns the quadrature, which must agree for all fields in the `PetscDS`

Not Collective

Input Parameter:
- `prob` - The `PetscDS` object

Output Parameter:
- `q` - The quadrature object

Level: intermediate

-seealso: `PetscDS`, `PetscQuadrature`, `PetscDSSetImplicit()`, `PetscDSSetDiscretization()`, `PetscDSAddDiscretization()`, `PetscDSGetNumFields()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSGetQuadrature"))
"""
function PetscDSGetQuadrature(petsclib::PetscLibType, prob::PetscDS)
    error("PetscDSGetQuadrature: no generated method for these argument types")
end

@for_petsc function PetscDSGetQuadrature(petsclib::$UnionPetscLib, prob::PetscDS )
	q_ = Ref{PetscQuadrature}()

    @chk ccall(
               (:PetscDSGetQuadrature, $petsc_library),
               PetscErrorCode,
               (PetscDS, Ptr{PetscQuadrature}),
               prob, q_,
              )

	q = q_[]

	return q
end 

"""
	f0::Ptr{Cvoid},f1::Ptr{Cvoid} = PetscDSGetRHSResidual(petsclib::PetscLibType,ds::PetscDS, f::PetscInt) 
Get the pointwise RHS residual function for explicit timestepping for a given test field

Not Collective

Input Parameters:
- `ds` - The `PetscDS`
- `f`  - The test field number

Output Parameters:
- `f0` - integrand for the test function term, see `PetscPointFn`
- `f1` - integrand for the test function gradient term, see `PetscPointFn`

Level: intermediate

-seealso: `PetscPointFn`, `PetscDS`, `PetscDSSetRHSResidual()`

# External Links
$(_doc_external("DT/PetscDSGetRHSResidual"))
"""
function PetscDSGetRHSResidual(petsclib::PetscLibType, ds::PetscDS, f::Integer)
    error("PetscDSGetRHSResidual: no generated method for these argument types")
end

@for_petsc function PetscDSGetRHSResidual(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt )
	f0_ = Ref{Ptr{Cvoid}}()
	f1_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:PetscDSGetRHSResidual, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}),
               ds, f, f0_, f1_,
              )

	f0 = f0_[]
	f1 = f1_[]

	return f0,f1
end 

"""
	f0::Ptr{Cvoid},f1::Ptr{Cvoid} = PetscDSGetResidual(petsclib::PetscLibType,ds::PetscDS, f::PetscInt) 
Get the pointwise residual function for a given test field

Not Collective

Input Parameters:
- `ds` - The `PetscDS`
- `f`  - The test field number

Output Parameters:
- `f0` - integrand for the test function term, see `PetscPointFn`
- `f1` - integrand for the test function gradient term, see `PetscPointFn`

Level: intermediate

-seealso: `PetscPointFn`, `PetscDS`, `PetscDSSetResidual()`

# External Links
$(_doc_external("DT/PetscDSGetResidual"))
"""
function PetscDSGetResidual(petsclib::PetscLibType, ds::PetscDS, f::Integer)
    error("PetscDSGetResidual: no generated method for these argument types")
end

@for_petsc function PetscDSGetResidual(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt )
	f0_ = Ref{Ptr{Cvoid}}()
	f1_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:PetscDSGetResidual, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}),
               ds, f, f0_, f1_,
              )

	f0 = f0_[]
	f1 = f1_[]

	return f0,f1
end 

"""
	r::Ptr{Cvoid} = PetscDSGetRiemannSolver(petsclib::PetscLibType,ds::PetscDS, f::PetscInt) 
Returns the Riemann solver for the given field

Not Collective

Input Parameters:
- `ds` - The `PetscDS` object
- `f`  - The field number

Output Parameter:
- `r` - Riemann solver, see `PetscRiemannFn`

Level: intermediate

-seealso: `PetscDS`, `PetscRiemannFn`, `PetscDSSetRiemannSolver()`

# External Links
$(_doc_external("DT/PetscDSGetRiemannSolver"))
"""
function PetscDSGetRiemannSolver(petsclib::PetscLibType, ds::PetscDS, f::Integer)
    error("PetscDSGetRiemannSolver: no generated method for these argument types")
end

@for_petsc function PetscDSGetRiemannSolver(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt )
	r_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:PetscDSGetRiemannSolver, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{Ptr{Cvoid}}),
               ds, f, r_,
              )

	r = r_[]

	return r
end 

"""
	dim::PetscInt = PetscDSGetSpatialDimension(petsclib::PetscLibType,prob::PetscDS) 
Returns the spatial dimension of the `PetscDS`, meaning the topological dimension of the discretizations

Not Collective

Input Parameter:
- `prob` - The `PetscDS` object

Output Parameter:
- `dim` - The spatial dimension

Level: beginner

-seealso: `PetscDS`, `PetscDSGetCoordinateDimension()`, `PetscDSGetNumFields()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSGetSpatialDimension"))
"""
function PetscDSGetSpatialDimension(petsclib::PetscLibType, prob::PetscDS)
    error("PetscDSGetSpatialDimension: no generated method for these argument types")
end

@for_petsc function PetscDSGetSpatialDimension(petsclib::$UnionPetscLib, prob::PetscDS )
	dim_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDSGetSpatialDimension, $petsc_library),
               PetscErrorCode,
               (PetscDS, Ptr{$PetscInt}),
               prob, dim_,
              )

	dim = dim_[]

	return dim
end 

"""
	T::Ptr{PetscTabulation} = PetscDSGetTabulation(petsclib::PetscLibType,prob::PetscDS) 
Return the basis tabulation at quadrature points for the volume discretization

Not Collective

Input Parameter:
- `prob` - The `PetscDS` object

Output Parameter:
- `T` - The basis function and derivatives tabulation at quadrature points for each field, see `PetscTabulation` for its details

Level: intermediate

-seealso: `PetscDS`, `PetscTabulation`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSGetTabulation"))
"""
function PetscDSGetTabulation(petsclib::PetscLibType, prob::PetscDS)
    error("PetscDSGetTabulation: no generated method for these argument types")
end

@for_petsc function PetscDSGetTabulation(petsclib::$UnionPetscLib, prob::PetscDS )
	T_ = Ref{Ptr{PetscTabulation}}()

    @chk ccall(
               (:PetscDSGetTabulation, $petsc_library),
               PetscErrorCode,
               (PetscDS, Ptr{Ptr{PetscTabulation}}),
               prob, T_,
              )

	T = T_[]

	return T
end 

"""
	Nc::PetscInt = PetscDSGetTotalComponents(petsclib::PetscLibType,prob::PetscDS) 
Returns the total number of components in this system

Not Collective

Input Parameter:
- `prob` - The `PetscDS` object

Output Parameter:
- `Nc` - The total number of components

Level: beginner

-seealso: `PetscDS`, `PetscDSGetNumFields()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSGetTotalComponents"))
"""
function PetscDSGetTotalComponents(petsclib::PetscLibType, prob::PetscDS)
    error("PetscDSGetTotalComponents: no generated method for these argument types")
end

@for_petsc function PetscDSGetTotalComponents(petsclib::$UnionPetscLib, prob::PetscDS )
	Nc_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDSGetTotalComponents, $petsc_library),
               PetscErrorCode,
               (PetscDS, Ptr{$PetscInt}),
               prob, Nc_,
              )

	Nc = Nc_[]

	return Nc
end 

"""
	dim::PetscInt = PetscDSGetTotalDimension(petsclib::PetscLibType,prob::PetscDS) 
Returns the total size of the approximation space for this system

Not Collective

Input Parameter:
- `prob` - The `PetscDS` object

Output Parameter:
- `dim` - The total problem dimension

Level: beginner

-seealso: `PetscDS`, `PetscDSGetNumFields()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSGetTotalDimension"))
"""
function PetscDSGetTotalDimension(petsclib::PetscLibType, prob::PetscDS)
    error("PetscDSGetTotalDimension: no generated method for these argument types")
end

@for_petsc function PetscDSGetTotalDimension(petsclib::$UnionPetscLib, prob::PetscDS )
	dim_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDSGetTotalDimension, $petsc_library),
               PetscErrorCode,
               (PetscDS, Ptr{$PetscInt}),
               prob, dim_,
              )

	dim = dim_[]

	return dim
end 

"""
	name::PetscDSType = PetscDSGetType(petsclib::PetscLibType,prob::PetscDS) 
Gets the `PetscDSType` name (as a string) from the `PetscDS`

Not Collective; No Fortran Support

Input Parameter:
- `prob` - The `PetscDS`

Output Parameter:
- `name` - The `PetscDSType` name

Level: intermediate

-seealso: `PetscDSType`, `PetscDS`, `PetscDSSetType()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSGetType"))
"""
function PetscDSGetType(petsclib::PetscLibType, prob::PetscDS)
    error("PetscDSGetType: no generated method for these argument types")
end

@for_petsc function PetscDSGetType(petsclib::$UnionPetscLib, prob::PetscDS )
	name_ = Ref{PetscDSType}()

    @chk ccall(
               (:PetscDSGetType, $petsc_library),
               PetscErrorCode,
               (PetscDS, Ptr{PetscDSType}),
               prob, name_,
              )

	name = name_[] == C_NULL ? "" : unsafe_string(name_[])

	return name
end 

"""
	update::Ptr{Cvoid} = PetscDSGetUpdate(petsclib::PetscLibType,ds::PetscDS, f::PetscInt) 
Get the pointwise update function for a given field

Not Collective

Input Parameters:
- `ds` - The `PetscDS`
- `f`  - The field number

Output Parameter:
- `update` - update function, see `PetscPointFn`

Level: intermediate

-seealso: `PetscDS`, `PetscPointFn`, `PetscDSSetUpdate()`, `PetscDSSetResidual()`

# External Links
$(_doc_external("DT/PetscDSGetUpdate"))
"""
function PetscDSGetUpdate(petsclib::PetscLibType, ds::PetscDS, f::Integer)
    error("PetscDSGetUpdate: no generated method for these argument types")
end

@for_petsc function PetscDSGetUpdate(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt )
	update_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:PetscDSGetUpdate, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{Ptr{Cvoid}}),
               ds, f, update_,
              )

	update = update_[]

	return update
end 

"""
	ub::Ptr{Cvoid},ctx::Ptr{Cvoid} = PetscDSGetUpperBound(petsclib::PetscLibType,ds::PetscDS, f::PetscInt) 
Get the pointwise upper bound function for a given field

Not Collective

Input Parameters:
- `ds` - The `PetscDS`
- `f`  - The field number

Output Parameters:
- `ub`  - upper bound function for the field, see `PetscPointBoundFn`
- `ctx` - upper bound context that was set with `PetscDSSetUpperBound()`

Level: intermediate

-seealso: `PetscDS`, `PetscPointBoundFn`, `PetscDSSetUpperBound()`, `PetscDSGetLowerBound()`, `PetscDSGetExactSolution()`

# External Links
$(_doc_external("DT/PetscDSGetUpperBound"))
"""
function PetscDSGetUpperBound(petsclib::PetscLibType, ds::PetscDS, f::Integer)
    error("PetscDSGetUpperBound: no generated method for these argument types")
end

@for_petsc function PetscDSGetUpperBound(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt )
	ub_ = Ref{Ptr{Cvoid}}()
	ctx_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:PetscDSGetUpperBound, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}),
               ds, f, ub_, ctx_,
              )

	ub = ub_[]
	ctx = ctx_[]

	return ub,ctx
end 

"""
	wf::PetscWeakForm = PetscDSGetWeakForm(petsclib::PetscLibType,ds::PetscDS) 
Returns the weak form object from within the `PetscDS`

Not Collective

Input Parameter:
- `ds` - The `PetscDS` object

Output Parameter:
- `wf` - The weak form object

Level: beginner

-seealso: `PetscWeakForm`, `PetscDSSetWeakForm()`, `PetscDSGetNumFields()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSGetWeakForm"))
"""
function PetscDSGetWeakForm(petsclib::PetscLibType, ds::PetscDS)
    error("PetscDSGetWeakForm: no generated method for these argument types")
end

@for_petsc function PetscDSGetWeakForm(petsclib::$UnionPetscLib, ds::PetscDS )
	wf_ = Ref{PetscWeakForm}()

    @chk ccall(
               (:PetscDSGetWeakForm, $petsc_library),
               PetscErrorCode,
               (PetscDS, Ptr{PetscWeakForm}),
               ds, wf_,
              )

	wf = wf_[]

	return wf
end 

"""
	f0::Ptr{PetscScalar},f1::Ptr{PetscScalar},g0::Ptr{PetscScalar},g1::Ptr{PetscScalar},g2::Ptr{PetscScalar},g3::Ptr{PetscScalar} = PetscDSGetWeakFormArrays(petsclib::PetscLibType,prob::PetscDS) 

# External Links
$(_doc_external("DT/PetscDSGetWeakFormArrays"))
"""
function PetscDSGetWeakFormArrays(petsclib::PetscLibType, prob::PetscDS)
    error("PetscDSGetWeakFormArrays: no generated method for these argument types")
end

@for_petsc function PetscDSGetWeakFormArrays(petsclib::$UnionPetscLib, prob::PetscDS )
	f0_ = Ref{Ptr{$PetscScalar}}()
	f1_ = Ref{Ptr{$PetscScalar}}()
	g0_ = Ref{Ptr{$PetscScalar}}()
	g1_ = Ref{Ptr{$PetscScalar}}()
	g2_ = Ref{Ptr{$PetscScalar}}()
	g3_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:PetscDSGetWeakFormArrays, $petsc_library),
               PetscErrorCode,
               (PetscDS, Ptr{Ptr{$PetscScalar}}, Ptr{Ptr{$PetscScalar}}, Ptr{Ptr{$PetscScalar}}, Ptr{Ptr{$PetscScalar}}, Ptr{Ptr{$PetscScalar}}, Ptr{Ptr{$PetscScalar}}),
               prob, f0_, f1_, g0_, g1_, g2_, g3_,
              )

	f0 = f0_[]
	f1 = f1_[]
	g0 = g0_[]
	g1 = g1_[]
	g2 = g2_[]
	g3 = g3_[]

	return f0,f1,g0,g1,g2,g3
end 

"""
	x::Ptr{PetscReal},basisReal::Ptr{PetscScalar},basisDerReal::Ptr{PetscScalar},testReal::Ptr{PetscScalar},testDerReal::Ptr{PetscScalar} = PetscDSGetWorkspace(petsclib::PetscLibType,prob::PetscDS) 
Get scratch storage used during discretization computations.

Not Collective

Input Parameter:
- `prob` - the `PetscDS`

Output Parameters:
- `x`            - array for real-valued quadrature point coordinates, or `NULL` if not needed
- `basisReal`    - array for the real-valued basis function values, or `NULL` if not needed
- `basisDerReal` - array for the real-valued basis function derivatives, or `NULL` if not needed
- `testReal`     - array for the real-valued test function values, or `NULL` if not needed
- `testDerReal`  - array for the real-valued test function derivatives, or `NULL` if not needed

Level: developer

-seealso: `PetscDS`, `PetscDSGetEvaluationArrays()`, `PetscDSGetWeakFormArrays()`

# External Links
$(_doc_external("DT/PetscDSGetWorkspace"))
"""
function PetscDSGetWorkspace(petsclib::PetscLibType, prob::PetscDS)
    error("PetscDSGetWorkspace: no generated method for these argument types")
end

@for_petsc function PetscDSGetWorkspace(petsclib::$UnionPetscLib, prob::PetscDS )
	x_ = Ref{Ptr{$PetscReal}}()
	basisReal_ = Ref{Ptr{$PetscScalar}}()
	basisDerReal_ = Ref{Ptr{$PetscScalar}}()
	testReal_ = Ref{Ptr{$PetscScalar}}()
	testDerReal_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:PetscDSGetWorkspace, $petsc_library),
               PetscErrorCode,
               (PetscDS, Ptr{Ptr{$PetscReal}}, Ptr{Ptr{$PetscScalar}}, Ptr{Ptr{$PetscScalar}}, Ptr{Ptr{$PetscScalar}}, Ptr{Ptr{$PetscScalar}}),
               prob, x_, basisReal_, basisDerReal_, testReal_, testDerReal_,
              )

	x = x_[]
	basisReal = basisReal_[]
	basisDerReal = basisDerReal_[]
	testReal = testReal_[]
	testDerReal = testDerReal_[]

	return x,basisReal,basisDerReal,testReal,testDerReal
end 

"""
	hasBdJac::PetscBool = PetscDSHasBdJacobian(petsclib::PetscLibType,ds::PetscDS) 
Indicates that boundary Jacobian functions have been set

Not Collective

Input Parameter:
- `ds` - The `PetscDS`

Output Parameter:
- `hasBdJac` - flag that pointwise function for the boundary Jacobian has been set

Level: intermediate

-seealso: `PetscDS`, `PetscDSHasJacobian()`, `PetscDSSetBdJacobian()`, `PetscDSGetBdJacobian()`

# External Links
$(_doc_external("DT/PetscDSHasBdJacobian"))
"""
function PetscDSHasBdJacobian(petsclib::PetscLibType, ds::PetscDS)
    error("PetscDSHasBdJacobian: no generated method for these argument types")
end

@for_petsc function PetscDSHasBdJacobian(petsclib::$UnionPetscLib, ds::PetscDS )
	hasBdJac_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscDSHasBdJacobian, $petsc_library),
               PetscErrorCode,
               (PetscDS, Ptr{PetscBool}),
               ds, hasBdJac_,
              )

	hasBdJac = hasBdJac_[]

	return hasBdJac
end 

"""
	hasBdJacPre::PetscBool = PetscDSHasBdJacobianPreconditioner(petsclib::PetscLibType,ds::PetscDS) 
Signals that boundary Jacobian preconditioner functions have been set with `PetscDSSetBdJacobianPreconditioner()`

Not Collective

Input Parameter:
- `ds` - The `PetscDS`

Output Parameter:
- `hasBdJacPre` - flag that pointwise function for the boundary Jacobian matrix to construct the preconditioner has been set

Level: intermediate

-seealso: `PetscDS`, `PetscDSHasJacobian()`, `PetscDSSetBdJacobian()`, `PetscDSGetBdJacobian()`

# External Links
$(_doc_external("DT/PetscDSHasBdJacobianPreconditioner"))
"""
function PetscDSHasBdJacobianPreconditioner(petsclib::PetscLibType, ds::PetscDS)
    error("PetscDSHasBdJacobianPreconditioner: no generated method for these argument types")
end

@for_petsc function PetscDSHasBdJacobianPreconditioner(petsclib::$UnionPetscLib, ds::PetscDS )
	hasBdJacPre_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscDSHasBdJacobianPreconditioner, $petsc_library),
               PetscErrorCode,
               (PetscDS, Ptr{PetscBool}),
               ds, hasBdJacPre_,
              )

	hasBdJacPre = hasBdJacPre_[]

	return hasBdJacPre
end 

"""
	hasDynJac::PetscBool = PetscDSHasDynamicJacobian(petsclib::PetscLibType,ds::PetscDS) 
Signals that a dynamic Jacobian, dF/du_t, has been set

Not Collective

Input Parameter:
- `ds` - The `PetscDS`

Output Parameter:
- `hasDynJac` - flag that pointwise function for dynamic Jacobian has been set

Level: intermediate

-seealso: `PetscDS`, `PetscDSGetDynamicJacobian()`, `PetscDSSetDynamicJacobian()`, `PetscDSGetJacobian()`

# External Links
$(_doc_external("DT/PetscDSHasDynamicJacobian"))
"""
function PetscDSHasDynamicJacobian(petsclib::PetscLibType, ds::PetscDS)
    error("PetscDSHasDynamicJacobian: no generated method for these argument types")
end

@for_petsc function PetscDSHasDynamicJacobian(petsclib::$UnionPetscLib, ds::PetscDS )
	hasDynJac_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscDSHasDynamicJacobian, $petsc_library),
               PetscErrorCode,
               (PetscDS, Ptr{PetscBool}),
               ds, hasDynJac_,
              )

	hasDynJac = hasDynJac_[]

	return hasDynJac
end 

"""
	hasJac::PetscBool = PetscDSHasJacobian(petsclib::PetscLibType,ds::PetscDS) 
Checks that the Jacobian functions have been set

Not Collective

Input Parameter:
- `ds` - The `PetscDS`

Output Parameter:
- `hasJac` - flag that indicates the pointwise function for the Jacobian has been set

Level: intermediate

-seealso: `PetscDS`, `PetscDSGetJacobianPreconditioner()`, `PetscDSSetJacobianPreconditioner()`, `PetscDSGetJacobian()`

# External Links
$(_doc_external("DT/PetscDSHasJacobian"))
"""
function PetscDSHasJacobian(petsclib::PetscLibType, ds::PetscDS)
    error("PetscDSHasJacobian: no generated method for these argument types")
end

@for_petsc function PetscDSHasJacobian(petsclib::$UnionPetscLib, ds::PetscDS )
	hasJac_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscDSHasJacobian, $petsc_library),
               PetscErrorCode,
               (PetscDS, Ptr{PetscBool}),
               ds, hasJac_,
              )

	hasJac = hasJac_[]

	return hasJac
end 

"""
	hasJacPre::PetscBool = PetscDSHasJacobianPreconditioner(petsclib::PetscLibType,ds::PetscDS) 
Checks if a Jacobian matrix for constructing a preconditioner has been set

Not Collective

Input Parameter:
- `ds` - The `PetscDS`

Output Parameter:
- `hasJacPre` - the flag

Level: intermediate

-seealso: `PetscDS`, `PetscDSGetJacobianPreconditioner()`, `PetscDSSetJacobianPreconditioner()`, `PetscDSGetJacobian()`

# External Links
$(_doc_external("DT/PetscDSHasJacobianPreconditioner"))
"""
function PetscDSHasJacobianPreconditioner(petsclib::PetscLibType, ds::PetscDS)
    error("PetscDSHasJacobianPreconditioner: no generated method for these argument types")
end

@for_petsc function PetscDSHasJacobianPreconditioner(petsclib::$UnionPetscLib, ds::PetscDS )
	hasJacPre_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscDSHasJacobianPreconditioner, $petsc_library),
               PetscErrorCode,
               (PetscDS, Ptr{PetscBool}),
               ds, hasJacPre_,
              )

	hasJacPre = hasJacPre_[]

	return hasJacPre
end 

"""
	PetscDSInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `PetscDS` package. It is called
from `PetscDLLibraryRegister()` when using dynamic libraries, and on the first call to `PetscDSCreate()`
when using static libraries.

Level: developer

-seealso: `PetscInitialize()`

# External Links
$(_doc_external("DM/PetscDSInitializePackage"))
"""
function PetscDSInitializePackage(petsclib::PetscLibType)
    error("PetscDSInitializePackage: no generated method for these argument types")
end

@for_petsc function PetscDSInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscDSInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	isCohesive::PetscBool = PetscDSIsCohesive(petsclib::PetscLibType,ds::PetscDS) 
Returns the flag indicating that this `PetscDS` is for a cohesive cell

Not Collective

Input Parameter:
- `ds` - The `PetscDS` object

Output Parameter:
- `isCohesive` - The flag

Level: developer

-seealso: `PetscDS`, `PetscDSGetNumCohesive()`, `PetscDSGetCohesive()`, `PetscDSSetCohesive()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSIsCohesive"))
"""
function PetscDSIsCohesive(petsclib::PetscLibType, ds::PetscDS)
    error("PetscDSIsCohesive: no generated method for these argument types")
end

@for_petsc function PetscDSIsCohesive(petsclib::$UnionPetscLib, ds::PetscDS )
	isCohesive_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscDSIsCohesive, $petsc_library),
               PetscErrorCode,
               (PetscDS, Ptr{PetscBool}),
               ds, isCohesive_,
              )

	isCohesive = isCohesive_[]

	return isCohesive
end 

"""
	qperm::PetscInt = PetscDSPermuteQuadPoint(petsclib::PetscLibType,ds::PetscDS, ornt::PetscInt, field::PetscInt, q::PetscInt) 
Permute a quadrature point index according to a cell orientation.

Not Collective

Input Parameters:
- `ds`    - the `PetscDS`
- `ornt`  - the cell orientation, in `[-Na, Na)` where `Na` is half the number of arrangements for the cell type
- `field` - the field number whose quadrature is used
- `q`     - the input quadrature point index in `[0, Nq)`

Output Parameter:
- `qperm` - the permuted quadrature point index

Level: developer

-seealso: `PetscDS`, `PetscQuadrature`, `PetscQuadratureComputePermutations()`, `DMPolytopeTypeGetNumArrangements()`

# External Links
$(_doc_external("DT/PetscDSPermuteQuadPoint"))
"""
function PetscDSPermuteQuadPoint(petsclib::PetscLibType, ds::PetscDS, ornt::Integer, field::Integer, q::Integer)
    error("PetscDSPermuteQuadPoint: no generated method for these argument types")
end

@for_petsc function PetscDSPermuteQuadPoint(petsclib::$UnionPetscLib, ds::PetscDS, ornt::$PetscInt, field::$PetscInt, q::$PetscInt )
	qperm_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDSPermuteQuadPoint, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}),
               ds, ornt, field, q, qperm_,
              )

	qperm = qperm_[]

	return qperm
end 

"""
	PetscDSRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Adds a new `PetscDS` implementation

Not Collective; No Fortran Support

Input Parameters:
- `sname`    - The name of a new user-defined creation routine
- `function` - The creation routine itself

-seealso: `PetscDSType`, `PetscDS`, `PetscDSRegisterAll()`, `PetscDSRegisterDestroy()`

# External Links
$(_doc_external("DT/PetscDSRegister"))
"""
function PetscDSRegister(petsclib::PetscLibType, sname::String, fnc::external)
    error("PetscDSRegister: no generated method for these argument types")
end

@for_petsc function PetscDSRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:PetscDSRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	PetscDSSelectDiscretizations(petsclib::PetscLibType,prob::PetscDS, numFields::PetscInt, fields::Vector{PetscInt}, minDegree::PetscInt, maxDegree::PetscInt, newprob::PetscDS) 
Copy discretizations to the new `PetscDS` with different field layout

Not Collective

Input Parameters:
- `prob`      - The `PetscDS` object
- `numFields` - Number of new fields
- `fields`    - Old field number for each new field
- `minDegree` - Minimum degree for a discretization, or `PETSC_DETERMINE` for no limit
- `maxDegree` - Maximum degree for a discretization, or `PETSC_DETERMINE` for no limit

Output Parameter:
- `newprob` - The `PetscDS` copy

Level: intermediate

-seealso: `PetscDS`, `PetscDSSelectEquations()`, `PetscDSCopyBoundary()`, `PetscDSSetResidual()`, `PetscDSSetJacobian()`, `PetscDSSetRiemannSolver()`, `PetscDSSetBdResidual()`, `PetscDSSetBdJacobian()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSSelectDiscretizations"))
"""
function PetscDSSelectDiscretizations(petsclib::PetscLibType, prob::PetscDS, numFields::Integer, fields::AbstractVector{<:Number}, minDegree::Integer, maxDegree::Integer, newprob::PetscDS)
    error("PetscDSSelectDiscretizations: no generated method for these argument types")
end

@for_petsc function PetscDSSelectDiscretizations(petsclib::$UnionPetscLib, prob::PetscDS, numFields::$PetscInt, fields::Vector{$PetscInt}, minDegree::$PetscInt, maxDegree::$PetscInt, newprob::PetscDS )

    @chk ccall(
               (:PetscDSSelectDiscretizations, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{$PetscInt}, $PetscInt, $PetscInt, PetscDS),
               prob, numFields, fields, minDegree, maxDegree, newprob,
              )


	return nothing
end 

"""
	PetscDSSelectEquations(petsclib::PetscLibType,prob::PetscDS, numFields::PetscInt, fields::Vector{PetscInt}, newprob::PetscDS) 
Copy pointwise function pointers to the new `PetscDS` with different field layout

Not Collective

Input Parameters:
- `prob`      - The `PetscDS` object
- `numFields` - Number of new fields
- `fields`    - Old field number for each new field

Output Parameter:
- `newprob` - The `PetscDS` copy

Level: intermediate

-seealso: `PetscDS`, `PetscDSSelectDiscretizations()`, `PetscDSCopyBoundary()`, `PetscDSSetResidual()`, `PetscDSSetJacobian()`, `PetscDSSetRiemannSolver()`, `PetscDSSetBdResidual()`, `PetscDSSetBdJacobian()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSSelectEquations"))
"""
function PetscDSSelectEquations(petsclib::PetscLibType, prob::PetscDS, numFields::Integer, fields::AbstractVector{<:Number}, newprob::PetscDS)
    error("PetscDSSelectEquations: no generated method for these argument types")
end

@for_petsc function PetscDSSelectEquations(petsclib::$UnionPetscLib, prob::PetscDS, numFields::$PetscInt, fields::Vector{$PetscInt}, newprob::PetscDS )

    @chk ccall(
               (:PetscDSSelectEquations, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{$PetscInt}, PetscDS),
               prob, numFields, fields, newprob,
              )


	return nothing
end 

"""
	PetscDSSetBdJacobian(petsclib::PetscLibType,ds::PetscDS, f::PetscInt, g::PetscInt, g0::Ptr{Cvoid}, g1::Ptr{Cvoid}, g2::Ptr{Cvoid}, g3::Ptr{Cvoid}) 
Set the pointwise boundary Jacobian function for given test and basis field

Not Collective

Input Parameters:
- `ds` - The PetscDS
- `f`  - The test field number
- `g`  - The field number
- `g0` - integrand for the test and basis function term, see `PetscBdPointJacFn`
- `g1` - integrand for the test function and basis function gradient term, see `PetscBdPointJacFn`
- `g2` - integrand for the test function gradient and basis function term, see `PetscBdPointJacFn`
- `g3` - integrand for the test function gradient and basis function gradient term, see `PetscBdPointJacFn`

Level: intermediate

-seealso: `PetscDS`, `PetscBdPointJacFn`, `PetscDSGetBdJacobian()`

# External Links
$(_doc_external("DT/PetscDSSetBdJacobian"))
"""
function PetscDSSetBdJacobian(petsclib::PetscLibType, ds::PetscDS, f::Integer, g::Integer, g0::Ptr{Cvoid}, g1::Ptr{Cvoid}, g2::Ptr{Cvoid}, g3::Ptr{Cvoid})
    error("PetscDSSetBdJacobian: no generated method for these argument types")
end

@for_petsc function PetscDSSetBdJacobian(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt, g::$PetscInt, g0::Ptr{Cvoid}, g1::Ptr{Cvoid}, g2::Ptr{Cvoid}, g3::Ptr{Cvoid} )

    @chk ccall(
               (:PetscDSSetBdJacobian, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, $PetscInt, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
               ds, f, g, g0, g1, g2, g3,
              )


	return nothing
end 

"""
	PetscDSSetBdJacobianPreconditioner(petsclib::PetscLibType,ds::PetscDS, f::PetscInt, g::PetscInt, g0::Ptr{Cvoid}, g1::Ptr{Cvoid}, g2::Ptr{Cvoid}, g3::Ptr{Cvoid}) 
Set the pointwise boundary Jacobian preconditioner function for given test and basis field that constructs the
matrix used to construct the preconditioner

Not Collective; No Fortran Support

Input Parameters:
- `ds` - The `PetscDS`
- `f`  - The test field number
- `g`  - The field number
- `g0` - integrand for the test and basis function term, see `PetscBdPointJacFn`
- `g1` - integrand for the test function and basis function gradient term, see `PetscBdPointJacFn`
- `g2` - integrand for the test function gradient and basis function term, see `PetscBdPointJacFn`
- `g3` - integrand for the test function gradient and basis function gradient term, see `PetscBdPointJacFn`

Level: intermediate

-seealso: `PetscDS`, `PetscBdPointJacFn`, `PetscDSGetBdJacobianPreconditioner()`

# External Links
$(_doc_external("DT/PetscDSSetBdJacobianPreconditioner"))
"""
function PetscDSSetBdJacobianPreconditioner(petsclib::PetscLibType, ds::PetscDS, f::Integer, g::Integer, g0::Ptr{Cvoid}, g1::Ptr{Cvoid}, g2::Ptr{Cvoid}, g3::Ptr{Cvoid})
    error("PetscDSSetBdJacobianPreconditioner: no generated method for these argument types")
end

@for_petsc function PetscDSSetBdJacobianPreconditioner(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt, g::$PetscInt, g0::Ptr{Cvoid}, g1::Ptr{Cvoid}, g2::Ptr{Cvoid}, g3::Ptr{Cvoid} )

    @chk ccall(
               (:PetscDSSetBdJacobianPreconditioner, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, $PetscInt, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
               ds, f, g, g0, g1, g2, g3,
              )


	return nothing
end 

"""
	PetscDSSetBdResidual(petsclib::PetscLibType,ds::PetscDS, f::PetscInt, f0::Ptr{Cvoid}, f1::Ptr{Cvoid}) 
Get the pointwise boundary residual function for a given test field

Not Collective

Input Parameters:
- `ds` - The `PetscDS`
- `f`  - The test field number
- `f0` - boundary integrand for the test function term, see `PetscBdPointFn`
- `f1` - boundary integrand for the test function gradient term, see `PetscBdPointFn`

Level: intermediate

-seealso: `PetscDS`, `PetscBdPointFn`, `PetscDSGetBdResidual()`

# External Links
$(_doc_external("DT/PetscDSSetBdResidual"))
"""
function PetscDSSetBdResidual(petsclib::PetscLibType, ds::PetscDS, f::Integer, f0::Ptr{Cvoid}, f1::Ptr{Cvoid})
    error("PetscDSSetBdResidual: no generated method for these argument types")
end

@for_petsc function PetscDSSetBdResidual(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt, f0::Ptr{Cvoid}, f1::Ptr{Cvoid} )

    @chk ccall(
               (:PetscDSSetBdResidual, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{Cvoid}, Ptr{Cvoid}),
               ds, f, f0, f1,
              )


	return nothing
end 

"""
	PetscDSSetCellParameters(petsclib::PetscLibType,ds::PetscDS, volume::PetscReal) 
Set the parameters for a particular cell

Not Collective

Input Parameters:
- `ds`     - The `PetscDS` object
- `volume` - The cell volume

Level: intermediate

-seealso: `PetscDS`, `PetscDSSetConstants()`, `PetscDSGetConstants()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSSetCellParameters"))
"""
function PetscDSSetCellParameters(petsclib::PetscLibType, ds::PetscDS, volume::Real)
    error("PetscDSSetCellParameters: no generated method for these argument types")
end

@for_petsc function PetscDSSetCellParameters(petsclib::$UnionPetscLib, ds::PetscDS, volume::$PetscReal )

    @chk ccall(
               (:PetscDSSetCellParameters, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscReal),
               ds, volume,
              )


	return nothing
end 

"""
	PetscDSSetCohesive(petsclib::PetscLibType,ds::PetscDS, f::PetscInt, isCohesive::PetscBool) 
Set the flag indicating that a field is cohesive, meaning it is defined on the interior of a cohesive cell

Not Collective

Input Parameters:
- `ds`         - The `PetscDS` object
- `f`          - The field index
- `isCohesive` - The flag for a cohesive field

Level: developer

-seealso: `PetscDS`, `PetscDSGetCohesive()`, `PetscDSIsCohesive()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSSetCohesive"))
"""
function PetscDSSetCohesive(petsclib::PetscLibType, ds::PetscDS, f::Integer, isCohesive::PetscBool)
    error("PetscDSSetCohesive: no generated method for these argument types")
end

@for_petsc function PetscDSSetCohesive(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt, isCohesive::PetscBool )

    @chk ccall(
               (:PetscDSSetCohesive, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, PetscBool),
               ds, f, isCohesive,
              )


	return nothing
end 

"""
	PetscDSSetConstants(petsclib::PetscLibType,ds::PetscDS, numConstants::PetscInt, constants::Vector{PetscScalar}) 
Set the array of constants passed to point functions from a `PetscDS`

Not Collective

Input Parameters:
- `ds`           - The `PetscDS` object
- `numConstants` - The number of constants
- `constants`    - The array of constants, `NULL` if there are none

Level: intermediate

-seealso: `PetscDS`, `PetscDSGetConstants()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSSetConstants"))
"""
function PetscDSSetConstants(petsclib::PetscLibType, ds::PetscDS, numConstants::Integer, constants::AbstractVector{<:Number})
    error("PetscDSSetConstants: no generated method for these argument types")
end

@for_petsc function PetscDSSetConstants(petsclib::$UnionPetscLib, ds::PetscDS, numConstants::$PetscInt, constants::Vector{$PetscScalar} )

    @chk ccall(
               (:PetscDSSetConstants, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{$PetscScalar}),
               ds, numConstants, constants,
              )


	return nothing
end 

"""
	PetscDSSetContext(petsclib::PetscLibType,ds::PetscDS, f::PetscInt, ctx::Ptr{Cvoid}) 
Sets the context that is passed back to some of the pointwise function callbacks used by this `PetscDS`

Not Collective

Input Parameters:
- `ds`  - The `PetscDS`
- `f`   - The field number
- `ctx` - the context

Level: intermediate

-seealso: `PetscDS`, `PetscPointFn`, `PetscDSGetContext()`

# External Links
$(_doc_external("DT/PetscDSSetContext"))
"""
function PetscDSSetContext(petsclib::PetscLibType, ds::PetscDS, f::Integer, ctx::Ptr{Cvoid})
    error("PetscDSSetContext: no generated method for these argument types")
end

@for_petsc function PetscDSSetContext(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:PetscDSSetContext, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{Cvoid}),
               ds, f, ctx,
              )


	return nothing
end 

"""
	PetscDSSetCoordinateDimension(petsclib::PetscLibType,prob::PetscDS, dimEmbed::PetscInt) 
Set the coordinate dimension of the `PetscDS`, meaning the dimension of the space into which the discretiaztions are embedded

Logically Collective

Input Parameters:
- `prob`     - The `PetscDS` object
- `dimEmbed` - The coordinate dimension

Level: beginner

-seealso: `PetscDS`, `PetscDSGetCoordinateDimension()`, `PetscDSGetSpatialDimension()`, `PetscDSGetNumFields()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSSetCoordinateDimension"))
"""
function PetscDSSetCoordinateDimension(petsclib::PetscLibType, prob::PetscDS, dimEmbed::Integer)
    error("PetscDSSetCoordinateDimension: no generated method for these argument types")
end

@for_petsc function PetscDSSetCoordinateDimension(petsclib::$UnionPetscLib, prob::PetscDS, dimEmbed::$PetscInt )

    @chk ccall(
               (:PetscDSSetCoordinateDimension, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt),
               prob, dimEmbed,
              )


	return nothing
end 

"""
	PetscDSSetDiscretization(petsclib::PetscLibType,prob::PetscDS, f::PetscInt, disc) 
Sets the discretization object for the given field

Not Collective

Input Parameters:
- `prob` - The `PetscDS` object
- `f`    - The field number
- `disc` - The discretization object, this can be a `PetscFE` or a `PetscFV`

Level: beginner

-seealso: `PetscDS`, `PetscFE`, `PetscFV`, `PetscDSGetDiscretization()`, `PetscDSAddDiscretization()`, `PetscDSGetNumFields()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSSetDiscretization"))
"""
function PetscDSSetDiscretization(petsclib::PetscLibType, prob::PetscDS, f::Integer, disc)
    error("PetscDSSetDiscretization: no generated method for these argument types")
end

@for_petsc function PetscDSSetDiscretization(petsclib::$UnionPetscLib, prob::PetscDS, f::$PetscInt, disc )

    @chk ccall(
               (:PetscDSSetDiscretization, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, PetscObject),
               prob, f, disc,
              )


	return nothing
end 

"""
	PetscDSSetDynamicJacobian(petsclib::PetscLibType,ds::PetscDS, f::PetscInt, g::PetscInt, g0::Ptr{Cvoid}, g1::Ptr{Cvoid}, g2::Ptr{Cvoid}, g3::Ptr{Cvoid}) 
Set the pointwise dynamic Jacobian, dF/du_t, function for given test and basis fields

Not Collective

Input Parameters:
- `ds` - The `PetscDS`
- `f`  - The test field number
- `g`  - The field number
- `g0` - integrand for the test and basis function term, see `PetscPointJacFn`
- `g1` - integrand for the test function and basis function gradient term, see `PetscPointJacFn`
- `g2` - integrand for the test function gradient and basis function term, see `PetscPointJacFn`
- `g3` - integrand for the test function gradient and basis function gradient term, see `PetscPointJacFn`

Level: intermediate

-seealso: `PetscDS`, `PetscDSGetDynamicJacobian()`, `PetscDSGetJacobian()`, `PetscPointJacFn`

# External Links
$(_doc_external("DT/PetscDSSetDynamicJacobian"))
"""
function PetscDSSetDynamicJacobian(petsclib::PetscLibType, ds::PetscDS, f::Integer, g::Integer, g0::Ptr{Cvoid}, g1::Ptr{Cvoid}, g2::Ptr{Cvoid}, g3::Ptr{Cvoid})
    error("PetscDSSetDynamicJacobian: no generated method for these argument types")
end

@for_petsc function PetscDSSetDynamicJacobian(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt, g::$PetscInt, g0::Ptr{Cvoid}, g1::Ptr{Cvoid}, g2::Ptr{Cvoid}, g3::Ptr{Cvoid} )

    @chk ccall(
               (:PetscDSSetDynamicJacobian, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, $PetscInt, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
               ds, f, g, g0, g1, g2, g3,
              )


	return nothing
end 

"""
	PetscDSSetExactSolution(petsclib::PetscLibType,prob::PetscDS, f::PetscInt, sol::Ptr{Cvoid}, ctx::Ptr{Cvoid}) 
Set the pointwise exact solution function for a given test field

Not Collective

Input Parameters:
- `prob` - The `PetscDS`
- `f`    - The test field number
- `sol`  - solution function for the test fields, see `PetscPointExactSolutionFn`
- `ctx`  - solution context or `NULL`

Level: intermediate

-seealso: `PetscDS`, `PetscPointExactSolutionFn`, `PetscDSGetExactSolution()`

# External Links
$(_doc_external("DT/PetscDSSetExactSolution"))
"""
function PetscDSSetExactSolution(petsclib::PetscLibType, prob::PetscDS, f::Integer, sol::Ptr{Cvoid}, ctx::Ptr{Cvoid})
    error("PetscDSSetExactSolution: no generated method for these argument types")
end

@for_petsc function PetscDSSetExactSolution(petsclib::$UnionPetscLib, prob::PetscDS, f::$PetscInt, sol::Ptr{Cvoid}, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:PetscDSSetExactSolution, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{Cvoid}, Ptr{Cvoid}),
               prob, f, sol, ctx,
              )


	return nothing
end 

"""
	PetscDSSetExactSolutionTimeDerivative(petsclib::PetscLibType,prob::PetscDS, f::PetscInt, sol::Ptr{Cvoid}, ctx::Ptr{Cvoid}) 
Set the pointwise time derivative of the exact solution function for a given test field

Not Collective

Input Parameters:
- `prob` - The `PetscDS`
- `f`    - The test field number
- `sol`  - time derivative of the solution function for the test fields, see `PetscPointExactSolutionFn`
- `ctx`  - the solution context or `NULL`

Level: intermediate

-seealso: `PetscDS`, `PetscPointExactSolutionFn`, `PetscDSGetExactSolutionTimeDerivative()`, `PetscDSSetExactSolution()`

# External Links
$(_doc_external("DT/PetscDSSetExactSolutionTimeDerivative"))
"""
function PetscDSSetExactSolutionTimeDerivative(petsclib::PetscLibType, prob::PetscDS, f::Integer, sol::Ptr{Cvoid}, ctx::Ptr{Cvoid})
    error("PetscDSSetExactSolutionTimeDerivative: no generated method for these argument types")
end

@for_petsc function PetscDSSetExactSolutionTimeDerivative(petsclib::$UnionPetscLib, prob::PetscDS, f::$PetscInt, sol::Ptr{Cvoid}, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:PetscDSSetExactSolutionTimeDerivative, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{Cvoid}, Ptr{Cvoid}),
               prob, f, sol, ctx,
              )


	return nothing
end 

"""
	PetscDSSetForceQuad(petsclib::PetscLibType,ds::PetscDS, forceQuad::PetscBool) 
Set the flag to force matching quadratures among the field discretizations

Logically collective on ds

Input Parameters:
- `ds`        - The `PetscDS` object
- `forceQuad` - The flag

Level: intermediate

-seealso: `PetscDS`, `PetscDSGetForceQuad()`, `PetscDSGetDiscretization()`, `PetscDSGetNumFields()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSSetForceQuad"))
"""
function PetscDSSetForceQuad(petsclib::PetscLibType, ds::PetscDS, forceQuad::PetscBool)
    error("PetscDSSetForceQuad: no generated method for these argument types")
end

@for_petsc function PetscDSSetForceQuad(petsclib::$UnionPetscLib, ds::PetscDS, forceQuad::PetscBool )

    @chk ccall(
               (:PetscDSSetForceQuad, $petsc_library),
               PetscErrorCode,
               (PetscDS, PetscBool),
               ds, forceQuad,
              )


	return nothing
end 

"""
	PetscDSSetFromOptions(petsclib::PetscLibType,prob::PetscDS) 
sets parameters in a `PetscDS` from the options database

Collective

Input Parameter:
- `prob` - the `PetscDS` object to set options for

Options Database Keys:
- `-petscds_type type`            - Set the `PetscDS` type
- `-petscds_view`                 - View the `PetscDS`
- `-petscds_jac_pre (true|false)` - Turn formation of a separate Jacobian preconditioner on or off
- `-bc_NAME ids`                  - comma separated list of label ids for the boundary condition NAME
- `-bc_NAME_comp comps`           - comma separated list of field components to constrain for the boundary condition NAME

Level: intermediate

-seealso: `PetscDS`, `PetscDSView()`

# External Links
$(_doc_external("DT/PetscDSSetFromOptions"))
"""
function PetscDSSetFromOptions(petsclib::PetscLibType, prob::PetscDS)
    error("PetscDSSetFromOptions: no generated method for these argument types")
end

@for_petsc function PetscDSSetFromOptions(petsclib::$UnionPetscLib, prob::PetscDS )

    @chk ccall(
               (:PetscDSSetFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscDS,),
               prob,
              )


	return nothing
end 

"""
	PetscDSSetImplicit(petsclib::PetscLibType,prob::PetscDS, f::PetscInt, implicit::PetscBool) 
Set the flag for implicit solve for this field. This is just a guide for `TSARKIMEX`

Not Collective

Input Parameters:
- `prob`     - The `PetscDS` object
- `f`        - The field number
- `implicit` - The flag indicating what kind of solve to use for this field

Level: developer

-seealso: `TSARKIMEX`, `PetscDSGetImplicit()`, `PetscDSSetDiscretization()`, `PetscDSAddDiscretization()`, `PetscDSGetNumFields()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSSetImplicit"))
"""
function PetscDSSetImplicit(petsclib::PetscLibType, prob::PetscDS, f::Integer, implicit::PetscBool)
    error("PetscDSSetImplicit: no generated method for these argument types")
end

@for_petsc function PetscDSSetImplicit(petsclib::$UnionPetscLib, prob::PetscDS, f::$PetscInt, implicit::PetscBool )

    @chk ccall(
               (:PetscDSSetImplicit, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, PetscBool),
               prob, f, implicit,
              )


	return nothing
end 

"""
	PetscDSSetIntegrationParameters(petsclib::PetscLibType,ds::PetscDS, fieldI::PetscInt, fieldJ::PetscInt) 
Set the parameters for a particular integration

Not Collective

Input Parameters:
- `ds`     - The `PetscDS` object
- `fieldI` - The test field for a given point function, or `PETSC_DETERMINE`
- `fieldJ` - The basis field for a given point function, or `PETSC_DETERMINE`

Level: intermediate

-seealso: `PetscDS`, `PetscDSSetConstants()`, `PetscDSGetConstants()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSSetIntegrationParameters"))
"""
function PetscDSSetIntegrationParameters(petsclib::PetscLibType, ds::PetscDS, fieldI::Integer, fieldJ::Integer)
    error("PetscDSSetIntegrationParameters: no generated method for these argument types")
end

@for_petsc function PetscDSSetIntegrationParameters(petsclib::$UnionPetscLib, ds::PetscDS, fieldI::$PetscInt, fieldJ::$PetscInt )

    @chk ccall(
               (:PetscDSSetIntegrationParameters, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, $PetscInt),
               ds, fieldI, fieldJ,
              )


	return nothing
end 

"""
	PetscDSSetJacobian(petsclib::PetscLibType,ds::PetscDS, f::PetscInt, g::PetscInt, g0::Ptr{Cvoid}, g1::Ptr{Cvoid}, g2::Ptr{Cvoid}, g3::Ptr{Cvoid}) 
Set the pointwise Jacobian function for given test and basis fields

Not Collective

Input Parameters:
- `ds` - The `PetscDS`
- `f`  - The test field number
- `g`  - The field number
- `g0` - integrand for the test and basis function term, see `PetscPointJacFn`
- `g1` - integrand for the test function and basis function gradient term, see `PetscPointJacFn`
- `g2` - integrand for the test function gradient and basis function term, see `PetscPointJacFn`
- `g3` - integrand for the test function gradient and basis function gradient term, see `PetscPointJacFn`

Level: intermediate

-seealso: `PetscDS`, `PetscDSGetJacobian()`, `PetscPointJacFn`

# External Links
$(_doc_external("DT/PetscDSSetJacobian"))
"""
function PetscDSSetJacobian(petsclib::PetscLibType, ds::PetscDS, f::Integer, g::Integer, g0::Ptr{Cvoid}, g1::Ptr{Cvoid}, g2::Ptr{Cvoid}, g3::Ptr{Cvoid})
    error("PetscDSSetJacobian: no generated method for these argument types")
end

@for_petsc function PetscDSSetJacobian(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt, g::$PetscInt, g0::Ptr{Cvoid}, g1::Ptr{Cvoid}, g2::Ptr{Cvoid}, g3::Ptr{Cvoid} )

    @chk ccall(
               (:PetscDSSetJacobian, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, $PetscInt, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
               ds, f, g, g0, g1, g2, g3,
              )


	return nothing
end 

"""
	PetscDSSetJacobianPreconditioner(petsclib::PetscLibType,ds::PetscDS, f::PetscInt, g::PetscInt, g0::Ptr{Cvoid}, g1::Ptr{Cvoid}, g2::Ptr{Cvoid}, g3::Ptr{Cvoid}) 
Set the pointwise Jacobian function for given test and basis fields that constructs the matrix used
to compute the preconditioner. If this is missing, the system matrix is used to build the preconditioner.

Not Collective

Input Parameters:
- `ds` - The `PetscDS`
- `f`  - The test field number
- `g`  - The field number
- `g0` - integrand for the test and basis function term, see `PetscPointJacFn`
- `g1` - integrand for the test function and basis function gradient term, see `PetscPointJacFn`
- `g2` - integrand for the test function gradient and basis function term, see `PetscPointJacFn`
- `g3` - integrand for the test function gradient and basis function gradient term, see `PetscPointJacFn`

Level: intermediate

-seealso: `PetscDS`, `PetscDSGetJacobianPreconditioner()`, `PetscDSSetJacobian()`, `PetscPointJacFn`

# External Links
$(_doc_external("DT/PetscDSSetJacobianPreconditioner"))
"""
function PetscDSSetJacobianPreconditioner(petsclib::PetscLibType, ds::PetscDS, f::Integer, g::Integer, g0::Ptr{Cvoid}, g1::Ptr{Cvoid}, g2::Ptr{Cvoid}, g3::Ptr{Cvoid})
    error("PetscDSSetJacobianPreconditioner: no generated method for these argument types")
end

@for_petsc function PetscDSSetJacobianPreconditioner(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt, g::$PetscInt, g0::Ptr{Cvoid}, g1::Ptr{Cvoid}, g2::Ptr{Cvoid}, g3::Ptr{Cvoid} )

    @chk ccall(
               (:PetscDSSetJacobianPreconditioner, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, $PetscInt, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
               ds, f, g, g0, g1, g2, g3,
              )


	return nothing
end 

"""
	PetscDSSetJetDegree(petsclib::PetscLibType,ds::PetscDS, f::PetscInt, k::PetscInt) 
Set the highest derivative for this field equation, or the k

Not Collective

Input Parameters:
- `ds` - The `PetscDS` object
- `f`  - The field number
- `k`  - The highest derivative we need to tabulate

Level: developer

-seealso: `PetscDS`, `PetscDSGetJetDegree()`, `PetscDSSetDiscretization()`, `PetscDSAddDiscretization()`, `PetscDSGetNumFields()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSSetJetDegree"))
"""
function PetscDSSetJetDegree(petsclib::PetscLibType, ds::PetscDS, f::Integer, k::Integer)
    error("PetscDSSetJetDegree: no generated method for these argument types")
end

@for_petsc function PetscDSSetJetDegree(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt, k::$PetscInt )

    @chk ccall(
               (:PetscDSSetJetDegree, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, $PetscInt),
               ds, f, k,
              )


	return nothing
end 

"""
	PetscDSSetLowerBound(petsclib::PetscLibType,ds::PetscDS, f::PetscInt, lb::Ptr{Cvoid}, ctx::Ptr{Cvoid}) 
Set the pointwise lower bound function for a given field

Not Collective

Input Parameters:
- `ds`  - The `PetscDS`
- `f`   - The field number
- `lb`  - lower bound function for the test fields, see `PetscPointBoundFn`
- `ctx` - lower bound context or `NULL` which will be passed to `lb`

Level: intermediate

-seealso: `PetscDS`, `PetscPointBoundFn`, `PetscDSGetLowerBound()`, `PetscDSGetUpperBound()`, `PetscDSGetExactSolution()`

# External Links
$(_doc_external("DT/PetscDSSetLowerBound"))
"""
function PetscDSSetLowerBound(petsclib::PetscLibType, ds::PetscDS, f::Integer, lb::Ptr{Cvoid}, ctx::Ptr{Cvoid})
    error("PetscDSSetLowerBound: no generated method for these argument types")
end

@for_petsc function PetscDSSetLowerBound(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt, lb::Ptr{Cvoid}, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:PetscDSSetLowerBound, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{Cvoid}, Ptr{Cvoid}),
               ds, f, lb, ctx,
              )


	return nothing
end 

"""
	PetscDSSetObjective(petsclib::PetscLibType,ds::PetscDS, f::PetscInt, obj::Ptr{Cvoid}) 
Set the pointwise objective function for a given test field

Not Collective

Input Parameters:
- `ds`  - The `PetscDS`
- `f`   - The test field number
- `obj` - integrand for the test function term, see `PetscPointFn`

Level: intermediate

-seealso: `PetscPointFn`, `PetscDS`, `PetscDSGetObjective()`, `PetscDSSetResidual()`

# External Links
$(_doc_external("DT/PetscDSSetObjective"))
"""
function PetscDSSetObjective(petsclib::PetscLibType, ds::PetscDS, f::Integer, obj::Ptr{Cvoid})
    error("PetscDSSetObjective: no generated method for these argument types")
end

@for_petsc function PetscDSSetObjective(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt, obj::Ptr{Cvoid} )

    @chk ccall(
               (:PetscDSSetObjective, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{Cvoid}),
               ds, f, obj,
              )


	return nothing
end 

"""
	PetscDSSetRHSResidual(petsclib::PetscLibType,ds::PetscDS, f::PetscInt, f0::Ptr{Cvoid}, f1::Ptr{Cvoid}) 
Set the pointwise residual function for explicit timestepping for a given test field

Not Collective

Input Parameters:
- `ds` - The `PetscDS`
- `f`  - The test field number
- `f0` - integrand for the test function term, see `PetscPointFn`
- `f1` - integrand for the test function gradient term, see `PetscPointFn`

Level: intermediate

-seealso: `PetscDS`, `PetscDSGetResidual()`

# External Links
$(_doc_external("DT/PetscDSSetRHSResidual"))
"""
function PetscDSSetRHSResidual(petsclib::PetscLibType, ds::PetscDS, f::Integer, f0::Ptr{Cvoid}, f1::Ptr{Cvoid})
    error("PetscDSSetRHSResidual: no generated method for these argument types")
end

@for_petsc function PetscDSSetRHSResidual(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt, f0::Ptr{Cvoid}, f1::Ptr{Cvoid} )

    @chk ccall(
               (:PetscDSSetRHSResidual, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{Cvoid}, Ptr{Cvoid}),
               ds, f, f0, f1,
              )


	return nothing
end 

"""
	PetscDSSetResidual(petsclib::PetscLibType,ds::PetscDS, f::PetscInt, f0::Ptr{Cvoid}, f1::Ptr{Cvoid}) 
Set the pointwise residual function for a given test field

Not Collective

Input Parameters:
- `ds` - The `PetscDS`
- `f`  - The test field number
- `f0` - integrand for the test function term, see `PetscPointFn`
- `f1` - integrand for the test function gradient term, see `PetscPointFn`

Level: intermediate

-seealso: `PetscPointFn`, `PetscDS`, `PetscDSGetResidual()`

# External Links
$(_doc_external("DT/PetscDSSetResidual"))
"""
function PetscDSSetResidual(petsclib::PetscLibType, ds::PetscDS, f::Integer, f0::Ptr{Cvoid}, f1::Ptr{Cvoid})
    error("PetscDSSetResidual: no generated method for these argument types")
end

@for_petsc function PetscDSSetResidual(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt, f0::Ptr{Cvoid}, f1::Ptr{Cvoid} )

    @chk ccall(
               (:PetscDSSetResidual, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{Cvoid}, Ptr{Cvoid}),
               ds, f, f0, f1,
              )


	return nothing
end 

"""
	PetscDSSetRiemannSolver(petsclib::PetscLibType,ds::PetscDS, f::PetscInt, r::Ptr{Cvoid}) 
Sets the Riemann solver for the given field

Not Collective

Input Parameters:
- `ds` - The `PetscDS` object
- `f`  - The field number
- `r`  - Riemann solver, see `PetscRiemannFn`

Level: intermediate

-seealso: `PetscDS`, `PetscRiemannFn`, `PetscDSGetRiemannSolver()`

# External Links
$(_doc_external("DT/PetscDSSetRiemannSolver"))
"""
function PetscDSSetRiemannSolver(petsclib::PetscLibType, ds::PetscDS, f::Integer, r::Ptr{Cvoid})
    error("PetscDSSetRiemannSolver: no generated method for these argument types")
end

@for_petsc function PetscDSSetRiemannSolver(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt, r::Ptr{Cvoid} )

    @chk ccall(
               (:PetscDSSetRiemannSolver, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{Cvoid}),
               ds, f, r,
              )


	return nothing
end 

"""
	PetscDSSetType(petsclib::PetscLibType,prob::PetscDS, name::PetscDSType) 
Builds a particular `PetscDS`

Collective; No Fortran Support

Input Parameters:
- `prob` - The `PetscDS` object
- `name` - The `PetscDSType`

Options Database Key:
- `-petscds_type type` - Sets the PetscDS type; use -help for a list of available types

Level: intermediate

-seealso: `PetscDSType`, `PetscDS`, `PetscDSGetType()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSSetType"))
"""
function PetscDSSetType(petsclib::PetscLibType, prob::PetscDS, name::PetscDSType)
    error("PetscDSSetType: no generated method for these argument types")
end

@for_petsc function PetscDSSetType(petsclib::$UnionPetscLib, prob::PetscDS, name::PetscDSType )

    @chk ccall(
               (:PetscDSSetType, $petsc_library),
               PetscErrorCode,
               (PetscDS, PetscDSType),
               prob, name,
              )


	return nothing
end 

"""
	PetscDSSetUp(petsclib::PetscLibType,prob::PetscDS) 
Construct data structures for the `PetscDS`

Collective

Input Parameter:
- `prob` - the `PetscDS` object to setup

Level: developer

-seealso: `PetscDS`, `PetscDSView()`, `PetscDSDestroy()`

# External Links
$(_doc_external("DT/PetscDSSetUp"))
"""
function PetscDSSetUp(petsclib::PetscLibType, prob::PetscDS)
    error("PetscDSSetUp: no generated method for these argument types")
end

@for_petsc function PetscDSSetUp(petsclib::$UnionPetscLib, prob::PetscDS )

    @chk ccall(
               (:PetscDSSetUp, $petsc_library),
               PetscErrorCode,
               (PetscDS,),
               prob,
              )


	return nothing
end 

"""
	PetscDSSetUpdate(petsclib::PetscLibType,ds::PetscDS, f::PetscInt, update::Ptr{Cvoid}) 
Set the pointwise update function for a given field

Not Collective

Input Parameters:
- `ds`     - The `PetscDS`
- `f`      - The field number
- `update` - update function, see `PetscPointFn`

Level: intermediate

-seealso: `PetscDS`, `PetscPointFn`, `PetscDSGetResidual()`

# External Links
$(_doc_external("DT/PetscDSSetUpdate"))
"""
function PetscDSSetUpdate(petsclib::PetscLibType, ds::PetscDS, f::Integer, update::Ptr{Cvoid})
    error("PetscDSSetUpdate: no generated method for these argument types")
end

@for_petsc function PetscDSSetUpdate(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt, update::Ptr{Cvoid} )

    @chk ccall(
               (:PetscDSSetUpdate, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{Cvoid}),
               ds, f, update,
              )


	return nothing
end 

"""
	PetscDSSetUpperBound(petsclib::PetscLibType,ds::PetscDS, f::PetscInt, ub::Ptr{Cvoid}, ctx::Ptr{Cvoid}) 
Set the pointwise upper bound function for a given field

Not Collective

Input Parameters:
- `ds`  - The `PetscDS`
- `f`   - The field number
- `ub`  - upper bound function for the test fields, see `PetscPointBoundFn`
- `ctx` - context or `NULL` that will be passed to `ub`

Level: intermediate

-seealso: `PetscDS`, `PetscPointBoundFn`, `PetscDSGetUpperBound()`, `PetscDSGetLowerBound()`, `PetscDSGetExactSolution()`

# External Links
$(_doc_external("DT/PetscDSSetUpperBound"))
"""
function PetscDSSetUpperBound(petsclib::PetscLibType, ds::PetscDS, f::Integer, ub::Ptr{Cvoid}, ctx::Ptr{Cvoid})
    error("PetscDSSetUpperBound: no generated method for these argument types")
end

@for_petsc function PetscDSSetUpperBound(petsclib::$UnionPetscLib, ds::PetscDS, f::$PetscInt, ub::Ptr{Cvoid}, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:PetscDSSetUpperBound, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, Ptr{Cvoid}, Ptr{Cvoid}),
               ds, f, ub, ctx,
              )


	return nothing
end 

"""
	PetscDSSetWeakForm(petsclib::PetscLibType,ds::PetscDS, wf::PetscWeakForm) 
Sets the weak form object to be used by the `PetscDS`

Not Collective

Input Parameters:
- `ds` - The `PetscDS` object
- `wf` - The weak form object

Level: beginner

-seealso: `PetscWeakForm`, `PetscDSGetWeakForm()`, `PetscDSGetNumFields()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSSetWeakForm"))
"""
function PetscDSSetWeakForm(petsclib::PetscLibType, ds::PetscDS, wf::PetscWeakForm)
    error("PetscDSSetWeakForm: no generated method for these argument types")
end

@for_petsc function PetscDSSetWeakForm(petsclib::$UnionPetscLib, ds::PetscDS, wf::PetscWeakForm )

    @chk ccall(
               (:PetscDSSetWeakForm, $petsc_library),
               PetscErrorCode,
               (PetscDS, PetscWeakForm),
               ds, wf,
              )


	return nothing
end 

"""
	PetscDSUpdateBoundary(petsclib::PetscLibType,ds::PetscDS, bd::PetscInt, type::DMBoundaryConditionType, name::String, label::DMLabel, Nv::PetscInt, values::Vector{PetscInt}, field::PetscInt, Nc::PetscInt, comps::Vector{PetscInt}, bcFunc::Ptr{Cvoid}, bcFunc_t::Ptr{Cvoid}, ctx::Ptr{Cvoid}) 
Change a boundary condition for the model.

Input Parameters:
- `ds`       - The `PetscDS` object
- `bd`       - The boundary condition number
- `type`     - The type of condition, e.g. `DM_BC_ESSENTIAL`/`DM_BC_ESSENTIAL_FIELD` (Dirichlet), or `DM_BC_NATURAL` (Neumann)
- `name`     - The boundary condition name
- `label`    - The label defining constrained points
- `Nv`       - The number of `DMLabel` ids for constrained points
- `values`   - An array of ids for constrained points
- `field`    - The field to constrain
- `Nc`       - The number of constrained field components
- `comps`    - An array of constrained component numbers
- `bcFunc`   - A pointwise function giving boundary values
- `bcFunc_t` - A pointwise function giving the time derivative of the boundary values, or `NULL`
- `ctx`      - An optional application context for `bcFunc`

Level: developer

-seealso: `PetscDS`, `PetscWeakForm`, `DMBoundaryConditionType`, `PetscDSAddBoundary()`, `PetscDSGetBoundary()`, `PetscDSGetNumBoundary()`, `DMLabel`

# External Links
$(_doc_external("DT/PetscDSUpdateBoundary"))
"""
function PetscDSUpdateBoundary(petsclib::PetscLibType, ds::PetscDS, bd::Integer, type::DMBoundaryConditionType, name::String, label::DMLabel, Nv::Integer, values::AbstractVector{<:Number}, field::Integer, Nc::Integer, comps::AbstractVector{<:Number}, bcFunc::Ptr{Cvoid}, bcFunc_t::Ptr{Cvoid}, ctx::Ptr{Cvoid})
    error("PetscDSUpdateBoundary: no generated method for these argument types")
end

@for_petsc function PetscDSUpdateBoundary(petsclib::$UnionPetscLib, ds::PetscDS, bd::$PetscInt, type::DMBoundaryConditionType, name::String, label::DMLabel, Nv::$PetscInt, values::Vector{$PetscInt}, field::$PetscInt, Nc::$PetscInt, comps::Vector{$PetscInt}, bcFunc::Ptr{Cvoid}, bcFunc_t::Ptr{Cvoid}, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:PetscDSUpdateBoundary, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt, DMBoundaryConditionType, Ptr{Cchar}, DMLabel, $PetscInt, Ptr{$PetscInt}, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
               ds, bd, type, name, label, Nv, values, field, Nc, comps, bcFunc, bcFunc_t, ctx,
              )


	return nothing
end 

"""
	PetscDSUpdateBoundaryLabels(petsclib::PetscLibType,ds::PetscDS, dm::AbstractPetscDM) 
Update `DMLabel` in each boundary condition using the label name and the input `DM`

Not Collective

Input Parameters:
- `ds` - The source `PetscDS` object
- `dm` - The `DM` holding labels

Level: intermediate

-seealso: `PetscDS`, `DMBoundary`, `DM`, `PetscDSCopyBoundary()`, `PetscDSCreate()`, `DMGetLabel()`

# External Links
$(_doc_external("DT/PetscDSUpdateBoundaryLabels"))
"""
function PetscDSUpdateBoundaryLabels(petsclib::PetscLibType, ds::PetscDS, dm::AbstractPetscDM)
    error("PetscDSUpdateBoundaryLabels: no generated method for these argument types")
end

@for_petsc function PetscDSUpdateBoundaryLabels(petsclib::$UnionPetscLib, ds::PetscDS, dm::AbstractPetscDM )

    @chk ccall(
               (:PetscDSUpdateBoundaryLabels, $petsc_library),
               PetscErrorCode,
               (PetscDS, CDM),
               ds, dm,
              )


	return nothing
end 

"""
	PetscDSUseJacobianPreconditioner(petsclib::PetscLibType,prob::PetscDS, useJacPre::PetscBool) 
Set whether to construct a Jacobian preconditioner

Not Collective

Input Parameters:
- `prob`      - The `PetscDS`
- `useJacPre` - flag that enables construction of a Jacobian preconditioner

Level: intermediate

-seealso: `PetscDS`, `PetscDSGetJacobianPreconditioner()`, `PetscDSSetJacobianPreconditioner()`, `PetscDSGetJacobian()`

# External Links
$(_doc_external("DT/PetscDSUseJacobianPreconditioner"))
"""
function PetscDSUseJacobianPreconditioner(petsclib::PetscLibType, prob::PetscDS, useJacPre::PetscBool)
    error("PetscDSUseJacobianPreconditioner: no generated method for these argument types")
end

@for_petsc function PetscDSUseJacobianPreconditioner(petsclib::$UnionPetscLib, prob::PetscDS, useJacPre::PetscBool )

    @chk ccall(
               (:PetscDSUseJacobianPreconditioner, $petsc_library),
               PetscErrorCode,
               (PetscDS, PetscBool),
               prob, useJacPre,
              )


	return nothing
end 

"""
	PetscDSView(petsclib::PetscLibType,prob::PetscDS, v::PetscViewer) 
Views a `PetscDS`

Collective

Input Parameters:
- `prob` - the `PetscDS` object to view
- `v`    - the viewer

Level: developer

-seealso: `PetscDSType`, `PetscDS`, `PetscViewer`, `PetscDSDestroy()`, `PetscDSViewFromOptions()`

# External Links
$(_doc_external("DT/PetscDSView"))
"""
function PetscDSView(petsclib::PetscLibType, prob::PetscDS, v::PetscViewer)
    error("PetscDSView: no generated method for these argument types")
end

@for_petsc function PetscDSView(petsclib::$UnionPetscLib, prob::PetscDS, v::PetscViewer )

    @chk ccall(
               (:PetscDSView, $petsc_library),
               PetscErrorCode,
               (PetscDS, PetscViewer),
               prob, v,
              )


	return nothing
end 

"""
	PetscDSViewFromOptions(petsclib::PetscLibType,A::PetscDS, obj, name::String) 
View a `PetscDS` based on values in the options database

Collective

Input Parameters:
- `A`    - the `PetscDS` object
- `obj`  - Optional object that provides the options prefix used in the search of the options database
- `name` - command line option

Options Database Key:
- `-name [viewertype][:...]` - option name and values. See `PetscObjectViewFromOptions()` for the possible arguments

Level: intermediate

-seealso: `PetscDSType`, `PetscDS`, `PetscDSView()`, `PetscObjectViewFromOptions()`, `PetscDSCreate()`

# External Links
$(_doc_external("DT/PetscDSViewFromOptions"))
"""
function PetscDSViewFromOptions(petsclib::PetscLibType, A::PetscDS, obj, name::String)
    error("PetscDSViewFromOptions: no generated method for these argument types")
end

@for_petsc function PetscDSViewFromOptions(petsclib::$UnionPetscLib, A::PetscDS, obj, name::String )

    @chk ccall(
               (:PetscDSViewFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscDS, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

