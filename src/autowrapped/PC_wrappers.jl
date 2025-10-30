# autodefined type arguments for class ------
mutable struct _n_PCRiCchardsonConvergedReason end
const PCRiCchardsonConvergedReason = Ptr{_n_PCRiCchardsonConvergedReason}

mutable struct PCModifySubMatricesFn end

mutable struct PCMGCoarseSpaceConstructorFn end

mutable struct PCShellPSolveFn end

# -------------------------------------------------------
"""
	PCSetType(petsclib::PetscLibType,pc::PC, type::PCType) 
Builds `PC` for a particular preconditioner type

Collective

Input Parameters:
- `pc`   - the preconditioner context
- `type` - a known method, see `PCType` for possible values

Options Database Key:
- `-pc_type <type>` - Sets `PC` type

-seealso: [](ch_ksp), `KSPSetType()`, `PCType`, `PCRegister()`, `PCCreate()`, `KSPGetPC()`

# External Links
$(_doc_external("Ksp/PCSetType"))
"""
function PCSetType(petsclib::PetscLibType, pc::PC, type::PCType) end

@for_petsc function PCSetType(petsclib::$UnionPetscLib, pc::PC, type::PCType )

    @chk ccall(
               (:PCSetType, $petsc_library),
               PetscErrorCode,
               (PC, PCType),
               pc, type,
              )


	return nothing
end 

"""
	type::PCType = PCGetType(petsclib::PetscLibType,pc::PC) 
Gets the `PCType` (as a string) from the `PC`
context.

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `type` - name of preconditioner method

Level: intermediate

-seealso: [](ch_ksp), `PC`, `PCType`, `PCSetType()`

# External Links
$(_doc_external("Ksp/PCGetType"))
"""
function PCGetType(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCGetType(petsclib::$UnionPetscLib, pc::PC )
	type_ = Ref{PCType}()

    @chk ccall(
               (:PCGetType, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PCType}),
               pc, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	PCSetFromOptions(petsclib::PetscLibType,pc::PC) 
Sets `PC` options from the options database.

Collective

Input Parameter:
- `pc` - the preconditioner context

Options Database Key:
- `-pc_type` - name of type, for example `bjacobi`

Level: advanced

-seealso: [](ch_ksp), `PC`, `PCSetType()`, `PCType`, `KSPSetFromOptions()`

# External Links
$(_doc_external("Ksp/PCSetFromOptions"))
"""
function PCSetFromOptions(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCSetFromOptions(petsclib::$UnionPetscLib, pc::PC )

    @chk ccall(
               (:PCSetFromOptions, $petsc_library),
               PetscErrorCode,
               (PC,),
               pc,
              )


	return nothing
end 

"""
	PCSetDM(petsclib::PetscLibType,pc::PC, dm::PetscDM) 
Sets the `DM` that may be used by some preconditioners

Logically Collective

Input Parameters:
- `pc` - the preconditioner context
- `dm` - the `DM`, can be `NULL` to remove any current `DM`

Level: intermediate

-seealso: [](ch_ksp), `PC`, `DM`, `PCGetDM()`, `KSPSetDM()`, `KSPGetDM()`, `SNESSetDM()`, `TSSetDM()`

# External Links
$(_doc_external("Ksp/PCSetDM"))
"""
function PCSetDM(petsclib::PetscLibType, pc::PC, dm::PetscDM) end

@for_petsc function PCSetDM(petsclib::$UnionPetscLib, pc::PC, dm::PetscDM )

    @chk ccall(
               (:PCSetDM, $petsc_library),
               PetscErrorCode,
               (PC, CDM),
               pc, dm,
              )


	return nothing
end 

"""
	PCGetDM(petsclib::PetscLibType,pc::PC, dm::PetscDM) 
Gets the `DM` that may be used by some preconditioners

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `dm` - the `DM`

Level: intermediate

-seealso: [](ch_ksp), `PC`, `DM`, `PCSetDM()`, `KSPSetDM()`, `KSPGetDM()`

# External Links
$(_doc_external("Ksp/PCGetDM"))
"""
function PCGetDM(petsclib::PetscLibType, pc::PC, dm::PetscDM) end

@for_petsc function PCGetDM(petsclib::$UnionPetscLib, pc::PC, dm::PetscDM )
	dm_ = Ref(dm.ptr)

    @chk ccall(
               (:PCGetDM, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{CDM}),
               pc, dm_,
              )

	dm.ptr = C_NULL

	return nothing
end 

"""
	PCSetApplicationContext(petsclib::PetscLibType,pc::PC, ctx::Cvoid) 
Sets the optional user

Logically Collective

Input Parameters:
- `pc`  - the `PC` context
- `ctx` - optional user context

Level: advanced

-seealso: [](ch_ksp), `PC`, `PCGetApplicationContext()`, `KSPSetApplicationContext()`, `KSPGetApplicationContext()`, `PetscObjectCompose()`

# External Links
$(_doc_external("Ksp/PCSetApplicationContext"))
"""
function PCSetApplicationContext(petsclib::PetscLibType, pc::PC, ctx::Cvoid) end

@for_petsc function PCSetApplicationContext(petsclib::$UnionPetscLib, pc::PC, ctx::Cvoid )

    @chk ccall(
               (:PCSetApplicationContext, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{Cvoid}),
               pc, ctx,
              )


	return nothing
end 

"""
	PCGetApplicationContext(petsclib::PetscLibType,pc::PC, ctx::PeCtx) 
Gets the user

Not Collective

Input Parameter:
- `pc` - `PC` context

Output Parameter:
- `ctx` - user context

Level: intermediate

-seealso: [](ch_ksp), `PC`, `PCSetApplicationContext()`, `KSPSetApplicationContext()`, `KSPGetApplicationContext()`

# External Links
$(_doc_external("Ksp/PCGetApplicationContext"))
"""
function PCGetApplicationContext(petsclib::PetscLibType, pc::PC, ctx::PeCtx) end

@for_petsc function PCGetApplicationContext(petsclib::$UnionPetscLib, pc::PC, ctx::PeCtx )

    @chk ccall(
               (:PCGetApplicationContext, $petsc_library),
               PetscErrorCode,
               (PC, PeCtx),
               pc, ctx,
              )


	return nothing
end 

"""
	PCReset(petsclib::PetscLibType,pc::PC) 
Resets a `PC` context to the state it was in before `PCSetUp()` was called, and removes any allocated `Vec` and `Mat` from its data structure

Collective

Input Parameter:
- `pc` - the `PC` preconditioner context

Level: developer

-seealso: [](ch_ksp), `PC`, `PCCreate()`, `PCSetUp()`

# External Links
$(_doc_external("Ksp/PCReset"))
"""
function PCReset(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCReset(petsclib::$UnionPetscLib, pc::PC )

    @chk ccall(
               (:PCReset, $petsc_library),
               PetscErrorCode,
               (PC,),
               pc,
              )


	return nothing
end 

"""
	PCDestroy(petsclib::PetscLibType,pc::PC) 
Destroys `PC` context that was created with `PCCreate()`.

Collective

Input Parameter:
- `pc` - the `PC` preconditioner context

Level: developer

-seealso: [](ch_ksp), `PC`, `PCCreate()`, `PCSetUp()`

# External Links
$(_doc_external("Ksp/PCDestroy"))
"""
function PCDestroy(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCDestroy(petsclib::$UnionPetscLib, pc::PC )

    @chk ccall(
               (:PCDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PC},),
               pc,
              )


	return nothing
end 

"""
	flag::PetscBool = PCGetDiagonalScale(petsclib::PetscLibType,pc::PC) 
Indicates if the preconditioner applies an additional left and right
scaling as needed by certain time-stepping codes.

Logically Collective

Input Parameter:
- `pc` - the `PC` preconditioner context

Output Parameter:
- `flag` - `PETSC_TRUE` if it applies the scaling

Level: developer

-seealso: [](ch_ksp), `PC`, `PCCreate()`, `PCSetUp()`, `PCDiagonalScaleLeft()`, `PCDiagonalScaleRight()`, `PCSetDiagonalScale()`

# External Links
$(_doc_external("Ksp/PCGetDiagonalScale"))
"""
function PCGetDiagonalScale(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCGetDiagonalScale(petsclib::$UnionPetscLib, pc::PC )
	flag_ = Ref{PetscBool}()

    @chk ccall(
               (:PCGetDiagonalScale, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PetscBool}),
               pc, flag_,
              )

	flag = flag_[]

	return flag
end 

"""
	PCSetDiagonalScale(petsclib::PetscLibType,pc::PC, s::PetscVec) 
Indicates the left scaling to use to apply an additional left and right
scaling as needed by certain time-stepping codes.

Logically Collective

Input Parameters:
- `pc` - the `PC` preconditioner context
- `s`  - scaling vector

Level: intermediate

-seealso: [](ch_ksp), `PCCreate()`, `PCSetUp()`, `PCDiagonalScaleLeft()`, `PCDiagonalScaleRight()`, `PCGetDiagonalScale()`

# External Links
$(_doc_external("Ksp/PCSetDiagonalScale"))
"""
function PCSetDiagonalScale(petsclib::PetscLibType, pc::PC, s::PetscVec) end

@for_petsc function PCSetDiagonalScale(petsclib::$UnionPetscLib, pc::PC, s::PetscVec )

    @chk ccall(
               (:PCSetDiagonalScale, $petsc_library),
               PetscErrorCode,
               (PC, CVec),
               pc, s,
              )


	return nothing
end 

"""
	PCDiagonalScaleLeft(petsclib::PetscLibType,pc::PC, in::PetscVec, out::PetscVec) 
Scales a vector by the left scaling as needed by certain time

Logically Collective

Input Parameters:
- `pc`  - the `PC` preconditioner context
- `in`  - input vector
- `out` - scaled vector (maybe the same as in)

Level: intermediate

-seealso: [](ch_ksp), `PCCreate()`, `PCSetUp()`, `PCSetDiagonalScale()`, `PCDiagonalScaleRight()`, `MatDiagonalScale()`

# External Links
$(_doc_external("Ksp/PCDiagonalScaleLeft"))
"""
function PCDiagonalScaleLeft(petsclib::PetscLibType, pc::PC, in::PetscVec, out::PetscVec) end

@for_petsc function PCDiagonalScaleLeft(petsclib::$UnionPetscLib, pc::PC, in::PetscVec, out::PetscVec )

    @chk ccall(
               (:PCDiagonalScaleLeft, $petsc_library),
               PetscErrorCode,
               (PC, CVec, CVec),
               pc, in, out,
              )


	return nothing
end 

"""
	PCDiagonalScaleRight(petsclib::PetscLibType,pc::PC, in::PetscVec, out::PetscVec) 
Scales a vector by the right scaling as needed by certain time

Logically Collective

Input Parameters:
- `pc`  - the `PC` preconditioner context
- `in`  - input vector
- `out` - scaled vector (maybe the same as in)

Level: intermediate

-seealso: [](ch_ksp), `PCCreate()`, `PCSetUp()`, `PCDiagonalScaleLeft()`, `PCSetDiagonalScale()`, `MatDiagonalScale()`

# External Links
$(_doc_external("Ksp/PCDiagonalScaleRight"))
"""
function PCDiagonalScaleRight(petsclib::PetscLibType, pc::PC, in::PetscVec, out::PetscVec) end

@for_petsc function PCDiagonalScaleRight(petsclib::$UnionPetscLib, pc::PC, in::PetscVec, out::PetscVec )

    @chk ccall(
               (:PCDiagonalScaleRight, $petsc_library),
               PetscErrorCode,
               (PC, CVec, CVec),
               pc, in, out,
              )


	return nothing
end 

"""
	PCSetUseAmat(petsclib::PetscLibType,pc::PC, flg::PetscBool) 
Sets a flag to indicate that when the preconditioner needs to apply (part of) the
operator during the preconditioning process it applies the Amat provided to `TSSetRHSJacobian()`,
`TSSetIJacobian()`, `SNESSetJacobian()`, `KSPSetOperators()` or `PCSetOperators()` not the Pmat.

Logically Collective

Input Parameters:
- `pc`  - the `PC` preconditioner context
- `flg` - `PETSC_TRUE` to use the Amat, `PETSC_FALSE` to use the Pmat (default is false)

Options Database Key:
- `-pc_use_amat <true,false>` - use the amat argument to `KSPSetOperators()` or `PCSetOperators()` to apply the operator

Level: intermediate

-seealso: [](ch_ksp), `PC`, `PCGetUseAmat()`, `PCBJACOBI`, `PCMG`, `PCFIELDSPLIT`, `PCCOMPOSITE`,
`KSPSetOperators()`, `PCSetOperators()`

# External Links
$(_doc_external("Ksp/PCSetUseAmat"))
"""
function PCSetUseAmat(petsclib::PetscLibType, pc::PC, flg::PetscBool) end

@for_petsc function PCSetUseAmat(petsclib::$UnionPetscLib, pc::PC, flg::PetscBool )

    @chk ccall(
               (:PCSetUseAmat, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, flg,
              )


	return nothing
end 

"""
	PCSetErrorIfFailure(petsclib::PetscLibType,pc::PC, flg::PetscBool) 
Causes `PC` to generate an error if a floating point exception, for example a zero pivot, is detected.

Logically Collective

Input Parameters:
- `pc`  - iterative context obtained from `PCCreate()`
- `flg` - `PETSC_TRUE` indicates you want the error generated

Level: advanced

-seealso: [](ch_ksp), `PC`, `KSPSetErrorIfNotConverged()`, `PCGetInitialGuessNonzero()`, `PCSetInitialGuessKnoll()`, `PCGetInitialGuessKnoll()`

# External Links
$(_doc_external("Ksp/PCSetErrorIfFailure"))
"""
function PCSetErrorIfFailure(petsclib::PetscLibType, pc::PC, flg::PetscBool) end

@for_petsc function PCSetErrorIfFailure(petsclib::$UnionPetscLib, pc::PC, flg::PetscBool )

    @chk ccall(
               (:PCSetErrorIfFailure, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, flg,
              )


	return nothing
end 

"""
	flg::PetscBool = PCGetUseAmat(petsclib::PetscLibType,pc::PC) 
Gets a flag to indicate that when the preconditioner needs to apply (part of) the
operator during the preconditioning process it applies the Amat provided to `TSSetRHSJacobian()`,
`TSSetIJacobian()`, `SNESSetJacobian()`, `KSPSetOperators()` or `PCSetOperators()` not the Pmat.

Logically Collective

Input Parameter:
- `pc` - the `PC` preconditioner context

Output Parameter:
- `flg` - `PETSC_TRUE` to use the Amat, `PETSC_FALSE` to use the Pmat (default is false)

Level: intermediate

-seealso: [](ch_ksp), `PC`, `PCSetUseAmat()`, `PCBJACOBI`, `PCMG`, `PCFIELDSPLIT`, `PCCOMPOSITE`

# External Links
$(_doc_external("Ksp/PCGetUseAmat"))
"""
function PCGetUseAmat(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCGetUseAmat(petsclib::$UnionPetscLib, pc::PC )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PCGetUseAmat, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PetscBool}),
               pc, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	PCSetKSPNestLevel(petsclib::PetscLibType,pc::PC, level::PetscInt) 
sets the amount of nesting the `KSP` that contains this `PC` has

Collective

Input Parameters:
- `pc`    - the `PC`
- `level` - the nest level

Level: developer

-seealso: [](ch_ksp), `KSPSetUp()`, `KSPSolve()`, `KSPDestroy()`, `KSP`, `KSPGMRES`, `KSPType`, `KSPGetNestLevel()`, `PCGetKSPNestLevel()`, `KSPSetNestLevel()`

# External Links
$(_doc_external("Ksp/PCSetKSPNestLevel"))
"""
function PCSetKSPNestLevel(petsclib::PetscLibType, pc::PC, level::PetscInt) end

@for_petsc function PCSetKSPNestLevel(petsclib::$UnionPetscLib, pc::PC, level::$PetscInt )

    @chk ccall(
               (:PCSetKSPNestLevel, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, level,
              )


	return nothing
end 

"""
	level::PetscInt = PCGetKSPNestLevel(petsclib::PetscLibType,pc::PC) 
gets the amount of nesting the `KSP` that contains this `PC` has

Not Collective

Input Parameter:
- `pc` - the `PC`

Output Parameter:
- `level` - the nest level

Level: developer

-seealso: [](ch_ksp), `KSPSetUp()`, `KSPSolve()`, `KSPDestroy()`, `KSP`, `KSPGMRES`, `KSPType`, `KSPSetNestLevel()`, `PCSetKSPNestLevel()`, `KSPGetNestLevel()`

# External Links
$(_doc_external("Ksp/PCGetKSPNestLevel"))
"""
function PCGetKSPNestLevel(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCGetKSPNestLevel(petsclib::$UnionPetscLib, pc::PC )
	level_ = Ref{$PetscInt}()

    @chk ccall(
               (:PCGetKSPNestLevel, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{$PetscInt}),
               pc, level_,
              )

	level = level_[]

	return level
end 

"""
	newpc::PC = PCCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Creates a preconditioner context, `PC`

Collective

Input Parameter:
- `comm` - MPI communicator

Output Parameter:
- `newpc` - location to put the `PC` preconditioner context

Level: developer

-seealso: [](ch_ksp), `PC`, `PCType`, `PCSetType`, `PCSetUp()`, `PCApply()`, `PCDestroy()`, `KSP`, `KSPGetPC()`

# External Links
$(_doc_external("Ksp/PCCreate"))
"""
function PCCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function PCCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	newpc_ = Ref{PC}()

    @chk ccall(
               (:PCCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PC}),
               comm, newpc_,
              )

	newpc = newpc_[]

	return newpc
end 

"""
	PCApply(petsclib::PetscLibType,pc::PC, x::PetscVec, y::PetscVec) 
Applies the preconditioner to a vector.

Collective

Input Parameters:
- `pc` - the `PC` preconditioner context
- `x`  - input vector

Output Parameter:
- `y` - output vector

Level: developer

-seealso: [](ch_ksp), `PC`, `PCApplyTranspose()`, `PCApplyBAorAB()`

# External Links
$(_doc_external("Ksp/PCApply"))
"""
function PCApply(petsclib::PetscLibType, pc::PC, x::PetscVec, y::PetscVec) end

@for_petsc function PCApply(petsclib::$UnionPetscLib, pc::PC, x::PetscVec, y::PetscVec )

    @chk ccall(
               (:PCApply, $petsc_library),
               PetscErrorCode,
               (PC, CVec, CVec),
               pc, x, y,
              )


	return nothing
end 

"""
	PCMatApply(petsclib::PetscLibType,pc::PC, X::PetscMat, Y::PetscMat) 
Applies the preconditioner to multiple vectors stored as a `MATDENSE`. Like `PCApply()`, `Y` and `X` must be different matrices.

Collective

Input Parameters:
- `pc` - the `PC` preconditioner context
- `X`  - block of input vectors

Output Parameter:
- `Y` - block of output vectors

Level: developer

-seealso: [](ch_ksp), `PC`, `PCApply()`, `KSPMatSolve()`

# External Links
$(_doc_external("Ksp/PCMatApply"))
"""
function PCMatApply(petsclib::PetscLibType, pc::PC, X::PetscMat, Y::PetscMat) end

@for_petsc function PCMatApply(petsclib::$UnionPetscLib, pc::PC, X::PetscMat, Y::PetscMat )

    @chk ccall(
               (:PCMatApply, $petsc_library),
               PetscErrorCode,
               (PC, CMat, CMat),
               pc, X, Y,
              )


	return nothing
end 

"""
	PCMatApplyTranspose(petsclib::PetscLibType,pc::PC, X::PetscMat, Y::PetscMat) 
Applies the transpose of preconditioner to multiple vectors stored as a `MATDENSE`. Like `PCApplyTranspose()`, `Y` and `X` must be different matrices.

Collective

Input Parameters:
- `pc` - the `PC` preconditioner context
- `X`  - block of input vectors

Output Parameter:
- `Y` - block of output vectors

Level: developer

-seealso: [](ch_ksp), `PC`, `PCApplyTranspose()`, `KSPMatSolveTranspose()`

# External Links
$(_doc_external("Ksp/PCMatApplyTranspose"))
"""
function PCMatApplyTranspose(petsclib::PetscLibType, pc::PC, X::PetscMat, Y::PetscMat) end

@for_petsc function PCMatApplyTranspose(petsclib::$UnionPetscLib, pc::PC, X::PetscMat, Y::PetscMat )

    @chk ccall(
               (:PCMatApplyTranspose, $petsc_library),
               PetscErrorCode,
               (PC, CMat, CMat),
               pc, X, Y,
              )


	return nothing
end 

"""
	PCApplySymmetricLeft(petsclib::PetscLibType,pc::PC, x::PetscVec, y::PetscVec) 
Applies the left part of a symmetric preconditioner to a vector.

Collective

Input Parameters:
- `pc` - the `PC` preconditioner context
- `x`  - input vector

Output Parameter:
- `y` - output vector

Level: developer

-seealso: [](ch_ksp), `PC`, `PCApply()`, `PCApplySymmetricRight()`

# External Links
$(_doc_external("Ksp/PCApplySymmetricLeft"))
"""
function PCApplySymmetricLeft(petsclib::PetscLibType, pc::PC, x::PetscVec, y::PetscVec) end

@for_petsc function PCApplySymmetricLeft(petsclib::$UnionPetscLib, pc::PC, x::PetscVec, y::PetscVec )

    @chk ccall(
               (:PCApplySymmetricLeft, $petsc_library),
               PetscErrorCode,
               (PC, CVec, CVec),
               pc, x, y,
              )


	return nothing
end 

"""
	PCApplySymmetricRight(petsclib::PetscLibType,pc::PC, x::PetscVec, y::PetscVec) 
Applies the right part of a symmetric preconditioner to a vector.

Collective

Input Parameters:
- `pc` - the `PC` preconditioner context
- `x`  - input vector

Output Parameter:
- `y` - output vector

Level: developer

-seealso: [](ch_ksp), `PC`, `PCApply()`, `PCApplySymmetricLeft()`

# External Links
$(_doc_external("Ksp/PCApplySymmetricRight"))
"""
function PCApplySymmetricRight(petsclib::PetscLibType, pc::PC, x::PetscVec, y::PetscVec) end

@for_petsc function PCApplySymmetricRight(petsclib::$UnionPetscLib, pc::PC, x::PetscVec, y::PetscVec )

    @chk ccall(
               (:PCApplySymmetricRight, $petsc_library),
               PetscErrorCode,
               (PC, CVec, CVec),
               pc, x, y,
              )


	return nothing
end 

"""
	PCApplyTranspose(petsclib::PetscLibType,pc::PC, x::PetscVec, y::PetscVec) 
Applies the transpose of preconditioner to a vector.

Collective

Input Parameters:
- `pc` - the `PC` preconditioner context
- `x`  - input vector

Output Parameter:
- `y` - output vector

Level: developer

-seealso: [](ch_ksp), `PC`, `PCApply()`, `PCApplyBAorAB()`, `PCApplyBAorABTranspose()`, `PCApplyTransposeExists()`

# External Links
$(_doc_external("Ksp/PCApplyTranspose"))
"""
function PCApplyTranspose(petsclib::PetscLibType, pc::PC, x::PetscVec, y::PetscVec) end

@for_petsc function PCApplyTranspose(petsclib::$UnionPetscLib, pc::PC, x::PetscVec, y::PetscVec )

    @chk ccall(
               (:PCApplyTranspose, $petsc_library),
               PetscErrorCode,
               (PC, CVec, CVec),
               pc, x, y,
              )


	return nothing
end 

"""
	flg::PetscBool = PCApplyTransposeExists(petsclib::PetscLibType,pc::PC) 
Test whether the preconditioner has a transpose apply operation

Collective

Input Parameter:
- `pc` - the `PC` preconditioner context

Output Parameter:
- `flg` - `PETSC_TRUE` if a transpose operation is defined

Level: developer

-seealso: [](ch_ksp), `PC`, `PCApplyTranspose()`

# External Links
$(_doc_external("Ksp/PCApplyTransposeExists"))
"""
function PCApplyTransposeExists(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCApplyTransposeExists(petsclib::$UnionPetscLib, pc::PC )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PCApplyTransposeExists, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PetscBool}),
               pc, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	PCApplyBAorAB(petsclib::PetscLibType,pc::PC, side::PCSide, x::PetscVec, y::PetscVec, work::PetscVec) 
Applies the preconditioner and operator to a vector. y = B*A*x  or  y = A*B*x.

Collective

Input Parameters:
- `pc`   - the `PC` preconditioner context
- `side` - indicates the preconditioner side, one of `PC_LEFT`, `PC_RIGHT`, or `PC_SYMMETRIC`
- `x`    - input vector
- `work` - work vector

Output Parameter:
- `y` - output vector

Level: developer

-seealso: [](ch_ksp), `PC`, `PCApply()`, `PCApplyTranspose()`, `PCApplyBAorABTranspose()`

# External Links
$(_doc_external("Ksp/PCApplyBAorAB"))
"""
function PCApplyBAorAB(petsclib::PetscLibType, pc::PC, side::PCSide, x::PetscVec, y::PetscVec, work::PetscVec) end

@for_petsc function PCApplyBAorAB(petsclib::$UnionPetscLib, pc::PC, side::PCSide, x::PetscVec, y::PetscVec, work::PetscVec )

    @chk ccall(
               (:PCApplyBAorAB, $petsc_library),
               PetscErrorCode,
               (PC, PCSide, CVec, CVec, CVec),
               pc, side, x, y, work,
              )


	return nothing
end 

"""
	PCApplyBAorABTranspose(petsclib::PetscLibType,pc::PC, side::PCSide, x::PetscVec, y::PetscVec, work::PetscVec) 
Applies the transpose of the preconditioner
and operator to a vector. That is, applies B^T * A^T with left preconditioning,
NOT (B*A)^T = A^T*B^T.

Collective

Input Parameters:
- `pc`   - the `PC` preconditioner context
- `side` - indicates the preconditioner side, one of `PC_LEFT`, `PC_RIGHT`, or `PC_SYMMETRIC`
- `x`    - input vector
- `work` - work vector

Output Parameter:
- `y` - output vector

Level: developer

-seealso: [](ch_ksp), `PC`, `PCApply()`, `PCApplyTranspose()`, `PCApplyBAorAB()`

# External Links
$(_doc_external("Ksp/PCApplyBAorABTranspose"))
"""
function PCApplyBAorABTranspose(petsclib::PetscLibType, pc::PC, side::PCSide, x::PetscVec, y::PetscVec, work::PetscVec) end

@for_petsc function PCApplyBAorABTranspose(petsclib::$UnionPetscLib, pc::PC, side::PCSide, x::PetscVec, y::PetscVec, work::PetscVec )

    @chk ccall(
               (:PCApplyBAorABTranspose, $petsc_library),
               PetscErrorCode,
               (PC, PCSide, CVec, CVec, CVec),
               pc, side, x, y, work,
              )


	return nothing
end 

"""
	exists::PetscBool = PCApplyRichardsonExists(petsclib::PetscLibType,pc::PC) 
Determines whether a particular preconditioner has a
built-in fast application of Richardson's method.

Not Collective

Input Parameter:
- `pc` - the preconditioner

Output Parameter:
- `exists` - `PETSC_TRUE` or `PETSC_FALSE`

Level: developer

-seealso: [](ch_ksp), `PC`, `KSPRICHARDSON`, `PCApplyRichardson()`

# External Links
$(_doc_external("Ksp/PCApplyRichardsonExists"))
"""
function PCApplyRichardsonExists(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCApplyRichardsonExists(petsclib::$UnionPetscLib, pc::PC )
	exists_ = Ref{PetscBool}()

    @chk ccall(
               (:PCApplyRichardsonExists, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PetscBool}),
               pc, exists_,
              )

	exists = exists_[]

	return exists
end 

"""
	outits::PetscInt = PCApplyRichardson(petsclib::PetscLibType,pc::PC, b::PetscVec, y::PetscVec, w::PetscVec, rtol::PetscReal, abstol::PetscReal, dtol::PetscReal, its::PetscInt, guesszero::PetscBool, reason::PCRiCchardsonConvergedReason) 
Applies several steps of Richardson iteration with
the particular preconditioner. This routine is usually used by the
Krylov solvers and not the application code directly.

Collective

Input Parameters:
- `pc`        - the `PC` preconditioner context
- `b`         - the right-hand side
- `w`         - one work vector
- `rtol`      - relative decrease in residual norm convergence criteria
- `abstol`    - absolute residual norm convergence criteria
- `dtol`      - divergence residual norm increase criteria
- `its`       - the number of iterations to apply.
- `guesszero` - if the input x contains nonzero initial guess

Output Parameters:
- `outits` - number of iterations actually used (for SOR this always equals its)
- `reason` - the reason the apply terminated
- `y`      - the solution (also contains initial guess if guesszero is `PETSC_FALSE`

Level: developer

-seealso: [](ch_ksp), `PC`, `PCApplyRichardsonExists()`

# External Links
$(_doc_external("Ksp/PCApplyRichardson"))
"""
function PCApplyRichardson(petsclib::PetscLibType, pc::PC, b::PetscVec, y::PetscVec, w::PetscVec, rtol::PetscReal, abstol::PetscReal, dtol::PetscReal, its::PetscInt, guesszero::PetscBool, reason::PCRiCchardsonConvergedReason) end

@for_petsc function PCApplyRichardson(petsclib::$UnionPetscLib, pc::PC, b::PetscVec, y::PetscVec, w::PetscVec, rtol::$PetscReal, abstol::$PetscReal, dtol::$PetscReal, its::$PetscInt, guesszero::PetscBool, reason::PCRiCchardsonConvergedReason )
	outits_ = Ref{$PetscInt}()

    @chk ccall(
               (:PCApplyRichardson, $petsc_library),
               PetscErrorCode,
               (PC, CVec, CVec, CVec, $PetscReal, $PetscReal, $PetscReal, $PetscInt, PetscBool, Ptr{$PetscInt}, Ptr{PCRiCchardsonConvergedReason}),
               pc, b, y, w, rtol, abstol, dtol, its, guesszero, outits_, reason,
              )

	outits = outits_[]

	return outits
end 

"""
	PCSetFailedReason(petsclib::PetscLibType,pc::PC, reason::PCFailedReason) 
Sets the reason a `PCSetUp()` failed or `PC_NOERROR` if it did not fail

Logically Collective

Input Parameters:
- `pc`     - the `PC` preconditioner context
- `reason` - the reason it failed

Level: advanced

-seealso: [](ch_ksp), `PC`, `PCCreate()`, `PCApply()`, `PCDestroy()`, `PCFailedReason`

# External Links
$(_doc_external("Ksp/PCSetFailedReason"))
"""
function PCSetFailedReason(petsclib::PetscLibType, pc::PC, reason::PCFailedReason) end

@for_petsc function PCSetFailedReason(petsclib::$UnionPetscLib, pc::PC, reason::PCFailedReason )

    @chk ccall(
               (:PCSetFailedReason, $petsc_library),
               PetscErrorCode,
               (PC, PCFailedReason),
               pc, reason,
              )


	return nothing
end 

"""
	PCGetFailedReason(petsclib::PetscLibType,pc::PC, reason::PCFailedReason) 
Gets the reason a `PCSetUp()` failed or `PC_NOERROR` if it did not fail

Not Collective

Input Parameter:
- `pc` - the `PC` preconditioner context

Output Parameter:
- `reason` - the reason it failed

Level: advanced

-seealso: [](ch_ksp), `PC`, `PCCreate()`, `PCApply()`, `PCDestroy()`, `PCSetFailedReason()`, `PCFailedReason`

# External Links
$(_doc_external("Ksp/PCGetFailedReason"))
"""
function PCGetFailedReason(petsclib::PetscLibType, pc::PC, reason::PCFailedReason) end

@for_petsc function PCGetFailedReason(petsclib::$UnionPetscLib, pc::PC, reason::PCFailedReason )

    @chk ccall(
               (:PCGetFailedReason, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PCFailedReason}),
               pc, reason,
              )


	return nothing
end 

"""
	PCReduceFailedReason(petsclib::PetscLibType,pc::PC) 
Reduce the failed reason among the MPI processes that share the `PC`

Collective

Input Parameter:
- `pc` - the `PC` preconditioner context

Level: advanced

-seealso: [](ch_ksp), `PC`, `PCCreate()`, `PCApply()`, `PCDestroy()`, `PCGetFailedReason()`, `PCSetFailedReason()`, `PCFailedReason`

# External Links
$(_doc_external("Ksp/PCReduceFailedReason"))
"""
function PCReduceFailedReason(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCReduceFailedReason(petsclib::$UnionPetscLib, pc::PC )

    @chk ccall(
               (:PCReduceFailedReason, $petsc_library),
               PetscErrorCode,
               (PC,),
               pc,
              )


	return nothing
end 

"""
	PCSetUp(petsclib::PetscLibType,pc::PC) 
Prepares for the use of a preconditioner. Performs all the one
can be used with `PCApply()`

Collective

Input Parameter:
- `pc` - the `PC` preconditioner context

Level: developer

-seealso: [](ch_ksp), `PC`, `PCCreate()`, `PCApply()`, `PCDestroy()`, `KSPSetUp()`, `PCSetUpOnBlocks()`

# External Links
$(_doc_external("Ksp/PCSetUp"))
"""
function PCSetUp(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCSetUp(petsclib::$UnionPetscLib, pc::PC )

    @chk ccall(
               (:PCSetUp, $petsc_library),
               PetscErrorCode,
               (PC,),
               pc,
              )


	return nothing
end 

"""
	PCSetUpOnBlocks(petsclib::PetscLibType,pc::PC) 
Sets up the preconditioner for each block in
the block Jacobi, overlapping Schwarz, and fieldsplit methods.

Collective

Input Parameter:
- `pc` - the `PC` preconditioner context

Level: developer

-seealso: [](ch_ksp), `PC`, `PCSetUp()`, `PCCreate()`, `PCApply()`, `PCDestroy()`

# External Links
$(_doc_external("Ksp/PCSetUpOnBlocks"))
"""
function PCSetUpOnBlocks(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCSetUpOnBlocks(petsclib::$UnionPetscLib, pc::PC )

    @chk ccall(
               (:PCSetUpOnBlocks, $petsc_library),
               PetscErrorCode,
               (PC,),
               pc,
              )


	return nothing
end 

"""
	PCSetModifySubMatrices(petsclib::PetscLibType,pc::PC, func::PCModifySubMatricesFn, ctx::Cvoid) 
Sets a user
submatrices that arise within certain subdomain-based preconditioners such as `PCASM`

Logically Collective

Input Parameters:
- `pc`   - the `PC` preconditioner context
- `func` - routine for modifying the submatrices, see `PCModifySubMatricesFn`
- `ctx`  - optional user-defined context (may be `NULL`)

Level: advanced

-seealso: [](ch_ksp), `PC`, `PCModifySubMatricesFn`, `PCBJACOBI`, `PCASM`, `PCModifySubMatrices()`

# External Links
$(_doc_external("Ksp/PCSetModifySubMatrices"))
"""
function PCSetModifySubMatrices(petsclib::PetscLibType, pc::PC, func::PCModifySubMatricesFn, ctx::Cvoid) end

@for_petsc function PCSetModifySubMatrices(petsclib::$UnionPetscLib, pc::PC, func::PCModifySubMatricesFn, ctx::Cvoid )

    @chk ccall(
               (:PCSetModifySubMatrices, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PCModifySubMatricesFn}, Ptr{Cvoid}),
               pc, func, ctx,
              )


	return nothing
end 

"""
	PCModifySubMatrices(petsclib::PetscLibType,pc::PC, nsub::PetscInt, row::Vector{IS}, col::Vector{IS}, submat::Vector{PetscMat}, ctx::Cvoid) 
Calls an optional user
certain preconditioners if one has been set with `PCSetModifySubMatrices()`.

Collective

Input Parameters:
- `pc`     - the `PC` preconditioner context
- `nsub`   - the number of local submatrices
- `row`    - an array of index sets that contain the global row numbers
that comprise each local submatrix
- `col`    - an array of index sets that contain the global column numbers
that comprise each local submatrix
- `submat` - array of local submatrices
- `ctx`    - optional user-defined context for private data for the
user-defined routine (may be `NULL`)

Output Parameter:
- `submat` - array of local submatrices (the entries of which may
have been modified)

Level: developer

-seealso: [](ch_ksp), `PC`, `PCModifySubMatricesFn`, `PCSetModifySubMatrices()`

# External Links
$(_doc_external("Ksp/PCModifySubMatrices"))
"""
function PCModifySubMatrices(petsclib::PetscLibType, pc::PC, nsub::PetscInt, row::Vector{IS}, col::Vector{IS}, submat::Vector{PetscMat}, ctx::Cvoid) end

@for_petsc function PCModifySubMatrices(petsclib::$UnionPetscLib, pc::PC, nsub::$PetscInt, row::Vector{IS}, col::Vector{IS}, submat::Vector{PetscMat}, ctx::Cvoid )

    @chk ccall(
               (:PCModifySubMatrices, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, Ptr{IS}, Ptr{IS}, Ptr{CMat}, Ptr{Cvoid}),
               pc, nsub, row, col, submat, ctx,
              )


	return nothing
end 

"""
	PCSetOperators(petsclib::PetscLibType,pc::PC, Amat::PetscMat, Pmat::PetscMat) 
Sets the matrix associated with the linear system and
a (possibly) different one associated with the preconditioner.

Logically Collective

Input Parameters:
- `pc`   - the `PC` preconditioner context
- `Amat` - the matrix that defines the linear system
- `Pmat` - the matrix to be used in constructing the preconditioner, usually the same as Amat.

Level: intermediate

-seealso: [](ch_ksp), `PC`, `PCGetOperators()`, `MatZeroEntries()`

# External Links
$(_doc_external("Ksp/PCSetOperators"))
"""
function PCSetOperators(petsclib::PetscLibType, pc::PC, Amat::PetscMat, Pmat::PetscMat) end

@for_petsc function PCSetOperators(petsclib::$UnionPetscLib, pc::PC, Amat::PetscMat, Pmat::PetscMat )

    @chk ccall(
               (:PCSetOperators, $petsc_library),
               PetscErrorCode,
               (PC, CMat, CMat),
               pc, Amat, Pmat,
              )


	return nothing
end 

"""
	PCSetReusePreconditioner(petsclib::PetscLibType,pc::PC, flag::PetscBool) 
reuse the current preconditioner even if the operator in the preconditioner `PC` has changed.

Logically Collective

Input Parameters:
- `pc`   - the `PC` preconditioner context
- `flag` - `PETSC_TRUE` do not compute a new preconditioner, `PETSC_FALSE` do compute a new preconditioner

Level: intermediate

-seealso: [](ch_ksp), `PC`, `PCGetOperators()`, `MatZeroEntries()`, `PCGetReusePreconditioner()`, `KSPSetReusePreconditioner()`

# External Links
$(_doc_external("Ksp/PCSetReusePreconditioner"))
"""
function PCSetReusePreconditioner(petsclib::PetscLibType, pc::PC, flag::PetscBool) end

@for_petsc function PCSetReusePreconditioner(petsclib::$UnionPetscLib, pc::PC, flag::PetscBool )

    @chk ccall(
               (:PCSetReusePreconditioner, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, flag,
              )


	return nothing
end 

"""
	flag::PetscBool = PCGetReusePreconditioner(petsclib::PetscLibType,pc::PC) 
Determines if the `PC` reuses the current preconditioner even if the operator in the preconditioner has changed.

Not Collective

Input Parameter:
- `pc` - the `PC` preconditioner context

Output Parameter:
- `flag` - `PETSC_TRUE` do not compute a new preconditioner, `PETSC_FALSE` do compute a new preconditioner

Level: intermediate

-seealso: [](ch_ksp), `PC`, `PCGetOperators()`, `MatZeroEntries()`, `PCSetReusePreconditioner()`

# External Links
$(_doc_external("Ksp/PCGetReusePreconditioner"))
"""
function PCGetReusePreconditioner(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCGetReusePreconditioner(petsclib::$UnionPetscLib, pc::PC )
	flag_ = Ref{PetscBool}()

    @chk ccall(
               (:PCGetReusePreconditioner, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PetscBool}),
               pc, flag_,
              )

	flag = flag_[]

	return flag
end 

"""
	PCGetOperators(petsclib::PetscLibType,pc::PC, Amat::PetscMat, Pmat::PetscMat) 
Gets the matrix associated with the linear system and
possibly a different one which is used to construct the preconditioner.

Not Collective, though parallel `Mat`s are returned if `pc` is parallel

Input Parameter:
- `pc` - the `PC` preconditioner context

Output Parameters:
- `Amat` - the matrix defining the linear system
- `Pmat` - the matrix from which the preconditioner is constructed, usually the same as Amat.

Level: intermediate

-seealso: [](ch_ksp), `PC`, `PCSetOperators()`, `KSPGetOperators()`, `KSPSetOperators()`, `PCGetOperatorsSet()`

# External Links
$(_doc_external("Ksp/PCGetOperators"))
"""
function PCGetOperators(petsclib::PetscLibType, pc::PC, Amat::PetscMat, Pmat::PetscMat) end

@for_petsc function PCGetOperators(petsclib::$UnionPetscLib, pc::PC, Amat::PetscMat, Pmat::PetscMat )
	Amat_ = Ref(Amat.ptr)
	Pmat_ = Ref(Pmat.ptr)

    @chk ccall(
               (:PCGetOperators, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{CMat}, Ptr{CMat}),
               pc, Amat_, Pmat_,
              )

	Amat.ptr = C_NULL
	Pmat.ptr = C_NULL

	return nothing
end 

"""
	mat::PetscBool,pmat::PetscBool = PCGetOperatorsSet(petsclib::PetscLibType,pc::PC) 
Determines if the matrix associated with the linear system and
possibly a different one associated with the preconditioner have been set in the `PC`.

Not Collective, though the results on all processes should be the same

Input Parameter:
- `pc` - the `PC` preconditioner context

Output Parameters:
- `mat`  - the matrix associated with the linear system was set
- `pmat` - matrix associated with the preconditioner was set, usually the same

Level: intermediate

-seealso: [](ch_ksp), `PC`, `PCSetOperators()`, `KSPGetOperators()`, `KSPSetOperators()`, `PCGetOperators()`

# External Links
$(_doc_external("Ksp/PCGetOperatorsSet"))
"""
function PCGetOperatorsSet(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCGetOperatorsSet(petsclib::$UnionPetscLib, pc::PC )
	mat_ = Ref{PetscBool}()
	pmat_ = Ref{PetscBool}()

    @chk ccall(
               (:PCGetOperatorsSet, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PetscBool}, Ptr{PetscBool}),
               pc, mat_, pmat_,
              )

	mat = mat_[]
	pmat = pmat_[]

	return mat,pmat
end 

"""
	PCFactorGetMatrix(petsclib::PetscLibType,pc::PC, mat::PetscMat) 
Gets the factored matrix from the
preconditioner context.  This routine is valid only for the `PCLU`,
`PCILU`, `PCCHOLESKY`, and `PCICC` methods.

Not Collective though `mat` is parallel if `pc` is parallel

Input Parameter:
- `pc` - the `PC` preconditioner context

Output Parameters:
- `mat` - the factored matrix

Level: advanced

-seealso: [](ch_ksp), `PC`, `PCLU`, `PCILU`, `PCCHOLESKY`, `PCICC`

# External Links
$(_doc_external("Ksp/PCFactorGetMatrix"))
"""
function PCFactorGetMatrix(petsclib::PetscLibType, pc::PC, mat::PetscMat) end

@for_petsc function PCFactorGetMatrix(petsclib::$UnionPetscLib, pc::PC, mat::PetscMat )
	mat_ = Ref(mat.ptr)

    @chk ccall(
               (:PCFactorGetMatrix, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{CMat}),
               pc, mat_,
              )

	mat.ptr = C_NULL

	return nothing
end 

"""
	PCSetOptionsPrefix(petsclib::PetscLibType,pc::PC, prefix::String) 
Sets the prefix used for searching for all
`PC` options in the database.

Logically Collective

Input Parameters:
- `pc`     - the `PC` preconditioner context
- `prefix` - the prefix string to prepend to all `PC` option requests

-seealso: [](ch_ksp), `PC`, `PCSetFromOptions`, `PCAppendOptionsPrefix()`, `PCGetOptionsPrefix()`

# External Links
$(_doc_external("Ksp/PCSetOptionsPrefix"))
"""
function PCSetOptionsPrefix(petsclib::PetscLibType, pc::PC, prefix::String) end

@for_petsc function PCSetOptionsPrefix(petsclib::$UnionPetscLib, pc::PC, prefix::String )

    @chk ccall(
               (:PCSetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{Cchar}),
               pc, prefix,
              )


	return nothing
end 

"""
	PCAppendOptionsPrefix(petsclib::PetscLibType,pc::PC, prefix::String) 
Appends to the prefix used for searching for all
`PC` options in the database.

Logically Collective

Input Parameters:
- `pc`     - the `PC` preconditioner context
- `prefix` - the prefix string to prepend to all `PC` option requests

-seealso: [](ch_ksp), `PC`, `PCSetFromOptions`, `PCSetOptionsPrefix()`, `PCGetOptionsPrefix()`

# External Links
$(_doc_external("Ksp/PCAppendOptionsPrefix"))
"""
function PCAppendOptionsPrefix(petsclib::PetscLibType, pc::PC, prefix::String) end

@for_petsc function PCAppendOptionsPrefix(petsclib::$UnionPetscLib, pc::PC, prefix::String )

    @chk ccall(
               (:PCAppendOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{Cchar}),
               pc, prefix,
              )


	return nothing
end 

"""
	PCGetOptionsPrefix(petsclib::PetscLibType,pc::PC, prefix::String) 
Gets the prefix used for searching for all
PC options in the database.

Not Collective

Input Parameter:
- `pc` - the `PC` preconditioner context

Output Parameter:
- `prefix` - pointer to the prefix string used, is returned

Level: advanced

-seealso: [](ch_ksp), `PC`, `PCSetFromOptions`, `PCSetOptionsPrefix()`, `PCAppendOptionsPrefix()`

# External Links
$(_doc_external("Ksp/PCGetOptionsPrefix"))
"""
function PCGetOptionsPrefix(petsclib::PetscLibType, pc::PC, prefix::String) end

@for_petsc function PCGetOptionsPrefix(petsclib::$UnionPetscLib, pc::PC, prefix::String )
	prefix_ = Ref(pointer(prefix))

    @chk ccall(
               (:PCGetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{Ptr{Cchar}}),
               pc, prefix_,
              )


	return nothing
end 

"""
	PCPreSolve(petsclib::PetscLibType,pc::PC, ksp::PetscKSP) 
Optional pre
the iterative solve itself. Used in conjunction with `PCPostSolve()`

Collective

Input Parameters:
- `pc`  - the `PC` preconditioner context
- `ksp` - the Krylov subspace context

Level: developer

-seealso: [](ch_ksp), `PC`, `PCPostSolve()`, `KSP`, `PCSetPostSetUp()`, `KSPSetPreSolve()`, `KSPSetPostSolve()`

# External Links
$(_doc_external("Ksp/PCPreSolve"))
"""
function PCPreSolve(petsclib::PetscLibType, pc::PC, ksp::PetscKSP) end

@for_petsc function PCPreSolve(petsclib::$UnionPetscLib, pc::PC, ksp::PetscKSP )

    @chk ccall(
               (:PCPreSolve, $petsc_library),
               PetscErrorCode,
               (PC, CKSP),
               pc, ksp,
              )


	return nothing
end 

"""
	PCSetPostSetUp(petsclib::PetscLibType,pc::PC, postsetup::external) 
Sets function called at the end of `PCSetUp()` to adjust the computed preconditioner

Logically Collective

Input Parameters:
- `pc`        - the preconditioner object
- `postsetup` - the function to call after `PCSetUp()`

Calling sequence of `postsetup`:
- `pc` - the `PC` context

Level: developer

-seealso: [](ch_ksp), `PC`, `PCSetUp()`

# External Links
$(_doc_external("Ksp/PCSetPostSetUp"))
"""
function PCSetPostSetUp(petsclib::PetscLibType, pc::PC, postsetup::external) end

@for_petsc function PCSetPostSetUp(petsclib::$UnionPetscLib, pc::PC, postsetup::external )

    @chk ccall(
               (:PCSetPostSetUp, $petsc_library),
               PetscErrorCode,
               (PC, external),
               pc, postsetup,
              )


	return nothing
end 

"""
	PCPostSolve(petsclib::PetscLibType,pc::PC, ksp::PetscKSP) 
Optional post
preconditioner-specific actions that must be performed after
the iterative solve itself.

Collective

Input Parameters:
- `pc`  - the `PC` preconditioner context
- `ksp` - the `KSP` Krylov subspace context

-seealso: [](ch_ksp), `PC`, `KSPSetPostSolve()`, `KSPSetPreSolve()`, `PCPreSolve()`, `KSPSolve()`

# External Links
$(_doc_external("Ksp/PCPostSolve"))
"""
function PCPostSolve(petsclib::PetscLibType, pc::PC, ksp::PetscKSP) end

@for_petsc function PCPostSolve(petsclib::$UnionPetscLib, pc::PC, ksp::PetscKSP )

    @chk ccall(
               (:PCPostSolve, $petsc_library),
               PetscErrorCode,
               (PC, CKSP),
               pc, ksp,
              )


	return nothing
end 

"""
	PCLoad(petsclib::PetscLibType,newdm::PC, viewer::PetscViewer) 
Loads a `PC` that has been stored in binary  with `PCView()`.

Collective

Input Parameters:
- `newdm`  - the newly loaded `PC`, this needs to have been created with `PCCreate()` or
some related function before a call to `PCLoad()`.
- `viewer` - binary file viewer `PETSCVIEWERBINARY`, obtained from `PetscViewerBinaryOpen()`

Level: intermediate

-seealso: [](ch_ksp), `PC`, `PetscViewerBinaryOpen()`, `PCView()`, `MatLoad()`, `VecLoad()`, `PETSCVIEWERBINARY`

# External Links
$(_doc_external("Ksp/PCLoad"))
"""
function PCLoad(petsclib::PetscLibType, newdm::PC, viewer::PetscViewer) end

@for_petsc function PCLoad(petsclib::$UnionPetscLib, newdm::PC, viewer::PetscViewer )

    @chk ccall(
               (:PCLoad, $petsc_library),
               PetscErrorCode,
               (PC, PetscViewer),
               newdm, viewer,
              )


	return nothing
end 

"""
	PCViewFromOptions(petsclib::PetscLibType,A::PC, obj::PetscObject, name::String) 
View (print or provide information about) the `PC`, based on options in the options database

Collective

Input Parameters:
- `A`    - the `PC` context
- `obj`  - Optional object that provides the options prefix
- `name` - command line option name

Level: developer

-seealso: [](ch_ksp), `PC`, `PCView`, `PetscObjectViewFromOptions()`, `PCCreate()`

# External Links
$(_doc_external("Ksp/PCViewFromOptions"))
"""
function PCViewFromOptions(petsclib::PetscLibType, A::PC, obj::PetscObject, name::String) end

@for_petsc function PCViewFromOptions(petsclib::$UnionPetscLib, A::PC, obj::PetscObject, name::String )

    @chk ccall(
               (:PCViewFromOptions, $petsc_library),
               PetscErrorCode,
               (PC, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

"""
	PCView(petsclib::PetscLibType,pc::PC, viewer::PetscViewer) 
Prints information about the `PC`

Collective

Input Parameters:
- `pc`     - the `PC` preconditioner context
- `viewer` - optional `PetscViewer` visualization context

Level: intermediate

-seealso: [](ch_ksp), `PC`, `PetscViewer`, `PetscViewerType`, `KSPView()`, `PetscViewerASCIIOpen()`

# External Links
$(_doc_external("Ksp/PCView"))
"""
function PCView(petsclib::PetscLibType, pc::PC, viewer::PetscViewer) end

@for_petsc function PCView(petsclib::$UnionPetscLib, pc::PC, viewer::PetscViewer )

    @chk ccall(
               (:PCView, $petsc_library),
               PetscErrorCode,
               (PC, PetscViewer),
               pc, viewer,
              )


	return nothing
end 

"""
	PCRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Adds a method (`PCType`) to the PETSc preconditioner package.

Not collective. No Fortran Support

Input Parameters:
- `sname`    - name of a new user-defined solver
- `function` - routine to create the method context which will be stored in a `PC` when `PCSetType()` is called

-seealso: [](ch_ksp), `PC`, `PCType`, `PCRegisterAll()`, `PCSetType()`, `PCShellSetContext()`, `PCShellSetApply()`, `PCSHELL`

# External Links
$(_doc_external("Ksp/PCRegister"))
"""
function PCRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function PCRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:PCRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	PCComputeOperator(petsclib::PetscLibType,pc::PC, mattype::MatType, mat::PetscMat) 
Computes the explicit preconditioned operator as a matrix `Mat`.

Collective

Input Parameters:
- `pc`      - the `PC` preconditioner object
- `mattype` - the `MatType` to be used for the operator

Output Parameter:
- `mat` - the explicit preconditioned operator

Level: advanced

-seealso: [](ch_ksp), `PC`, `KSPComputeOperator()`, `MatType`

# External Links
$(_doc_external("Ksp/PCComputeOperator"))
"""
function PCComputeOperator(petsclib::PetscLibType, pc::PC, mattype::MatType, mat::PetscMat) end

@for_petsc function PCComputeOperator(petsclib::$UnionPetscLib, pc::PC, mattype::MatType, mat::PetscMat )
	mat_ = Ref(mat.ptr)

    @chk ccall(
               (:PCComputeOperator, $petsc_library),
               PetscErrorCode,
               (PC, MatType, Ptr{CMat}),
               pc, mattype, mat_,
              )

	mat.ptr = C_NULL

	return nothing
end 

"""
	PCSetCoordinates(petsclib::PetscLibType,pc::PC, dim::PetscInt, nloc::PetscInt, coords::Vector{PetscReal}) 
sets the coordinates of all the nodes (degrees of freedom in the vector) on the local process

Collective

Input Parameters:
- `pc`     - the `PC` preconditioner context
- `dim`    - the dimension of the coordinates 1, 2, or 3
- `nloc`   - the blocked size of the coordinates array
- `coords` - the coordinates array

Level: intermediate

-seealso: [](ch_ksp), `PC`, `MatSetNearNullSpace()`

# External Links
$(_doc_external("Ksp/PCSetCoordinates"))
"""
function PCSetCoordinates(petsclib::PetscLibType, pc::PC, dim::PetscInt, nloc::PetscInt, coords::Vector{PetscReal}) end

@for_petsc function PCSetCoordinates(petsclib::$UnionPetscLib, pc::PC, dim::$PetscInt, nloc::$PetscInt, coords::Vector{$PetscReal} )

    @chk ccall(
               (:PCSetCoordinates, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, $PetscInt, Ptr{$PetscReal}),
               pc, dim, nloc, coords,
              )


	return nothing
end 

"""
	num_levels::PetscInt = PCGetInterpolations(petsclib::PetscLibType,pc::PC, interpolations::Vector{PetscMat}) 
Gets interpolation matrices for all levels (except level 0)

Logically Collective

Input Parameter:
- `pc` - the precondition context

Output Parameters:
- `num_levels`     - the number of levels
- `interpolations` - the interpolation matrices (size of `num_levels`-1)

Level: advanced

-seealso: [](ch_ksp), `PC`, `PCMG`, `PCMGGetRestriction()`, `PCMGSetInterpolation()`, `PCMGGetInterpolation()`, `PCGetCoarseOperators()`

# External Links
$(_doc_external("Ksp/PCGetInterpolations"))
"""
function PCGetInterpolations(petsclib::PetscLibType, pc::PC, interpolations::Vector{PetscMat}) end

@for_petsc function PCGetInterpolations(petsclib::$UnionPetscLib, pc::PC, interpolations::Vector{PetscMat} )
	num_levels_ = Ref{$PetscInt}()
	interpolations_ = Ref(pointer(interpolations))

    @chk ccall(
               (:PCGetInterpolations, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{$PetscInt}, Ptr{Ptr{CMat}}),
               pc, num_levels_, interpolations_,
              )

	num_levels = num_levels_[]

	return num_levels
end 

"""
	num_levels::PetscInt = PCGetCoarseOperators(petsclib::PetscLibType,pc::PC, coarseOperators::Vector{PetscMat}) 
Gets coarse operator matrices for all levels (except the finest level)

Logically Collective

Input Parameter:
- `pc` - the precondition context

Output Parameters:
- `num_levels`      - the number of levels
- `coarseOperators` - the coarse operator matrices (size of `num_levels`-1)

Level: advanced

-seealso: [](ch_ksp), `PC`, `PCMG`, `PCMGGetRestriction()`, `PCMGSetInterpolation()`, `PCMGGetRScale()`, `PCMGGetInterpolation()`, `PCGetInterpolations()`

# External Links
$(_doc_external("Ksp/PCGetCoarseOperators"))
"""
function PCGetCoarseOperators(petsclib::PetscLibType, pc::PC, coarseOperators::Vector{PetscMat}) end

@for_petsc function PCGetCoarseOperators(petsclib::$UnionPetscLib, pc::PC, coarseOperators::Vector{PetscMat} )
	num_levels_ = Ref{$PetscInt}()
	coarseOperators_ = Ref(pointer(coarseOperators))

    @chk ccall(
               (:PCGetCoarseOperators, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{$PetscInt}, Ptr{Ptr{CMat}}),
               pc, num_levels_, coarseOperators_,
              )

	num_levels = num_levels_[]

	return num_levels
end 

"""
	PCMatSetApplyOperation(petsclib::PetscLibType,pc::PC, matop::MatOperation) 
Set which matrix operation of the matrix implements `PCApply()` for `PCMAT`.

Logically collective

Input Parameters:
- `pc`    - An instance of `PCMAT`
- `matop` - The selected `MatOperation`

Level: intermediate

-seealso: [](ch_ksp), `PCMAT`, `PCMatGetApplyOperation()`, `PCApply()`, `MatOperation`

# External Links
$(_doc_external("Ksp/PCMatSetApplyOperation"))
"""
function PCMatSetApplyOperation(petsclib::PetscLibType, pc::PC, matop::MatOperation) end

@for_petsc function PCMatSetApplyOperation(petsclib::$UnionPetscLib, pc::PC, matop::MatOperation )

    @chk ccall(
               (:PCMatSetApplyOperation, $petsc_library),
               PetscErrorCode,
               (PC, MatOperation),
               pc, matop,
              )


	return nothing
end 

"""
	PCMatGetApplyOperation(petsclib::PetscLibType,pc::PC, matop::MatOperation) 
Get which matrix operation of the matrix implements `PCApply()` for `PCMAT`.

Logically collective

Input Parameter:
- `pc` - An instance of `PCMAT`

Output Parameter:
- `matop` - The `MatOperation`

Level: intermediate

-seealso: [](ch_ksp), `PCMAT`, `PCMatSetApplyOperation()`, `PCApply()`, `MatOperation`

# External Links
$(_doc_external("Ksp/PCMatGetApplyOperation"))
"""
function PCMatGetApplyOperation(petsclib::PetscLibType, pc::PC, matop::MatOperation) end

@for_petsc function PCMatGetApplyOperation(petsclib::$UnionPetscLib, pc::PC, matop::MatOperation )

    @chk ccall(
               (:PCMatGetApplyOperation, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{MatOperation}),
               pc, matop,
              )


	return nothing
end 

"""
	PCASMSetLocalSubdomains(petsclib::PetscLibType,pc::PC, n::PetscInt, is::Vector{IS}, is_loc::Vector{IS}) 
Sets the local subdomains (for this processor only) for the additive Schwarz preconditioner `PCASM`.

Collective

Input Parameters:
- `pc`       - the preconditioner context
- `n`        - the number of subdomains for this processor (default value = 1)
- `is`       - the index set that defines the subdomains for this processor (or `NULL` for PETSc to determine subdomains)
the values of the `is` array are copied so you can free the array (not the `IS` in the array) after this call
- `is_local` - the index sets that define the local part of the subdomains for this processor, not used unless `PCASMType` is `PC_ASM_RESTRICT`
(or `NULL` to not provide these). The values of the `is_local` array are copied so you can free the array
(not the `IS` in the array) after this call

Options Database Key:
- `-pc_asm_local_blocks <blks>` - Sets number of local blocks

Level: advanced

-seealso: [](ch_ksp), `PCASM`, `PCASMSetTotalSubdomains()`, `PCASMSetOverlap()`, `PCASMGetSubKSP()`,
`PCASMCreateSubdomains2D()`, `PCASMGetLocalSubdomains()`, `PCASMType`, `PCASMSetType()`, `PCGASM`

# External Links
$(_doc_external("Ksp/PCASMSetLocalSubdomains"))
"""
function PCASMSetLocalSubdomains(petsclib::PetscLibType, pc::PC, n::PetscInt, is::Vector{IS}, is_loc::Vector{IS}) end

@for_petsc function PCASMSetLocalSubdomains(petsclib::$UnionPetscLib, pc::PC, n::$PetscInt, is::Vector{IS}, is_loc::Vector{IS} )

    @chk ccall(
               (:PCASMSetLocalSubdomains, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, Ptr{IS}, Ptr{IS}),
               pc, n, is, is_loc,
              )


	return nothing
end 

"""
	PCASMSetTotalSubdomains(petsclib::PetscLibType,pc::PC, N::PetscInt, is::Vector{IS}, is_loc::Vector{IS}) 
Sets the subdomains for all processors for the
additive Schwarz preconditioner, `PCASM`.

Collective, all MPI ranks must pass in the same array of `IS`

Input Parameters:
- `pc`       - the preconditioner context
- `N`        - the number of subdomains for all processors
- `is`       - the index sets that define the subdomains for all processors (or `NULL` to ask PETSc to determine the subdomains)
the values of the `is` array are copied so you can free the array (not the `IS` in the array) after this call
- `is_local` - the index sets that define the local part of the subdomains for this processor (or `NULL` to not provide this information)
The values of the `is_local` array are copied so you can free the array (not the `IS` in the array) after this call

Options Database Key:
- `-pc_asm_blocks <blks>` - Sets total blocks

Level: advanced

-seealso: [](ch_ksp), `PCASM`, `PCASMSetLocalSubdomains()`, `PCASMSetOverlap()`, `PCASMGetSubKSP()`,
`PCASMCreateSubdomains2D()`, `PCGASM`

# External Links
$(_doc_external("Ksp/PCASMSetTotalSubdomains"))
"""
function PCASMSetTotalSubdomains(petsclib::PetscLibType, pc::PC, N::PetscInt, is::Vector{IS}, is_loc::Vector{IS}) end

@for_petsc function PCASMSetTotalSubdomains(petsclib::$UnionPetscLib, pc::PC, N::$PetscInt, is::Vector{IS}, is_loc::Vector{IS} )

    @chk ccall(
               (:PCASMSetTotalSubdomains, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, Ptr{IS}, Ptr{IS}),
               pc, N, is, is_loc,
              )


	return nothing
end 

"""
	PCASMSetOverlap(petsclib::PetscLibType,pc::PC, ovl::PetscInt) 
Sets the overlap between a pair of subdomains for the
additive Schwarz preconditioner, `PCASM`.

Logically Collective

Input Parameters:
- `pc`  - the preconditioner context
- `ovl` - the amount of overlap between subdomains (ovl >= 0, default value = 1)

Options Database Key:
- `-pc_asm_overlap <ovl>` - Sets overlap

Level: intermediate

-seealso: [](ch_ksp), `PCASM`, `PCASMSetTotalSubdomains()`, `PCASMSetLocalSubdomains()`, `PCASMGetSubKSP()`,
`PCASMCreateSubdomains2D()`, `PCASMGetLocalSubdomains()`, `MatIncreaseOverlap()`, `PCGASM`

# External Links
$(_doc_external("Ksp/PCASMSetOverlap"))
"""
function PCASMSetOverlap(petsclib::PetscLibType, pc::PC, ovl::PetscInt) end

@for_petsc function PCASMSetOverlap(petsclib::$UnionPetscLib, pc::PC, ovl::$PetscInt )

    @chk ccall(
               (:PCASMSetOverlap, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, ovl,
              )


	return nothing
end 

"""
	PCASMSetType(petsclib::PetscLibType,pc::PC, type::PCASMType) 
Sets the type of restriction and interpolation used
for local problems in the additive Schwarz method, `PCASM`.

Logically Collective

Input Parameters:
- `pc`   - the preconditioner context
- `type` - variant of `PCASM`, one of
-seealso: [](ch_ksp), `PCASM`, `PCASMSetTotalSubdomains()`, `PCASMGetSubKSP()`,
`PCASMCreateSubdomains2D()`, `PCASMType`, `PCASMSetLocalType()`, `PCASMGetLocalType()`, `PCGASM`

# External Links
$(_doc_external("Ksp/PCASMSetType"))
"""
function PCASMSetType(petsclib::PetscLibType, pc::PC, type::PCASMType) end

@for_petsc function PCASMSetType(petsclib::$UnionPetscLib, pc::PC, type::PCASMType )

    @chk ccall(
               (:PCASMSetType, $petsc_library),
               PetscErrorCode,
               (PC, PCASMType),
               pc, type,
              )


	return nothing
end 

"""
	type::PCASMType = PCASMGetType(petsclib::PetscLibType,pc::PC) 
Gets the type of restriction and interpolation used
for local problems in the additive Schwarz method, `PCASM`.

Logically Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `type` - variant of `PCASM`, one of
-seealso: [](ch_ksp), `PCASM`, `PCASMSetTotalSubdomains()`, `PCASMGetSubKSP()`, `PCGASM`,
`PCASMCreateSubdomains2D()`, `PCASMType`, `PCASMSetType()`, `PCASMSetLocalType()`, `PCASMGetLocalType()`

# External Links
$(_doc_external("Ksp/PCASMGetType"))
"""
function PCASMGetType(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCASMGetType(petsclib::$UnionPetscLib, pc::PC )
	type_ = Ref{PCASMType}()

    @chk ccall(
               (:PCASMGetType, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PCASMType}),
               pc, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	PCASMSetLocalType(petsclib::PetscLibType,pc::PC, type::PCCompositeType) 
Sets the type of composition used for local problems in the additive Schwarz method, `PCASM`.

Logically Collective

Input Parameters:
- `pc`   - the preconditioner context
- `type` - type of composition, one of
-seealso: [](ch_ksp), `PCASM`, `PCASMSetType()`, `PCASMGetType()`, `PCASMGetLocalType()`, `PCASMType`, `PCCompositeType`

# External Links
$(_doc_external("Ksp/PCASMSetLocalType"))
"""
function PCASMSetLocalType(petsclib::PetscLibType, pc::PC, type::PCCompositeType) end

@for_petsc function PCASMSetLocalType(petsclib::$UnionPetscLib, pc::PC, type::PCCompositeType )

    @chk ccall(
               (:PCASMSetLocalType, $petsc_library),
               PetscErrorCode,
               (PC, PCCompositeType),
               pc, type,
              )


	return nothing
end 

"""
	type::PCCompositeType = PCASMGetLocalType(petsclib::PetscLibType,pc::PC) 
Gets the type of composition used for local problems in the additive Schwarz method, `PCASM`.

Logically Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `type` - type of composition, one of
-seealso: [](ch_ksp), `PCASM`, `PCASMSetType()`, `PCASMGetType()`, `PCASMSetLocalType()`, `PCASMType`, `PCCompositeType`

# External Links
$(_doc_external("Ksp/PCASMGetLocalType"))
"""
function PCASMGetLocalType(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCASMGetLocalType(petsclib::$UnionPetscLib, pc::PC )
	type_ = Ref{PCCompositeType}()

    @chk ccall(
               (:PCASMGetLocalType, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PCCompositeType}),
               pc, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	PCASMSetSortIndices(petsclib::PetscLibType,pc::PC, doSort::PetscBool) 
Determines whether subdomain indices are sorted.

Logically Collective

Input Parameters:
- `pc`     - the preconditioner context
- `doSort` - sort the subdomain indices

Level: intermediate

-seealso: [](ch_ksp), `PCASM`, `PCASMSetLocalSubdomains()`, `PCASMSetTotalSubdomains()`, `PCASMGetSubKSP()`,
`PCASMCreateSubdomains2D()`

# External Links
$(_doc_external("Ksp/PCASMSetSortIndices"))
"""
function PCASMSetSortIndices(petsclib::PetscLibType, pc::PC, doSort::PetscBool) end

@for_petsc function PCASMSetSortIndices(petsclib::$UnionPetscLib, pc::PC, doSort::PetscBool )

    @chk ccall(
               (:PCASMSetSortIndices, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, doSort,
              )


	return nothing
end 

"""
	n_loc::PetscInt,first_loc::PetscInt = PCASMGetSubKSP(petsclib::PetscLibType,pc::PC, ksp::Vector{PetscKSP}) 
Gets the local `KSP` contexts for all blocks on
this processor.

Collective iff first_local is requested

Input Parameter:
- `pc` - the preconditioner context

Output Parameters:
- `n_local`     - the number of blocks on this processor or `NULL`
- `first_local` - the global number of the first block on this processor or `NULL`, all processors must request or all must pass `NULL`
- `ksp`         - the array of `KSP` contexts

Level: advanced

-seealso: [](ch_ksp), `PCASM`, `PCASMSetTotalSubdomains()`, `PCASMSetOverlap()`,
`PCASMCreateSubdomains2D()`,

# External Links
$(_doc_external("Ksp/PCASMGetSubKSP"))
"""
function PCASMGetSubKSP(petsclib::PetscLibType, pc::PC, ksp::Vector{PetscKSP}) end

@for_petsc function PCASMGetSubKSP(petsclib::$UnionPetscLib, pc::PC, ksp::Vector{PetscKSP} )
	n_loc_ = Ref{$PetscInt}()
	first_loc_ = Ref{$PetscInt}()
	ksp_ = Ref(pointer(ksp))

    @chk ccall(
               (:PCASMGetSubKSP, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{CKSP}}),
               pc, n_loc_, first_loc_, ksp_,
              )

	n_loc = n_loc_[]
	first_loc = first_loc_[]

	return n_loc,first_loc
end 

"""
	outis::Vector{IS} = PCASMCreateSubdomains(petsclib::PetscLibType,A::PetscMat, n::PetscInt) 
Creates the index sets for the overlapping Schwarz
preconditioner, `PCASM`,  for any problem on a general grid.

Collective

Input Parameters:
- `A` - The global matrix operator
- `n` - the number of local blocks

Output Parameter:
- `outis` - the array of index sets defining the subdomains

Level: advanced

-seealso: [](ch_ksp), `PCASM`, `PCASMSetLocalSubdomains()`, `PCASMDestroySubdomains()`

# External Links
$(_doc_external("Ksp/PCASMCreateSubdomains"))
"""
function PCASMCreateSubdomains(petsclib::PetscLibType, A::PetscMat, n::PetscInt) end

@for_petsc function PCASMCreateSubdomains(petsclib::$UnionPetscLib, A::PetscMat, n::$PetscInt )
	outis_ = Ref{Ptr{IS}}()

    @chk ccall(
               (:PCASMCreateSubdomains, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{Ptr{IS}}),
               A, n, outis_,
              )

	outis = unsafe_wrap(Array, outis_[], VecGetLocalSize(petsclib, x); own = false)

	return outis
end 

"""
	PCASMDestroySubdomains(petsclib::PetscLibType,n::PetscInt, is::Vector{IS}, is_loc::Vector{IS}) 
Destroys the index sets created with
`PCASMCreateSubdomains()`. Should be called after setting subdomains with `PCASMSetLocalSubdomains()`.

Collective

Input Parameters:
- `n`        - the number of index sets
- `is`       - the array of index sets
- `is_local` - the array of local index sets, can be `NULL`

Level: advanced

-seealso: [](ch_ksp), `PCASM`, `PCASMCreateSubdomains()`, `PCASMSetLocalSubdomains()`

# External Links
$(_doc_external("Ksp/PCASMDestroySubdomains"))
"""
function PCASMDestroySubdomains(petsclib::PetscLibType, n::PetscInt, is::Vector{IS}, is_loc::Vector{IS}) end

@for_petsc function PCASMDestroySubdomains(petsclib::$UnionPetscLib, n::$PetscInt, is::Vector{IS}, is_loc::Vector{IS} )
	is_ = Ref(pointer(is))
	is_loc_ = Ref(pointer(is_loc))

    @chk ccall(
               (:PCASMDestroySubdomains, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{Ptr{IS}}, Ptr{Ptr{IS}}),
               n, is_, is_loc_,
              )


	return nothing
end 

"""
	Nsub::PetscInt,is::Vector{IS},is_loc::Vector{IS} = PCASMCreateSubdomains2D(petsclib::PetscLibType,m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, dof::PetscInt, overlap::PetscInt) 
Creates the index sets for the overlapping Schwarz
preconditioner, `PCASM`, for a two-dimensional problem on a regular grid.

Not Collective

Input Parameters:
- `m`       - the number of mesh points in the x direction
- `n`       - the number of mesh points in the y direction
- `M`       - the number of subdomains in the x direction
- `N`       - the number of subdomains in the y direction
- `dof`     - degrees of freedom per node
- `overlap` - overlap in mesh lines

Output Parameters:
- `Nsub`     - the number of subdomains created
- `is`       - array of index sets defining overlapping (if overlap > 0) subdomains
- `is_local` - array of index sets defining non-overlapping subdomains

Level: advanced

-seealso: [](ch_ksp), `PCASM`, `PCASMSetTotalSubdomains()`, `PCASMSetLocalSubdomains()`, `PCASMGetSubKSP()`,
`PCASMSetOverlap()`

# External Links
$(_doc_external("Ksp/PCASMCreateSubdomains2D"))
"""
function PCASMCreateSubdomains2D(petsclib::PetscLibType, m::PetscInt, n::PetscInt, M::PetscInt, N::PetscInt, dof::PetscInt, overlap::PetscInt) end

@for_petsc function PCASMCreateSubdomains2D(petsclib::$UnionPetscLib, m::$PetscInt, n::$PetscInt, M::$PetscInt, N::$PetscInt, dof::$PetscInt, overlap::$PetscInt )
	Nsub_ = Ref{$PetscInt}()
	is_ = Ref{Ptr{IS}}()
	is_loc_ = Ref{Ptr{IS}}()

    @chk ccall(
               (:PCASMCreateSubdomains2D, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{Ptr{IS}}, Ptr{Ptr{IS}}),
               m, n, M, N, dof, overlap, Nsub_, is_, is_loc_,
              )

	Nsub = Nsub_[]
	is = unsafe_wrap(Array, is_[], VecGetLocalSize(petsclib, x); own = false)
	is_loc = unsafe_wrap(Array, is_loc_[], VecGetLocalSize(petsclib, x); own = false)

	return Nsub,is,is_loc
end 

"""
	n::PetscInt = PCASMGetLocalSubdomains(petsclib::PetscLibType,pc::PC, is::Vector{IS}, is_loc::Vector{IS}) 
Gets the local subdomains (for this processor
only) for the additive Schwarz preconditioner, `PCASM`.

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameters:
- `n`        - if requested, the number of subdomains for this processor (default value = 1)
- `is`       - if requested, the index sets that define the subdomains for this processor
- `is_local` - if requested, the index sets that define the local part of the subdomains for this processor (can be `NULL`)

Level: advanced

-seealso: [](ch_ksp), `PCASM`, `PCASMSetTotalSubdomains()`, `PCASMSetOverlap()`, `PCASMGetSubKSP()`,
`PCASMCreateSubdomains2D()`, `PCASMSetLocalSubdomains()`, `PCASMGetLocalSubmatrices()`

# External Links
$(_doc_external("Ksp/PCASMGetLocalSubdomains"))
"""
function PCASMGetLocalSubdomains(petsclib::PetscLibType, pc::PC, is::Vector{IS}, is_loc::Vector{IS}) end

@for_petsc function PCASMGetLocalSubdomains(petsclib::$UnionPetscLib, pc::PC, is::Vector{IS}, is_loc::Vector{IS} )
	n_ = Ref{$PetscInt}()
	is_ = Ref(pointer(is))
	is_loc_ = Ref(pointer(is_loc))

    @chk ccall(
               (:PCASMGetLocalSubdomains, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{$PetscInt}, Ptr{Ptr{IS}}, Ptr{Ptr{IS}}),
               pc, n_, is_, is_loc_,
              )

	n = n_[]

	return n
end 

"""
	n::PetscInt = PCASMGetLocalSubmatrices(petsclib::PetscLibType,pc::PC, mat::Vector{PetscMat}) 
Gets the local submatrices (for this processor
only) for the additive Schwarz preconditioner, `PCASM`.

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameters:
- `n`   - if requested, the number of matrices for this processor (default value = 1)
- `mat` - if requested, the matrices

Level: advanced

-seealso: [](ch_ksp), `PCASM`, `PCASMSetTotalSubdomains()`, `PCASMSetOverlap()`, `PCASMGetSubKSP()`,
`PCASMCreateSubdomains2D()`, `PCASMSetLocalSubdomains()`, `PCASMGetLocalSubdomains()`, `PCSetModifySubMatrices()`

# External Links
$(_doc_external("Ksp/PCASMGetLocalSubmatrices"))
"""
function PCASMGetLocalSubmatrices(petsclib::PetscLibType, pc::PC, mat::Vector{PetscMat}) end

@for_petsc function PCASMGetLocalSubmatrices(petsclib::$UnionPetscLib, pc::PC, mat::Vector{PetscMat} )
	n_ = Ref{$PetscInt}()
	mat_ = Ref(pointer(mat))

    @chk ccall(
               (:PCASMGetLocalSubmatrices, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{$PetscInt}, Ptr{Ptr{CMat}}),
               pc, n_, mat_,
              )

	n = n_[]

	return n
end 

"""
	PCASMSetDMSubdomains(petsclib::PetscLibType,pc::PC, flg::PetscBool) 
Indicates whether to use `DMCreateDomainDecomposition()` to define the subdomains, whenever possible.

Logically Collective

Input Parameters:
- `pc`  - the preconditioner
- `flg` - boolean indicating whether to use subdomains defined by the `DM`

Options Database Key:
- `-pc_asm_dm_subdomains <bool>` - use subdomains defined by the `DM` with `DMCreateDomainDecomposition()`

Level: intermediate

-seealso: [](ch_ksp), `PCASM`, `PCASMGetDMSubdomains()`, `PCASMSetTotalSubdomains()`, `PCASMSetOverlap()`
`PCASMCreateSubdomains2D()`, `PCASMSetLocalSubdomains()`, `PCASMGetLocalSubdomains()`

# External Links
$(_doc_external("Ksp/PCASMSetDMSubdomains"))
"""
function PCASMSetDMSubdomains(petsclib::PetscLibType, pc::PC, flg::PetscBool) end

@for_petsc function PCASMSetDMSubdomains(petsclib::$UnionPetscLib, pc::PC, flg::PetscBool )

    @chk ccall(
               (:PCASMSetDMSubdomains, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, flg,
              )


	return nothing
end 

"""
	flg::PetscBool = PCASMGetDMSubdomains(petsclib::PetscLibType,pc::PC) 
Returns flag indicating whether to use `DMCreateDomainDecomposition()` to define the subdomains, whenever possible.

Not Collective

Input Parameter:
- `pc` - the preconditioner

Output Parameter:
- `flg` - boolean indicating whether to use subdomains defined by the `DM`

Level: intermediate

-seealso: [](ch_ksp), `PCASM`, `PCASMSetDMSubdomains()`, `PCASMSetTotalSubdomains()`, `PCASMSetOverlap()`
`PCASMCreateSubdomains2D()`, `PCASMSetLocalSubdomains()`, `PCASMGetLocalSubdomains()`

# External Links
$(_doc_external("Ksp/PCASMGetDMSubdomains"))
"""
function PCASMGetDMSubdomains(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCASMGetDMSubdomains(petsclib::$UnionPetscLib, pc::PC )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PCASMGetDMSubdomains, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PetscBool}),
               pc, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	sub_mat_type::MatType = PCASMGetSubMatType(petsclib::PetscLibType,pc::PC) 
Gets the matrix type used for `PCASM` subsolves, as a string.

Not Collective

Input Parameter:
- `pc` - the `PC`

Output Parameter:
- `sub_mat_type` - name of matrix type

Level: advanced

-seealso: [](ch_ksp), `PCASM`, `PCASMSetSubMatType()`, `PCSetType()`, `VecSetType()`, `MatType`, `Mat`

# External Links
$(_doc_external("Ksp/PCASMGetSubMatType"))
"""
function PCASMGetSubMatType(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCASMGetSubMatType(petsclib::$UnionPetscLib, pc::PC )
	sub_mat_type_ = Ref{MatType}()

    @chk ccall(
               (:PCASMGetSubMatType, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{MatType}),
               pc, sub_mat_type_,
              )

	sub_mat_type = unsafe_string(sub_mat_type_[])

	return sub_mat_type
end 

"""
	PCASMSetSubMatType(petsclib::PetscLibType,pc::PC, sub_mat_type::MatType) 
Set the type of matrix used for `PCASM` subsolves

Collective

Input Parameters:
- `pc`           - the `PC` object
- `sub_mat_type` - the `MatType`

Options Database Key:
- `-pc_asm_sub_mat_type  <sub_mat_type>` - Sets the matrix type used for subsolves, for example, seqaijviennacl.
If you specify a base name like aijviennacl, the corresponding sequential type is assumed.

-seealso: [](ch_ksp), `PCASM`, `PCASMGetSubMatType()`, `PCSetType()`, `VecSetType()`, `MatType`, `Mat`

# External Links
$(_doc_external("Ksp/PCASMSetSubMatType"))
"""
function PCASMSetSubMatType(petsclib::PetscLibType, pc::PC, sub_mat_type::MatType) end

@for_petsc function PCASMSetSubMatType(petsclib::$UnionPetscLib, pc::PC, sub_mat_type::MatType )

    @chk ccall(
               (:PCASMSetSubMatType, $petsc_library),
               PetscErrorCode,
               (PC, MatType),
               pc, sub_mat_type,
              )


	return nothing
end 

"""
	PCMPIServerBegin(petsclib::PetscLibType) 
starts a server that runs on the `rank != 0` MPI processes waiting to process requests for
parallel `KSP` solves and management of parallel `KSP` objects.

Logically Collective on all MPI processes except rank 0

Options Database Keys:
- `-mpi_linear_solver_server`                   - causes the PETSc program to start in MPI linear solver server mode where only the first MPI rank runs user code
- `-mpi_linear_solver_server_view`              - displays information about all the linear systems solved by the MPI linear solver server at the conclusion of the program
- `-mpi_linear_solver_server_use_shared_memory` - use shared memory when communicating matrices and vectors to server processes (default where supported)

Level: developer

-seealso: [](sec_pcmpi), `PCMPIServerEnd()`, `PCMPI`, `KSPCheckPCMPI()`

# External Links
$(_doc_external("Ksp/PCMPIServerBegin"))
"""
function PCMPIServerBegin(petsclib::PetscLibType) end

@for_petsc function PCMPIServerBegin(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PCMPIServerBegin, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PCMPIServerEnd(petsclib::PetscLibType) 
ends a server that runs on the rank != 0 MPI processes waiting to process requests for
parallel KSP solves and management of parallel `KSP` objects.

Logically Collective on all MPI ranks except 0

Level: developer

-seealso: [](sec_pcmpi), `PCMPIServerBegin()`, `PCMPI`, `KSPCheckPCMPI()`

# External Links
$(_doc_external("Ksp/PCMPIServerEnd"))
"""
function PCMPIServerEnd(petsclib::PetscLibType) end

@for_petsc function PCMPIServerEnd(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PCMPIServerEnd, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PCMPIGetKSP(petsclib::PetscLibType,pc::PC, innerksp::PetscKSP) 
Gets the `KSP` created by the `PCMPI`

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `innerksp` - the inner `KSP`

Level: advanced

-seealso: [](ch_ksp), `KSP`, `PCMPI`, `PCREDISTRIBUTE`

# External Links
$(_doc_external("Ksp/PCMPIGetKSP"))
"""
function PCMPIGetKSP(petsclib::PetscLibType, pc::PC, innerksp::PetscKSP) end

@for_petsc function PCMPIGetKSP(petsclib::$UnionPetscLib, pc::PC, innerksp::PetscKSP )
	innerksp_ = Ref(innerksp.ptr)

    @chk ccall(
               (:PCMPIGetKSP, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{CKSP}),
               pc, innerksp_,
              )

	innerksp.ptr = C_NULL

	return nothing
end 

"""
	PCPythonSetType(petsclib::PetscLibType,pc::PC, pyname::String) 
Initialize a `PC` object implemented in Python, a `PCPYTHON`.

Collective

Input Parameters:
- `pc`  - the preconditioner (`PC`) context.
- `pyname`  - full dotted Python name [package].module[.{class|function}]

Options Database Key:
- `-pc_python_type <pyname>`  - python class

Level: intermediate

-seealso: [](ch_ksp), `PC`, `PCSHELL`, `PCCreate()`, `PCSetType()`, `PCPYTHON`, `PetscPythonInitialize()`

# External Links
$(_doc_external("Ksp/PCPythonSetType"))
"""
function PCPythonSetType(petsclib::PetscLibType, pc::PC, pyname::String) end

@for_petsc function PCPythonSetType(petsclib::$UnionPetscLib, pc::PC, pyname::String )

    @chk ccall(
               (:PCPythonSetType, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{Cchar}),
               pc, pyname,
              )


	return nothing
end 

"""
	pyname::String = PCPythonGetType(petsclib::PetscLibType,pc::PC) 
Get the type of a `PC` object implemented in Python, a `PCPYTHON`.

Not Collective

Input Parameter:
- `pc`  - the preconditioner (`PC`) context.

Output Parameter:
- `pyname`  - full dotted Python name [package].module[.{class|function}]

Level: intermediate

-seealso: [](ch_ksp), `PC`, `PCSHELL`, `PCCreate()`, `PCSetType()`, `PCPYTHON`, `PetscPythonInitialize()`, `PCPythonSetType()`

# External Links
$(_doc_external("Ksp/PCPythonGetType"))
"""
function PCPythonGetType(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCPythonGetType(petsclib::$UnionPetscLib, pc::PC )
	pyname_ = Ref{Ptr{Cchar}}()

    @chk ccall(
               (:PCPythonGetType, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{Ptr{Cchar}}),
               pc, pyname_,
              )

	pyname = unsafe_wrap(Array, pyname_[], VecGetLocalSize(petsclib, x); own = false)

	return pyname
end 

"""
	PCKSPSetKSP(petsclib::PetscLibType,pc::PC, ksp::PetscKSP) 
Sets the `KSP` context for a `PCKSP`.

Collective

Input Parameters:
- `pc`  - the preconditioner context
- `ksp` - the `KSP` solver

Level: advanced

-seealso: [](ch_ksp), `PCKSP`, `PCKSPGetKSP()`

# External Links
$(_doc_external("Ksp/PCKSPSetKSP"))
"""
function PCKSPSetKSP(petsclib::PetscLibType, pc::PC, ksp::PetscKSP) end

@for_petsc function PCKSPSetKSP(petsclib::$UnionPetscLib, pc::PC, ksp::PetscKSP )

    @chk ccall(
               (:PCKSPSetKSP, $petsc_library),
               PetscErrorCode,
               (PC, CKSP),
               pc, ksp,
              )


	return nothing
end 

"""
	PCKSPGetKSP(petsclib::PetscLibType,pc::PC, ksp::PetscKSP) 
Gets the `KSP` context for a `PCKSP`.

Not Collective but ksp returned is parallel if pc was parallel

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `ksp` - the `KSP` solver

-seealso: [](ch_ksp), `PCKSP`, `PCKSPSetKSP()`

# External Links
$(_doc_external("Ksp/PCKSPGetKSP"))
"""
function PCKSPGetKSP(petsclib::PetscLibType, pc::PC, ksp::PetscKSP) end

@for_petsc function PCKSPGetKSP(petsclib::$UnionPetscLib, pc::PC, ksp::PetscKSP )
	ksp_ = Ref(ksp.ptr)

    @chk ccall(
               (:PCKSPGetKSP, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{CKSP}),
               pc, ksp_,
              )

	ksp.ptr = C_NULL

	return nothing
end 

"""
	PCMGResidualDefault(petsclib::PetscLibType,mat::PetscMat, b::PetscVec, x::PetscVec, r::PetscVec) 
Default routine to calculate the residual.

Collective

Input Parameters:
- `mat` - the matrix
- `b`   - the right-hand side
- `x`   - the approximate solution

Output Parameter:
- `r` - location to store the residual

Level: developer

-seealso: [](ch_ksp), `PCMG`, `PCMGSetResidual()`, `PCMGSetMatResidual()`

# External Links
$(_doc_external("Ksp/PCMGResidualDefault"))
"""
function PCMGResidualDefault(petsclib::PetscLibType, mat::PetscMat, b::PetscVec, x::PetscVec, r::PetscVec) end

@for_petsc function PCMGResidualDefault(petsclib::$UnionPetscLib, mat::PetscMat, b::PetscVec, x::PetscVec, r::PetscVec )

    @chk ccall(
               (:PCMGResidualDefault, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec, CVec),
               mat, b, x, r,
              )


	return nothing
end 

"""
	PCMGResidualTransposeDefault(petsclib::PetscLibType,mat::PetscMat, b::PetscVec, x::PetscVec, r::PetscVec) 
Default routine to calculate the residual of the transposed linear system

Collective

Input Parameters:
- `mat` - the matrix
- `b`   - the right-hand side
- `x`   - the approximate solution

Output Parameter:
- `r` - location to store the residual

Level: developer

-seealso: [](ch_ksp), `PCMG`, `PCMGSetResidualTranspose()`, `PCMGMatResidualTransposeDefault()`

# External Links
$(_doc_external("Ksp/PCMGResidualTransposeDefault"))
"""
function PCMGResidualTransposeDefault(petsclib::PetscLibType, mat::PetscMat, b::PetscVec, x::PetscVec, r::PetscVec) end

@for_petsc function PCMGResidualTransposeDefault(petsclib::$UnionPetscLib, mat::PetscMat, b::PetscVec, x::PetscVec, r::PetscVec )

    @chk ccall(
               (:PCMGResidualTransposeDefault, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec, CVec),
               mat, b, x, r,
              )


	return nothing
end 

"""
	PCMGMatResidualDefault(petsclib::PetscLibType,mat::PetscMat, b::PetscMat, x::PetscMat, r::PetscMat) 
Default routine to calculate the residual.

Collective

Input Parameters:
- `mat` - the matrix
- `b`   - the right-hand side
- `x`   - the approximate solution

Output Parameter:
- `r` - location to store the residual

Level: developer

-seealso: [](ch_ksp), `PCMG`, `PCMGSetMatResidual()`, `PCMGResidualDefault()`

# External Links
$(_doc_external("Ksp/PCMGMatResidualDefault"))
"""
function PCMGMatResidualDefault(petsclib::PetscLibType, mat::PetscMat, b::PetscMat, x::PetscMat, r::PetscMat) end

@for_petsc function PCMGMatResidualDefault(petsclib::$UnionPetscLib, mat::PetscMat, b::PetscMat, x::PetscMat, r::PetscMat )

    @chk ccall(
               (:PCMGMatResidualDefault, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat, CMat),
               mat, b, x, r,
              )


	return nothing
end 

"""
	PCMGMatResidualTransposeDefault(petsclib::PetscLibType,mat::PetscMat, b::PetscMat, x::PetscMat, r::PetscMat) 
Default routine to calculate the residual of the transposed linear system

Collective

Input Parameters:
- `mat` - the matrix
- `b`   - the right-hand side
- `x`   - the approximate solution

Output Parameter:
- `r` - location to store the residual

Level: developer

-seealso: [](ch_ksp), `PCMG`, `PCMGSetMatResidualTranspose()`

# External Links
$(_doc_external("Ksp/PCMGMatResidualTransposeDefault"))
"""
function PCMGMatResidualTransposeDefault(petsclib::PetscLibType, mat::PetscMat, b::PetscMat, x::PetscMat, r::PetscMat) end

@for_petsc function PCMGMatResidualTransposeDefault(petsclib::$UnionPetscLib, mat::PetscMat, b::PetscMat, x::PetscMat, r::PetscMat )

    @chk ccall(
               (:PCMGMatResidualTransposeDefault, $petsc_library),
               PetscErrorCode,
               (CMat, CMat, CMat, CMat),
               mat, b, x, r,
              )


	return nothing
end 

"""
	PCMGGetCoarseSolve(petsclib::PetscLibType,pc::PC, ksp::PetscKSP) 
Gets the solver context to be used on the coarse grid.

Not Collective

Input Parameter:
- `pc` - the multigrid context

Output Parameter:
- `ksp` - the coarse grid solver context

Level: advanced

-seealso: [](ch_ksp), `PCMG`, `PCMGGetSmootherUp()`, `PCMGGetSmootherDown()`, `PCMGGetSmoother()`

# External Links
$(_doc_external("Ksp/PCMGGetCoarseSolve"))
"""
function PCMGGetCoarseSolve(petsclib::PetscLibType, pc::PC, ksp::PetscKSP) end

@for_petsc function PCMGGetCoarseSolve(petsclib::$UnionPetscLib, pc::PC, ksp::PetscKSP )
	ksp_ = Ref(ksp.ptr)

    @chk ccall(
               (:PCMGGetCoarseSolve, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{CKSP}),
               pc, ksp_,
              )

	ksp.ptr = C_NULL

	return nothing
end 

"""
	PCMGSetResidual(petsclib::PetscLibType,pc::PC, l::PetscInt, residual::external, mat::PetscMat) 
Sets the function to be used to calculate the residual on the lth level.

Logically Collective

Input Parameters:
- `pc`       - the multigrid context
- `l`        - the level (0 is coarsest) to supply
- `residual` - function used to form residual, if none is provided the previously provide one is used, if no
previous one were provided then a default is used
- `mat`      - matrix associated with residual

Level: advanced

-seealso: [](ch_ksp), `PCMG`, `PCMGResidualDefault()`

# External Links
$(_doc_external("Ksp/PCMGSetResidual"))
"""
function PCMGSetResidual(petsclib::PetscLibType, pc::PC, l::PetscInt, residual::external, mat::PetscMat) end

@for_petsc function PCMGSetResidual(petsclib::$UnionPetscLib, pc::PC, l::$PetscInt, residual::external, mat::PetscMat )

    @chk ccall(
               (:PCMGSetResidual, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, external, CMat),
               pc, l, residual, mat,
              )


	return nothing
end 

"""
	PCMGSetResidualTranspose(petsclib::PetscLibType,pc::PC, l::PetscInt, residualt::external, mat::PetscMat) 
Sets the function to be used to calculate the residual of the transposed linear system
on the lth level.

Logically Collective

Input Parameters:
- `pc`        - the multigrid context
- `l`         - the level (0 is coarsest) to supply
- `residualt` - function used to form transpose of residual, if none is provided the previously provide one is used, if no
previous one were provided then a default is used
- `mat`       - matrix associated with residual

Level: advanced

-seealso: [](ch_ksp), `PCMG`, `PCMGResidualTransposeDefault()`

# External Links
$(_doc_external("Ksp/PCMGSetResidualTranspose"))
"""
function PCMGSetResidualTranspose(petsclib::PetscLibType, pc::PC, l::PetscInt, residualt::external, mat::PetscMat) end

@for_petsc function PCMGSetResidualTranspose(petsclib::$UnionPetscLib, pc::PC, l::$PetscInt, residualt::external, mat::PetscMat )

    @chk ccall(
               (:PCMGSetResidualTranspose, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, external, CMat),
               pc, l, residualt, mat,
              )


	return nothing
end 

"""
	PCMGSetInterpolation(petsclib::PetscLibType,pc::PC, l::PetscInt, mat::PetscMat) 
Sets the function to be used to calculate the
interpolation from l-1 to the lth level

Logically Collective

Input Parameters:
- `pc`  - the multigrid context
- `mat` - the interpolation operator
- `l`   - the level (0 is coarsest) to supply [do not supply 0]

Level: advanced

-seealso: [](ch_ksp), `PCMG`, `PCMGSetRestriction()`

# External Links
$(_doc_external("Ksp/PCMGSetInterpolation"))
"""
function PCMGSetInterpolation(petsclib::PetscLibType, pc::PC, l::PetscInt, mat::PetscMat) end

@for_petsc function PCMGSetInterpolation(petsclib::$UnionPetscLib, pc::PC, l::$PetscInt, mat::PetscMat )

    @chk ccall(
               (:PCMGSetInterpolation, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, CMat),
               pc, l, mat,
              )


	return nothing
end 

"""
	PCMGSetOperators(petsclib::PetscLibType,pc::PC, l::PetscInt, Amat::PetscMat, Pmat::PetscMat) 
Sets operator and matrix from which to construct a preconditioner for lth level

Logically Collective

Input Parameters:
- `pc`   - the multigrid context
- `Amat` - the operator
- `Pmat` - the preconditioning operator
- `l`    - the level (0 is the coarsest) to supply

Level: advanced

-seealso: [](ch_ksp), `PCMG`, `PCMGSetGalerkin()`, `PCMGSetRestriction()`, `PCMGSetInterpolation()`

# External Links
$(_doc_external("Ksp/PCMGSetOperators"))
"""
function PCMGSetOperators(petsclib::PetscLibType, pc::PC, l::PetscInt, Amat::PetscMat, Pmat::PetscMat) end

@for_petsc function PCMGSetOperators(petsclib::$UnionPetscLib, pc::PC, l::$PetscInt, Amat::PetscMat, Pmat::PetscMat )

    @chk ccall(
               (:PCMGSetOperators, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, CMat, CMat),
               pc, l, Amat, Pmat,
              )


	return nothing
end 

"""
	PCMGGetInterpolation(petsclib::PetscLibType,pc::PC, l::PetscInt, mat::PetscMat) 
Gets the function to be used to calculate the
interpolation from l-1 to the lth level

Logically Collective

Input Parameters:
- `pc` - the multigrid context
- `l`  - the level (0 is coarsest) to supply [Do not supply 0]

Output Parameter:
- `mat` - the interpolation matrix, can be `NULL`

Level: advanced

-seealso: [](ch_ksp), `PCMG`, `PCMGGetRestriction()`, `PCMGSetInterpolation()`, `PCMGGetRScale()`

# External Links
$(_doc_external("Ksp/PCMGGetInterpolation"))
"""
function PCMGGetInterpolation(petsclib::PetscLibType, pc::PC, l::PetscInt, mat::PetscMat) end

@for_petsc function PCMGGetInterpolation(petsclib::$UnionPetscLib, pc::PC, l::$PetscInt, mat::PetscMat )
	mat_ = Ref(mat.ptr)

    @chk ccall(
               (:PCMGGetInterpolation, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, Ptr{CMat}),
               pc, l, mat_,
              )

	mat.ptr = C_NULL

	return nothing
end 

"""
	PCMGSetRestriction(petsclib::PetscLibType,pc::PC, l::PetscInt, mat::PetscMat) 
Sets the function to be used to restrict dual vectors
from level l to l-1.

Logically Collective

Input Parameters:
- `pc`  - the multigrid context
- `l`   - the level (0 is coarsest) to supply [Do not supply 0]
- `mat` - the restriction matrix

Level: advanced

-seealso: [](ch_ksp), `PCMG`, `PCMGSetInterpolation()`

# External Links
$(_doc_external("Ksp/PCMGSetRestriction"))
"""
function PCMGSetRestriction(petsclib::PetscLibType, pc::PC, l::PetscInt, mat::PetscMat) end

@for_petsc function PCMGSetRestriction(petsclib::$UnionPetscLib, pc::PC, l::$PetscInt, mat::PetscMat )

    @chk ccall(
               (:PCMGSetRestriction, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, CMat),
               pc, l, mat,
              )


	return nothing
end 

"""
	PCMGGetRestriction(petsclib::PetscLibType,pc::PC, l::PetscInt, mat::PetscMat) 
Gets the function to be used to restrict dual (i.e. residual) vectors
from level l to l-1.

Logically Collective

Input Parameters:
- `pc` - the multigrid context
- `l`  - the level (0 is coarsest) to supply [Do not supply 0]

Output Parameter:
- `mat` - the restriction matrix

Level: advanced

-seealso: [](ch_ksp), `PCMG`, `PCMGGetInterpolation()`, `PCMGSetRestriction()`, `PCMGGetRScale()`, `PCMGGetInjection()`

# External Links
$(_doc_external("Ksp/PCMGGetRestriction"))
"""
function PCMGGetRestriction(petsclib::PetscLibType, pc::PC, l::PetscInt, mat::PetscMat) end

@for_petsc function PCMGGetRestriction(petsclib::$UnionPetscLib, pc::PC, l::$PetscInt, mat::PetscMat )
	mat_ = Ref(mat.ptr)

    @chk ccall(
               (:PCMGGetRestriction, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, Ptr{CMat}),
               pc, l, mat_,
              )

	mat.ptr = C_NULL

	return nothing
end 

"""
	PCMGSetRScale(petsclib::PetscLibType,pc::PC, l::PetscInt, rscale::PetscVec) 
Sets the pointwise scaling for the restriction operator from level l to l

Logically Collective

Input Parameters:
- `pc`     - the multigrid context
- `l`      - the level (0 is coarsest) to supply [Do not supply 0]
- `rscale` - the scaling

Level: advanced

-seealso: [](ch_ksp), `PCMG`, `PCMGSetInterpolation()`, `PCMGSetRestriction()`, `PCMGGetRScale()`, `PCMGSetInjection()`

# External Links
$(_doc_external("Ksp/PCMGSetRScale"))
"""
function PCMGSetRScale(petsclib::PetscLibType, pc::PC, l::PetscInt, rscale::PetscVec) end

@for_petsc function PCMGSetRScale(petsclib::$UnionPetscLib, pc::PC, l::$PetscInt, rscale::PetscVec )

    @chk ccall(
               (:PCMGSetRScale, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, CVec),
               pc, l, rscale,
              )


	return nothing
end 

"""
	PCMGGetRScale(petsclib::PetscLibType,pc::PC, l::PetscInt, rscale::PetscVec) 
Gets the pointwise scaling for the restriction operator from level l to l

Collective

Input Parameters:
- `pc`     - the multigrid context
- `rscale` - the scaling
- `l`      - the level (0 is coarsest) to supply [Do not supply 0]

Level: advanced

-seealso: [](ch_ksp), `PCMG`, `PCMGSetInterpolation()`, `PCMGGetRestriction()`, `PCMGGetInjection()`

# External Links
$(_doc_external("Ksp/PCMGGetRScale"))
"""
function PCMGGetRScale(petsclib::PetscLibType, pc::PC, l::PetscInt, rscale::PetscVec) end

@for_petsc function PCMGGetRScale(petsclib::$UnionPetscLib, pc::PC, l::$PetscInt, rscale::PetscVec )
	rscale_ = Ref(rscale.ptr)

    @chk ccall(
               (:PCMGGetRScale, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, Ptr{CVec}),
               pc, l, rscale_,
              )

	rscale.ptr = C_NULL

	return nothing
end 

"""
	PCMGSetInjection(petsclib::PetscLibType,pc::PC, l::PetscInt, mat::PetscMat) 
Sets the function to be used to inject primal (i.e. solution) vectors
from level l to l-1.

Logically Collective

Input Parameters:
- `pc`  - the multigrid context
- `l`   - the level (0 is coarsest) to supply [Do not supply 0]
- `mat` - the injection matrix

Level: advanced

-seealso: [](ch_ksp), `PCMG`, `PCMGSetRestriction()`

# External Links
$(_doc_external("Ksp/PCMGSetInjection"))
"""
function PCMGSetInjection(petsclib::PetscLibType, pc::PC, l::PetscInt, mat::PetscMat) end

@for_petsc function PCMGSetInjection(petsclib::$UnionPetscLib, pc::PC, l::$PetscInt, mat::PetscMat )

    @chk ccall(
               (:PCMGSetInjection, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, CMat),
               pc, l, mat,
              )


	return nothing
end 

"""
	PCMGGetInjection(petsclib::PetscLibType,pc::PC, l::PetscInt, mat::PetscMat) 
Gets the function to be used to inject primal vectors (i.e. solutions)
from level l to l-1.

Logically Collective

Input Parameters:
- `pc` - the multigrid context
- `l`  - the level (0 is coarsest) to supply [Do not supply 0]

Output Parameter:
- `mat` - the restriction matrix (may be `NULL` if no injection is available).

Level: advanced

-seealso: [](ch_ksp), `PCMG`, `PCMGSetInjection()`, `PCMGetGetRestriction()`

# External Links
$(_doc_external("Ksp/PCMGGetInjection"))
"""
function PCMGGetInjection(petsclib::PetscLibType, pc::PC, l::PetscInt, mat::PetscMat) end

@for_petsc function PCMGGetInjection(petsclib::$UnionPetscLib, pc::PC, l::$PetscInt, mat::PetscMat )
	mat_ = Ref(mat.ptr)

    @chk ccall(
               (:PCMGGetInjection, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, Ptr{CMat}),
               pc, l, mat_,
              )

	mat.ptr = C_NULL

	return nothing
end 

"""
	PCMGGetSmoother(petsclib::PetscLibType,pc::PC, l::PetscInt, ksp::PetscKSP) 
Gets the `KSP` context to be used as smoother for
both pre- and post-smoothing.  Call both `PCMGGetSmootherUp()` and
`PCMGGetSmootherDown()` to use different functions for pre- and
post-smoothing.

Not Collective, ksp returned is parallel if pc is

Input Parameters:
- `pc` - the multigrid context
- `l`  - the level (0 is coarsest) to supply

Output Parameter:
- `ksp` - the smoother

-seealso: [](ch_ksp), `PCMG`, `PCMGGetSmootherUp()`, `PCMGGetSmootherDown()`, `PCMGGetCoarseSolve()`

# External Links
$(_doc_external("Ksp/PCMGGetSmoother"))
"""
function PCMGGetSmoother(petsclib::PetscLibType, pc::PC, l::PetscInt, ksp::PetscKSP) end

@for_petsc function PCMGGetSmoother(petsclib::$UnionPetscLib, pc::PC, l::$PetscInt, ksp::PetscKSP )
	ksp_ = Ref(ksp.ptr)

    @chk ccall(
               (:PCMGGetSmoother, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, Ptr{CKSP}),
               pc, l, ksp_,
              )

	ksp.ptr = C_NULL

	return nothing
end 

"""
	PCMGGetSmootherUp(petsclib::PetscLibType,pc::PC, l::PetscInt, ksp::PetscKSP) 
Gets the KSP context to be used as smoother after
coarse grid correction (post-smoother).

Not Collective, ksp returned is parallel if pc is

Input Parameters:
- `pc` - the multigrid context
- `l`  - the level (0 is coarsest) to supply

Output Parameter:
- `ksp` - the smoother

Level: advanced

-seealso: [](ch_ksp), `PCMG`, `PCMGGetSmootherDown()`

# External Links
$(_doc_external("Ksp/PCMGGetSmootherUp"))
"""
function PCMGGetSmootherUp(petsclib::PetscLibType, pc::PC, l::PetscInt, ksp::PetscKSP) end

@for_petsc function PCMGGetSmootherUp(petsclib::$UnionPetscLib, pc::PC, l::$PetscInt, ksp::PetscKSP )
	ksp_ = Ref(ksp.ptr)

    @chk ccall(
               (:PCMGGetSmootherUp, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, Ptr{CKSP}),
               pc, l, ksp_,
              )

	ksp.ptr = C_NULL

	return nothing
end 

"""
	PCMGGetSmootherDown(petsclib::PetscLibType,pc::PC, l::PetscInt, ksp::PetscKSP) 
Gets the `KSP` context to be used as smoother before
coarse grid correction (pre-smoother).

Not Collective, ksp returned is parallel if pc is

Input Parameters:
- `pc` - the multigrid context
- `l`  - the level (0 is coarsest) to supply

Output Parameter:
- `ksp` - the smoother

Level: advanced

-seealso: [](ch_ksp), `PCMG`, `PCMGGetSmootherUp()`, `PCMGGetSmoother()`

# External Links
$(_doc_external("Ksp/PCMGGetSmootherDown"))
"""
function PCMGGetSmootherDown(petsclib::PetscLibType, pc::PC, l::PetscInt, ksp::PetscKSP) end

@for_petsc function PCMGGetSmootherDown(petsclib::$UnionPetscLib, pc::PC, l::$PetscInt, ksp::PetscKSP )
	ksp_ = Ref(ksp.ptr)

    @chk ccall(
               (:PCMGGetSmootherDown, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, Ptr{CKSP}),
               pc, l, ksp_,
              )

	ksp.ptr = C_NULL

	return nothing
end 

"""
	PCMGSetCycleTypeOnLevel(petsclib::PetscLibType,pc::PC, l::PetscInt, c::PCMGCycleType) 
Sets the type of cycle (aka cycle index) to run on the specified level.

Logically Collective

Input Parameters:
- `pc` - the multigrid context
- `l`  - the level (0 is coarsest)
- `c`  - either `PC_MG_CYCLE_V` or `PC_MG_CYCLE_W`

Level: advanced

-seealso: [](ch_ksp), `PCMG`, `PCMGCycleType`, `PCMGSetCycleType()`

# External Links
$(_doc_external("Ksp/PCMGSetCycleTypeOnLevel"))
"""
function PCMGSetCycleTypeOnLevel(petsclib::PetscLibType, pc::PC, l::PetscInt, c::PCMGCycleType) end

@for_petsc function PCMGSetCycleTypeOnLevel(petsclib::$UnionPetscLib, pc::PC, l::$PetscInt, c::PCMGCycleType )

    @chk ccall(
               (:PCMGSetCycleTypeOnLevel, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, PCMGCycleType),
               pc, l, c,
              )


	return nothing
end 

"""
	PCMGSetRhs(petsclib::PetscLibType,pc::PC, l::PetscInt, c::PetscVec) 
Sets the vector to be used to store the right

Logically Collective

Input Parameters:
- `pc` - the multigrid context
- `l`  - the level (0 is coarsest) this is to be used for
- `c`  - the Vec

Level: advanced

-seealso: [](ch_ksp), `PCMG`, `PCMGSetX()`, `PCMGSetR()`

# External Links
$(_doc_external("Ksp/PCMGSetRhs"))
"""
function PCMGSetRhs(petsclib::PetscLibType, pc::PC, l::PetscInt, c::PetscVec) end

@for_petsc function PCMGSetRhs(petsclib::$UnionPetscLib, pc::PC, l::$PetscInt, c::PetscVec )

    @chk ccall(
               (:PCMGSetRhs, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, CVec),
               pc, l, c,
              )


	return nothing
end 

"""
	PCMGSetX(petsclib::PetscLibType,pc::PC, l::PetscInt, c::PetscVec) 
Sets the vector to be used to store the solution on a particular level.

Logically Collective

Input Parameters:
- `pc` - the multigrid context
- `l`  - the level (0 is coarsest) this is to be used for (do not supply the finest level)
- `c`  - the Vec

Level: advanced

-seealso: [](ch_ksp), `PCMG`, `PCMGSetRhs()`, `PCMGSetR()`

# External Links
$(_doc_external("Ksp/PCMGSetX"))
"""
function PCMGSetX(petsclib::PetscLibType, pc::PC, l::PetscInt, c::PetscVec) end

@for_petsc function PCMGSetX(petsclib::$UnionPetscLib, pc::PC, l::$PetscInt, c::PetscVec )

    @chk ccall(
               (:PCMGSetX, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, CVec),
               pc, l, c,
              )


	return nothing
end 

"""
	PCMGSetR(petsclib::PetscLibType,pc::PC, l::PetscInt, c::PetscVec) 
Sets the vector to be used to store the residual on a particular level.

Logically Collective

Input Parameters:
- `pc` - the multigrid context
- `l`  - the level (0 is coarsest) this is to be used for
- `c`  - the Vec

Level: advanced

-seealso: [](ch_ksp), `PCMG`, `PCMGSetRhs()`, `PCMGSetX()`

# External Links
$(_doc_external("Ksp/PCMGSetR"))
"""
function PCMGSetR(petsclib::PetscLibType, pc::PC, l::PetscInt, c::PetscVec) end

@for_petsc function PCMGSetR(petsclib::$UnionPetscLib, pc::PC, l::$PetscInt, c::PetscVec )

    @chk ccall(
               (:PCMGSetR, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, CVec),
               pc, l, c,
              )


	return nothing
end 

"""
	PCMGSetLevels(petsclib::PetscLibType,pc::PC, levels::PetscInt, comms::MPI_Comm) 
Sets the number of levels to use with `PCMG`.
Must be called before any other `PCMG` routine.

Logically Collective

Input Parameters:
- `pc`     - the preconditioner context
- `levels` - the number of levels
- `comms`  - optional communicators for each level; this is to allow solving the coarser problems
on smaller sets of processes. For processes that are not included in the computation
you must pass `MPI_COMM_NULL`. Use comms = `NULL` to specify that all processes
should participate in each level of problem.

Options Database Key:
- `-pc_mg_levels <levels>` - set the number of levels to use

Level: intermediate

-seealso: [](ch_ksp), `PCMGSetType()`, `PCMGGetLevels()`

# External Links
$(_doc_external("Ksp/PCMGSetLevels"))
"""
function PCMGSetLevels(petsclib::PetscLibType, pc::PC, levels::PetscInt, comms::MPI_Comm) end

@for_petsc function PCMGSetLevels(petsclib::$UnionPetscLib, pc::PC, levels::$PetscInt, comms::MPI_Comm )

    @chk ccall(
               (:PCMGSetLevels, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, Ptr{MPI_Comm}),
               pc, levels, comms,
              )


	return nothing
end 

"""
	levels::PetscInt = PCMGGetLevels(petsclib::PetscLibType,pc::PC) 
Gets the number of levels to use with `PCMG`.

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `levels` - the number of levels

Level: advanced

-seealso: [](ch_ksp), `PCMG`, `PCMGSetLevels()`

# External Links
$(_doc_external("Ksp/PCMGGetLevels"))
"""
function PCMGGetLevels(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCMGGetLevels(petsclib::$UnionPetscLib, pc::PC )
	levels_ = Ref{$PetscInt}()

    @chk ccall(
               (:PCMGGetLevels, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{$PetscInt}),
               pc, levels_,
              )

	levels = levels_[]

	return levels
end 

"""
	gc::PetscReal,oc::PetscReal = PCMGGetGridComplexity(petsclib::PetscLibType,pc::PC) 
compute operator and grid complexity of the `PCMG` hierarchy

Input Parameter:
- `pc` - the preconditioner context

Output Parameters:
- `gc` - grid complexity = sum_i(n_i) / n_0
- `oc` - operator complexity = sum_i(nnz_i) / nnz_0

Level: advanced

-seealso: [](ch_ksp), `PCMG`, `PCMGGetLevels()`, `PCMGSetLevels()`

# External Links
$(_doc_external("Ksp/PCMGGetGridComplexity"))
"""
function PCMGGetGridComplexity(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCMGGetGridComplexity(petsclib::$UnionPetscLib, pc::PC )
	gc_ = Ref{$PetscReal}()
	oc_ = Ref{$PetscReal}()

    @chk ccall(
               (:PCMGGetGridComplexity, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{$PetscReal}, Ptr{$PetscReal}),
               pc, gc_, oc_,
              )

	gc = gc_[]
	oc = oc_[]

	return gc,oc
end 

"""
	PCMGSetType(petsclib::PetscLibType,pc::PC, form::PCMGType) 
Determines the form of multigrid to use, either
multiplicative, additive, full, or the Kaskade algorithm.

Logically Collective

Input Parameters:
- `pc`   - the preconditioner context
- `form` - multigrid form, one of `PC_MG_MULTIPLICATIVE`, `PC_MG_ADDITIVE`, `PC_MG_FULL`, `PC_MG_KASKADE`

Options Database Key:
- `-pc_mg_type <form>` - Sets <form>, one of multiplicative, additive, full, kaskade

Level: advanced

-seealso: [](ch_ksp), `PCMGType`, `PCMG`, `PCMGGetLevels()`, `PCMGSetLevels()`, `PCMGGetType()`, `PCMGCycleType`

# External Links
$(_doc_external("Ksp/PCMGSetType"))
"""
function PCMGSetType(petsclib::PetscLibType, pc::PC, form::PCMGType) end

@for_petsc function PCMGSetType(petsclib::$UnionPetscLib, pc::PC, form::PCMGType )

    @chk ccall(
               (:PCMGSetType, $petsc_library),
               PetscErrorCode,
               (PC, PCMGType),
               pc, form,
              )


	return nothing
end 

"""
	type::PCMGType = PCMGGetType(petsclib::PetscLibType,pc::PC) 
Finds the form of multigrid the `PCMG` is using  multiplicative, additive, full, or the Kaskade algorithm.

Logically Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `type` - one of `PC_MG_MULTIPLICATIVE`, `PC_MG_ADDITIVE`, `PC_MG_FULL`, `PC_MG_KASKADE`, `PCMGCycleType`

Level: advanced

-seealso: [](ch_ksp), `PCMGType`, `PCMG`, `PCMGGetLevels()`, `PCMGSetLevels()`, `PCMGSetType()`

# External Links
$(_doc_external("Ksp/PCMGGetType"))
"""
function PCMGGetType(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCMGGetType(petsclib::$UnionPetscLib, pc::PC )
	type_ = Ref{PCMGType}()

    @chk ccall(
               (:PCMGGetType, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PCMGType}),
               pc, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	PCMGSetCycleType(petsclib::PetscLibType,pc::PC, n::PCMGCycleType) 
Sets the type cycles to use.  Use `PCMGSetCycleTypeOnLevel()` for more
complicated cycling.

Logically Collective

Input Parameters:
- `pc` - the multigrid context
- `n`  - either `PC_MG_CYCLE_V` or `PC_MG_CYCLE_W`

Options Database Key:
- `-pc_mg_cycle_type <v,w>` - provide the cycle desired

Level: advanced

-seealso: [](ch_ksp), `PCMG`, `PCMGSetCycleTypeOnLevel()`, `PCMGType`, `PCMGCycleType`

# External Links
$(_doc_external("Ksp/PCMGSetCycleType"))
"""
function PCMGSetCycleType(petsclib::PetscLibType, pc::PC, n::PCMGCycleType) end

@for_petsc function PCMGSetCycleType(petsclib::$UnionPetscLib, pc::PC, n::PCMGCycleType )

    @chk ccall(
               (:PCMGSetCycleType, $petsc_library),
               PetscErrorCode,
               (PC, PCMGCycleType),
               pc, n,
              )


	return nothing
end 

"""
	PCMGMultiplicativeSetCycles(petsclib::PetscLibType,pc::PC, n::PetscInt) 
Sets the number of cycles to use for each preconditioner step
of multigrid when `PCMGType` is `PC_MG_MULTIPLICATIVE`

Logically Collective

Input Parameters:
- `pc` - the multigrid context
- `n`  - number of cycles (default is 1)

Options Database Key:
- `-pc_mg_multiplicative_cycles n` - set the number of cycles

Level: advanced

-seealso: [](ch_ksp), `PCMGSetCycleTypeOnLevel()`, `PCMGSetCycleType()`, `PCMGCycleType`, `PCMGType`

# External Links
$(_doc_external("Ksp/PCMGMultiplicativeSetCycles"))
"""
function PCMGMultiplicativeSetCycles(petsclib::PetscLibType, pc::PC, n::PetscInt) end

@for_petsc function PCMGMultiplicativeSetCycles(petsclib::$UnionPetscLib, pc::PC, n::$PetscInt )

    @chk ccall(
               (:PCMGMultiplicativeSetCycles, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, n,
              )


	return nothing
end 

"""
	PCMGSetGalerkin(petsclib::PetscLibType,pc::PC, use::PCMGGalerkinType) 
Causes the coarser grid matrices to be computed from the
finest grid via the Galerkin process: A_i-1 = r_i * A_i * p_i

Logically Collective

Input Parameters:
- `pc`  - the multigrid context
- `use` - one of `PC_MG_GALERKIN_BOTH`, `PC_MG_GALERKIN_PMAT`, `PC_MG_GALERKIN_MAT`, or `PC_MG_GALERKIN_NONE`

Options Database Key:
- `-pc_mg_galerkin <both,pmat,mat,none>` - set the matrices to form via the Galerkin process

Level: intermediate

-seealso: [](ch_ksp), `PCMG`, `PCMGGetGalerkin()`, `PCMGGalerkinType`

# External Links
$(_doc_external("Ksp/PCMGSetGalerkin"))
"""
function PCMGSetGalerkin(petsclib::PetscLibType, pc::PC, use::PCMGGalerkinType) end

@for_petsc function PCMGSetGalerkin(petsclib::$UnionPetscLib, pc::PC, use::PCMGGalerkinType )

    @chk ccall(
               (:PCMGSetGalerkin, $petsc_library),
               PetscErrorCode,
               (PC, PCMGGalerkinType),
               pc, use,
              )


	return nothing
end 

"""
	PCMGGetGalerkin(petsclib::PetscLibType,pc::PC, galerkin::PCMGGalerkinType) 
Checks if Galerkin multigrid is being used, i.e. A_i

Not Collective

Input Parameter:
- `pc` - the multigrid context

Output Parameter:
- `galerkin` - one of `PC_MG_GALERKIN_BOTH`,`PC_MG_GALERKIN_PMAT`,`PC_MG_GALERKIN_MAT`, `PC_MG_GALERKIN_NONE`, or `PC_MG_GALERKIN_EXTERNAL`

Level: intermediate

-seealso: [](ch_ksp), `PCMG`, `PCMGSetGalerkin()`, `PCMGGalerkinType`

# External Links
$(_doc_external("Ksp/PCMGGetGalerkin"))
"""
function PCMGGetGalerkin(petsclib::PetscLibType, pc::PC, galerkin::PCMGGalerkinType) end

@for_petsc function PCMGGetGalerkin(petsclib::$UnionPetscLib, pc::PC, galerkin::PCMGGalerkinType )

    @chk ccall(
               (:PCMGGetGalerkin, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PCMGGalerkinType}),
               pc, galerkin,
              )


	return nothing
end 

"""
	PCMGSetAdaptCoarseSpaceType(petsclib::PetscLibType,pc::PC, ctype::PCMGCoarseSpaceType) 
Set the type of adaptive coarse space.

Adapts or creates the interpolator based upon a vector space which should be accurately captured by the next coarser mesh, and thus accurately interpolated.

Logically Collective

Input Parameters:
- `pc`    - the multigrid context
- `ctype` - the type of coarse space

Options Database Keys:
- `-pc_mg_adapt_interp_n <int>`             - The number of modes to use
- `-pc_mg_adapt_interp_coarse_space <type>` - The type of coarse space: none, `polynomial`, `harmonic`, `eigenvector`, `generalized_eigenvector`, `gdsw`

Level: intermediate

-seealso: [](ch_ksp), `PCMG`, `PCMGCoarseSpaceType`, `PCMGGetAdaptCoarseSpaceType()`, `PCMGSetGalerkin()`, `PCMGSetAdaptInterpolation()`, `DM`

# External Links
$(_doc_external("Ksp/PCMGSetAdaptCoarseSpaceType"))
"""
function PCMGSetAdaptCoarseSpaceType(petsclib::PetscLibType, pc::PC, ctype::PCMGCoarseSpaceType) end

@for_petsc function PCMGSetAdaptCoarseSpaceType(petsclib::$UnionPetscLib, pc::PC, ctype::PCMGCoarseSpaceType )

    @chk ccall(
               (:PCMGSetAdaptCoarseSpaceType, $petsc_library),
               PetscErrorCode,
               (PC, PCMGCoarseSpaceType),
               pc, ctype,
              )


	return nothing
end 

"""
	ctype::PCMGCoarseSpaceType = PCMGGetAdaptCoarseSpaceType(petsclib::PetscLibType,pc::PC) 
Get the type of adaptive coarse space.

Not Collective

Input Parameter:
- `pc` - the multigrid context

Output Parameter:
- `ctype` - the type of coarse space

Level: intermediate

-seealso: [](ch_ksp), `PCMG`, `PCMGCoarseSpaceType`, `PCMGSetAdaptCoarseSpaceType()`, `PCMGSetGalerkin()`, `PCMGSetAdaptInterpolation()`

# External Links
$(_doc_external("Ksp/PCMGGetAdaptCoarseSpaceType"))
"""
function PCMGGetAdaptCoarseSpaceType(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCMGGetAdaptCoarseSpaceType(petsclib::$UnionPetscLib, pc::PC )
	ctype_ = Ref{PCMGCoarseSpaceType}()

    @chk ccall(
               (:PCMGGetAdaptCoarseSpaceType, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PCMGCoarseSpaceType}),
               pc, ctype_,
              )

	ctype = unsafe_string(ctype_[])

	return ctype
end 

"""
	PCMGSetAdaptInterpolation(petsclib::PetscLibType,pc::PC, adapt::PetscBool) 
Adapt the interpolator based upon a vector space which should be accurately captured by the next coarser mesh, and thus accurately interpolated.

Logically Collective

Input Parameters:
- `pc`    - the multigrid context
- `adapt` - flag for adaptation of the interpolator

Options Database Keys:
- `-pc_mg_adapt_interp`                     - Turn on adaptation
- `-pc_mg_adapt_interp_n <int>`             - The number of modes to use, should be divisible by dimension
- `-pc_mg_adapt_interp_coarse_space <type>` - The type of coarse space: polynomial, harmonic, eigenvector, generalized_eigenvector

Level: intermediate

-seealso: [](ch_ksp), `PCMG`, `PCMGGetAdaptInterpolation()`, `PCMGSetGalerkin()`, `PCMGGetAdaptCoarseSpaceType()`, `PCMGSetAdaptCoarseSpaceType()`

# External Links
$(_doc_external("Ksp/PCMGSetAdaptInterpolation"))
"""
function PCMGSetAdaptInterpolation(petsclib::PetscLibType, pc::PC, adapt::PetscBool) end

@for_petsc function PCMGSetAdaptInterpolation(petsclib::$UnionPetscLib, pc::PC, adapt::PetscBool )

    @chk ccall(
               (:PCMGSetAdaptInterpolation, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, adapt,
              )


	return nothing
end 

"""
	adapt::PetscBool = PCMGGetAdaptInterpolation(petsclib::PetscLibType,pc::PC) 
Get the flag to adapt the interpolator based upon a vector space which should be accurately captured by the next coarser mesh,
and thus accurately interpolated.

Not Collective

Input Parameter:
- `pc` - the multigrid context

Output Parameter:
- `adapt` - flag for adaptation of the interpolator

Level: intermediate

-seealso: [](ch_ksp), `PCMG`, `PCMGSetAdaptInterpolation()`, `PCMGSetGalerkin()`, `PCMGGetAdaptCoarseSpaceType()`, `PCMGSetAdaptCoarseSpaceType()`

# External Links
$(_doc_external("Ksp/PCMGGetAdaptInterpolation"))
"""
function PCMGGetAdaptInterpolation(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCMGGetAdaptInterpolation(petsclib::$UnionPetscLib, pc::PC )
	adapt_ = Ref{PetscBool}()

    @chk ccall(
               (:PCMGGetAdaptInterpolation, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PetscBool}),
               pc, adapt_,
              )

	adapt = adapt_[]

	return adapt
end 

"""
	PCMGSetAdaptCR(petsclib::PetscLibType,pc::PC, cr::PetscBool) 
Monitor the coarse space quality using an auxiliary solve with compatible relaxation.

Logically Collective

Input Parameters:
- `pc` - the multigrid context
- `cr` - flag for compatible relaxation

Options Database Key:
- `-pc_mg_adapt_cr` - Turn on compatible relaxation

Level: intermediate

-seealso: [](ch_ksp), `PCMG`, `PCMGGetAdaptCR()`, `PCMGSetAdaptInterpolation()`, `PCMGSetGalerkin()`, `PCMGGetAdaptCoarseSpaceType()`, `PCMGSetAdaptCoarseSpaceType()`

# External Links
$(_doc_external("Ksp/PCMGSetAdaptCR"))
"""
function PCMGSetAdaptCR(petsclib::PetscLibType, pc::PC, cr::PetscBool) end

@for_petsc function PCMGSetAdaptCR(petsclib::$UnionPetscLib, pc::PC, cr::PetscBool )

    @chk ccall(
               (:PCMGSetAdaptCR, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, cr,
              )


	return nothing
end 

"""
	cr::PetscBool = PCMGGetAdaptCR(petsclib::PetscLibType,pc::PC) 
Get the flag to monitor coarse space quality using an auxiliary solve with compatible relaxation.

Not Collective

Input Parameter:
- `pc` - the multigrid context

Output Parameter:
- `cr` - flag for compatible relaxaion

Level: intermediate

-seealso: [](ch_ksp), `PCMGSetAdaptCR()`, `PCMGGetAdaptInterpolation()`, `PCMGSetGalerkin()`, `PCMGGetAdaptCoarseSpaceType()`, `PCMGSetAdaptCoarseSpaceType()`

# External Links
$(_doc_external("Ksp/PCMGGetAdaptCR"))
"""
function PCMGGetAdaptCR(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCMGGetAdaptCR(petsclib::$UnionPetscLib, pc::PC )
	cr_ = Ref{PetscBool}()

    @chk ccall(
               (:PCMGGetAdaptCR, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PetscBool}),
               pc, cr_,
              )

	cr = cr_[]

	return cr
end 

"""
	PCMGSetNumberSmooth(petsclib::PetscLibType,pc::PC, n::PetscInt) 
Sets the number of pre and post
on all levels.  Use `PCMGDistinctSmoothUp()` to create separate up and down smoothers if you want different numbers of
pre- and post-smoothing steps.

Logically Collective

Input Parameters:
- `pc` - the multigrid context
- `n`  - the number of smoothing steps

Options Database Key:
- `-mg_levels_ksp_max_it <n>` - Sets number of pre and post-smoothing steps

Level: advanced

-seealso: [](ch_ksp), `PCMG`, `PCMGSetDistinctSmoothUp()`

# External Links
$(_doc_external("Ksp/PCMGSetNumberSmooth"))
"""
function PCMGSetNumberSmooth(petsclib::PetscLibType, pc::PC, n::PetscInt) end

@for_petsc function PCMGSetNumberSmooth(petsclib::$UnionPetscLib, pc::PC, n::$PetscInt )

    @chk ccall(
               (:PCMGSetNumberSmooth, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, n,
              )


	return nothing
end 

"""
	PCMGSetDistinctSmoothUp(petsclib::PetscLibType,pc::PC) 
sets the up (post) smoother to be a separate `KSP` from the down (pre) smoother on all levels
and adds the suffix _up to the options name

Logically Collective

Input Parameter:
- `pc` - the preconditioner context

Options Database Key:
- `-pc_mg_distinct_smoothup <bool>` - use distinct smoothing objects

Level: advanced

-seealso: [](ch_ksp), `PCMG`, `PCMGSetNumberSmooth()`

# External Links
$(_doc_external("Ksp/PCMGSetDistinctSmoothUp"))
"""
function PCMGSetDistinctSmoothUp(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCMGSetDistinctSmoothUp(petsclib::$UnionPetscLib, pc::PC )

    @chk ccall(
               (:PCMGSetDistinctSmoothUp, $petsc_library),
               PetscErrorCode,
               (PC,),
               pc,
              )


	return nothing
end 

"""
	PCMGRegisterCoarseSpaceConstructor(petsclib::PetscLibType,name::String, fnc::PCMGCoarseSpaceConstructorFn) 
Adds a method to the `PCMG` package for coarse space construction.

Not Collective, No Fortran Support

Input Parameters:
- `name`     - name of the constructor
- `function` - constructor routine, see `PCMGCoarseSpaceConstructorFn`

Level: advanced

-seealso: [](ch_ksp), `PCMGCoarseSpaceConstructorFn`, `PCMG`, `PCMGGetCoarseSpaceConstructor()`, `PCRegister()`

# External Links
$(_doc_external("Ksp/PCMGRegisterCoarseSpaceConstructor"))
"""
function PCMGRegisterCoarseSpaceConstructor(petsclib::PetscLibType, name::String, fnc::PCMGCoarseSpaceConstructorFn) end

@for_petsc function PCMGRegisterCoarseSpaceConstructor(petsclib::$UnionPetscLib, name::String, fnc::PCMGCoarseSpaceConstructorFn )

    @chk ccall(
               (:PCMGRegisterCoarseSpaceConstructor, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{PCMGCoarseSpaceConstructorFn}),
               name, fnc,
              )


	return nothing
end 

"""
	PCMGGetCoarseSpaceConstructor(petsclib::PetscLibType,name::String, fnc::PCMGCoarseSpaceConstructorFn) 
Returns the given coarse space construction method.

Not Collective, No Fortran Support

Input Parameter:
- `name` - name of the constructor

Output Parameter:
- `function` - constructor routine

Level: advanced

-seealso: [](ch_ksp), `PCMGCoarseSpaceConstructorFn`, `PCMG`, `PCMGRegisterCoarseSpaceConstructor()`, `PCRegister()`

# External Links
$(_doc_external("Ksp/PCMGGetCoarseSpaceConstructor"))
"""
function PCMGGetCoarseSpaceConstructor(petsclib::PetscLibType, name::String, fnc::PCMGCoarseSpaceConstructorFn) end

@for_petsc function PCMGGetCoarseSpaceConstructor(petsclib::$UnionPetscLib, name::String, fnc::PCMGCoarseSpaceConstructorFn )

    @chk ccall(
               (:PCMGGetCoarseSpaceConstructor, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, PCMGCoarseSpaceConstructorFn),
               name, fnc,
              )


	return nothing
end 

"""
	PCHMGSetReuseInterpolation(petsclib::PetscLibType,pc::PC, reuse::PetscBool) 
Reuse the interpolation matrices in `PCHMG` after changing the matrices numerical values

Logically Collective

Input Parameters:
- `pc`    - the `PCHMG` context
- `reuse` - `PETSC_TRUE` indicates that `PCHMG` will reuse the interpolations

Options Database Key:
- `-pc_hmg_reuse_interpolation <true | false>` - Whether or not to reuse the interpolations. If true, it potentially save the compute time.

Level: beginner

-seealso: [](ch_ksp), `PCHMG`, `PCGAMG`, `PCHMGSetUseSubspaceCoarsening()`, `PCHMGSetCoarseningComponent()`, `PCHMGSetInnerPCType()`

# External Links
$(_doc_external("Ksp/PCHMGSetReuseInterpolation"))
"""
function PCHMGSetReuseInterpolation(petsclib::PetscLibType, pc::PC, reuse::PetscBool) end

@for_petsc function PCHMGSetReuseInterpolation(petsclib::$UnionPetscLib, pc::PC, reuse::PetscBool )

    @chk ccall(
               (:PCHMGSetReuseInterpolation, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, reuse,
              )


	return nothing
end 

"""
	PCHMGSetUseSubspaceCoarsening(petsclib::PetscLibType,pc::PC, subspace::PetscBool) 
Use subspace coarsening in `PCHMG`

Logically Collective

Input Parameters:
- `pc`       - the `PCHMG` context
- `subspace` - `PETSC_TRUE` indicates that `PCHMG` will use the subspace coarsening

Options Database Key:
- `-pc_hmg_use_subspace_coarsening  <true | false>` - Whether or not to use subspace coarsening (that is, coarsen a submatrix).

Level: beginner

-seealso: [](ch_ksp), `PCHMG`, `PCHMGSetReuseInterpolation()`, `PCHMGSetCoarseningComponent()`, `PCHMGSetInnerPCType()`

# External Links
$(_doc_external("Ksp/PCHMGSetUseSubspaceCoarsening"))
"""
function PCHMGSetUseSubspaceCoarsening(petsclib::PetscLibType, pc::PC, subspace::PetscBool) end

@for_petsc function PCHMGSetUseSubspaceCoarsening(petsclib::$UnionPetscLib, pc::PC, subspace::PetscBool )

    @chk ccall(
               (:PCHMGSetUseSubspaceCoarsening, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, subspace,
              )


	return nothing
end 

"""
	PCHMGSetInnerPCType(petsclib::PetscLibType,pc::PC, type::PCType) 
Set an inner `PC` type to be used in the `PCHMG` preconditioner. That is the method used to compute
the hierarchy of restriction operators.

Logically Collective

Input Parameters:
- `pc`   - the `PCHMG` context
- `type` - `PCHYPRE` or `PCGAMG` coarsening algorithm

Options Database Key:
- `-hmg_inner_pc_type <hypre, gamg>` - What method is used to coarsen matrix

Level: beginner

-seealso: [](ch_ksp), `PCHMG`, `PCType`, `PCHMGSetReuseInterpolation()`, `PCHMGSetUseSubspaceCoarsening()`, `PCHMGSetCoarseningComponent()`

# External Links
$(_doc_external("Ksp/PCHMGSetInnerPCType"))
"""
function PCHMGSetInnerPCType(petsclib::PetscLibType, pc::PC, type::PCType) end

@for_petsc function PCHMGSetInnerPCType(petsclib::$UnionPetscLib, pc::PC, type::PCType )

    @chk ccall(
               (:PCHMGSetInnerPCType, $petsc_library),
               PetscErrorCode,
               (PC, PCType),
               pc, type,
              )


	return nothing
end 

"""
	PCHMGSetCoarseningComponent(petsclib::PetscLibType,pc::PC, component::PetscInt) 
Set which component of the PDE is used for the subspace

Logically Collective

Input Parameters:
- `pc`        - the `PCHMG` context
- `component` - which component `PC` will coarsen

Options Database Key:
- `-pc_hmg_coarsening_component <i>` - Which component is chosen for the subspace-based coarsening algorithm

Level: beginner

-seealso: [](ch_ksp), `PCHMG`, `PCType`, `PCGAMG`, `PCHMGSetReuseInterpolation()`, `PCHMGSetUseSubspaceCoarsening()`, `PCHMGSetInnerPCType()`

# External Links
$(_doc_external("Ksp/PCHMGSetCoarseningComponent"))
"""
function PCHMGSetCoarseningComponent(petsclib::PetscLibType, pc::PC, component::PetscInt) end

@for_petsc function PCHMGSetCoarseningComponent(petsclib::$UnionPetscLib, pc::PC, component::$PetscInt )

    @chk ccall(
               (:PCHMGSetCoarseningComponent, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, component,
              )


	return nothing
end 

"""
	PCHMGUseMatMAIJ(petsclib::PetscLibType,pc::PC, usematmaij::PetscBool) 
Set a flag that indicates if or not to use `MATMAIJ` for the interpolation matrices to save memory

Logically Collective

Input Parameters:
- `pc`         - the `PCHMG` context
- `usematmaij` - `PETSC_TRUE` (default) to use `MATMAIJ` for interpolations.

Options Database Key:
- `-pc_hmg_use_matmaij` - <true | false >

Level: beginner

-seealso: [](ch_ksp), `PCHMG`, `PCType`, `PCGAMG`

# External Links
$(_doc_external("Ksp/PCHMGUseMatMAIJ"))
"""
function PCHMGUseMatMAIJ(petsclib::PetscLibType, pc::PC, usematmaij::PetscBool) end

@for_petsc function PCHMGUseMatMAIJ(petsclib::$UnionPetscLib, pc::PC, usematmaij::PetscBool )

    @chk ccall(
               (:PCHMGUseMatMAIJ, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, usematmaij,
              )


	return nothing
end 

"""
	PCTelescopeGetKSP(petsclib::PetscLibType,pc::PC, subksp::PetscKSP) 
Gets the `KSP` created by the telescoping `PC`.

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `subksp` - the `KSP` defined on the smaller set of processes

Level: advanced

-seealso: [](ch_ksp), `PC`, `KSP`, `PCTELESCOPE`

# External Links
$(_doc_external("Ksp/PCTelescopeGetKSP"))
"""
function PCTelescopeGetKSP(petsclib::PetscLibType, pc::PC, subksp::PetscKSP) end

@for_petsc function PCTelescopeGetKSP(petsclib::$UnionPetscLib, pc::PC, subksp::PetscKSP )
	subksp_ = Ref(subksp.ptr)

    @chk ccall(
               (:PCTelescopeGetKSP, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{CKSP}),
               pc, subksp_,
              )

	subksp.ptr = C_NULL

	return nothing
end 

"""
	fact::PetscInt = PCTelescopeGetReductionFactor(petsclib::PetscLibType,pc::PC) 
Gets the factor by which the original number of MPI processes has been reduced by that was set by
`PCTelescopeSetReductionFactor()`

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `fact` - the reduction factor

Level: advanced

-seealso: [](ch_ksp), `PC`, `PCTELESCOPE`, `PCTelescopeSetReductionFactor()`

# External Links
$(_doc_external("Ksp/PCTelescopeGetReductionFactor"))
"""
function PCTelescopeGetReductionFactor(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCTelescopeGetReductionFactor(petsclib::$UnionPetscLib, pc::PC )
	fact_ = Ref{$PetscInt}()

    @chk ccall(
               (:PCTelescopeGetReductionFactor, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{$PetscInt}),
               pc, fact_,
              )

	fact = fact_[]

	return fact
end 

"""
	fact::PetscInt = PCTelescopeSetReductionFactor(petsclib::PetscLibType,pc::PC) 
Sets the factor by which the original number of MPI processes will been reduced by when
constructing the subcommunicator to be used with the `PCTELESCOPE`.

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `fact` - the reduction factor

Level: advanced

-seealso: [](ch_ksp), `PCTELESCOPE`, `PCTelescopeGetReductionFactor()`

# External Links
$(_doc_external("Ksp/PCTelescopeSetReductionFactor"))
"""
function PCTelescopeSetReductionFactor(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCTelescopeSetReductionFactor(petsclib::$UnionPetscLib, pc::PC )
	fact_ = Ref{$PetscInt}()

    @chk ccall(
               (:PCTelescopeSetReductionFactor, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, fact_,
              )

	fact = fact_[]

	return fact
end 

"""
	v::PetscBool = PCTelescopeGetIgnoreDM(petsclib::PetscLibType,pc::PC) 
Get the flag indicating if any `DM` attached to the `PC` will be used in constructing the `PC` on the
reduced number of MPI processes

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `v` - the flag

Level: advanced

-seealso: [](ch_ksp), `DM`, `PCTELESCOPE`, `PCTelescopeSetIgnoreDM()`

# External Links
$(_doc_external("Ksp/PCTelescopeGetIgnoreDM"))
"""
function PCTelescopeGetIgnoreDM(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCTelescopeGetIgnoreDM(petsclib::$UnionPetscLib, pc::PC )
	v_ = Ref{PetscBool}()

    @chk ccall(
               (:PCTelescopeGetIgnoreDM, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PetscBool}),
               pc, v_,
              )

	v = v_[]

	return v
end 

"""
	v::PetscBool = PCTelescopeSetIgnoreDM(petsclib::PetscLibType,pc::PC) 
Set a flag to ignore any `DM` attached to the `PC` when constructing the `PC` on the
reduced number of MPI processes

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `v` - Use `PETSC_TRUE` to ignore any `DM`

Level: advanced

-seealso: [](ch_ksp), `DM`, `PCTELESCOPE`, `PCTelescopeGetIgnoreDM()`

# External Links
$(_doc_external("Ksp/PCTelescopeSetIgnoreDM"))
"""
function PCTelescopeSetIgnoreDM(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCTelescopeSetIgnoreDM(petsclib::$UnionPetscLib, pc::PC )
	v_ = Ref{PetscBool}()

    @chk ccall(
               (:PCTelescopeSetIgnoreDM, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, v_,
              )

	v = v_[]

	return v
end 

"""
	v::PetscBool = PCTelescopeGetUseCoarseDM(petsclib::PetscLibType,pc::PC) 
Get the flag indicating if the coarse `DM` attached to `DM` associated with the `PC` will be used in constructing
the `PC` on the reduced number of MPI processes

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `v` - the flag

Level: advanced

-seealso: [](ch_ksp), `DM`, `PCTELESCOPE`, `PCTelescopeSetIgnoreDM()`, `PCTelescopeSetUseCoarseDM()`

# External Links
$(_doc_external("Ksp/PCTelescopeGetUseCoarseDM"))
"""
function PCTelescopeGetUseCoarseDM(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCTelescopeGetUseCoarseDM(petsclib::$UnionPetscLib, pc::PC )
	v_ = Ref{PetscBool}()

    @chk ccall(
               (:PCTelescopeGetUseCoarseDM, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PetscBool}),
               pc, v_,
              )

	v = v_[]

	return v
end 

"""
	v::PetscBool = PCTelescopeSetUseCoarseDM(petsclib::PetscLibType,pc::PC) 
Set a flag to query the `DM` attached to the `PC` if it also has a coarse `DM` and utilize that `DM`
in constructing the `PC` on the reduced number of MPI processes

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `v` - Use `PETSC_FALSE` to ignore any coarse `DM`

Level: advanced

-seealso: [](ch_ksp), `DM`, `PCTELESCOPE`, `PCTelescopeSetIgnoreDM()`

# External Links
$(_doc_external("Ksp/PCTelescopeSetUseCoarseDM"))
"""
function PCTelescopeSetUseCoarseDM(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCTelescopeSetUseCoarseDM(petsclib::$UnionPetscLib, pc::PC )
	v_ = Ref{PetscBool}()

    @chk ccall(
               (:PCTelescopeSetUseCoarseDM, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, v_,
              )

	v = v_[]

	return v
end 

"""
	v::PetscBool = PCTelescopeGetIgnoreKSPComputeOperators(petsclib::PetscLibType,pc::PC) 
Get the flag indicating if `KSPComputeOperators()` will be used to construct
the matrix on the reduced number of MPI processes

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `v` - the flag

Level: advanced

-seealso: [](ch_ksp), `PCTELESCOPE`, `PCTelescopeSetIgnoreDM()`, `PCTelescopeSetUseCoarseDM()`, `PCTelescopeSetIgnoreKSPComputeOperators()`

# External Links
$(_doc_external("Ksp/PCTelescopeGetIgnoreKSPComputeOperators"))
"""
function PCTelescopeGetIgnoreKSPComputeOperators(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCTelescopeGetIgnoreKSPComputeOperators(petsclib::$UnionPetscLib, pc::PC )
	v_ = Ref{PetscBool}()

    @chk ccall(
               (:PCTelescopeGetIgnoreKSPComputeOperators, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PetscBool}),
               pc, v_,
              )

	v = v_[]

	return v
end 

"""
	v::PetscBool = PCTelescopeSetIgnoreKSPComputeOperators(petsclib::PetscLibType,pc::PC) 
Set a flag to have `PCTELESCOPE` ignore the function provided to `KSPComputeOperators()` in
constructint the matrix on the reduced number of MPI processes

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `v` - Use `PETSC_TRUE` to ignore the function (if defined) set via `KSPSetComputeOperators()` on `pc`

Level: advanced

-seealso: [](ch_ksp), `PCTELESCOPE`, `PCTelescopeSetIgnoreDM()`, `PCTelescopeSetUseCoarseDM()`, `PCTelescopeGetIgnoreKSPComputeOperators()`

# External Links
$(_doc_external("Ksp/PCTelescopeSetIgnoreKSPComputeOperators"))
"""
function PCTelescopeSetIgnoreKSPComputeOperators(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCTelescopeSetIgnoreKSPComputeOperators(petsclib::$UnionPetscLib, pc::PC )
	v_ = Ref{PetscBool}()

    @chk ccall(
               (:PCTelescopeSetIgnoreKSPComputeOperators, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, v_,
              )

	v = v_[]

	return v
end 

"""
	PCTelescopeGetDM(petsclib::PetscLibType,pc::PC, subdm::PetscDM) 
Get the re

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `subdm` - The re-partitioned `DM`

Level: advanced

-seealso: [](ch_ksp), `DM`, `PCTELESCOPE`, `PCTelescopeSetIgnoreDM()`, `PCTelescopeSetUseCoarseDM()`, `PCTelescopeGetIgnoreKSPComputeOperators()`

# External Links
$(_doc_external("Ksp/PCTelescopeGetDM"))
"""
function PCTelescopeGetDM(petsclib::PetscLibType, pc::PC, subdm::PetscDM) end

@for_petsc function PCTelescopeGetDM(petsclib::$UnionPetscLib, pc::PC, subdm::PetscDM )
	subdm_ = Ref(subdm.ptr)

    @chk ccall(
               (:PCTelescopeGetDM, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{CDM}),
               pc, subdm_,
              )

	subdm.ptr = C_NULL

	return nothing
end 

"""
	PCTelescopeSetSubcommType(petsclib::PetscLibType,pc::PC, subcommtype::PetscSubcommType) 
set subcommunicator type `PetscSubcommType` (interlaced or contiguous) to be used when
the subcommunicator is generated from the given `PC`

Logically Collective

Input Parameters:
- `pc`          - the preconditioner context
- `subcommtype` - the subcommunicator type (see `PetscSubcommType`)

Level: advanced

-seealso: [](ch_ksp), `PetscSubcommType`, `PetscSubcomm`, `PCTELESCOPE`, `PCTelescopeGetSubcommType()`

# External Links
$(_doc_external("Ksp/PCTelescopeSetSubcommType"))
"""
function PCTelescopeSetSubcommType(petsclib::PetscLibType, pc::PC, subcommtype::PetscSubcommType) end

@for_petsc function PCTelescopeSetSubcommType(petsclib::$UnionPetscLib, pc::PC, subcommtype::PetscSubcommType )

    @chk ccall(
               (:PCTelescopeSetSubcommType, $petsc_library),
               PetscErrorCode,
               (PC, PetscSubcommType),
               pc, subcommtype,
              )


	return nothing
end 

"""
	subcommtype::PetscSubcommType = PCTelescopeGetSubcommType(petsclib::PetscLibType,pc::PC) 
Get the subcommunicator type `PetscSubcommType` (interlaced or contiguous) set with `PCTelescopeSetSubcommType()`

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `subcommtype` - the subcommunicator type (see `PetscSubcommType`)

Level: advanced

-seealso: [](ch_ksp), `PetscSubcomm`, `PetscSubcommType`, `PCTELESCOPE`, `PCTelescopeSetSubcommType()`

# External Links
$(_doc_external("Ksp/PCTelescopeGetSubcommType"))
"""
function PCTelescopeGetSubcommType(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCTelescopeGetSubcommType(petsclib::$UnionPetscLib, pc::PC )
	subcommtype_ = Ref{PetscSubcommType}()

    @chk ccall(
               (:PCTelescopeGetSubcommType, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PetscSubcommType}),
               pc, subcommtype_,
              )

	subcommtype = unsafe_string(subcommtype_[])

	return subcommtype
end 

"""
	PCShellGetContext(petsclib::PetscLibType,pc::PC, ctx::Cvoid) 
Returns the user

Not Collective

Input Parameter:
- `pc` - of type `PCSHELL`

Output Parameter:
- `ctx` - the user provided context

Level: advanced

-seealso: [](ch_ksp), `PC`, `PCSHELL`, `PCShellSetContext()`, `PCShellSetApply()`, `PCShellSetDestroy()`

# External Links
$(_doc_external("Ksp/PCShellGetContext"))
"""
function PCShellGetContext(petsclib::PetscLibType, pc::PC, ctx::Cvoid) end

@for_petsc function PCShellGetContext(petsclib::$UnionPetscLib, pc::PC, ctx::Cvoid )

    @chk ccall(
               (:PCShellGetContext, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{Cvoid}),
               pc, ctx,
              )


	return nothing
end 

"""
	PCShellSetContext(petsclib::PetscLibType,pc::PC, ctx::Cvoid) 
sets the context for a shell `PC` that can be accessed with `PCShellGetContext()`

Logically Collective

Input Parameters:
- `pc`  - the `PC` of type `PCSHELL`
- `ctx` - the context

Level: advanced

-seealso: [](ch_ksp), `PC`, `PCShellGetContext()`, `PCSHELL`, `PCShellSetApply()`, `PCShellSetDestroy()`

# External Links
$(_doc_external("Ksp/PCShellSetContext"))
"""
function PCShellSetContext(petsclib::PetscLibType, pc::PC, ctx::Cvoid) end

@for_petsc function PCShellSetContext(petsclib::$UnionPetscLib, pc::PC, ctx::Cvoid )

    @chk ccall(
               (:PCShellSetContext, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{Cvoid}),
               pc, ctx,
              )


	return nothing
end 

"""
	PCShellSetDestroy(petsclib::PetscLibType,pc::PC, destroy::external) 
Sets routine to use to destroy the user

Logically Collective

Input Parameters:
- `pc`      - the preconditioner context
- `destroy` - the application-provided destroy routine

Calling sequence of `destroy`:
- `pc` - the preconditioner

Level: intermediate

-seealso: [](ch_ksp), `PCSHELL`, `PCShellSetApply()`, `PCShellSetContext()`, `PCShellGetContext()`

# External Links
$(_doc_external("Ksp/PCShellSetDestroy"))
"""
function PCShellSetDestroy(petsclib::PetscLibType, pc::PC, destroy::external) end

@for_petsc function PCShellSetDestroy(petsclib::$UnionPetscLib, pc::PC, destroy::external )

    @chk ccall(
               (:PCShellSetDestroy, $petsc_library),
               PetscErrorCode,
               (PC, external),
               pc, destroy,
              )


	return nothing
end 

"""
	PCShellSetSetUp(petsclib::PetscLibType,pc::PC, setup::external) 
Sets routine to use to "setup" the preconditioner whenever the
matrix operator is changed.

Logically Collective

Input Parameters:
- `pc`    - the preconditioner context
- `setup` - the application-provided setup routine

Calling sequence of `setup`:
- `pc` - the preconditioner

Level: intermediate

-seealso: [](ch_ksp), `PCSHELL`, `PCShellSetApplyRichardson()`, `PCShellSetApply()`, `PCShellSetContext()`, , `PCShellGetContext()`

# External Links
$(_doc_external("Ksp/PCShellSetSetUp"))
"""
function PCShellSetSetUp(petsclib::PetscLibType, pc::PC, setup::external) end

@for_petsc function PCShellSetSetUp(petsclib::$UnionPetscLib, pc::PC, setup::external )

    @chk ccall(
               (:PCShellSetSetUp, $petsc_library),
               PetscErrorCode,
               (PC, external),
               pc, setup,
              )


	return nothing
end 

"""
	PCShellSetView(petsclib::PetscLibType,pc::PC, view::external) 
Sets routine to use as viewer of a `PCSHELL` shell preconditioner

Logically Collective

Input Parameters:
- `pc`   - the preconditioner context
- `view` - the application-provided view routine

Calling sequence of `view`:
- `pc` - the preconditioner
- `v`  - viewer

Level: advanced

-seealso: [](ch_ksp), `PC`, `PCSHELL`, `PCShellSetApplyRichardson()`, `PCShellSetSetUp()`, `PCShellSetApplyTranspose()`, `PCShellSetContext()`, `PCShellGetContext()`

# External Links
$(_doc_external("Ksp/PCShellSetView"))
"""
function PCShellSetView(petsclib::PetscLibType, pc::PC, view::external) end

@for_petsc function PCShellSetView(petsclib::$UnionPetscLib, pc::PC, view::external )

    @chk ccall(
               (:PCShellSetView, $petsc_library),
               PetscErrorCode,
               (PC, external),
               pc, view,
              )


	return nothing
end 

"""
	PCShellSetApply(petsclib::PetscLibType,pc::PC, apply::external) 
Sets routine to use as preconditioner.

Logically Collective

Input Parameters:
- `pc`    - the preconditioner context
- `apply` - the application-provided preconditioning routine

Calling sequence of `apply`:
- `pc`   - the preconditioner, get the application context with `PCShellGetContext()`
- `xin`  - input vector
- `xout` - output vector

Level: intermediate

-seealso: [](ch_ksp), `PCSHELL`, `PCShellSetApplyRichardson()`, `PCShellSetSetUp()`, `PCShellSetApplyTranspose()`, `PCShellSetContext()`, `PCShellSetApplyBA()`, `PCShellSetApplySymmetricRight()`, `PCShellSetApplySymmetricLeft()`, `PCShellGetContext()`

# External Links
$(_doc_external("Ksp/PCShellSetApply"))
"""
function PCShellSetApply(petsclib::PetscLibType, pc::PC, apply::external) end

@for_petsc function PCShellSetApply(petsclib::$UnionPetscLib, pc::PC, apply::external )

    @chk ccall(
               (:PCShellSetApply, $petsc_library),
               PetscErrorCode,
               (PC, external),
               pc, apply,
              )


	return nothing
end 

"""
	PCShellSetMatApply(petsclib::PetscLibType,pc::PC, matapply::external) 
Sets routine to use as preconditioner on a block of vectors.

Logically Collective

Input Parameters:
- `pc`       - the preconditioner context
- `matapply` - the application-provided preconditioning routine

Calling sequence of `matapply`:
- `pc`   - the preconditioner
- `Xin`  - input block of vectors represented as a dense `Mat`
- `Xout` - output block of vectors represented as a dense `Mat`

Level: advanced

-seealso: [](ch_ksp), `PCSHELL`, `PCShellSetApply()`, `PCShellSetContext()`, `PCShellGetContext()`

# External Links
$(_doc_external("Ksp/PCShellSetMatApply"))
"""
function PCShellSetMatApply(petsclib::PetscLibType, pc::PC, matapply::external) end

@for_petsc function PCShellSetMatApply(petsclib::$UnionPetscLib, pc::PC, matapply::external )

    @chk ccall(
               (:PCShellSetMatApply, $petsc_library),
               PetscErrorCode,
               (PC, external),
               pc, matapply,
              )


	return nothing
end 

"""
	PCShellSetApplySymmetricLeft(petsclib::PetscLibType,pc::PC, apply::external) 
Sets routine to use as left preconditioner (when the `PC_SYMMETRIC` is used).

Logically Collective

Input Parameters:
- `pc`    - the preconditioner context
- `apply` - the application-provided left preconditioning routine

Calling sequence of `apply`:
- `pc`   - the preconditioner
- `xin`  - input vector
- `xout` - output vector

Level: advanced

-seealso: [](ch_ksp), `PCSHELL`, `PCShellSetApply()`, `PCShellSetSetUp()`, `PCShellSetApplyTranspose()`, `PCShellSetContext()`

# External Links
$(_doc_external("Ksp/PCShellSetApplySymmetricLeft"))
"""
function PCShellSetApplySymmetricLeft(petsclib::PetscLibType, pc::PC, apply::external) end

@for_petsc function PCShellSetApplySymmetricLeft(petsclib::$UnionPetscLib, pc::PC, apply::external )

    @chk ccall(
               (:PCShellSetApplySymmetricLeft, $petsc_library),
               PetscErrorCode,
               (PC, external),
               pc, apply,
              )


	return nothing
end 

"""
	PCShellSetApplySymmetricRight(petsclib::PetscLibType,pc::PC, apply::external) 
Sets routine to use as right preconditioner (when the `PC_SYMMETRIC` is used).

Logically Collective

Input Parameters:
- `pc`    - the preconditioner context
- `apply` - the application-provided right preconditioning routine

Calling sequence of `apply`:
- `pc`   - the preconditioner
- `xin`  - input vector
- `xout` - output vector

Level: advanced

-seealso: [](ch_ksp), `PCSHELL`, `PCShellSetApply()`, `PCShellSetApplySymmetricLeft()`, `PCShellSetSetUp()`, `PCShellSetApplyTranspose()`, `PCShellSetContext()`, `PCShellGetContext()`

# External Links
$(_doc_external("Ksp/PCShellSetApplySymmetricRight"))
"""
function PCShellSetApplySymmetricRight(petsclib::PetscLibType, pc::PC, apply::external) end

@for_petsc function PCShellSetApplySymmetricRight(petsclib::$UnionPetscLib, pc::PC, apply::external )

    @chk ccall(
               (:PCShellSetApplySymmetricRight, $petsc_library),
               PetscErrorCode,
               (PC, external),
               pc, apply,
              )


	return nothing
end 

"""
	PCShellSetApplyBA(petsclib::PetscLibType,pc::PC, applyBA::external) 
Sets routine to use as the preconditioner times the operator.

Logically Collective

Input Parameters:
- `pc`      - the preconditioner context
- `applyBA` - the application-provided BA routine

Calling sequence of `applyBA`:
- `pc`   - the preconditioner
- `side` - `PC_LEFT`, `PC_RIGHT`, or `PC_SYMMETRIC`
- `xin`  - input vector
- `xout` - output vector
- `w`    - work vector

Level: intermediate

-seealso: [](ch_ksp), `PCSHELL`, `PCShellSetApplyRichardson()`, `PCShellSetSetUp()`, `PCShellSetApplyTranspose()`, `PCShellSetContext()`, `PCShellSetApply()`, `PCShellGetContext()`, `PCSide`

# External Links
$(_doc_external("Ksp/PCShellSetApplyBA"))
"""
function PCShellSetApplyBA(petsclib::PetscLibType, pc::PC, applyBA::external) end

@for_petsc function PCShellSetApplyBA(petsclib::$UnionPetscLib, pc::PC, applyBA::external )

    @chk ccall(
               (:PCShellSetApplyBA, $petsc_library),
               PetscErrorCode,
               (PC, external),
               pc, applyBA,
              )


	return nothing
end 

"""
	PCShellSetApplyTranspose(petsclib::PetscLibType,pc::PC, applytranspose::external) 
Sets routine to use as preconditioner transpose.

Logically Collective

Input Parameters:
- `pc`             - the preconditioner context
- `applytranspose` - the application-provided preconditioning transpose routine

Calling sequence of `applytranspose`:
- `pc`   - the preconditioner
- `xin`  - input vector
- `xout` - output vector

Level: intermediate

-seealso: [](ch_ksp), `PCSHELL`, `PCShellSetApplyRichardson()`, `PCShellSetSetUp()`, `PCShellSetApply()`, `PCShellSetContext()`, `PCShellSetApplyBA()`, `PCShellGetContext()`

# External Links
$(_doc_external("Ksp/PCShellSetApplyTranspose"))
"""
function PCShellSetApplyTranspose(petsclib::PetscLibType, pc::PC, applytranspose::external) end

@for_petsc function PCShellSetApplyTranspose(petsclib::$UnionPetscLib, pc::PC, applytranspose::external )

    @chk ccall(
               (:PCShellSetApplyTranspose, $petsc_library),
               PetscErrorCode,
               (PC, external),
               pc, applytranspose,
              )


	return nothing
end 

"""
	PCShellSetMatApplyTranspose(petsclib::PetscLibType,pc::PC, matapplytranspose::external) 
Sets routine to use as preconditioner transpose.

Logically Collective

Input Parameters:
- `pc`                - the preconditioner context
- `matapplytranspose` - the application-provided preconditioning transpose routine

Calling sequence of `matapplytranspose`:
- `pc`   - the preconditioner
- `xin`  - input matrix
- `xout` - output matrix

Level: intermediate

-seealso: [](ch_ksp), `PCSHELL`, `PCShellSetApplyRichardson()`, `PCShellSetSetUp()`, `PCShellSetApply()`, `PCShellSetContext()`, `PCShellSetApplyBA()`, `PCShellGetContext()`

# External Links
$(_doc_external("Ksp/PCShellSetMatApplyTranspose"))
"""
function PCShellSetMatApplyTranspose(petsclib::PetscLibType, pc::PC, matapplytranspose::external) end

@for_petsc function PCShellSetMatApplyTranspose(petsclib::$UnionPetscLib, pc::PC, matapplytranspose::external )

    @chk ccall(
               (:PCShellSetMatApplyTranspose, $petsc_library),
               PetscErrorCode,
               (PC, external),
               pc, matapplytranspose,
              )


	return nothing
end 

"""
	PCShellSetPreSolve(petsclib::PetscLibType,pc::PC, presolve::PCShellPSolveFn) 
Sets routine to apply to the operators/vectors before a `KSPSolve()` is
applied. This usually does something like scale the linear system in some application
specific way.

Logically Collective

Input Parameters:
- `pc`       - the preconditioner context
- `presolve` - the application-provided presolve routine, see `PCShellPSolveFn`

Level: advanced

-seealso: [](ch_ksp), `PCSHELL`, `PCShellPSolveFn`, `PCShellSetApplyRichardson()`, `PCShellSetSetUp()`, `PCShellSetApplyTranspose()`, `PCShellSetPostSolve()`, `PCShellSetContext()`, `PCShellGetContext()`

# External Links
$(_doc_external("Ksp/PCShellSetPreSolve"))
"""
function PCShellSetPreSolve(petsclib::PetscLibType, pc::PC, presolve::PCShellPSolveFn) end

@for_petsc function PCShellSetPreSolve(petsclib::$UnionPetscLib, pc::PC, presolve::PCShellPSolveFn )

    @chk ccall(
               (:PCShellSetPreSolve, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PCShellPSolveFn}),
               pc, presolve,
              )


	return nothing
end 

"""
	PCShellSetPostSolve(petsclib::PetscLibType,pc::PC, postsolve::PCShellPSolveFn) 
Sets routine to apply to the operators/vectors after a `KSPSolve()` is
applied. This usually does something like scale the linear system in some application
specific way.

Logically Collective

Input Parameters:
- `pc`        - the preconditioner context
- `postsolve` - the application-provided postsolve routine, see `PCShellPSolveFn`

Level: advanced

-seealso: [](ch_ksp), `PCSHELL`, `PCShellPSolveFn`, `PCShellSetApplyRichardson()`, `PCShellSetSetUp()`, `PCShellSetApplyTranspose()`, `PCShellSetPreSolve()`, `PCShellSetContext()`, `PCShellGetContext()`

# External Links
$(_doc_external("Ksp/PCShellSetPostSolve"))
"""
function PCShellSetPostSolve(petsclib::PetscLibType, pc::PC, postsolve::PCShellPSolveFn) end

@for_petsc function PCShellSetPostSolve(petsclib::$UnionPetscLib, pc::PC, postsolve::PCShellPSolveFn )

    @chk ccall(
               (:PCShellSetPostSolve, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PCShellPSolveFn}),
               pc, postsolve,
              )


	return nothing
end 

"""
	PCShellSetName(petsclib::PetscLibType,pc::PC, name::String) 
Sets an optional name to associate with a `PCSHELL`
preconditioner.

Not Collective

Input Parameters:
- `pc`   - the preconditioner context
- `name` - character string describing shell preconditioner

Level: intermediate

-seealso: [](ch_ksp), `PCSHELL`, `PCShellGetName()`, `PetscObjectSetName()`, `PetscObjectGetName()`

# External Links
$(_doc_external("Ksp/PCShellSetName"))
"""
function PCShellSetName(petsclib::PetscLibType, pc::PC, name::String) end

@for_petsc function PCShellSetName(petsclib::$UnionPetscLib, pc::PC, name::String )

    @chk ccall(
               (:PCShellSetName, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{Cchar}),
               pc, name,
              )


	return nothing
end 

"""
	PCShellGetName(petsclib::PetscLibType,pc::PC, name::String) 
Gets an optional name that the user has set for a `PCSHELL` with `PCShellSetName()`
preconditioner.

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `name` - character string describing shell preconditioner (you should not free this)

Level: intermediate

-seealso: [](ch_ksp), `PCSHELL`, `PCShellSetName()`, `PetscObjectSetName()`, `PetscObjectGetName()`

# External Links
$(_doc_external("Ksp/PCShellGetName"))
"""
function PCShellGetName(petsclib::PetscLibType, pc::PC, name::String) end

@for_petsc function PCShellGetName(petsclib::$UnionPetscLib, pc::PC, name::String )
	name_ = Ref(pointer(name))

    @chk ccall(
               (:PCShellGetName, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{Ptr{Cchar}}),
               pc, name_,
              )


	return nothing
end 

"""
	PCShellSetApplyRichardson(petsclib::PetscLibType,pc::PC, apply::external) 
Sets routine to use as preconditioner
in Richardson iteration.

Logically Collective

Input Parameters:
- `pc`    - the preconditioner context
- `apply` - the application-provided preconditioning routine

Calling sequence of `apply`:
- `pc`               - the preconditioner
- `b`                - right-hand side
- `x`                - current iterate
- `r`                - work space
- `rtol`             - relative tolerance of residual norm to stop at
- `abstol`           - absolute tolerance of residual norm to stop at
- `dtol`             - if residual norm increases by this factor than return
- `maxits`           - number of iterations to run
- `zeroinitialguess` - `PETSC_TRUE` if `x` is known to be initially zero
- `its`              - returns the number of iterations used
- `reason`           - returns the reason the iteration has converged

Level: advanced

-seealso: [](ch_ksp), `PCSHELL`, `PCShellSetApply()`, `PCShellSetContext()`, `PCRichardsonConvergedReason()`, `PCShellGetContext()`, `KSPRICHARDSON`

# External Links
$(_doc_external("Ksp/PCShellSetApplyRichardson"))
"""
function PCShellSetApplyRichardson(petsclib::PetscLibType, pc::PC, apply::external) end

@for_petsc function PCShellSetApplyRichardson(petsclib::$UnionPetscLib, pc::PC, apply::external )

    @chk ccall(
               (:PCShellSetApplyRichardson, $petsc_library),
               PetscErrorCode,
               (PC, external),
               pc, apply,
              )


	return nothing
end 

"""
	PCISSetUseStiffnessScaling(petsclib::PetscLibType,pc::PC, use::PetscBool) 
Tells `PCIS` to construct partition of unity using
the local matrices' diagonal entries

Logically Collective

Input Parameters:
- `pc`  - the preconditioning context
- `use` - whether or not it should use matrix diagonal to build partition of unity.

Level: intermediate

-seealso: [](ch_ksp), `PCBDDC`, `PCNN`, `PCISSetSubdomainDiagonalScaling()`, `PCISScatterArrayNToVecB()`,
`PCISSetSubdomainScalingFactor()`,
`PCISReset()`, `PCISInitialize()`, `PCISApplyInvSchur()`, `PCISApplySchur()`

# External Links
$(_doc_external("Ksp/PCISSetUseStiffnessScaling"))
"""
function PCISSetUseStiffnessScaling(petsclib::PetscLibType, pc::PC, use::PetscBool) end

@for_petsc function PCISSetUseStiffnessScaling(petsclib::$UnionPetscLib, pc::PC, use::PetscBool )

    @chk ccall(
               (:PCISSetUseStiffnessScaling, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, use,
              )


	return nothing
end 

"""
	PCISSetSubdomainDiagonalScaling(petsclib::PetscLibType,pc::PC, scaling_factors::PetscVec) 
Set diagonal scaling for `PCIS`.

Logically Collective

Input Parameters:
- `pc`              - the preconditioning context
- `scaling_factors` - scaling factors for the subdomain

Level: intermediate

-seealso: [](ch_ksp), `PCBDDC`, `PCNN`, `PCISScatterArrayNToVecB()`,
`PCISSetSubdomainScalingFactor()`, `PCISSetUseStiffnessScaling()`,
`PCISReset()`, `PCISInitialize()`, `PCISApplyInvSchur()`, `PCISApplySchur()`

# External Links
$(_doc_external("Ksp/PCISSetSubdomainDiagonalScaling"))
"""
function PCISSetSubdomainDiagonalScaling(petsclib::PetscLibType, pc::PC, scaling_factors::PetscVec) end

@for_petsc function PCISSetSubdomainDiagonalScaling(petsclib::$UnionPetscLib, pc::PC, scaling_factors::PetscVec )

    @chk ccall(
               (:PCISSetSubdomainDiagonalScaling, $petsc_library),
               PetscErrorCode,
               (PC, CVec),
               pc, scaling_factors,
              )


	return nothing
end 

"""
	PCISSetSubdomainScalingFactor(petsclib::PetscLibType,pc::PC, scal::PetscScalar) 
Set scaling factor for `PCIS`.

Not Collective

Input Parameters:
- `pc`   - the preconditioning context
- `scal` - scaling factor for the subdomain

Level: intermediate

-seealso: [](ch_ksp), `PCBDDC`, `PCNN`, `PCISScatterArrayNToVecB()`,
`PCISSetSubdomainDiagonalScaling()`, `PCISSetUseStiffnessScaling()`,
`PCISReset()`, `PCISInitialize()`, `PCISApplyInvSchur()`, `PCISApplySchur()`

# External Links
$(_doc_external("Ksp/PCISSetSubdomainScalingFactor"))
"""
function PCISSetSubdomainScalingFactor(petsclib::PetscLibType, pc::PC, scal::PetscScalar) end

@for_petsc function PCISSetSubdomainScalingFactor(petsclib::$UnionPetscLib, pc::PC, scal::$PetscScalar )

    @chk ccall(
               (:PCISSetSubdomainScalingFactor, $petsc_library),
               PetscErrorCode,
               (PC, $PetscScalar),
               pc, scal,
              )


	return nothing
end 

"""
	PCISSetUp(petsclib::PetscLibType,pc::PC, computematrices::PetscBool, computesolvers::PetscBool) 
sets up the `PC_IS` portion of `PCNN` and `PCBDDC` preconditioner context as part of their setup process

Input Parameters:
- `pc`              - the `PC` object, must be of type `PCNN` or `PCBDDC`
- `computematrices` - Extract the blocks `A_II`, `A_BI`, `A_IB` and `A_BB` from the matrix
- `computesolvers`  - Create the `KSP` for the local Dirichlet and Neumann problems

Level: advanced

-seealso: [](ch_ksp), `PCBDDC`, `PCNN`, `PCISSetUseStiffnessScaling()`, `PCISSetSubdomainDiagonalScaling()`, `PCISScatterArrayNToVecB()`,
`PCISSetSubdomainScalingFactor()`,
`PCISReset()`, `PCISApplySchur()`, `PCISApplyInvSchur()`

# External Links
$(_doc_external("Ksp/PCISSetUp"))
"""
function PCISSetUp(petsclib::PetscLibType, pc::PC, computematrices::PetscBool, computesolvers::PetscBool) end

@for_petsc function PCISSetUp(petsclib::$UnionPetscLib, pc::PC, computematrices::PetscBool, computesolvers::PetscBool )

    @chk ccall(
               (:PCISSetUp, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool, PetscBool),
               pc, computematrices, computesolvers,
              )


	return nothing
end 

"""
	PCISReset(petsclib::PetscLibType,pc::PC) 
Removes all the `PC_IS` parts of the `PC` implementation data structure

Input Parameter:
- `pc` - the `PC` object, must be of type `PCNN` or `PCBDDC`

Level: advanced

-seealso: [](ch_ksp), `PCISSetUseStiffnessScaling()`, `PCISSetSubdomainDiagonalScaling()`, `PCISScatterArrayNToVecB()`, `PCISSetSubdomainScalingFactor()`,
`PCISInitialize()`, `PCISApplySchur()`, `PCISApplyInvSchur()`

# External Links
$(_doc_external("Ksp/PCISReset"))
"""
function PCISReset(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCISReset(petsclib::$UnionPetscLib, pc::PC )

    @chk ccall(
               (:PCISReset, $petsc_library),
               PetscErrorCode,
               (PC,),
               pc,
              )


	return nothing
end 

"""
	PCISInitialize(petsclib::PetscLibType,pc::PC) 
initializes the `PC_IS` portion of `PCNN` and `PCBDDC` preconditioner context

Input Parameter:
- `pc` - the `PC` object, must be of type `PCNN` or `PCBDDC`

Level: advanced

-seealso: [](ch_ksp), `PCBDDC`, `PCNN`, `PCISSetUseStiffnessScaling()`, `PCISSetSubdomainDiagonalScaling()`, `PCISScatterArrayNToVecB()`,
`PCISSetSubdomainScalingFactor()`,
`PCISReset()`, `PCISApplySchur()`, `PCISApplyInvSchur()`

# External Links
$(_doc_external("Ksp/PCISInitialize"))
"""
function PCISInitialize(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCISInitialize(petsclib::$UnionPetscLib, pc::PC )

    @chk ccall(
               (:PCISInitialize, $petsc_library),
               PetscErrorCode,
               (PC,),
               pc,
              )


	return nothing
end 

"""
	PCISApplySchur(petsclib::PetscLibType,pc::PC, v::PetscVec, vec1_B::PetscVec, vec2_B::PetscVec, vec1_D::PetscVec, vec2_D::PetscVec) 
applies the Schur complement arising from the `MATIS` inside the `PCNN` preconditioner

Input Parameters:
- `pc`     - preconditioner context
- `v`      - vector to which the Schur complement is to be applied (it is NOT modified inside this function, UNLESS vec2_B is null)
- `vec1_B` - location to store the result of Schur complement applied to chunk
- `vec2_B` - workspace or `NULL`, `v` is used as workspace in that case
- `vec1_D` - work space
- `vec2_D` - work space

Level: advanced

-seealso: [](ch_ksp), `PCBDDC`, `PCNN`, `PCISSetUseStiffnessScaling()`, `PCISSetSubdomainDiagonalScaling()`, `PCISScatterArrayNToVecB()`,
`PCISSetSubdomainScalingFactor()`, `PCISApplyInvSchur()`,
`PCISReset()`, `PCISInitialize()`

# External Links
$(_doc_external("Ksp/PCISApplySchur"))
"""
function PCISApplySchur(petsclib::PetscLibType, pc::PC, v::PetscVec, vec1_B::PetscVec, vec2_B::PetscVec, vec1_D::PetscVec, vec2_D::PetscVec) end

@for_petsc function PCISApplySchur(petsclib::$UnionPetscLib, pc::PC, v::PetscVec, vec1_B::PetscVec, vec2_B::PetscVec, vec1_D::PetscVec, vec2_D::PetscVec )

    @chk ccall(
               (:PCISApplySchur, $petsc_library),
               PetscErrorCode,
               (PC, CVec, CVec, CVec, CVec, CVec),
               pc, v, vec1_B, vec2_B, vec1_D, vec2_D,
              )


	return nothing
end 

"""
	PCISScatterArrayNToVecB(petsclib::PetscLibType,pc::PC, array_N::PetscScalar, v_B::PetscVec, imode::InsertMode, smode::ScatterMode) 
Scatters interface node values from a big array (of all local nodes, interior or interface,
including ghosts) into an interface vector, when in `SCATTER_FORWARD` mode, or vice-versa, when in `SCATTER_REVERSE`
mode.

Input Parameters:
- `pc`      - preconditioner context
- `array_N` - [when in `SCATTER_FORWARD` mode] Array to be scattered into the vector otherwise output array
- `imode`   - insert mode, `ADD_VALUES` or `INSERT_VALUES`
- `smode`   - scatter mode, `SCATTER_FORWARD` or `SCATTER_REVERSE` mode]
- `v_B`     - [when in `SCATTER_REVERSE` mode] Vector to be scattered into the array, otherwise output vector

Level: advanced

-seealso: [](ch_ksp), `PCBDDC`, `PCNN`, `PCISSetUseStiffnessScaling()`, `PCISSetSubdomainDiagonalScaling()`,
`PCISSetSubdomainScalingFactor()`, `PCISApplySchur()`, `PCISApplyInvSchur()`,
`PCISReset()`, `PCISInitialize()`, `InsertMode`

# External Links
$(_doc_external("Ksp/PCISScatterArrayNToVecB"))
"""
function PCISScatterArrayNToVecB(petsclib::PetscLibType, pc::PC, array_N::PetscScalar, v_B::PetscVec, imode::InsertMode, smode::ScatterMode) end

@for_petsc function PCISScatterArrayNToVecB(petsclib::$UnionPetscLib, pc::PC, array_N::$PetscScalar, v_B::PetscVec, imode::InsertMode, smode::ScatterMode )

    @chk ccall(
               (:PCISScatterArrayNToVecB, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{$PetscScalar}, CVec, InsertMode, ScatterMode),
               pc, array_N, v_B, imode, smode,
              )


	return nothing
end 

"""
	PCISApplyInvSchur(petsclib::PetscLibType,pc::PC, b::PetscVec, x::PetscVec, vec1_N::PetscVec, vec2_N::PetscVec) 
Solves the Neumann problem related to applying the inverse of the Schur complement.

Input Parameters:
- `pc`     - preconditioner context
- `b`      - vector of local interface nodes (including ghosts)
- `x`      - vector of local interface nodes (including ghosts); returns the application of the inverse of the Schur complement to `b`
- `vec1_N` - vector of local nodes (interior and interface, including ghosts); used as work space
- `vec2_N` - vector of local nodes (interior and interface, including ghosts); used as work space

Level: advanced

-seealso: [](ch_ksp), `PCBDDC`, `PCNN`, `PCISSetUseStiffnessScaling()`, `PCISSetSubdomainDiagonalScaling()`, `PCISScatterArrayNToVecB()`,
`PCISSetSubdomainScalingFactor()`,
`PCISReset()`, `PCISInitialize()`

# External Links
$(_doc_external("Ksp/PCISApplyInvSchur"))
"""
function PCISApplyInvSchur(petsclib::PetscLibType, pc::PC, b::PetscVec, x::PetscVec, vec1_N::PetscVec, vec2_N::PetscVec) end

@for_petsc function PCISApplyInvSchur(petsclib::$UnionPetscLib, pc::PC, b::PetscVec, x::PetscVec, vec1_N::PetscVec, vec2_N::PetscVec )

    @chk ccall(
               (:PCISApplyInvSchur, $petsc_library),
               PetscErrorCode,
               (PC, CVec, CVec, CVec, CVec),
               pc, b, x, vec1_N, vec2_N,
              )


	return nothing
end 

"""
	PCEisenstatSetOmega(petsclib::PetscLibType,pc::PC, omega::PetscReal) 
Sets the SSOR relaxation coefficient, omega,
to use with Eisenstat's trick (where omega = 1.0 by default)

Logically Collective

Input Parameters:
- `pc`    - the preconditioner context
- `omega` - relaxation coefficient (0 < omega < 2)

Options Database Key:
- `-pc_eisenstat_omega <omega>` - Sets omega

Level: intermediate

-seealso: [](ch_ksp), `PCSORSetOmega()`, `PCEISENSTAT`

# External Links
$(_doc_external("Ksp/PCEisenstatSetOmega"))
"""
function PCEisenstatSetOmega(petsclib::PetscLibType, pc::PC, omega::PetscReal) end

@for_petsc function PCEisenstatSetOmega(petsclib::$UnionPetscLib, pc::PC, omega::$PetscReal )

    @chk ccall(
               (:PCEisenstatSetOmega, $petsc_library),
               PetscErrorCode,
               (PC, $PetscReal),
               pc, omega,
              )


	return nothing
end 

"""
	PCEisenstatSetNoDiagonalScaling(petsclib::PetscLibType,pc::PC, flg::PetscBool) 
Causes the Eisenstat preconditioner, `PCEISENSTAT`
not to do additional diagonal preconditioning. For matrices with a constant
along the diagonal, this may save a small amount of work.

Logically Collective

Input Parameters:
- `pc`  - the preconditioner context
- `flg` - `PETSC_TRUE` turns off diagonal scaling inside the algorithm

Options Database Key:
- `-pc_eisenstat_no_diagonal_scaling` - Activates `PCEisenstatSetNoDiagonalScaling()`

Level: intermediate

-seealso: [](ch_ksp), `PCEisenstatSetOmega()`, `PCEISENSTAT`

# External Links
$(_doc_external("Ksp/PCEisenstatSetNoDiagonalScaling"))
"""
function PCEisenstatSetNoDiagonalScaling(petsclib::PetscLibType, pc::PC, flg::PetscBool) end

@for_petsc function PCEisenstatSetNoDiagonalScaling(petsclib::$UnionPetscLib, pc::PC, flg::PetscBool )

    @chk ccall(
               (:PCEisenstatSetNoDiagonalScaling, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, flg,
              )


	return nothing
end 

"""
	omega::PetscReal = PCEisenstatGetOmega(petsclib::PetscLibType,pc::PC) 
Gets the SSOR relaxation coefficient, omega,
to use with Eisenstat's trick (where omega = 1.0 by default).

Logically Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `omega` - relaxation coefficient (0 < omega < 2)

Options Database Key:
- `-pc_eisenstat_omega <omega>` - Sets omega

-seealso: [](ch_ksp), `PCEISENSTAT`, `PCSORGetOmega()`, `PCEisenstatSetOmega()`

# External Links
$(_doc_external("Ksp/PCEisenstatGetOmega"))
"""
function PCEisenstatGetOmega(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCEisenstatGetOmega(petsclib::$UnionPetscLib, pc::PC )
	omega_ = Ref{$PetscReal}()

    @chk ccall(
               (:PCEisenstatGetOmega, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{$PetscReal}),
               pc, omega_,
              )

	omega = omega_[]

	return omega
end 

"""
	flg::PetscBool = PCEisenstatGetNoDiagonalScaling(petsclib::PetscLibType,pc::PC) 
Tells if the Eisenstat preconditioner
not to do additional diagonal preconditioning. For matrices with a constant
along the diagonal, this may save a small amount of work.

Logically Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `flg` - `PETSC_TRUE` means there is no diagonal scaling applied

Options Database Key:
- `-pc_eisenstat_no_diagonal_scaling` - Activates `PCEisenstatSetNoDiagonalScaling()`

Level: intermediate

-seealso: , `PCEISENSTAT`, `PCEisenstatGetOmega()`

# External Links
$(_doc_external("Ksp/PCEisenstatGetNoDiagonalScaling"))
"""
function PCEisenstatGetNoDiagonalScaling(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCEisenstatGetNoDiagonalScaling(petsclib::$UnionPetscLib, pc::PC )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PCEisenstatGetNoDiagonalScaling, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PetscBool}),
               pc, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	PCJacobiSetUseAbs(petsclib::PetscLibType,pc::PC, flg::PetscBool) 
Causes the Jacobi preconditioner `PCJACOBI` to use the
absolute values of the diagonal divisors in the preconditioner

Logically Collective

Input Parameters:
- `pc`  - the preconditioner context
- `flg` - whether to use absolute values or not

Options Database Key:
- `-pc_jacobi_abs <bool>` - use absolute values

-seealso: [](ch_ksp), `PCJACOBI`, `PCJacobiSetType()`, `PCJacobiGetUseAbs()`

# External Links
$(_doc_external("Ksp/PCJacobiSetUseAbs"))
"""
function PCJacobiSetUseAbs(petsclib::PetscLibType, pc::PC, flg::PetscBool) end

@for_petsc function PCJacobiSetUseAbs(petsclib::$UnionPetscLib, pc::PC, flg::PetscBool )

    @chk ccall(
               (:PCJacobiSetUseAbs, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, flg,
              )


	return nothing
end 

"""
	flg::PetscBool = PCJacobiGetUseAbs(petsclib::PetscLibType,pc::PC) 
Determines if the Jacobi preconditioner `PCJACOBI` uses the
absolute values of the diagonal divisors in the preconditioner

Logically Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `flg` - whether to use absolute values or not

Level: intermediate

-seealso: [](ch_ksp), `PCJACOBI`, `PCJacobiSetType()`, `PCJacobiSetUseAbs()`, `PCJacobiGetType()`

# External Links
$(_doc_external("Ksp/PCJacobiGetUseAbs"))
"""
function PCJacobiGetUseAbs(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCJacobiGetUseAbs(petsclib::$UnionPetscLib, pc::PC )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PCJacobiGetUseAbs, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PetscBool}),
               pc, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	PCJacobiSetRowl1Scale(petsclib::PetscLibType,pc::PC, scale::PetscReal) 
Set scaling of off
Remark 6.1 in "Multigrid Smoothers for Ultraparallel Computing", Baker et al, with 0.5 scaling

Logically Collective

Input Parameters:
- `pc`    - the preconditioner context
- `scale` - scaling

Options Database Key:
- `-pc_jacobi_rowl1_scale <real>` - use absolute values

Level: intermediate

-seealso: [](ch_ksp), `PCJACOBI`, `PCJacobiSetType()`, `PCJacobiGetRowl1Scale()`

# External Links
$(_doc_external("Ksp/PCJacobiSetRowl1Scale"))
"""
function PCJacobiSetRowl1Scale(petsclib::PetscLibType, pc::PC, scale::PetscReal) end

@for_petsc function PCJacobiSetRowl1Scale(petsclib::$UnionPetscLib, pc::PC, scale::$PetscReal )

    @chk ccall(
               (:PCJacobiSetRowl1Scale, $petsc_library),
               PetscErrorCode,
               (PC, $PetscReal),
               pc, scale,
              )


	return nothing
end 

"""
	scale::PetscReal = PCJacobiGetRowl1Scale(petsclib::PetscLibType,pc::PC) 
Get scaling of off

Logically Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `scale` - scaling

Level: intermediate

-seealso: [](ch_ksp), `PCJACOBI`, `PCJacobiSetType()`, `PCJacobiSetRowl1Scale()`, `PCJacobiGetType()`

# External Links
$(_doc_external("Ksp/PCJacobiGetRowl1Scale"))
"""
function PCJacobiGetRowl1Scale(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCJacobiGetRowl1Scale(petsclib::$UnionPetscLib, pc::PC )
	scale_ = Ref{$PetscReal}()

    @chk ccall(
               (:PCJacobiGetRowl1Scale, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{$PetscReal}),
               pc, scale_,
              )

	scale = scale_[]

	return scale
end 

"""
	PCJacobiSetFixDiagonal(petsclib::PetscLibType,pc::PC, flg::PetscBool) 
Check for zero values on the diagonal and replace them with 1.0

Logically Collective

Input Parameters:
- `pc`  - the preconditioner context
- `flg` - the boolean flag

Options Database Key:
- `-pc_jacobi_fixdiagonal <bool>` - check for zero values on the diagonal

-seealso: [](ch_ksp), `PCJACOBI`, `PCJacobiSetType()`, `PCJacobiGetFixDiagonal()`, `PCJacobiSetUseAbs()`

# External Links
$(_doc_external("Ksp/PCJacobiSetFixDiagonal"))
"""
function PCJacobiSetFixDiagonal(petsclib::PetscLibType, pc::PC, flg::PetscBool) end

@for_petsc function PCJacobiSetFixDiagonal(petsclib::$UnionPetscLib, pc::PC, flg::PetscBool )

    @chk ccall(
               (:PCJacobiSetFixDiagonal, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, flg,
              )


	return nothing
end 

"""
	flg::PetscBool = PCJacobiGetFixDiagonal(petsclib::PetscLibType,pc::PC) 
Determines if the Jacobi preconditioner `PCJACOBI` checks for zero diagonal terms

Logically Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `flg` - the boolean flag

Options Database Key:
- `-pc_jacobi_fixdiagonal <bool>` - Fix 0 terms on diagonal by using 1

Level: intermediate

-seealso: [](ch_ksp), `PCJACOBI`, `PCJacobiSetType()`, `PCJacobiSetFixDiagonal()`

# External Links
$(_doc_external("Ksp/PCJacobiGetFixDiagonal"))
"""
function PCJacobiGetFixDiagonal(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCJacobiGetFixDiagonal(petsclib::$UnionPetscLib, pc::PC )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PCJacobiGetFixDiagonal, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PetscBool}),
               pc, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	PCJacobiGetDiagonal(petsclib::PetscLibType,pc::PC, diagonal::PetscVec, diagonal_sqrt::PetscVec) 
Returns copy of the diagonal and/or diagonal squareroot `Vec`

Logically Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameters:
- `diagonal`      - Copy of `Vec` of the inverted diagonal
- `diagonal_sqrt` - Copy of `Vec` of the inverted square root diagonal

Level: developer

-seealso: [](ch_ksp), `PCJACOBI`, `PCJacobiSetType()`

# External Links
$(_doc_external("Ksp/PCJacobiGetDiagonal"))
"""
function PCJacobiGetDiagonal(petsclib::PetscLibType, pc::PC, diagonal::PetscVec, diagonal_sqrt::PetscVec) end

@for_petsc function PCJacobiGetDiagonal(petsclib::$UnionPetscLib, pc::PC, diagonal::PetscVec, diagonal_sqrt::PetscVec )

    @chk ccall(
               (:PCJacobiGetDiagonal, $petsc_library),
               PetscErrorCode,
               (PC, CVec, CVec),
               pc, diagonal, diagonal_sqrt,
              )


	return nothing
end 

"""
	PCJacobiSetType(petsclib::PetscLibType,pc::PC, type::PCJacobiType) 
Causes the Jacobi preconditioner to use either the diagonal, the maximum entry in each row,
of the sum of rows entries for the diagonal preconditioner

Logically Collective

Input Parameters:
- `pc`   - the preconditioner context
- `type` - `PC_JACOBI_DIAGONAL`, `PC_JACOBI_ROWL1`, `PC_JACOBI_ROWMAX`, `PC_JACOBI_ROWSUM`

Options Database Key:
- `-pc_jacobi_type <diagonal,rowl1,rowmax,rowsum>` - the type of diagonal matrix to use for Jacobi

Level: intermediate

-seealso: [](ch_ksp), `PCJACOBI`, `PCJacobiSetUseAbs()`, `PCJacobiGetType()`

# External Links
$(_doc_external("Ksp/PCJacobiSetType"))
"""
function PCJacobiSetType(petsclib::PetscLibType, pc::PC, type::PCJacobiType) end

@for_petsc function PCJacobiSetType(petsclib::$UnionPetscLib, pc::PC, type::PCJacobiType )

    @chk ccall(
               (:PCJacobiSetType, $petsc_library),
               PetscErrorCode,
               (PC, PCJacobiType),
               pc, type,
              )


	return nothing
end 

"""
	type::PCJacobiType = PCJacobiGetType(petsclib::PetscLibType,pc::PC) 
Gets how the diagonal matrix is produced for the preconditioner

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `type` - `PC_JACOBI_DIAGONAL`, `PC_JACOBI_ROWL1`, `PC_JACOBI_ROWMAX`, `PC_JACOBI_ROWSUM`

Level: intermediate

-seealso: [](ch_ksp), `PCJACOBI`, `PCJacobiSetUseAbs()`, `PCJacobiSetType()`

# External Links
$(_doc_external("Ksp/PCJacobiGetType"))
"""
function PCJacobiGetType(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCJacobiGetType(petsclib::$UnionPetscLib, pc::PC )
	type_ = Ref{PCJacobiType}()

    @chk ccall(
               (:PCJacobiGetType, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PCJacobiType}),
               pc, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	PCGalerkinSetRestriction(petsclib::PetscLibType,pc::PC, R::PetscMat) 
Sets the restriction operator for the `PCGALERKIN` preconditioner

Logically Collective

Input Parameters:
- `pc` - the preconditioner context
- `R`  - the restriction operator

Level: intermediate

-seealso: [](ch_ksp), `PC`, `PCCreate()`, `PCSetType()`, `PCType`, `PCGALERKIN`,
`PCGalerkinSetInterpolation()`, `PCGalerkinGetKSP()`

# External Links
$(_doc_external("Ksp/PCGalerkinSetRestriction"))
"""
function PCGalerkinSetRestriction(petsclib::PetscLibType, pc::PC, R::PetscMat) end

@for_petsc function PCGalerkinSetRestriction(petsclib::$UnionPetscLib, pc::PC, R::PetscMat )

    @chk ccall(
               (:PCGalerkinSetRestriction, $petsc_library),
               PetscErrorCode,
               (PC, CMat),
               pc, R,
              )


	return nothing
end 

"""
	PCGalerkinSetInterpolation(petsclib::PetscLibType,pc::PC, P::PetscMat) 
Sets the interpolation operator for the `PCGALERKIN` preconditioner

Logically Collective

Input Parameters:
- `pc` - the preconditioner context
- `P`  - the interpolation operator

Level: intermediate

-seealso: [](ch_ksp), `PC`, `PCCreate()`, `PCSetType()`, `PCType`, `PCGALERKIN`,
`PCGalerkinSetRestriction()`, `PCGalerkinGetKSP()`

# External Links
$(_doc_external("Ksp/PCGalerkinSetInterpolation"))
"""
function PCGalerkinSetInterpolation(petsclib::PetscLibType, pc::PC, P::PetscMat) end

@for_petsc function PCGalerkinSetInterpolation(petsclib::$UnionPetscLib, pc::PC, P::PetscMat )

    @chk ccall(
               (:PCGalerkinSetInterpolation, $petsc_library),
               PetscErrorCode,
               (PC, CMat),
               pc, P,
              )


	return nothing
end 

"""
	PCGalerkinSetComputeSubmatrix(petsclib::PetscLibType,pc::PC, computeAsub::external, ctx::Cvoid) 
Provide a routine that will be called to compute the Galerkin submatrix

Logically Collective

Input Parameters:
- `pc`          - the preconditioner context
- `computeAsub` - routine that computes the submatrix from the global matrix
- `ctx`         - context used by the routine, or `NULL`

Calling sequence of `computeAsub`:
- `pc`  - the `PCGALERKIN` preconditioner
- `A`   - the matrix in the `PCGALERKIN`
- `Ap`  - the computed submatrix from any previous computation, if `NULL` it has not previously been computed
- `cAp` - the submatrix computed by this routine
- `ctx` - optional user-defined function context

Level: intermediate

-seealso: [](ch_ksp), `PC`, `PCCreate()`, `PCSetType()`, `PCType`, `PCGALERKIN`,
`PCGalerkinSetRestriction()`, `PCGalerkinSetInterpolation()`, `PCGalerkinGetKSP()`

# External Links
$(_doc_external("Ksp/PCGalerkinSetComputeSubmatrix"))
"""
function PCGalerkinSetComputeSubmatrix(petsclib::PetscLibType, pc::PC, computeAsub::external, ctx::Cvoid) end

@for_petsc function PCGalerkinSetComputeSubmatrix(petsclib::$UnionPetscLib, pc::PC, computeAsub::external, ctx::Cvoid )

    @chk ccall(
               (:PCGalerkinSetComputeSubmatrix, $petsc_library),
               PetscErrorCode,
               (PC, external, Ptr{Cvoid}),
               pc, computeAsub, ctx,
              )


	return nothing
end 

"""
	PCGalerkinGetKSP(petsclib::PetscLibType,pc::PC, ksp::PetscKSP) 
Gets the `KSP` object in the `PCGALERKIN`

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `ksp` - the `KSP` object

Level: intermediate

-seealso: [](ch_ksp), `PC`, `PCCreate()`, `PCSetType()`, `PCType`, `PCGALERKIN`,
`PCGalerkinSetRestriction()`, `PCGalerkinSetInterpolation()`, `PCGalerkinSetComputeSubmatrix()`

# External Links
$(_doc_external("Ksp/PCGalerkinGetKSP"))
"""
function PCGalerkinGetKSP(petsclib::PetscLibType, pc::PC, ksp::PetscKSP) end

@for_petsc function PCGalerkinGetKSP(petsclib::$UnionPetscLib, pc::PC, ksp::PetscKSP )
	ksp_ = Ref(ksp.ptr)

    @chk ccall(
               (:PCGalerkinGetKSP, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{CKSP}),
               pc, ksp_,
              )

	ksp.ptr = C_NULL

	return nothing
end 

"""
	PCPARMSSetGlobal(petsclib::PetscLibType,pc::PC, type::PCPARMSGlobalType) 
Sets the global preconditioner to be used in `PCPARMS`.

Collective

Input Parameters:
- `pc`   - the preconditioner context
- `type` - the global preconditioner type, one of
-seealso: [](ch_ksp), `PCPARMS`, `PCPARMSSetLocal()`

# External Links
$(_doc_external("Ksp/PCPARMSSetGlobal"))
"""
function PCPARMSSetGlobal(petsclib::PetscLibType, pc::PC, type::PCPARMSGlobalType) end

@for_petsc function PCPARMSSetGlobal(petsclib::$UnionPetscLib, pc::PC, type::PCPARMSGlobalType )

    @chk ccall(
               (:PCPARMSSetGlobal, $petsc_library),
               PetscErrorCode,
               (PC, PCPARMSGlobalType),
               pc, type,
              )


	return nothing
end 

"""
	PCPARMSSetLocal(petsclib::PetscLibType,pc::PC, type::PCPARMSLocalType) 
Sets the local preconditioner to be used in `PCPARMS`.

Collective

Input Parameters:
- `pc`   - the preconditioner context
- `type` - the local preconditioner type, one of
-seealso: [](ch_ksp), `PCPARMS`, `PCPARMSSetGlobal()`, `PCPARMSSetNonsymPerm()`


# External Links
$(_doc_external("Ksp/PCPARMSSetLocal"))
"""
function PCPARMSSetLocal(petsclib::PetscLibType, pc::PC, type::PCPARMSLocalType) end

@for_petsc function PCPARMSSetLocal(petsclib::$UnionPetscLib, pc::PC, type::PCPARMSLocalType )

    @chk ccall(
               (:PCPARMSSetLocal, $petsc_library),
               PetscErrorCode,
               (PC, PCPARMSLocalType),
               pc, type,
              )


	return nothing
end 

"""
	PCPARMSSetSolveTolerances(petsclib::PetscLibType,pc::PC, tol::PetscReal, maxits::PetscInt) 
Sets the convergence tolerance and the maximum iterations for the
inner GMRES solver, when the Schur global preconditioner is used.

Collective

Input Parameters:
- `pc`     - the preconditioner context
- `tol`    - the convergence tolerance
- `maxits` - the maximum number of iterations to use

Options Database Keys:
- `-pc_parms_solve_tol` - set the tolerance for local solve
- `-pc_parms_max_it`    - set the maximum number of inner iterations

Level: intermediate

-seealso: [](ch_ksp), `PCPARMS`, `PCPARMSSetSolveRestart()`

# External Links
$(_doc_external("Ksp/PCPARMSSetSolveTolerances"))
"""
function PCPARMSSetSolveTolerances(petsclib::PetscLibType, pc::PC, tol::PetscReal, maxits::PetscInt) end

@for_petsc function PCPARMSSetSolveTolerances(petsclib::$UnionPetscLib, pc::PC, tol::$PetscReal, maxits::$PetscInt )

    @chk ccall(
               (:PCPARMSSetSolveTolerances, $petsc_library),
               PetscErrorCode,
               (PC, $PetscReal, $PetscInt),
               pc, tol, maxits,
              )


	return nothing
end 

"""
	PCPARMSSetSolveRestart(petsclib::PetscLibType,pc::PC, restart::PetscInt) 
Sets the number of iterations at which the
inner GMRES solver restarts.

Collective

Input Parameters:
- `pc`      - the preconditioner context
- `restart` - maximum dimension of the Krylov subspace

Options Database Key:
- `-pc_parms_max_dim` - sets the inner Krylov dimension

Level: intermediate

-seealso: [](ch_ksp), `PCPARMS`, `PCPARMSSetSolveTolerances()`

# External Links
$(_doc_external("Ksp/PCPARMSSetSolveRestart"))
"""
function PCPARMSSetSolveRestart(petsclib::PetscLibType, pc::PC, restart::PetscInt) end

@for_petsc function PCPARMSSetSolveRestart(petsclib::$UnionPetscLib, pc::PC, restart::$PetscInt )

    @chk ccall(
               (:PCPARMSSetSolveRestart, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, restart,
              )


	return nothing
end 

"""
	PCPARMSSetNonsymPerm(petsclib::PetscLibType,pc::PC, nonsym::PetscBool) 
Sets the type of permutation for the ARMS preconditioner: the standard
symmetric ARMS or the non-symmetric ARMS (ARMS-ddPQ).

Collective

Input Parameters:
- `pc`     - the preconditioner context
- `nonsym` - `PETSC_TRUE` indicates the non-symmetric ARMS is used;
`PETSC_FALSE` indicates the symmetric ARMS is used

Options Database Key:
- `-pc_parms_nonsymmetric_perm` - sets the use of nonsymmetric permutation

Level: intermediate

-seealso: [](ch_ksp), `PCPARMS`

# External Links
$(_doc_external("Ksp/PCPARMSSetNonsymPerm"))
"""
function PCPARMSSetNonsymPerm(petsclib::PetscLibType, pc::PC, nonsym::PetscBool) end

@for_petsc function PCPARMSSetNonsymPerm(petsclib::$UnionPetscLib, pc::PC, nonsym::PetscBool )

    @chk ccall(
               (:PCPARMSSetNonsymPerm, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, nonsym,
              )


	return nothing
end 

"""
	PCPARMSSetFill(petsclib::PetscLibType,pc::PC, lfil0::PetscInt, lfil1::PetscInt, lfil2::PetscInt) 
Sets the fill
Consider the original matrix A = [B F; E C] and the approximate version
M = [LB 0; E/UB I]*[UB LB F; 0 S].

Collective

Input Parameters:
- `pc`    - the preconditioner context
- `lfil0` - the level of fill-in kept in LB, UB, E/UB and LB F
- `lfil1` - the level of fill-in kept in S
- `lfil2` - the level of fill-in kept in the L and U parts of the LU factorization of S

Options Database Keys:
- `-pc_parms_lfil_ilu_arms` - set the amount of fill-in for ilut, iluk and arms
- `-pc_parms_lfil_schur`    - set the amount of fill-in for schur
- `-pc_parms_lfil_ilut_L_U` - set the amount of fill-in for ILUT L and U

Level: intermediate

-seealso: [](ch_ksp), `PCPARMS`

# External Links
$(_doc_external("Ksp/PCPARMSSetFill"))
"""
function PCPARMSSetFill(petsclib::PetscLibType, pc::PC, lfil0::PetscInt, lfil1::PetscInt, lfil2::PetscInt) end

@for_petsc function PCPARMSSetFill(petsclib::$UnionPetscLib, pc::PC, lfil0::$PetscInt, lfil1::$PetscInt, lfil2::$PetscInt )

    @chk ccall(
               (:PCPARMSSetFill, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, $PetscInt, $PetscInt),
               pc, lfil0, lfil1, lfil2,
              )


	return nothing
end 

"""
	PCExoticSetType(petsclib::PetscLibType,pc::PC, type::PCExoticType) 
Sets the type of coarse grid interpolation to use

Logically Collective

Input Parameters:
- `pc`   - the preconditioner context
- `type` - either `PC_EXOTIC_FACE` or `PC_EXOTIC_WIREBASKET` (defaults to face)

Options Database Keys:
- `-pc_exotic_type <face,wirebasket>` - use a coarse grid point for each face, or edge and vertex

-seealso: [](ch_ksp), `PCEXOTIC`, `PCExoticType()`

# External Links
$(_doc_external("Ksp/PCExoticSetType"))
"""
function PCExoticSetType(petsclib::PetscLibType, pc::PC, type::PCExoticType) end

@for_petsc function PCExoticSetType(petsclib::$UnionPetscLib, pc::PC, type::PCExoticType )

    @chk ccall(
               (:PCExoticSetType, $petsc_library),
               PetscErrorCode,
               (PC, PCExoticType),
               pc, type,
              )


	return nothing
end 

"""
	PCFactorSetUpMatSolverType(petsclib::PetscLibType,pc::PC) 
Can be called after `KSPSetOperators()` or `PCSetOperators()`, causes `MatGetFactor()` to be called so then one may
set the options for that particular factorization object.

Input Parameter:
- `pc` - the preconditioner context

-seealso: [](ch_ksp), `PCCHOLESKY`, `PCLU`, `PCFactorSetMatSolverType()`, `PCFactorGetMatrix()`

# External Links
$(_doc_external("Ksp/PCFactorSetUpMatSolverType"))
"""
function PCFactorSetUpMatSolverType(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCFactorSetUpMatSolverType(petsclib::$UnionPetscLib, pc::PC )

    @chk ccall(
               (:PCFactorSetUpMatSolverType, $petsc_library),
               PetscErrorCode,
               (PC,),
               pc,
              )


	return nothing
end 

"""
	PCFactorSetZeroPivot(petsclib::PetscLibType,pc::PC, zero::PetscReal) 
Sets the size at which smaller pivots are declared to be zero

Logically Collective

Input Parameters:
- `pc`   - the preconditioner context
- `zero` - all pivots smaller than this will be considered zero

Options Database Key:
- `-pc_factor_zeropivot <zero>` - Sets tolerance for what is considered a zero pivot

Level: intermediate

-seealso: [](ch_ksp), `PCCHOLESKY`, `PCLU`, `PCFactorSetShiftType()`, `PCFactorSetShiftAmount()`

# External Links
$(_doc_external("Ksp/PCFactorSetZeroPivot"))
"""
function PCFactorSetZeroPivot(petsclib::PetscLibType, pc::PC, zero::PetscReal) end

@for_petsc function PCFactorSetZeroPivot(petsclib::$UnionPetscLib, pc::PC, zero::$PetscReal )

    @chk ccall(
               (:PCFactorSetZeroPivot, $petsc_library),
               PetscErrorCode,
               (PC, $PetscReal),
               pc, zero,
              )


	return nothing
end 

"""
	PCFactorSetShiftType(petsclib::PetscLibType,pc::PC, shifttype::MatFactorShiftType) 
adds a particular type of quantity to the diagonal of the matrix during
numerical factorization, thus the matrix has nonzero pivots

Logically Collective

Input Parameters:
- `pc`        - the preconditioner context
- `shifttype` - type of shift; one of `MAT_SHIFT_NONE`, `MAT_SHIFT_NONZERO`, `MAT_SHIFT_POSITIVE_DEFINITE`, `MAT_SHIFT_INBLOCKS`

Options Database Key:
- `-pc_factor_shift_type <shifttype>` - Sets shift type; use '-help' for a list of available types

Level: intermediate

-seealso: [](ch_ksp), `PCCHOLESKY`, `PCLU`, `PCFactorSetZeroPivot()`, `PCFactorSetShiftAmount()`

# External Links
$(_doc_external("Ksp/PCFactorSetShiftType"))
"""
function PCFactorSetShiftType(petsclib::PetscLibType, pc::PC, shifttype::MatFactorShiftType) end

@for_petsc function PCFactorSetShiftType(petsclib::$UnionPetscLib, pc::PC, shifttype::MatFactorShiftType )

    @chk ccall(
               (:PCFactorSetShiftType, $petsc_library),
               PetscErrorCode,
               (PC, MatFactorShiftType),
               pc, shifttype,
              )


	return nothing
end 

"""
	PCFactorSetShiftAmount(petsclib::PetscLibType,pc::PC, shiftamount::PetscReal) 
adds a quantity to the diagonal of the matrix during
numerical factorization, thus the matrix has nonzero pivots

Logically Collective

Input Parameters:
- `pc`          - the preconditioner context
- `shiftamount` - amount of shift or `PETSC_DECIDE` for the default

Options Database Key:
- `-pc_factor_shift_amount <shiftamount>` - Sets shift amount or -1 for the default

Level: intermediate

-seealso: [](ch_ksp), `PCCHOLESKY`, `PCLU`, `PCFactorSetZeroPivot()`, `PCFactorSetShiftType()`

# External Links
$(_doc_external("Ksp/PCFactorSetShiftAmount"))
"""
function PCFactorSetShiftAmount(petsclib::PetscLibType, pc::PC, shiftamount::PetscReal) end

@for_petsc function PCFactorSetShiftAmount(petsclib::$UnionPetscLib, pc::PC, shiftamount::$PetscReal )

    @chk ccall(
               (:PCFactorSetShiftAmount, $petsc_library),
               PetscErrorCode,
               (PC, $PetscReal),
               pc, shiftamount,
              )


	return nothing
end 

"""
	PCFactorSetDropTolerance(petsclib::PetscLibType,pc::PC, dt::PetscReal, dtcol::PetscReal, maxrowcount::PetscInt) 
The preconditioner will use an `PCILU`
based on a drop tolerance.

Logically Collective

Input Parameters:
- `pc`          - the preconditioner context
- `dt`          - the drop tolerance, try from 1.e-10 to .1
- `dtcol`       - tolerance for column pivot, good values [0.1 to 0.01]
- `maxrowcount` - the max number of nonzeros allowed in a row, best value
depends on the number of nonzeros in row of original matrix

Options Database Key:
- `-pc_factor_drop_tolerance <dt,dtcol,maxrowcount>` - Sets drop tolerance

Level: intermediate

-seealso: [](ch_ksp), `PCILU`

# External Links
$(_doc_external("Ksp/PCFactorSetDropTolerance"))
"""
function PCFactorSetDropTolerance(petsclib::PetscLibType, pc::PC, dt::PetscReal, dtcol::PetscReal, maxrowcount::PetscInt) end

@for_petsc function PCFactorSetDropTolerance(petsclib::$UnionPetscLib, pc::PC, dt::$PetscReal, dtcol::$PetscReal, maxrowcount::$PetscInt )

    @chk ccall(
               (:PCFactorSetDropTolerance, $petsc_library),
               PetscErrorCode,
               (PC, $PetscReal, $PetscReal, $PetscInt),
               pc, dt, dtcol, maxrowcount,
              )


	return nothing
end 

"""
	pivot::PetscReal = PCFactorGetZeroPivot(petsclib::PetscLibType,pc::PC) 
Gets the tolerance used to define a zero privot

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `pivot` - the tolerance

Level: intermediate

-seealso: [](ch_ksp), `PCLU`, `PCCHOLESKY`, `PCFactorSetZeroPivot()`

# External Links
$(_doc_external("Ksp/PCFactorGetZeroPivot"))
"""
function PCFactorGetZeroPivot(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCFactorGetZeroPivot(petsclib::$UnionPetscLib, pc::PC )
	pivot_ = Ref{$PetscReal}()

    @chk ccall(
               (:PCFactorGetZeroPivot, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{$PetscReal}),
               pc, pivot_,
              )

	pivot = pivot_[]

	return pivot
end 

"""
	shift::PetscReal = PCFactorGetShiftAmount(petsclib::PetscLibType,pc::PC) 
Gets the tolerance used to define a zero privot

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `shift` - how much to shift the diagonal entry

Level: intermediate

-seealso: [](ch_ksp), `PCLU`, `PCCHOLESKY`, `PCFactorSetShiftAmount()`, `PCFactorSetShiftType()`, `PCFactorGetShiftType()`

# External Links
$(_doc_external("Ksp/PCFactorGetShiftAmount"))
"""
function PCFactorGetShiftAmount(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCFactorGetShiftAmount(petsclib::$UnionPetscLib, pc::PC )
	shift_ = Ref{$PetscReal}()

    @chk ccall(
               (:PCFactorGetShiftAmount, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{$PetscReal}),
               pc, shift_,
              )

	shift = shift_[]

	return shift
end 

"""
	type::MatFactorShiftType = PCFactorGetShiftType(petsclib::PetscLibType,pc::PC) 
Gets the type of shift, if any, done when a zero pivot is detected

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `type` - one of `MAT_SHIFT_NONE`, `MAT_SHIFT_NONZERO`, `MAT_SHIFT_POSITIVE_DEFINITE`, or `MAT_SHIFT_INBLOCKS`

Level: intermediate

-seealso: [](ch_ksp), `PCLU`, `PCCHOLESKY`, `PCFactorSetShiftType()`, `MatFactorShiftType`, `PCFactorSetShiftAmount()`, `PCFactorGetShiftAmount()`

# External Links
$(_doc_external("Ksp/PCFactorGetShiftType"))
"""
function PCFactorGetShiftType(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCFactorGetShiftType(petsclib::$UnionPetscLib, pc::PC )
	type_ = Ref{MatFactorShiftType}()

    @chk ccall(
               (:PCFactorGetShiftType, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{MatFactorShiftType}),
               pc, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	levels::PetscInt = PCFactorGetLevels(petsclib::PetscLibType,pc::PC) 
Gets the number of levels of fill to use.

Logically Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `levels` - number of levels of fill

Level: intermediate

-seealso: [](ch_ksp), `PCILU`, `PCICC`, `PCFactorSetLevels()`

# External Links
$(_doc_external("Ksp/PCFactorGetLevels"))
"""
function PCFactorGetLevels(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCFactorGetLevels(petsclib::$UnionPetscLib, pc::PC )
	levels_ = Ref{$PetscInt}()

    @chk ccall(
               (:PCFactorGetLevels, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{$PetscInt}),
               pc, levels_,
              )

	levels = levels_[]

	return levels
end 

"""
	PCFactorSetLevels(petsclib::PetscLibType,pc::PC, levels::PetscInt) 
Sets the number of levels of fill to use.

Logically Collective

Input Parameters:
- `pc`     - the preconditioner context
- `levels` - number of levels of fill

Options Database Key:
- `-pc_factor_levels <levels>` - Sets fill level

Level: intermediate

-seealso: [](ch_ksp), `PCILU`, `PCICC`, `PCFactorGetLevels()`

# External Links
$(_doc_external("Ksp/PCFactorSetLevels"))
"""
function PCFactorSetLevels(petsclib::PetscLibType, pc::PC, levels::PetscInt) end

@for_petsc function PCFactorSetLevels(petsclib::$UnionPetscLib, pc::PC, levels::$PetscInt )

    @chk ccall(
               (:PCFactorSetLevels, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, levels,
              )


	return nothing
end 

"""
	PCFactorSetAllowDiagonalFill(petsclib::PetscLibType,pc::PC, flg::PetscBool) 
Causes all diagonal matrix entries to be
treated as level 0 fill even if there is no non-zero location.

Logically Collective

Input Parameters:
- `pc`  - the preconditioner context
- `flg` - `PETSC_TRUE` to turn on, `PETSC_FALSE` to turn off

Options Database Key:
- `-pc_factor_diagonal_fill <bool>` - allow the diagonal fill

-seealso: [](ch_ksp), `PCILU`, `PCICC`, `PCFactorGetAllowDiagonalFill()`

# External Links
$(_doc_external("Ksp/PCFactorSetAllowDiagonalFill"))
"""
function PCFactorSetAllowDiagonalFill(petsclib::PetscLibType, pc::PC, flg::PetscBool) end

@for_petsc function PCFactorSetAllowDiagonalFill(petsclib::$UnionPetscLib, pc::PC, flg::PetscBool )

    @chk ccall(
               (:PCFactorSetAllowDiagonalFill, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, flg,
              )


	return nothing
end 

"""
	flg::PetscBool = PCFactorGetAllowDiagonalFill(petsclib::PetscLibType,pc::PC) 
Determines if all diagonal matrix entries are
treated as level 0 fill even if there is no non-zero location.

Logically Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `flg` - `PETSC_TRUE` to turn on, `PETSC_FALSE` to turn off

-seealso: [](ch_ksp), `PCILU`, `PCICC`, `PCFactorSetAllowDiagonalFill()`

# External Links
$(_doc_external("Ksp/PCFactorGetAllowDiagonalFill"))
"""
function PCFactorGetAllowDiagonalFill(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCFactorGetAllowDiagonalFill(petsclib::$UnionPetscLib, pc::PC )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PCFactorGetAllowDiagonalFill, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PetscBool}),
               pc, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	PCFactorReorderForNonzeroDiagonal(petsclib::PetscLibType,pc::PC, rtol::PetscReal) 
reorders rows/columns of matrix to remove zeros from diagonal

Logically Collective

Input Parameters:
- `pc`   - the preconditioner context
- `rtol` - diagonal entries smaller than this in absolute value are considered zero

Options Database Key:
- `-pc_factor_nonzeros_along_diagonal <tol>` - perform the reordering with the given tolerance

Level: intermediate

-seealso: [](ch_ksp), `PCILU`, `PCICC`, `PCFactorSetFill()`, `PCFactorSetShiftAmount()`, `PCFactorSetZeroPivot()`, `MatReorderForNonzeroDiagonal()`

# External Links
$(_doc_external("Ksp/PCFactorReorderForNonzeroDiagonal"))
"""
function PCFactorReorderForNonzeroDiagonal(petsclib::PetscLibType, pc::PC, rtol::PetscReal) end

@for_petsc function PCFactorReorderForNonzeroDiagonal(petsclib::$UnionPetscLib, pc::PC, rtol::$PetscReal )

    @chk ccall(
               (:PCFactorReorderForNonzeroDiagonal, $petsc_library),
               PetscErrorCode,
               (PC, $PetscReal),
               pc, rtol,
              )


	return nothing
end 

"""
	PCFactorSetMatSolverType(petsclib::PetscLibType,pc::PC, stype::MatSolverType) 
sets the solver package that is used to perform the factorization

Logically Collective

Input Parameters:
- `pc`    - the preconditioner context
- `stype` - for example, `MATSOLVERSUPERLU`, `MATSOLVERSUPERLU_DIST`, `MATSOLVERMUMPS`

Options Database Key:
- `-pc_factor_mat_solver_type <stype>` - petsc, superlu, superlu_dist, mumps, cusparse

Level: intermediate

-seealso: [](ch_ksp), `PCLU`, `PCCHOLESKY`, `MatGetFactor()`, `MatSolverType`, `PCFactorGetMatSolverType()`, `MatSolverTypeRegister()`,
`MatInitializePackage()`, `MATSOLVERSUPERLU`, `MATSOLVERSUPERLU_DIST`, `MATSOLVERMUMPS`, `MatSolverTypeGet()`

# External Links
$(_doc_external("Ksp/PCFactorSetMatSolverType"))
"""
function PCFactorSetMatSolverType(petsclib::PetscLibType, pc::PC, stype::MatSolverType) end

@for_petsc function PCFactorSetMatSolverType(petsclib::$UnionPetscLib, pc::PC, stype::MatSolverType )

    @chk ccall(
               (:PCFactorSetMatSolverType, $petsc_library),
               PetscErrorCode,
               (PC, MatSolverType),
               pc, stype,
              )


	return nothing
end 

"""
	stype::MatSolverType = PCFactorGetMatSolverType(petsclib::PetscLibType,pc::PC) 
gets the solver package that is used to perform the factorization

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `stype` - for example, `MATSOLVERSUPERLU`, `MATSOLVERSUPERLU_DIST`, `MATSOLVERMUMPS`

Level: intermediate

-seealso: [](ch_ksp), `PCLU`, `PCCHOLESKY`, `MatGetFactor()`, `MatSolverType`, `MATSOLVERSUPERLU`,
`MATSOLVERSUPERLU_DIST`, `MATSOLVERMUMPS`

# External Links
$(_doc_external("Ksp/PCFactorGetMatSolverType"))
"""
function PCFactorGetMatSolverType(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCFactorGetMatSolverType(petsclib::$UnionPetscLib, pc::PC )
	stype_ = Ref{MatSolverType}()

    @chk ccall(
               (:PCFactorGetMatSolverType, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{MatSolverType}),
               pc, stype_,
              )

	stype = unsafe_string(stype_[])

	return stype
end 

"""
	PCFactorSetFill(petsclib::PetscLibType,pc::PC, fill::PetscReal) 
Indicate the amount of fill you expect in the factored matrix,
fill = number nonzeros in factor/number nonzeros in original matrix.

Not Collective, each process can expect a different amount of fill

Input Parameters:
- `pc`   - the preconditioner context
- `fill` - amount of expected fill

Options Database Key:
- `-pc_factor_fill <fill>` - Sets fill amount

Level: intermediate

-seealso: [](ch_ksp), `PCLU`, `PCCHOLESKY`, `PCILU`, `PCICC`, `PCFactorSetReuseFill()`

# External Links
$(_doc_external("Ksp/PCFactorSetFill"))
"""
function PCFactorSetFill(petsclib::PetscLibType, pc::PC, fill::PetscReal) end

@for_petsc function PCFactorSetFill(petsclib::$UnionPetscLib, pc::PC, fill::$PetscReal )

    @chk ccall(
               (:PCFactorSetFill, $petsc_library),
               PetscErrorCode,
               (PC, $PetscReal),
               pc, fill,
              )


	return nothing
end 

"""
	PCFactorSetUseInPlace(petsclib::PetscLibType,pc::PC, flg::PetscBool) 
Tells the preconditioner to do an in

Logically Collective

Input Parameters:
- `pc`  - the preconditioner context
- `flg` - `PETSC_TRUE` to enable, `PETSC_FALSE` to disable

Options Database Key:
- `-pc_factor_in_place <true,false>` - Activate/deactivate in-place factorization

-seealso: [](ch_ksp), `PC`, `Mat`, `PCLU`, `PCCHOLESKY`, `PCILU`, `PCICC`, `PCFactorGetUseInPlace()`

# External Links
$(_doc_external("Ksp/PCFactorSetUseInPlace"))
"""
function PCFactorSetUseInPlace(petsclib::PetscLibType, pc::PC, flg::PetscBool) end

@for_petsc function PCFactorSetUseInPlace(petsclib::$UnionPetscLib, pc::PC, flg::PetscBool )

    @chk ccall(
               (:PCFactorSetUseInPlace, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, flg,
              )


	return nothing
end 

"""
	flg::PetscBool = PCFactorGetUseInPlace(petsclib::PetscLibType,pc::PC) 
Determines if an in

Logically Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `flg` - `PETSC_TRUE` to enable, `PETSC_FALSE` to disable

Level: intermediate

-seealso: [](ch_ksp), `PCLU`, `PCCHOLESKY`, `PCILU`, `PCICC`, `PCFactorSetUseInPlace()`

# External Links
$(_doc_external("Ksp/PCFactorGetUseInPlace"))
"""
function PCFactorGetUseInPlace(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCFactorGetUseInPlace(petsclib::$UnionPetscLib, pc::PC )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PCFactorGetUseInPlace, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PetscBool}),
               pc, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	PCFactorSetMatOrderingType(petsclib::PetscLibType,pc::PC, ordering::MatOrderingType) 
Sets the ordering routine (to reduce fill) to
be used in the `PCLU`, `PCCHOLESKY`, `PCILU`,  or `PCICC` preconditioners

Logically Collective

Input Parameters:
- `pc`       - the preconditioner context
- `ordering` - the matrix ordering name, for example, `MATORDERINGND` or `MATORDERINGRCM`

Options Database Key:
- `-pc_factor_mat_ordering_type <nd,rcm,...,external>` - Sets ordering routine

Level: intermediate

-seealso: [](ch_ksp), `PCLU`, `PCCHOLESKY`, `PCILU`, `PCICC`, `MatOrderingType`, `MATORDERINGEXTERNAL`, `MATORDERINGND`, `MATORDERINGRCM`

# External Links
$(_doc_external("Ksp/PCFactorSetMatOrderingType"))
"""
function PCFactorSetMatOrderingType(petsclib::PetscLibType, pc::PC, ordering::MatOrderingType) end

@for_petsc function PCFactorSetMatOrderingType(petsclib::$UnionPetscLib, pc::PC, ordering::MatOrderingType )

    @chk ccall(
               (:PCFactorSetMatOrderingType, $petsc_library),
               PetscErrorCode,
               (PC, MatOrderingType),
               pc, ordering,
              )


	return nothing
end 

"""
	PCFactorSetColumnPivot(petsclib::PetscLibType,pc::PC, dtcol::PetscReal) 
Determines when column pivoting is done during matrix factorization.
For PETSc dense matrices column pivoting is always done, for PETSc sparse matrices
it is never done. For the MATLAB and `MATSOLVERSUPERLU` factorization this is used.

Logically Collective

Input Parameters:
- `pc`    - the preconditioner context
- `dtcol` - 0.0 implies no pivoting, 1.0 complete pivoting (slower, requires more memory but more stable)

Options Database Key:
- `-pc_factor_pivoting <dtcol>` - perform the pivoting with the given tolerance

Level: intermediate

-seealso: [](ch_ksp), `PCLU`, `PCCHOLESKY`, `PCILU`, `PCICC`, `PCILUSetMatOrdering()`, `PCFactorSetPivotInBlocks()`

# External Links
$(_doc_external("Ksp/PCFactorSetColumnPivot"))
"""
function PCFactorSetColumnPivot(petsclib::PetscLibType, pc::PC, dtcol::PetscReal) end

@for_petsc function PCFactorSetColumnPivot(petsclib::$UnionPetscLib, pc::PC, dtcol::$PetscReal )

    @chk ccall(
               (:PCFactorSetColumnPivot, $petsc_library),
               PetscErrorCode,
               (PC, $PetscReal),
               pc, dtcol,
              )


	return nothing
end 

"""
	PCFactorSetPivotInBlocks(petsclib::PetscLibType,pc::PC, pivot::PetscBool) 
Determines if pivoting is done while factoring each block
with `MATBAIJ` or `MATSBAIJ` matrices

Logically Collective

Input Parameters:
- `pc`    - the preconditioner context
- `pivot` - `PETSC_TRUE` or `PETSC_FALSE`

Options Database Key:
- `-pc_factor_pivot_in_blocks <true,false>` - Pivot inside matrix dense blocks for `MATBAIJ` and `MATSBAIJ`

Level: intermediate

-seealso: [](ch_ksp), `PCLU`, `PCCHOLESKY`, `PCILU`, `PCICC`, `PCILUSetMatOrdering()`, `PCFactorSetColumnPivot()`

# External Links
$(_doc_external("Ksp/PCFactorSetPivotInBlocks"))
"""
function PCFactorSetPivotInBlocks(petsclib::PetscLibType, pc::PC, pivot::PetscBool) end

@for_petsc function PCFactorSetPivotInBlocks(petsclib::$UnionPetscLib, pc::PC, pivot::PetscBool )

    @chk ccall(
               (:PCFactorSetPivotInBlocks, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, pivot,
              )


	return nothing
end 

"""
	PCFactorSetReuseFill(petsclib::PetscLibType,pc::PC, flag::PetscBool) 
When matrices with different nonzero structure are factored,
this causes later ones to use the fill ratio computed in the initial factorization.

Logically Collective

Input Parameters:
- `pc`   - the preconditioner context
- `flag` - `PETSC_TRUE` to reuse else `PETSC_FALSE`

Options Database Key:
- `-pc_factor_reuse_fill` - Activates `PCFactorSetReuseFill()`

Level: intermediate

-seealso: [](ch_ksp), `PCLU`, `PCCHOLESKY`, `PCILU`, `PCICC`, `PCFactorSetReuseOrdering()`, `PCFactorSetFill()`

# External Links
$(_doc_external("Ksp/PCFactorSetReuseFill"))
"""
function PCFactorSetReuseFill(petsclib::PetscLibType, pc::PC, flag::PetscBool) end

@for_petsc function PCFactorSetReuseFill(petsclib::$UnionPetscLib, pc::PC, flag::PetscBool )

    @chk ccall(
               (:PCFactorSetReuseFill, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, flag,
              )


	return nothing
end 

"""
	PCFactorSetReuseOrdering(petsclib::PetscLibType,pc::PC, flag::PetscBool) 
When similar matrices are factored, this
causes the ordering computed in the first factor to be used for all
following factors.

Logically Collective

Input Parameters:
- `pc`   - the preconditioner context
- `flag` - `PETSC_TRUE` to reuse else `PETSC_FALSE`

Options Database Key:
- `-pc_factor_reuse_ordering` - Activate `PCFactorSetReuseOrdering()`

Level: intermediate

-seealso: [](ch_ksp), `PCLU`, `PCCHOLESKY`, `PCFactorSetReuseFill()`

# External Links
$(_doc_external("Ksp/PCFactorSetReuseOrdering"))
"""
function PCFactorSetReuseOrdering(petsclib::PetscLibType, pc::PC, flag::PetscBool) end

@for_petsc function PCFactorSetReuseOrdering(petsclib::$UnionPetscLib, pc::PC, flag::PetscBool )

    @chk ccall(
               (:PCFactorSetReuseOrdering, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, flag,
              )


	return nothing
end 

"""
	PCCompositeSetType(petsclib::PetscLibType,pc::PC, type::PCCompositeType) 
Sets the type of composite preconditioner.

Logically Collective

Input Parameters:
- `pc`   - the preconditioner context
- `type` - `PC_COMPOSITE_ADDITIVE` (default), `PC_COMPOSITE_MULTIPLICATIVE`, `PC_COMPOSITE_SPECIAL`

Options Database Key:
- `-pc_composite_type <type: one of multiplicative, additive, special>` - Sets composite preconditioner type

Level: advanced

-seealso: [](ch_ksp), `PCCOMPOSITE`, `PC_COMPOSITE_ADDITIVE`, `PC_COMPOSITE_MULTIPLICATIVE`, `PC_COMPOSITE_SPECIAL`, `PCCompositeType`,
`PCCompositeGetType()`

# External Links
$(_doc_external("Ksp/PCCompositeSetType"))
"""
function PCCompositeSetType(petsclib::PetscLibType, pc::PC, type::PCCompositeType) end

@for_petsc function PCCompositeSetType(petsclib::$UnionPetscLib, pc::PC, type::PCCompositeType )

    @chk ccall(
               (:PCCompositeSetType, $petsc_library),
               PetscErrorCode,
               (PC, PCCompositeType),
               pc, type,
              )


	return nothing
end 

"""
	type::PCCompositeType = PCCompositeGetType(petsclib::PetscLibType,pc::PC) 
Gets the type of composite preconditioner.

Logically Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `type` - `PC_COMPOSITE_ADDITIVE` (default), `PC_COMPOSITE_MULTIPLICATIVE`, `PC_COMPOSITE_SPECIAL`

Level: advanced

-seealso: [](ch_ksp), `PCCOMPOSITE`, `PC_COMPOSITE_ADDITIVE`, `PC_COMPOSITE_MULTIPLICATIVE`, `PC_COMPOSITE_SPECIAL`, `PCCompositeType`,
`PCCompositeSetType()`

# External Links
$(_doc_external("Ksp/PCCompositeGetType"))
"""
function PCCompositeGetType(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCCompositeGetType(petsclib::$UnionPetscLib, pc::PC )
	type_ = Ref{PCCompositeType}()

    @chk ccall(
               (:PCCompositeGetType, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PCCompositeType}),
               pc, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	PCCompositeSpecialSetAlpha(petsclib::PetscLibType,pc::PC, alpha::PetscScalar) 
Sets alpha for the special composite preconditioner, `PC_COMPOSITE_SPECIAL`,
for \alpha I + R + S

Logically Collective

Input Parameters:
- `pc`    - the preconditioner context
- `alpha` - scale on identity

Level: developer

-seealso: [](ch_ksp), `PCCOMPOSITE`, `PC_COMPOSITE_ADDITIVE`, `PC_COMPOSITE_MULTIPLICATIVE`, `PC_COMPOSITE_SPECIAL`, `PCCompositeType`,
`PCCompositeSetType()`, `PCCompositeGetType()`

# External Links
$(_doc_external("Ksp/PCCompositeSpecialSetAlpha"))
"""
function PCCompositeSpecialSetAlpha(petsclib::PetscLibType, pc::PC, alpha::PetscScalar) end

@for_petsc function PCCompositeSpecialSetAlpha(petsclib::$UnionPetscLib, pc::PC, alpha::$PetscScalar )

    @chk ccall(
               (:PCCompositeSpecialSetAlpha, $petsc_library),
               PetscErrorCode,
               (PC, $PetscScalar),
               pc, alpha,
              )


	return nothing
end 

"""
	PCCompositeSpecialSetAlphaMat(petsclib::PetscLibType,pc::PC, alpha_mat::PetscMat) 

# External Links
$(_doc_external("Ksp/PCCompositeSpecialSetAlphaMat"))
"""
function PCCompositeSpecialSetAlphaMat(petsclib::PetscLibType, pc::PC, alpha_mat::PetscMat) end

@for_petsc function PCCompositeSpecialSetAlphaMat(petsclib::$UnionPetscLib, pc::PC, alpha_mat::PetscMat )

    @chk ccall(
               (:PCCompositeSpecialSetAlphaMat, $petsc_library),
               PetscErrorCode,
               (PC, CMat),
               pc, alpha_mat,
              )


	return nothing
end 

"""
	PCCompositeAddPCType(petsclib::PetscLibType,pc::PC, type::PCType) 
Adds another `PC` of the given type to the composite `PC`.

Collective

Input Parameters:
- `pc`   - the preconditioner context
- `type` - the type of the new preconditioner

Level: intermediate

-seealso: [](ch_ksp), `PCCOMPOSITE`, `PCCompositeAddPC()`, `PCCompositeGetNumberPC()`

# External Links
$(_doc_external("Ksp/PCCompositeAddPCType"))
"""
function PCCompositeAddPCType(petsclib::PetscLibType, pc::PC, type::PCType) end

@for_petsc function PCCompositeAddPCType(petsclib::$UnionPetscLib, pc::PC, type::PCType )

    @chk ccall(
               (:PCCompositeAddPCType, $petsc_library),
               PetscErrorCode,
               (PC, PCType),
               pc, type,
              )


	return nothing
end 

"""
	PCCompositeAddPC(petsclib::PetscLibType,pc::PC, subpc::PC) 
Adds another `PC` to the composite `PC`.

Collective

Input Parameters:
- `pc`    - the preconditioner context
- `subpc` - the new preconditioner

Level: intermediate

-seealso: [](ch_ksp), `PCCOMPOSITE`, `PCCompositeAddPCType()`, `PCCompositeGetNumberPC()`

# External Links
$(_doc_external("Ksp/PCCompositeAddPC"))
"""
function PCCompositeAddPC(petsclib::PetscLibType, pc::PC, subpc::PC) end

@for_petsc function PCCompositeAddPC(petsclib::$UnionPetscLib, pc::PC, subpc::PC )

    @chk ccall(
               (:PCCompositeAddPC, $petsc_library),
               PetscErrorCode,
               (PC, PC),
               pc, subpc,
              )


	return nothing
end 

"""
	num::PetscInt = PCCompositeGetNumberPC(petsclib::PetscLibType,pc::PC) 
Gets the number of `PC` objects in the composite `PC`.

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `num` - the number of sub pcs

Level: developer

-seealso: [](ch_ksp), `PCCOMPOSITE`, `PCCompositeGetPC()`, `PCCompositeAddPC()`, `PCCompositeAddPCType()`

# External Links
$(_doc_external("Ksp/PCCompositeGetNumberPC"))
"""
function PCCompositeGetNumberPC(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCCompositeGetNumberPC(petsclib::$UnionPetscLib, pc::PC )
	num_ = Ref{$PetscInt}()

    @chk ccall(
               (:PCCompositeGetNumberPC, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{$PetscInt}),
               pc, num_,
              )

	num = num_[]

	return num
end 

"""
	PCCompositeGetPC(petsclib::PetscLibType,pc::PC, n::PetscInt, subpc::PC) 
Gets one of the `PC` objects in the composite `PC`.

Not Collective

Input Parameters:
- `pc` - the preconditioner context
- `n`  - the number of the pc requested

Output Parameter:
- `subpc` - the PC requested

Level: intermediate

-seealso: [](ch_ksp), `PCCOMPOSITE`, `PCCompositeAddPCType()`, `PCCompositeGetNumberPC()`, `PCSetOperators()`

# External Links
$(_doc_external("Ksp/PCCompositeGetPC"))
"""
function PCCompositeGetPC(petsclib::PetscLibType, pc::PC, n::PetscInt, subpc::PC) end

@for_petsc function PCCompositeGetPC(petsclib::$UnionPetscLib, pc::PC, n::$PetscInt, subpc::PC )

    @chk ccall(
               (:PCCompositeGetPC, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, Ptr{PC}),
               pc, n, subpc,
              )


	return nothing
end 

"""
	PCBDDCSetDiscreteGradient(petsclib::PetscLibType,pc::PC, G::PetscMat, order::PetscInt, field::PetscInt, global_::PetscBool, conforming::PetscBool) 
Sets the discrete gradient to be used by the `PCBDDC` preconditioner

Collective

Input Parameters:
- `pc`         - the preconditioning context
- `G`          - the discrete gradient matrix (in `MATAIJ` format)
- `order`      - the order of the Nedelec space (1 for the lowest order)
- `field`      - the field id of the Nedelec dofs (not used if the fields have not been specified)
- `global`     - the type of global ordering for the rows of `G`
- `conforming` - whether the mesh is conforming or not

Level: advanced

-seealso: [](ch_ksp), `PCBDDC`, `PCBDDCSetDofsSplitting()`, `PCBDDCSetDofsSplittingLocal()`, `MATAIJ`, `PCBDDCSetDivergenceMat()`

# External Links
$(_doc_external("Ksp/PCBDDCSetDiscreteGradient"))
"""
function PCBDDCSetDiscreteGradient(petsclib::PetscLibType, pc::PC, G::PetscMat, order::PetscInt, field::PetscInt, global_::PetscBool, conforming::PetscBool) end

@for_petsc function PCBDDCSetDiscreteGradient(petsclib::$UnionPetscLib, pc::PC, G::PetscMat, order::$PetscInt, field::$PetscInt, global_::PetscBool, conforming::PetscBool )

    @chk ccall(
               (:PCBDDCSetDiscreteGradient, $petsc_library),
               PetscErrorCode,
               (PC, CMat, $PetscInt, $PetscInt, PetscBool, PetscBool),
               pc, G, order, field, global_, conforming,
              )


	return nothing
end 

"""
	PCBDDCSetDivergenceMat(petsclib::PetscLibType,pc::PC, divudotp::PetscMat, trans::PetscBool, vl2l::IS) 
Sets the linear operator representing .. for the `PCBDDC` preconditioner

Collective

Input Parameters:
- `pc`       - the preconditioning context
- `divudotp` - the matrix (must be of type `MATIS`)
- `trans`    - if `PETSC_FALSE` (resp. `PETSC_TRUE`), then pressures are in the test (trial) space and velocities are in the trial (test) space.
- `vl2l`     - optional index set describing the local (wrt the local matrix in `divudotp`) to local (wrt the local matrix
in the matrix used to construct the preconditioner) map for the velocities

Level: advanced

-seealso: [](ch_ksp), `PCBDDC`, `PCBDDCSetDiscreteGradient()`

# External Links
$(_doc_external("Ksp/PCBDDCSetDivergenceMat"))
"""
function PCBDDCSetDivergenceMat(petsclib::PetscLibType, pc::PC, divudotp::PetscMat, trans::PetscBool, vl2l::IS) end

@for_petsc function PCBDDCSetDivergenceMat(petsclib::$UnionPetscLib, pc::PC, divudotp::PetscMat, trans::PetscBool, vl2l::IS )

    @chk ccall(
               (:PCBDDCSetDivergenceMat, $petsc_library),
               PetscErrorCode,
               (PC, CMat, PetscBool, IS),
               pc, divudotp, trans, vl2l,
              )


	return nothing
end 

"""
	PCBDDCSetChangeOfBasisMat(petsclib::PetscLibType,pc::PC, change::PetscMat, interior::PetscBool) 
Set user defined change of basis for dofs

Collective

Input Parameters:
- `pc`       - the preconditioning context
- `change`   - the change of basis matrix
- `interior` - whether or not the change of basis modifies interior dofs

Level: intermediate

-seealso: [](ch_ksp), `PCBDDC`

# External Links
$(_doc_external("Ksp/PCBDDCSetChangeOfBasisMat"))
"""
function PCBDDCSetChangeOfBasisMat(petsclib::PetscLibType, pc::PC, change::PetscMat, interior::PetscBool) end

@for_petsc function PCBDDCSetChangeOfBasisMat(petsclib::$UnionPetscLib, pc::PC, change::PetscMat, interior::PetscBool )

    @chk ccall(
               (:PCBDDCSetChangeOfBasisMat, $petsc_library),
               PetscErrorCode,
               (PC, CMat, PetscBool),
               pc, change, interior,
              )


	return nothing
end 

"""
	PCBDDCSetPrimalVerticesIS(petsclib::PetscLibType,pc::PC, PrimalVertices::IS) 
Set additional user defined primal vertices in `PCBDDC`

Collective

Input Parameters:
- `pc`             - the preconditioning context
- `PrimalVertices` - index set of primal vertices in global numbering (can be empty)

Level: intermediate

-seealso: [](ch_ksp), `PCBDDC`, `PCBDDCGetPrimalVerticesIS()`, `PCBDDCSetPrimalVerticesLocalIS()`, `PCBDDCGetPrimalVerticesLocalIS()`

# External Links
$(_doc_external("Ksp/PCBDDCSetPrimalVerticesIS"))
"""
function PCBDDCSetPrimalVerticesIS(petsclib::PetscLibType, pc::PC, PrimalVertices::IS) end

@for_petsc function PCBDDCSetPrimalVerticesIS(petsclib::$UnionPetscLib, pc::PC, PrimalVertices::IS )

    @chk ccall(
               (:PCBDDCSetPrimalVerticesIS, $petsc_library),
               PetscErrorCode,
               (PC, IS),
               pc, PrimalVertices,
              )


	return nothing
end 

"""
	PCBDDCGetPrimalVerticesIS(petsclib::PetscLibType,pc::PC, is::IS) 
Get user defined primal vertices set with `PCBDDCSetPrimalVerticesIS()`

Collective

Input Parameter:
- `pc` - the preconditioning context

Output Parameter:
- `is` - index set of primal vertices in global numbering (`NULL` if not set)

Level: intermediate

-seealso: [](ch_ksp), `PCBDDC`, `PCBDDCSetPrimalVerticesIS()`, `PCBDDCSetPrimalVerticesLocalIS()`, `PCBDDCGetPrimalVerticesLocalIS()`

# External Links
$(_doc_external("Ksp/PCBDDCGetPrimalVerticesIS"))
"""
function PCBDDCGetPrimalVerticesIS(petsclib::PetscLibType, pc::PC, is::IS) end

@for_petsc function PCBDDCGetPrimalVerticesIS(petsclib::$UnionPetscLib, pc::PC, is::IS )

    @chk ccall(
               (:PCBDDCGetPrimalVerticesIS, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{IS}),
               pc, is,
              )


	return nothing
end 

"""
	PCBDDCSetPrimalVerticesLocalIS(petsclib::PetscLibType,pc::PC, PrimalVertices::IS) 
Set additional user defined primal vertices in `PCBDDC`

Collective

Input Parameters:
- `pc`             - the preconditioning context
- `PrimalVertices` - index set of primal vertices in local numbering (can be empty)

Level: intermediate

-seealso: [](ch_ksp), `PCBDDC`, `PCBDDCSetPrimalVerticesIS()`, `PCBDDCGetPrimalVerticesIS()`, `PCBDDCGetPrimalVerticesLocalIS()`

# External Links
$(_doc_external("Ksp/PCBDDCSetPrimalVerticesLocalIS"))
"""
function PCBDDCSetPrimalVerticesLocalIS(petsclib::PetscLibType, pc::PC, PrimalVertices::IS) end

@for_petsc function PCBDDCSetPrimalVerticesLocalIS(petsclib::$UnionPetscLib, pc::PC, PrimalVertices::IS )

    @chk ccall(
               (:PCBDDCSetPrimalVerticesLocalIS, $petsc_library),
               PetscErrorCode,
               (PC, IS),
               pc, PrimalVertices,
              )


	return nothing
end 

"""
	PCBDDCGetPrimalVerticesLocalIS(petsclib::PetscLibType,pc::PC, is::IS) 
Get user defined primal vertices set with `PCBDDCSetPrimalVerticesLocalIS()`

Collective

Input Parameter:
- `pc` - the preconditioning context

Output Parameter:
- `is` - index set of primal vertices in local numbering (`NULL` if not set)

Level: intermediate

-seealso: [](ch_ksp), `PCBDDC`, `PCBDDCSetPrimalVerticesIS()`, `PCBDDCGetPrimalVerticesIS()`, `PCBDDCSetPrimalVerticesLocalIS()`

# External Links
$(_doc_external("Ksp/PCBDDCGetPrimalVerticesLocalIS"))
"""
function PCBDDCGetPrimalVerticesLocalIS(petsclib::PetscLibType, pc::PC, is::IS) end

@for_petsc function PCBDDCGetPrimalVerticesLocalIS(petsclib::$UnionPetscLib, pc::PC, is::IS )

    @chk ccall(
               (:PCBDDCGetPrimalVerticesLocalIS, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{IS}),
               pc, is,
              )


	return nothing
end 

"""
	PCBDDCSetCoarseningRatio(petsclib::PetscLibType,pc::PC, k::PetscInt) 
Set coarsening ratio used in the multi

Logically Collective

Input Parameters:
- `pc` - the preconditioning context
- `k`  - coarsening ratio (H/h at the coarser level)

Options Database Key:
- `-pc_bddc_coarsening_ratio <int>` - Set the coarsening ratio used in multi-level coarsening

Level: intermediate

-seealso: [](ch_ksp), `PCBDDC`, `PCBDDCSetLevels()`

# External Links
$(_doc_external("Ksp/PCBDDCSetCoarseningRatio"))
"""
function PCBDDCSetCoarseningRatio(petsclib::PetscLibType, pc::PC, k::PetscInt) end

@for_petsc function PCBDDCSetCoarseningRatio(petsclib::$UnionPetscLib, pc::PC, k::$PetscInt )

    @chk ccall(
               (:PCBDDCSetCoarseningRatio, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, k,
              )


	return nothing
end 

"""
	PCBDDCSetLevels(petsclib::PetscLibType,pc::PC, levels::PetscInt) 
Sets the maximum number of additional levels allowed for multilevel `PCBDDC`

Logically Collective

Input Parameters:
- `pc`     - the preconditioning context
- `levels` - the maximum number of levels

Options Database Key:
- `-pc_bddc_levels <int>` - Set maximum number of levels for multilevel

Level: intermediate

-seealso: [](ch_ksp), `PCBDDC`, `PCBDDCSetCoarseningRatio()`

# External Links
$(_doc_external("Ksp/PCBDDCSetLevels"))
"""
function PCBDDCSetLevels(petsclib::PetscLibType, pc::PC, levels::PetscInt) end

@for_petsc function PCBDDCSetLevels(petsclib::$UnionPetscLib, pc::PC, levels::$PetscInt )

    @chk ccall(
               (:PCBDDCSetLevels, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, levels,
              )


	return nothing
end 

"""
	PCBDDCSetDirichletBoundaries(petsclib::PetscLibType,pc::PC, DirichletBoundaries::IS) 
Set the `IS` defining Dirichlet boundaries for the global problem.

Collective

Input Parameters:
- `pc`                  - the preconditioning context
- `DirichletBoundaries` - parallel `IS` defining the Dirichlet boundaries

Level: intermediate

-seealso: [](ch_ksp), `PCBDDC`, `PCBDDCSetDirichletBoundariesLocal()`, `MatZeroRows()`, `MatZeroRowsColumns()`

# External Links
$(_doc_external("Ksp/PCBDDCSetDirichletBoundaries"))
"""
function PCBDDCSetDirichletBoundaries(petsclib::PetscLibType, pc::PC, DirichletBoundaries::IS) end

@for_petsc function PCBDDCSetDirichletBoundaries(petsclib::$UnionPetscLib, pc::PC, DirichletBoundaries::IS )

    @chk ccall(
               (:PCBDDCSetDirichletBoundaries, $petsc_library),
               PetscErrorCode,
               (PC, IS),
               pc, DirichletBoundaries,
              )


	return nothing
end 

"""
	PCBDDCSetDirichletBoundariesLocal(petsclib::PetscLibType,pc::PC, DirichletBoundaries::IS) 
Set the `IS` defining Dirichlet boundaries for the global problem in local ordering.

Collective

Input Parameters:
- `pc`                  - the preconditioning context
- `DirichletBoundaries` - parallel `IS` defining the Dirichlet boundaries (in local ordering)

Level: intermediate

-seealso: [](ch_ksp), `PCBDDC`, `PCBDDCSetDirichletBoundaries()`, `MatZeroRows()`, `MatZeroRowsColumns()`

# External Links
$(_doc_external("Ksp/PCBDDCSetDirichletBoundariesLocal"))
"""
function PCBDDCSetDirichletBoundariesLocal(petsclib::PetscLibType, pc::PC, DirichletBoundaries::IS) end

@for_petsc function PCBDDCSetDirichletBoundariesLocal(petsclib::$UnionPetscLib, pc::PC, DirichletBoundaries::IS )

    @chk ccall(
               (:PCBDDCSetDirichletBoundariesLocal, $petsc_library),
               PetscErrorCode,
               (PC, IS),
               pc, DirichletBoundaries,
              )


	return nothing
end 

"""
	PCBDDCSetNeumannBoundaries(petsclib::PetscLibType,pc::PC, NeumannBoundaries::IS) 
Set the `IS` defining Neumann boundaries for the global problem.

Collective

Input Parameters:
- `pc`                - the preconditioning context
- `NeumannBoundaries` - parallel `IS` defining the Neumann boundaries

Level: intermediate

-seealso: [](ch_ksp), `PCBDDC`, `PCBDDCSetNeumannBoundariesLocal()`

# External Links
$(_doc_external("Ksp/PCBDDCSetNeumannBoundaries"))
"""
function PCBDDCSetNeumannBoundaries(petsclib::PetscLibType, pc::PC, NeumannBoundaries::IS) end

@for_petsc function PCBDDCSetNeumannBoundaries(petsclib::$UnionPetscLib, pc::PC, NeumannBoundaries::IS )

    @chk ccall(
               (:PCBDDCSetNeumannBoundaries, $petsc_library),
               PetscErrorCode,
               (PC, IS),
               pc, NeumannBoundaries,
              )


	return nothing
end 

"""
	PCBDDCSetNeumannBoundariesLocal(petsclib::PetscLibType,pc::PC, NeumannBoundaries::IS) 
Set the `IS` defining Neumann boundaries for the global problem in local ordering.

Collective

Input Parameters:
- `pc`                - the preconditioning context
- `NeumannBoundaries` - parallel `IS` defining the subdomain part of Neumann boundaries (in local ordering)

Level: intermediate

-seealso: [](ch_ksp), `PCBDDC`, `PCBDDCSetNeumannBoundaries()`, `PCBDDCGetDirichletBoundaries()`

# External Links
$(_doc_external("Ksp/PCBDDCSetNeumannBoundariesLocal"))
"""
function PCBDDCSetNeumannBoundariesLocal(petsclib::PetscLibType, pc::PC, NeumannBoundaries::IS) end

@for_petsc function PCBDDCSetNeumannBoundariesLocal(petsclib::$UnionPetscLib, pc::PC, NeumannBoundaries::IS )

    @chk ccall(
               (:PCBDDCSetNeumannBoundariesLocal, $petsc_library),
               PetscErrorCode,
               (PC, IS),
               pc, NeumannBoundaries,
              )


	return nothing
end 

"""
	PCBDDCGetDirichletBoundaries(petsclib::PetscLibType,pc::PC, DirichletBoundaries::IS) 
Get parallel `IS` for Dirichlet boundaries

Collective

Input Parameter:
- `pc` - the preconditioning context

Output Parameter:
- `DirichletBoundaries` - index set defining the Dirichlet boundaries

Level: intermediate

-seealso: [](ch_ksp), `PCBDDC`, `PCBDDCSetDirichletBoundaries()`

# External Links
$(_doc_external("Ksp/PCBDDCGetDirichletBoundaries"))
"""
function PCBDDCGetDirichletBoundaries(petsclib::PetscLibType, pc::PC, DirichletBoundaries::IS) end

@for_petsc function PCBDDCGetDirichletBoundaries(petsclib::$UnionPetscLib, pc::PC, DirichletBoundaries::IS )

    @chk ccall(
               (:PCBDDCGetDirichletBoundaries, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{IS}),
               pc, DirichletBoundaries,
              )


	return nothing
end 

"""
	PCBDDCGetDirichletBoundariesLocal(petsclib::PetscLibType,pc::PC, DirichletBoundaries::IS) 
Get parallel `IS` for Dirichlet boundaries (in local ordering)

Collective

Input Parameter:
- `pc` - the preconditioning context

Output Parameter:
- `DirichletBoundaries` - index set defining the subdomain part of Dirichlet boundaries

Level: intermediate

-seealso: [](ch_ksp), `PCBDDC`, `PCBDDCGetDirichletBoundaries()`, `PCBDDCSetDirichletBoundaries()`

# External Links
$(_doc_external("Ksp/PCBDDCGetDirichletBoundariesLocal"))
"""
function PCBDDCGetDirichletBoundariesLocal(petsclib::PetscLibType, pc::PC, DirichletBoundaries::IS) end

@for_petsc function PCBDDCGetDirichletBoundariesLocal(petsclib::$UnionPetscLib, pc::PC, DirichletBoundaries::IS )

    @chk ccall(
               (:PCBDDCGetDirichletBoundariesLocal, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{IS}),
               pc, DirichletBoundaries,
              )


	return nothing
end 

"""
	PCBDDCGetNeumannBoundaries(petsclib::PetscLibType,pc::PC, NeumannBoundaries::IS) 
Get parallel `IS` for Neumann boundaries

Not Collective

Input Parameter:
- `pc` - the preconditioning context

Output Parameter:
- `NeumannBoundaries` - index set defining the Neumann boundaries

Level: intermediate

-seealso: [](ch_ksp), `PCBDDC`, `PCBDDCSetNeumannBoundaries()`, `PCBDDCGetDirichletBoundaries()`, `PCBDDCSetDirichletBoundaries()`

# External Links
$(_doc_external("Ksp/PCBDDCGetNeumannBoundaries"))
"""
function PCBDDCGetNeumannBoundaries(petsclib::PetscLibType, pc::PC, NeumannBoundaries::IS) end

@for_petsc function PCBDDCGetNeumannBoundaries(petsclib::$UnionPetscLib, pc::PC, NeumannBoundaries::IS )

    @chk ccall(
               (:PCBDDCGetNeumannBoundaries, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{IS}),
               pc, NeumannBoundaries,
              )


	return nothing
end 

"""
	PCBDDCGetNeumannBoundariesLocal(petsclib::PetscLibType,pc::PC, NeumannBoundaries::IS) 
Get parallel `IS` for Neumann boundaries (in local ordering)

Not Collective

Input Parameter:
- `pc` - the preconditioning context

Output Parameter:
- `NeumannBoundaries` - index set defining the subdomain part of Neumann boundaries

Level: intermediate

-seealso: [](ch_ksp), `PCBDDC`, `PCBDDCSetNeumannBoundaries()`, `PCBDDCSetNeumannBoundariesLocal()`, `PCBDDCGetNeumannBoundaries()`

# External Links
$(_doc_external("Ksp/PCBDDCGetNeumannBoundariesLocal"))
"""
function PCBDDCGetNeumannBoundariesLocal(petsclib::PetscLibType, pc::PC, NeumannBoundaries::IS) end

@for_petsc function PCBDDCGetNeumannBoundariesLocal(petsclib::$UnionPetscLib, pc::PC, NeumannBoundaries::IS )

    @chk ccall(
               (:PCBDDCGetNeumannBoundariesLocal, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{IS}),
               pc, NeumannBoundaries,
              )


	return nothing
end 

"""
	PCBDDCSetLocalAdjacencyGraph(petsclib::PetscLibType,pc::PC, nvtxs::PetscInt, xadj::Vector{PetscInt}, adjncy::Vector{PetscInt}, copymode::PetscCopyMode) 
Set adjacency structure (CSR graph) of the local degrees of freedom.

Not collective

Input Parameters:
- `pc`       - the preconditioning context.
- `nvtxs`    - number of local vertices of the graph (i.e., the number of local dofs).
- `xadj`     - CSR format row pointers for the connectivity of the dofs
- `adjncy`   - CSR format column pointers for the connectivity of the dofs
- `copymode` - supported modes are `PETSC_COPY_VALUES`, `PETSC_USE_POINTER` or `PETSC_OWN_POINTER`.

Level: intermediate

-seealso: [](ch_ksp), `PCBDDC`, `PetscCopyMode`

# External Links
$(_doc_external("Ksp/PCBDDCSetLocalAdjacencyGraph"))
"""
function PCBDDCSetLocalAdjacencyGraph(petsclib::PetscLibType, pc::PC, nvtxs::PetscInt, xadj::Vector{PetscInt}, adjncy::Vector{PetscInt}, copymode::PetscCopyMode) end

@for_petsc function PCBDDCSetLocalAdjacencyGraph(petsclib::$UnionPetscLib, pc::PC, nvtxs::$PetscInt, xadj::Vector{$PetscInt}, adjncy::Vector{$PetscInt}, copymode::PetscCopyMode )

    @chk ccall(
               (:PCBDDCSetLocalAdjacencyGraph, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, PetscCopyMode),
               pc, nvtxs, xadj, adjncy, copymode,
              )


	return nothing
end 

"""
	PCBDDCSetDofsSplittingLocal(petsclib::PetscLibType,pc::PC, n_is::PetscInt, ISForDofs::Vector{IS}) 
Set the `IS` defining fields of the local subdomain matrix

Collective

Input Parameters:
- `pc`        - the preconditioning context
- `n_is`      - number of index sets defining the fields, must be the same on all MPI processes
- `ISForDofs` - array of `IS` describing the fields in local ordering

Level: intermediate

-seealso: [](ch_ksp), `PCBDDC`, `PCBDDCSetDofsSplitting()`

# External Links
$(_doc_external("Ksp/PCBDDCSetDofsSplittingLocal"))
"""
function PCBDDCSetDofsSplittingLocal(petsclib::PetscLibType, pc::PC, n_is::PetscInt, ISForDofs::Vector{IS}) end

@for_petsc function PCBDDCSetDofsSplittingLocal(petsclib::$UnionPetscLib, pc::PC, n_is::$PetscInt, ISForDofs::Vector{IS} )

    @chk ccall(
               (:PCBDDCSetDofsSplittingLocal, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, Ptr{IS}),
               pc, n_is, ISForDofs,
              )


	return nothing
end 

"""
	PCBDDCSetDofsSplitting(petsclib::PetscLibType,pc::PC, n_is::PetscInt, ISForDofs::Vector{IS}) 
Set the `IS` defining fields of the global matrix

Collective

Input Parameters:
- `pc`        - the preconditioning context
- `n_is`      - number of index sets defining the fields
- `ISForDofs` - array of `IS` describing the fields in global ordering

Level: intermediate

-seealso: [](ch_ksp), `PCBDDC`, `PCBDDCSetDofsSplittingLocal()`

# External Links
$(_doc_external("Ksp/PCBDDCSetDofsSplitting"))
"""
function PCBDDCSetDofsSplitting(petsclib::PetscLibType, pc::PC, n_is::PetscInt, ISForDofs::Vector{IS}) end

@for_petsc function PCBDDCSetDofsSplitting(petsclib::$UnionPetscLib, pc::PC, n_is::$PetscInt, ISForDofs::Vector{IS} )

    @chk ccall(
               (:PCBDDCSetDofsSplitting, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, Ptr{IS}),
               pc, n_is, ISForDofs,
              )


	return nothing
end 

"""
	PCBDDCMatFETIDPGetRHS(petsclib::PetscLibType,fetidp_mat::PetscMat, standard_rhs::PetscVec, fetidp_flux_rhs::PetscVec) 
Compute the right

Collective

Input Parameters:
- `fetidp_mat`   - the FETI-DP matrix object obtained by a call to `PCBDDCCreateFETIDPOperators()`
- `standard_rhs` - the right-hand side of the original linear system

Output Parameter:
- `fetidp_flux_rhs` - the right-hand side for the FETI-DP linear system

Level: developer

-seealso: [](ch_ksp), `PCBDDC`, `PCBDDCCreateFETIDPOperators()`, `PCBDDCMatFETIDPGetSolution()`

# External Links
$(_doc_external("Ksp/PCBDDCMatFETIDPGetRHS"))
"""
function PCBDDCMatFETIDPGetRHS(petsclib::PetscLibType, fetidp_mat::PetscMat, standard_rhs::PetscVec, fetidp_flux_rhs::PetscVec) end

@for_petsc function PCBDDCMatFETIDPGetRHS(petsclib::$UnionPetscLib, fetidp_mat::PetscMat, standard_rhs::PetscVec, fetidp_flux_rhs::PetscVec )

    @chk ccall(
               (:PCBDDCMatFETIDPGetRHS, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               fetidp_mat, standard_rhs, fetidp_flux_rhs,
              )


	return nothing
end 

"""
	PCBDDCMatFETIDPGetSolution(petsclib::PetscLibType,fetidp_mat::PetscMat, fetidp_flux_sol::PetscVec, standard_sol::PetscVec) 
Compute the physical solution using the solution of the FETI

Collective

Input Parameters:
- `fetidp_mat`      - the FETI-DP matrix obtained by a call to `PCBDDCCreateFETIDPOperators()`
- `fetidp_flux_sol` - the solution of the FETI-DP linear system`

Output Parameter:
- `standard_sol` - the solution defined on the physical domain

Level: developer

-seealso: [](ch_ksp), `PCBDDC`, `PCBDDCCreateFETIDPOperators()`, `PCBDDCMatFETIDPGetRHS()`

# External Links
$(_doc_external("Ksp/PCBDDCMatFETIDPGetSolution"))
"""
function PCBDDCMatFETIDPGetSolution(petsclib::PetscLibType, fetidp_mat::PetscMat, fetidp_flux_sol::PetscVec, standard_sol::PetscVec) end

@for_petsc function PCBDDCMatFETIDPGetSolution(petsclib::$UnionPetscLib, fetidp_mat::PetscMat, fetidp_flux_sol::PetscVec, standard_sol::PetscVec )

    @chk ccall(
               (:PCBDDCMatFETIDPGetSolution, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec),
               fetidp_mat, fetidp_flux_sol, standard_sol,
              )


	return nothing
end 

"""
	fetidp_mat::PetscMat,fetidp_pc::PC = PCBDDCCreateFETIDPOperators(petsclib::PetscLibType,pc::PC, fully_redundant::PetscBool, prefix::String) 
Create FETI

Collective

Input Parameters:
- `pc`              - the `PCBDDC` preconditioning context (setup should have been called before)
- `fully_redundant` - true for a fully redundant set of Lagrange multipliers
- `prefix`          - optional options database prefix for the objects to be created (can be `NULL`)

Output Parameters:
- `fetidp_mat` - shell FETI-DP matrix object
- `fetidp_pc`  - shell Dirichlet preconditioner for FETI-DP matrix

Level: developer

-seealso: [](ch_ksp), `KSPFETIDP`, `PCBDDC`, `PCBDDCMatFETIDPGetRHS()`, `PCBDDCMatFETIDPGetSolution()`

# External Links
$(_doc_external("Ksp/PCBDDCCreateFETIDPOperators"))
"""
function PCBDDCCreateFETIDPOperators(petsclib::PetscLibType, pc::PC, fully_redundant::PetscBool, prefix::String) end

@for_petsc function PCBDDCCreateFETIDPOperators(petsclib::$UnionPetscLib, pc::PC, fully_redundant::PetscBool, prefix::String )
	fetidp_mat_ = Ref{CMat}()
	fetidp_pc_ = Ref{PC}()

    @chk ccall(
               (:PCBDDCCreateFETIDPOperators, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool, Ptr{Cchar}, Ptr{CMat}, Ptr{PC}),
               pc, fully_redundant, prefix, fetidp_mat_, fetidp_pc_,
              )

	fetidp_mat = PetscMat(fetidp_mat_[], petsclib)
	fetidp_pc = fetidp_pc_[]

	return fetidp_mat,fetidp_pc
end 

"""
	PCBDDCInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `PCBDDC` package. It is called
from `PCInitializePackage()`.

Level: developer

-seealso: [](ch_ksp), `PetscInitialize()`, `PCBDDCFinalizePackage()`

# External Links
$(_doc_external("Ksp/PCBDDCInitializePackage"))
"""
function PCBDDCInitializePackage(petsclib::PetscLibType) end

@for_petsc function PCBDDCInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PCBDDCInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PCBDDCFinalizePackage(petsclib::PetscLibType) 
This function frees everything from the `PCBDDC` package. It is
called from `PetscFinalize()` automatically.

Level: developer

-seealso: [](ch_ksp), `PetscFinalize()`, `PCBDDCInitializePackage()`

# External Links
$(_doc_external("Ksp/PCBDDCFinalizePackage"))
"""
function PCBDDCFinalizePackage(petsclib::PetscLibType) end

@for_petsc function PCBDDCFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PCBDDCFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PCHPDDMSetAuxiliaryMat(petsclib::PetscLibType,pc::PC, is::IS, A::PetscMat, setup::external, ctx::Cvoid) 

# External Links
$(_doc_external("Ksp/PCHPDDMSetAuxiliaryMat"))
"""
function PCHPDDMSetAuxiliaryMat(petsclib::PetscLibType, pc::PC, is::IS, A::PetscMat, setup::external, ctx::Cvoid) end

@for_petsc function PCHPDDMSetAuxiliaryMat(petsclib::$UnionPetscLib, pc::PC, is::IS, A::PetscMat, setup::external, ctx::Cvoid )

    @chk ccall(
               (:PCHPDDMSetAuxiliaryMat, $petsc_library),
               PetscErrorCode,
               (PC, IS, CMat, external, Ptr{Cvoid}),
               pc, is, A, setup, ctx,
              )


	return nothing
end 

"""
	PCHPDDMHasNeumannMat(petsclib::PetscLibType,pc::PC, has::PetscBool) 

# External Links
$(_doc_external("Ksp/PCHPDDMHasNeumannMat"))
"""
function PCHPDDMHasNeumannMat(petsclib::PetscLibType, pc::PC, has::PetscBool) end

@for_petsc function PCHPDDMHasNeumannMat(petsclib::$UnionPetscLib, pc::PC, has::PetscBool )

    @chk ccall(
               (:PCHPDDMHasNeumannMat, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, has,
              )


	return nothing
end 

"""
	PCHPDDMSetRHSMat(petsclib::PetscLibType,pc::PC, B::PetscMat) 

# External Links
$(_doc_external("Ksp/PCHPDDMSetRHSMat"))
"""
function PCHPDDMSetRHSMat(petsclib::PetscLibType, pc::PC, B::PetscMat) end

@for_petsc function PCHPDDMSetRHSMat(petsclib::$UnionPetscLib, pc::PC, B::PetscMat )

    @chk ccall(
               (:PCHPDDMSetRHSMat, $petsc_library),
               PetscErrorCode,
               (PC, CMat),
               pc, B,
              )


	return nothing
end 

"""
	gc::PetscReal,oc::PetscReal = PCHPDDMGetComplexities(petsclib::PetscLibType,pc::PC) 

# External Links
$(_doc_external("Ksp/PCHPDDMGetComplexities"))
"""
function PCHPDDMGetComplexities(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCHPDDMGetComplexities(petsclib::$UnionPetscLib, pc::PC )
	gc_ = Ref{$PetscReal}()
	oc_ = Ref{$PetscReal}()

    @chk ccall(
               (:PCHPDDMGetComplexities, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{$PetscReal}, Ptr{$PetscReal}),
               pc, gc_, oc_,
              )

	gc = gc_[]
	oc = oc_[]

	return gc,oc
end 

"""
	PCHPDDMSetCoarseCorrectionType(petsclib::PetscLibType,pc::PC, type::PCHPDDMCoarseCorrectionType) 

# External Links
$(_doc_external("Ksp/PCHPDDMSetCoarseCorrectionType"))
"""
function PCHPDDMSetCoarseCorrectionType(petsclib::PetscLibType, pc::PC, type::PCHPDDMCoarseCorrectionType) end

@for_petsc function PCHPDDMSetCoarseCorrectionType(petsclib::$UnionPetscLib, pc::PC, type::PCHPDDMCoarseCorrectionType )

    @chk ccall(
               (:PCHPDDMSetCoarseCorrectionType, $petsc_library),
               PetscErrorCode,
               (PC, PCHPDDMCoarseCorrectionType),
               pc, type,
              )


	return nothing
end 

"""
	type::PCHPDDMCoarseCorrectionType = PCHPDDMGetCoarseCorrectionType(petsclib::PetscLibType,pc::PC) 

# External Links
$(_doc_external("Ksp/PCHPDDMGetCoarseCorrectionType"))
"""
function PCHPDDMGetCoarseCorrectionType(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCHPDDMGetCoarseCorrectionType(petsclib::$UnionPetscLib, pc::PC )
	type_ = Ref{PCHPDDMCoarseCorrectionType}()

    @chk ccall(
               (:PCHPDDMGetCoarseCorrectionType, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PCHPDDMCoarseCorrectionType}),
               pc, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	PCHPDDMSetSTShareSubKSP(petsclib::PetscLibType,pc::PC, share::PetscBool) 

# External Links
$(_doc_external("Ksp/PCHPDDMSetSTShareSubKSP"))
"""
function PCHPDDMSetSTShareSubKSP(petsclib::PetscLibType, pc::PC, share::PetscBool) end

@for_petsc function PCHPDDMSetSTShareSubKSP(petsclib::$UnionPetscLib, pc::PC, share::PetscBool )

    @chk ccall(
               (:PCHPDDMSetSTShareSubKSP, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, share,
              )


	return nothing
end 

"""
	share::PetscBool = PCHPDDMGetSTShareSubKSP(petsclib::PetscLibType,pc::PC) 

# External Links
$(_doc_external("Ksp/PCHPDDMGetSTShareSubKSP"))
"""
function PCHPDDMGetSTShareSubKSP(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCHPDDMGetSTShareSubKSP(petsclib::$UnionPetscLib, pc::PC )
	share_ = Ref{PetscBool}()

    @chk ccall(
               (:PCHPDDMGetSTShareSubKSP, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PetscBool}),
               pc, share_,
              )

	share = share_[]

	return share
end 

"""
	PCHPDDMSetDeflationMat(petsclib::PetscLibType,pc::PC, is::IS, U::PetscMat) 

# External Links
$(_doc_external("Ksp/PCHPDDMSetDeflationMat"))
"""
function PCHPDDMSetDeflationMat(petsclib::PetscLibType, pc::PC, is::IS, U::PetscMat) end

@for_petsc function PCHPDDMSetDeflationMat(petsclib::$UnionPetscLib, pc::PC, is::IS, U::PetscMat )

    @chk ccall(
               (:PCHPDDMSetDeflationMat, $petsc_library),
               PetscErrorCode,
               (PC, IS, CMat),
               pc, is, U,
              )


	return nothing
end 

"""
	PCHPDDMInitializePackage(petsclib::PetscLibType) 

# External Links
$(_doc_external("Ksp/PCHPDDMInitializePackage"))
"""
function PCHPDDMInitializePackage(petsclib::PetscLibType) end

@for_petsc function PCHPDDMInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PCHPDDMInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PCHPDDMFinalizePackage(petsclib::PetscLibType) 

# External Links
$(_doc_external("Ksp/PCHPDDMFinalizePackage"))
"""
function PCHPDDMFinalizePackage(petsclib::PetscLibType) end

@for_petsc function PCHPDDMFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PCHPDDMFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PCPatchSetSaveOperators(petsclib::PetscLibType,pc::PC, flg::PetscBool) 

# External Links
$(_doc_external("Ksp/PCPatchSetSaveOperators"))
"""
function PCPatchSetSaveOperators(petsclib::PetscLibType, pc::PC, flg::PetscBool) end

@for_petsc function PCPatchSetSaveOperators(petsclib::$UnionPetscLib, pc::PC, flg::PetscBool )

    @chk ccall(
               (:PCPatchSetSaveOperators, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, flg,
              )


	return nothing
end 

"""
	flg::PetscBool = PCPatchGetSaveOperators(petsclib::PetscLibType,pc::PC) 

# External Links
$(_doc_external("Ksp/PCPatchGetSaveOperators"))
"""
function PCPatchGetSaveOperators(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCPatchGetSaveOperators(petsclib::$UnionPetscLib, pc::PC )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PCPatchGetSaveOperators, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PetscBool}),
               pc, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	PCPatchSetPrecomputeElementTensors(petsclib::PetscLibType,pc::PC, flg::PetscBool) 

# External Links
$(_doc_external("Ksp/PCPatchSetPrecomputeElementTensors"))
"""
function PCPatchSetPrecomputeElementTensors(petsclib::PetscLibType, pc::PC, flg::PetscBool) end

@for_petsc function PCPatchSetPrecomputeElementTensors(petsclib::$UnionPetscLib, pc::PC, flg::PetscBool )

    @chk ccall(
               (:PCPatchSetPrecomputeElementTensors, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, flg,
              )


	return nothing
end 

"""
	flg::PetscBool = PCPatchGetPrecomputeElementTensors(petsclib::PetscLibType,pc::PC) 

# External Links
$(_doc_external("Ksp/PCPatchGetPrecomputeElementTensors"))
"""
function PCPatchGetPrecomputeElementTensors(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCPatchGetPrecomputeElementTensors(petsclib::$UnionPetscLib, pc::PC )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PCPatchGetPrecomputeElementTensors, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PetscBool}),
               pc, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	PCPatchSetPartitionOfUnity(petsclib::PetscLibType,pc::PC, flg::PetscBool) 

# External Links
$(_doc_external("Ksp/PCPatchSetPartitionOfUnity"))
"""
function PCPatchSetPartitionOfUnity(petsclib::PetscLibType, pc::PC, flg::PetscBool) end

@for_petsc function PCPatchSetPartitionOfUnity(petsclib::$UnionPetscLib, pc::PC, flg::PetscBool )

    @chk ccall(
               (:PCPatchSetPartitionOfUnity, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, flg,
              )


	return nothing
end 

"""
	flg::PetscBool = PCPatchGetPartitionOfUnity(petsclib::PetscLibType,pc::PC) 

# External Links
$(_doc_external("Ksp/PCPatchGetPartitionOfUnity"))
"""
function PCPatchGetPartitionOfUnity(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCPatchGetPartitionOfUnity(petsclib::$UnionPetscLib, pc::PC )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PCPatchGetPartitionOfUnity, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PetscBool}),
               pc, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	npatch::PetscInt = PCPatchGetSubKSP(petsclib::PetscLibType,pc::PC, ksp::PetscKSP) 

# External Links
$(_doc_external("Ksp/PCPatchGetSubKSP"))
"""
function PCPatchGetSubKSP(petsclib::PetscLibType, pc::PC, ksp::PetscKSP) end

@for_petsc function PCPatchGetSubKSP(petsclib::$UnionPetscLib, pc::PC, ksp::PetscKSP )
	npatch_ = Ref{$PetscInt}()

    @chk ccall(
               (:PCPatchGetSubKSP, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{$PetscInt}, CKSP),
               pc, npatch_, ksp,
              )

	npatch = npatch_[]

	return npatch
end 

"""
	PCPatchSetSubMatType(petsclib::PetscLibType,pc::PC, sub_mat_type::MatType) 

# External Links
$(_doc_external("Ksp/PCPatchSetSubMatType"))
"""
function PCPatchSetSubMatType(petsclib::PetscLibType, pc::PC, sub_mat_type::MatType) end

@for_petsc function PCPatchSetSubMatType(petsclib::$UnionPetscLib, pc::PC, sub_mat_type::MatType )

    @chk ccall(
               (:PCPatchSetSubMatType, $petsc_library),
               PetscErrorCode,
               (PC, MatType),
               pc, sub_mat_type,
              )


	return nothing
end 

"""
	sub_mat_type::MatType = PCPatchGetSubMatType(petsclib::PetscLibType,pc::PC) 

# External Links
$(_doc_external("Ksp/PCPatchGetSubMatType"))
"""
function PCPatchGetSubMatType(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCPatchGetSubMatType(petsclib::$UnionPetscLib, pc::PC )
	sub_mat_type_ = Ref{MatType}()

    @chk ccall(
               (:PCPatchGetSubMatType, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{MatType}),
               pc, sub_mat_type_,
              )

	sub_mat_type = unsafe_string(sub_mat_type_[])

	return sub_mat_type
end 

"""
	PCPatchSetCellNumbering(petsclib::PetscLibType,pc::PC, cellNumbering::PetscSection) 

# External Links
$(_doc_external("Ksp/PCPatchSetCellNumbering"))
"""
function PCPatchSetCellNumbering(petsclib::PetscLibType, pc::PC, cellNumbering::PetscSection) end

@for_petsc function PCPatchSetCellNumbering(petsclib::$UnionPetscLib, pc::PC, cellNumbering::PetscSection )

    @chk ccall(
               (:PCPatchSetCellNumbering, $petsc_library),
               PetscErrorCode,
               (PC, PetscSection),
               pc, cellNumbering,
              )


	return nothing
end 

"""
	PCPatchGetCellNumbering(petsclib::PetscLibType,pc::PC, cellNumbering::PetscSection) 

# External Links
$(_doc_external("Ksp/PCPatchGetCellNumbering"))
"""
function PCPatchGetCellNumbering(petsclib::PetscLibType, pc::PC, cellNumbering::PetscSection) end

@for_petsc function PCPatchGetCellNumbering(petsclib::$UnionPetscLib, pc::PC, cellNumbering::PetscSection )

    @chk ccall(
               (:PCPatchGetCellNumbering, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PetscSection}),
               pc, cellNumbering,
              )


	return nothing
end 

"""
	ctx::Cvoid = PCPatchSetConstructType(petsclib::PetscLibType,pc::PC, ctype::PCPatchConstructType, func::external) 

# External Links
$(_doc_external("Ksp/PCPatchSetConstructType"))
"""
function PCPatchSetConstructType(petsclib::PetscLibType, pc::PC, ctype::PCPatchConstructType, func::external) end

@for_petsc function PCPatchSetConstructType(petsclib::$UnionPetscLib, pc::PC, ctype::PCPatchConstructType, func::external )
	ctx_ = Ref{Cvoid}()

    @chk ccall(
               (:PCPatchSetConstructType, $petsc_library),
               PetscErrorCode,
               (PC, PCPatchConstructType, external, Ptr{Cvoid}),
               pc, ctype, func, ctx_,
              )

	ctx = ctx_[]

	return ctx
end 

"""
	bs::PetscInt,nodesPerCell::PetscInt,subspaceOffsets::PetscInt,ghostBcNodes::PetscInt,globalBcNodes::PetscInt = PCPatchSetDiscretisationInfo(petsclib::PetscLibType,pc::PC, nsubspaces::PetscInt, dms::PetscDM, cellNodeMap::PetscInt, numGhostBcs::PetscInt, numGlobalBcs::PetscInt) 

# External Links
$(_doc_external("Ksp/PCPatchSetDiscretisationInfo"))
"""
function PCPatchSetDiscretisationInfo(petsclib::PetscLibType, pc::PC, nsubspaces::PetscInt, dms::PetscDM, cellNodeMap::PetscInt, numGhostBcs::PetscInt, numGlobalBcs::PetscInt) end

@for_petsc function PCPatchSetDiscretisationInfo(petsclib::$UnionPetscLib, pc::PC, nsubspaces::$PetscInt, dms::PetscDM, cellNodeMap::$PetscInt, numGhostBcs::$PetscInt, numGlobalBcs::$PetscInt )
	dms_ = Ref(dms.ptr)
	bs_ = Ref{$PetscInt}()
	nodesPerCell_ = Ref{$PetscInt}()
	subspaceOffsets_ = Ref{$PetscInt}()
	ghostBcNodes_ = Ref{$PetscInt}()
	globalBcNodes_ = Ref{$PetscInt}()

    @chk ccall(
               (:PCPatchSetDiscretisationInfo, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, Ptr{CDM}, Ptr{$PetscInt}, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}, $PetscInt, Ptr{$PetscInt}),
               pc, nsubspaces, dms_, bs_, nodesPerCell_, cellNodeMap, subspaceOffsets_, numGhostBcs, ghostBcNodes_, numGlobalBcs, globalBcNodes_,
              )

	dms.ptr = C_NULL
	bs = bs_[]
	nodesPerCell = nodesPerCell_[]
	subspaceOffsets = subspaceOffsets_[]
	ghostBcNodes = ghostBcNodes_[]
	globalBcNodes = globalBcNodes_[]

	return bs,nodesPerCell,subspaceOffsets,ghostBcNodes,globalBcNodes
end 

"""
	PCPatchSetComputeFunction(petsclib::PetscLibType,pc::PC, func::external, ctx::Cvoid) 
Set the callback function used to compute patch residuals

Logically Collective

Input Parameters:
- `pc`   - The `PC`
- `func` - The callback function
- `ctx`  - The user context

Calling sequence of `func`:
- `pc`               - The `PC`
- `point`            - The point
- `x`                - The input solution (not used in linear problems)
- `f`                - The patch residual vector
- `cellIS`           - An array of the cell numbers
- `n`                - The size of `dofsArray`
- `dofsArray`        - The dofmap for the dofs to be solved for
- `dofsArrayWithAll` - The dofmap for all dofs on the patch
- `ctx`              - The user context

Level: advanced

-seealso: [](ch_ksp), `PCPatchSetComputeOperator()`, `PCPatchGetComputeOperator()`, `PCPatchSetDiscretisationInfo()`, `PCPatchSetComputeFunctionInteriorFacets()`

# External Links
$(_doc_external("Ksp/PCPatchSetComputeFunction"))
"""
function PCPatchSetComputeFunction(petsclib::PetscLibType, pc::PC, func::external, ctx::Cvoid) end

@for_petsc function PCPatchSetComputeFunction(petsclib::$UnionPetscLib, pc::PC, func::external, ctx::Cvoid )

    @chk ccall(
               (:PCPatchSetComputeFunction, $petsc_library),
               PetscErrorCode,
               (PC, external, Ptr{Cvoid}),
               pc, func, ctx,
              )


	return nothing
end 

"""
	PCPatchSetComputeFunctionInteriorFacets(petsclib::PetscLibType,pc::PC, func::external, ctx::Cvoid) 
Set the callback function used to compute facet integrals for patch residuals

Logically Collective

Input Parameters:
- `pc`   - The `PC`
- `func` - The callback function
- `ctx`  - The user context

Calling sequence of `func`:
- `pc`               - The `PC`
- `point`            - The point
- `x`                - The input solution (not used in linear problems)
- `f`                - The patch residual vector
- `facetIS`          - An array of the facet numbers
- `n`                - The size of `dofsArray`
- `dofsArray`        - The dofmap for the dofs to be solved for
- `dofsArrayWithAll` - The dofmap for all dofs on the patch
- `ctx`              - The user context

Level: advanced

-seealso: [](ch_ksp), `PCPatchSetComputeOperator()`, `PCPatchGetComputeOperator()`, `PCPatchSetDiscretisationInfo()`, `PCPatchSetComputeFunction()`

# External Links
$(_doc_external("Ksp/PCPatchSetComputeFunctionInteriorFacets"))
"""
function PCPatchSetComputeFunctionInteriorFacets(petsclib::PetscLibType, pc::PC, func::external, ctx::Cvoid) end

@for_petsc function PCPatchSetComputeFunctionInteriorFacets(petsclib::$UnionPetscLib, pc::PC, func::external, ctx::Cvoid )

    @chk ccall(
               (:PCPatchSetComputeFunctionInteriorFacets, $petsc_library),
               PetscErrorCode,
               (PC, external, Ptr{Cvoid}),
               pc, func, ctx,
              )


	return nothing
end 

"""
	PCPatchSetComputeOperator(petsclib::PetscLibType,pc::PC, func::external, ctx::Cvoid) 
Set the callback function used to compute patch matrices

Logically Collective

Input Parameters:
- `pc`   - The `PC`
- `func` - The callback function
- `ctx`  - The user context

Calling sequence of `func`:
- `pc`               - The `PC`
- `point`            - The point
- `x`                - The input solution (not used in linear problems)
- `mat`              - The patch matrix
- `facetIS`          - An array of the cell numbers
- `n`                - The size of `dofsArray`
- `dofsArray`        - The dofmap for the dofs to be solved for
- `dofsArrayWithAll` - The dofmap for all dofs on the patch
- `ctx`              - The user context

Level: advanced

-seealso: [](ch_ksp), `PCPatchGetComputeOperator()`, `PCPatchSetComputeFunction()`, `PCPatchSetDiscretisationInfo()`

# External Links
$(_doc_external("Ksp/PCPatchSetComputeOperator"))
"""
function PCPatchSetComputeOperator(petsclib::PetscLibType, pc::PC, func::external, ctx::Cvoid) end

@for_petsc function PCPatchSetComputeOperator(petsclib::$UnionPetscLib, pc::PC, func::external, ctx::Cvoid )

    @chk ccall(
               (:PCPatchSetComputeOperator, $petsc_library),
               PetscErrorCode,
               (PC, external, Ptr{Cvoid}),
               pc, func, ctx,
              )


	return nothing
end 

"""
	PCPatchSetComputeOperatorInteriorFacets(petsclib::PetscLibType,pc::PC, func::external, ctx::Cvoid) 
Set the callback function used to compute facet integrals for patch matrices

Logically Collective

Input Parameters:
- `pc`   - The `PC`
- `func` - The callback function
- `ctx`  - The user context

Calling sequence of `func`:
- `pc`               - The `PC`
- `point`            - The point
- `x`                - The input solution (not used in linear problems)
- `mat`              - The patch matrix
- `facetIS`          - An array of the facet numbers
- `n`                - The size of `dofsArray`
- `dofsArray`        - The dofmap for the dofs to be solved for
- `dofsArrayWithAll` - The dofmap for all dofs on the patch
- `ctx`              - The user context

Level: advanced

-seealso: [](ch_ksp), `PCPatchGetComputeOperator()`, `PCPatchSetComputeFunction()`, `PCPatchSetDiscretisationInfo()`

# External Links
$(_doc_external("Ksp/PCPatchSetComputeOperatorInteriorFacets"))
"""
function PCPatchSetComputeOperatorInteriorFacets(petsclib::PetscLibType, pc::PC, func::external, ctx::Cvoid) end

@for_petsc function PCPatchSetComputeOperatorInteriorFacets(petsclib::$UnionPetscLib, pc::PC, func::external, ctx::Cvoid )

    @chk ccall(
               (:PCPatchSetComputeOperatorInteriorFacets, $petsc_library),
               PetscErrorCode,
               (PC, external, Ptr{Cvoid}),
               pc, func, ctx,
              )


	return nothing
end 

"""
	n_loc::PetscInt,first_loc::PetscInt = PCBJacobiGetSubKSP(petsclib::PetscLibType,pc::PC, ksp::Vector{PetscKSP}) 
Gets the local `KSP` contexts for all blocks on
this processor.

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameters:
- `n_local`     - the number of blocks on this processor, or NULL
- `first_local` - the global number of the first block on this processor, or NULL
- `ksp`         - the array of KSP contexts

Level: advanced

-seealso: [](ch_ksp), `PCBJACOBI`, `PCASM`, `PCASMGetSubKSP()`

# External Links
$(_doc_external("Ksp/PCBJacobiGetSubKSP"))
"""
function PCBJacobiGetSubKSP(petsclib::PetscLibType, pc::PC, ksp::Vector{PetscKSP}) end

@for_petsc function PCBJacobiGetSubKSP(petsclib::$UnionPetscLib, pc::PC, ksp::Vector{PetscKSP} )
	n_loc_ = Ref{$PetscInt}()
	first_loc_ = Ref{$PetscInt}()
	ksp_ = Ref(pointer(ksp))

    @chk ccall(
               (:PCBJacobiGetSubKSP, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{CKSP}}),
               pc, n_loc_, first_loc_, ksp_,
              )

	n_loc = n_loc_[]
	first_loc = first_loc_[]

	return n_loc,first_loc
end 

"""
	PCBJacobiSetTotalBlocks(petsclib::PetscLibType,pc::PC, blocks::PetscInt, lens::Vector{PetscInt}) 
Sets the global number of blocks for the block
Jacobi preconditioner.

Collective

Input Parameters:
- `pc`     - the preconditioner context
- `blocks` - the number of blocks
- `lens`   - [optional] integer array containing the size of each block

Options Database Key:
- `-pc_bjacobi_blocks <blocks>` - Sets the number of global blocks

Level: intermediate

-seealso: [](ch_ksp), `PCBJACOBI`, `PCSetUseAmat()`, `PCBJacobiSetLocalBlocks()`

# External Links
$(_doc_external("Ksp/PCBJacobiSetTotalBlocks"))
"""
function PCBJacobiSetTotalBlocks(petsclib::PetscLibType, pc::PC, blocks::PetscInt, lens::Vector{PetscInt}) end

@for_petsc function PCBJacobiSetTotalBlocks(petsclib::$UnionPetscLib, pc::PC, blocks::$PetscInt, lens::Vector{$PetscInt} )

    @chk ccall(
               (:PCBJacobiSetTotalBlocks, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, Ptr{$PetscInt}),
               pc, blocks, lens,
              )


	return nothing
end 

"""
	blocks::PetscInt,lens::Vector{PetscInt} = PCBJacobiGetTotalBlocks(petsclib::PetscLibType,pc::PC) 
Gets the global number of blocks for the block
Jacobi, `PCBJACOBI`, preconditioner.

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameters:
- `blocks` - the number of blocks
- `lens`   - integer array containing the size of each block

Level: intermediate

-seealso: [](ch_ksp), `PCBJACOBI`, `PCSetUseAmat()`, `PCBJacobiGetLocalBlocks()`

# External Links
$(_doc_external("Ksp/PCBJacobiGetTotalBlocks"))
"""
function PCBJacobiGetTotalBlocks(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCBJacobiGetTotalBlocks(petsclib::$UnionPetscLib, pc::PC )
	blocks_ = Ref{$PetscInt}()
	lens_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:PCBJacobiGetTotalBlocks, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
               pc, blocks_, lens_,
              )

	blocks = blocks_[]
	lens = unsafe_wrap(Array, lens_[], VecGetLocalSize(petsclib, x); own = false)

	return blocks,lens
end 

"""
	PCBJacobiSetLocalBlocks(petsclib::PetscLibType,pc::PC, blocks::PetscInt, lens::Vector{PetscInt}) 
Sets the local number of blocks for the block
Jacobi, `PCBJACOBI`,  preconditioner.

Not Collective

Input Parameters:
- `pc`     - the preconditioner context
- `blocks` - the number of blocks
- `lens`   - [optional] integer array containing size of each block

Options Database Key:
- `-pc_bjacobi_local_blocks <blocks>` - Sets the number of local blocks

Level: intermediate

-seealso: [](ch_ksp), `PCBJACOBI`, `PCSetUseAmat()`, `PCBJacobiSetTotalBlocks()`

# External Links
$(_doc_external("Ksp/PCBJacobiSetLocalBlocks"))
"""
function PCBJacobiSetLocalBlocks(petsclib::PetscLibType, pc::PC, blocks::PetscInt, lens::Vector{PetscInt}) end

@for_petsc function PCBJacobiSetLocalBlocks(petsclib::$UnionPetscLib, pc::PC, blocks::$PetscInt, lens::Vector{$PetscInt} )

    @chk ccall(
               (:PCBJacobiSetLocalBlocks, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, Ptr{$PetscInt}),
               pc, blocks, lens,
              )


	return nothing
end 

"""
	PCBJacobiGetLocalBlocks(petsclib::PetscLibType,pc::PC, blocks::PetscInt, lens::Vector{PetscInt}) 
Gets the local number of blocks for the block
Jacobi, `PCBJACOBI`, preconditioner.

Not Collective

Input Parameters:
- `pc`     - the preconditioner context
- `blocks` - the number of blocks
- `lens`   - [optional] integer array containing size of each block

Level: intermediate

-seealso: [](ch_ksp), `PCBJACOBI`, `PCSetUseAmat()`, `PCBJacobiGetTotalBlocks()`

# External Links
$(_doc_external("Ksp/PCBJacobiGetLocalBlocks"))
"""
function PCBJacobiGetLocalBlocks(petsclib::PetscLibType, pc::PC, blocks::PetscInt, lens::Vector{PetscInt}) end

@for_petsc function PCBJacobiGetLocalBlocks(petsclib::$UnionPetscLib, pc::PC, blocks::$PetscInt, lens::Vector{$PetscInt} )
	lens_ = Ref(pointer(lens))

    @chk ccall(
               (:PCBJacobiGetLocalBlocks, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}),
               pc, blocks, lens_,
              )


	return nothing
end 

"""
	PCBJKOKKOSSetKSP(petsclib::PetscLibType,pc::PC, ksp::PetscKSP) 

# External Links
$(_doc_external("Ksp/PCBJKOKKOSSetKSP"))
"""
function PCBJKOKKOSSetKSP(petsclib::PetscLibType, pc::PC, ksp::PetscKSP) end

@for_petsc function PCBJKOKKOSSetKSP(petsclib::$UnionPetscLib, pc::PC, ksp::PetscKSP )

    @chk ccall(
               (:PCBJKOKKOSSetKSP, $petsc_library),
               PetscErrorCode,
               (PC, CKSP),
               pc, ksp,
              )


	return nothing
end 

"""
	PCBJKOKKOSGetKSP(petsclib::PetscLibType,pc::PC, ksp::PetscKSP) 

# External Links
$(_doc_external("Ksp/PCBJKOKKOSGetKSP"))
"""
function PCBJKOKKOSGetKSP(petsclib::PetscLibType, pc::PC, ksp::PetscKSP) end

@for_petsc function PCBJKOKKOSGetKSP(petsclib::$UnionPetscLib, pc::PC, ksp::PetscKSP )
	ksp_ = Ref(ksp.ptr)

    @chk ccall(
               (:PCBJKOKKOSGetKSP, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{CKSP}),
               pc, ksp_,
              )

	ksp.ptr = C_NULL

	return nothing
end 

"""
	PCSPAISetEpsilon(petsclib::PetscLibType,pc::PC, epsilon1::PetscReal) 


Input Parameters:
- `pc`       - the preconditioner
- `epsilon1` - the tolerance (default .4)

Level: intermediate

-seealso: [](ch_ksp), `PCSPAI`, `PCSetType()`

# External Links
$(_doc_external("Ksp/PCSPAISetEpsilon"))
"""
function PCSPAISetEpsilon(petsclib::PetscLibType, pc::PC, epsilon1::PetscReal) end

@for_petsc function PCSPAISetEpsilon(petsclib::$UnionPetscLib, pc::PC, epsilon1::$PetscReal )

    @chk ccall(
               (:PCSPAISetEpsilon, $petsc_library),
               PetscErrorCode,
               (PC, $PetscReal),
               pc, epsilon1,
              )


	return nothing
end 

"""
	PCSPAISetNBSteps(petsclib::PetscLibType,pc::PC, nbsteps1::PetscInt) 
set maximum number of improvement steps per row in
the `PCSPAI` preconditioner

Input Parameters:
- `pc`       - the preconditioner
- `nbsteps1` - number of steps (default 5)

-seealso: [](ch_ksp), `PCSPAI`, `PCSetType()`, `PCSPAISetMaxNew()`

# External Links
$(_doc_external("Ksp/PCSPAISetNBSteps"))
"""
function PCSPAISetNBSteps(petsclib::PetscLibType, pc::PC, nbsteps1::PetscInt) end

@for_petsc function PCSPAISetNBSteps(petsclib::$UnionPetscLib, pc::PC, nbsteps1::$PetscInt )

    @chk ccall(
               (:PCSPAISetNBSteps, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, nbsteps1,
              )


	return nothing
end 

"""
	PCSPAISetMax(petsclib::PetscLibType,pc::PC, max1::PetscInt) 
set the size of various working buffers in the `PCSPAI` preconditioner

Input Parameters:
- `pc`   - the preconditioner
- `max1` - size (default is 5000)

Level: intermediate

-seealso: [](ch_ksp), `PCSPAI`, `PCSetType()`

# External Links
$(_doc_external("Ksp/PCSPAISetMax"))
"""
function PCSPAISetMax(petsclib::PetscLibType, pc::PC, max1::PetscInt) end

@for_petsc function PCSPAISetMax(petsclib::$UnionPetscLib, pc::PC, max1::$PetscInt )

    @chk ccall(
               (:PCSPAISetMax, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, max1,
              )


	return nothing
end 

"""
	PCSPAISetMaxNew(petsclib::PetscLibType,pc::PC, maxnew1::PetscInt) 
set maximum number of new nonzero candidates per step in the `PCSPAI` preconditioner

Input Parameters:
- `pc`      - the preconditioner
- `maxnew1` - maximum number (default 5)

Level: intermediate

-seealso: [](ch_ksp), `PCSPAI`, `PCSetType()`, `PCSPAISetNBSteps()`

# External Links
$(_doc_external("Ksp/PCSPAISetMaxNew"))
"""
function PCSPAISetMaxNew(petsclib::PetscLibType, pc::PC, maxnew1::PetscInt) end

@for_petsc function PCSPAISetMaxNew(petsclib::$UnionPetscLib, pc::PC, maxnew1::$PetscInt )

    @chk ccall(
               (:PCSPAISetMaxNew, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, maxnew1,
              )


	return nothing
end 

"""
	PCSPAISetBlockSize(petsclib::PetscLibType,pc::PC, block_size1::PetscInt) 
set the block size for the `PCSPAI` preconditioner

Input Parameters:
- `pc`          - the preconditioner
- `block_size1` - block size (default 1)

Level: intermediate

-seealso: [](ch_ksp), `PCSPAI`, `PCSetType()`

# External Links
$(_doc_external("Ksp/PCSPAISetBlockSize"))
"""
function PCSPAISetBlockSize(petsclib::PetscLibType, pc::PC, block_size1::PetscInt) end

@for_petsc function PCSPAISetBlockSize(petsclib::$UnionPetscLib, pc::PC, block_size1::$PetscInt )

    @chk ccall(
               (:PCSPAISetBlockSize, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, block_size1,
              )


	return nothing
end 

"""
	PCSPAISetCacheSize(petsclib::PetscLibType,pc::PC, cache_size::PetscInt) 
specify cache size in the `PCSPAI` preconditioner

Input Parameters:
- `pc`         - the preconditioner
- `cache_size` - cache size {0,1,2,3,4,5} (default 5)

Level: intermediate

-seealso: [](ch_ksp), `PCSPAI`, `PCSetType()`

# External Links
$(_doc_external("Ksp/PCSPAISetCacheSize"))
"""
function PCSPAISetCacheSize(petsclib::PetscLibType, pc::PC, cache_size::PetscInt) end

@for_petsc function PCSPAISetCacheSize(petsclib::$UnionPetscLib, pc::PC, cache_size::$PetscInt )

    @chk ccall(
               (:PCSPAISetCacheSize, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, cache_size,
              )


	return nothing
end 

"""
	PCSPAISetVerbose(petsclib::PetscLibType,pc::PC, verbose::PetscInt) 
verbosity level for the `PCSPAI` preconditioner

Input Parameters:
- `pc`      - the preconditioner
- `verbose` - level (default 1)

Level: intermediate

-seealso: [](ch_ksp), `PCSPAI`, `PCSetType()`

# External Links
$(_doc_external("Ksp/PCSPAISetVerbose"))
"""
function PCSPAISetVerbose(petsclib::PetscLibType, pc::PC, verbose::PetscInt) end

@for_petsc function PCSPAISetVerbose(petsclib::$UnionPetscLib, pc::PC, verbose::$PetscInt )

    @chk ccall(
               (:PCSPAISetVerbose, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, verbose,
              )


	return nothing
end 

"""
	PCSPAISetSp(petsclib::PetscLibType,pc::PC, sp::PetscInt) 
specify a symmetric matrix sparsity pattern in the `PCSPAI` preconditioner

Input Parameters:
- `pc` - the preconditioner
- `sp` - 0 or 1

Level: intermediate

-seealso: [](ch_ksp), `PCSPAI`, `PCSetType()`

# External Links
$(_doc_external("Ksp/PCSPAISetSp"))
"""
function PCSPAISetSp(petsclib::PetscLibType, pc::PC, sp::PetscInt) end

@for_petsc function PCSPAISetSp(petsclib::$UnionPetscLib, pc::PC, sp::$PetscInt )

    @chk ccall(
               (:PCSPAISetSp, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, sp,
              )


	return nothing
end 

"""
	PCRedundantSetNumber(petsclib::PetscLibType,pc::PC, nredundant::PetscInt) 
Sets the number of redundant preconditioner contexts.

Logically Collective

Input Parameters:
- `pc`         - the preconditioner context
- `nredundant` - number of redundant preconditioner contexts; for example if you are using 64 MPI processes and
use an nredundant of 4 there will be 4 parallel solves each on 16 = 64/4 processes.

Level: advanced

-seealso: [](ch_ksp), `PCREDUNDANT`

# External Links
$(_doc_external("Ksp/PCRedundantSetNumber"))
"""
function PCRedundantSetNumber(petsclib::PetscLibType, pc::PC, nredundant::PetscInt) end

@for_petsc function PCRedundantSetNumber(petsclib::$UnionPetscLib, pc::PC, nredundant::$PetscInt )

    @chk ccall(
               (:PCRedundantSetNumber, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, nredundant,
              )


	return nothing
end 

"""
	PCRedundantSetScatter(petsclib::PetscLibType,pc::PC, in::VecScatter, out::VecScatter) 
Sets the scatter used to copy values into the
redundant local solve and the scatter to move them back into the global
vector.

Logically Collective

Input Parameters:
- `pc`  - the preconditioner context
- `in`  - the scatter to move the values in
- `out` - the scatter to move them out

Level: advanced

-seealso: [](ch_ksp), `PCREDUNDANT`

# External Links
$(_doc_external("Ksp/PCRedundantSetScatter"))
"""
function PCRedundantSetScatter(petsclib::PetscLibType, pc::PC, in::VecScatter, out::VecScatter) end

@for_petsc function PCRedundantSetScatter(petsclib::$UnionPetscLib, pc::PC, in::VecScatter, out::VecScatter )

    @chk ccall(
               (:PCRedundantSetScatter, $petsc_library),
               PetscErrorCode,
               (PC, VecScatter, VecScatter),
               pc, in, out,
              )


	return nothing
end 

"""
	PCRedundantGetKSP(petsclib::PetscLibType,pc::PC, innerksp::PetscKSP) 
Gets the less parallel `KSP` created by the redundant `PC`.

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `innerksp` - the `KSP` on the smaller set of processes

Level: advanced

-seealso: [](ch_ksp), `PCREDUNDANT`

# External Links
$(_doc_external("Ksp/PCRedundantGetKSP"))
"""
function PCRedundantGetKSP(petsclib::PetscLibType, pc::PC, innerksp::PetscKSP) end

@for_petsc function PCRedundantGetKSP(petsclib::$UnionPetscLib, pc::PC, innerksp::PetscKSP )
	innerksp_ = Ref(innerksp.ptr)

    @chk ccall(
               (:PCRedundantGetKSP, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{CKSP}),
               pc, innerksp_,
              )

	innerksp.ptr = C_NULL

	return nothing
end 

"""
	PCRedundantGetOperators(petsclib::PetscLibType,pc::PC, mat::PetscMat, pmat::PetscMat) 
gets the sequential linear system matrix and matrix used to construct the preconditioner

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameters:
- `mat`  - the matrix
- `pmat` - the (possibly different) matrix used to construct the preconditioner

Level: advanced

-seealso: [](ch_ksp), `PCREDUNDANT`

# External Links
$(_doc_external("Ksp/PCRedundantGetOperators"))
"""
function PCRedundantGetOperators(petsclib::PetscLibType, pc::PC, mat::PetscMat, pmat::PetscMat) end

@for_petsc function PCRedundantGetOperators(petsclib::$UnionPetscLib, pc::PC, mat::PetscMat, pmat::PetscMat )
	mat_ = Ref(mat.ptr)
	pmat_ = Ref(pmat.ptr)

    @chk ccall(
               (:PCRedundantGetOperators, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{CMat}, Ptr{CMat}),
               pc, mat_, pmat_,
              )

	mat.ptr = C_NULL
	pmat.ptr = C_NULL

	return nothing
end 

"""
	PCLMVMSetUpdateVec(petsclib::PetscLibType,pc::PC, X::PetscVec) 
Set the vector to be used as solution update for the internal LMVM matrix.

Input Parameters:
- `pc` - The preconditioner
- `X`  - Solution vector

Level: intermediate

-seealso: `MatLMVMUpdate()`, `PCLMVMSetMatLMVM()`

# External Links
$(_doc_external("Ksp/PCLMVMSetUpdateVec"))
"""
function PCLMVMSetUpdateVec(petsclib::PetscLibType, pc::PC, X::PetscVec) end

@for_petsc function PCLMVMSetUpdateVec(petsclib::$UnionPetscLib, pc::PC, X::PetscVec )

    @chk ccall(
               (:PCLMVMSetUpdateVec, $petsc_library),
               PetscErrorCode,
               (PC, CVec),
               pc, X,
              )


	return nothing
end 

"""
	PCLMVMSetMatLMVM(petsclib::PetscLibType,pc::PC, B::PetscMat) 
Replaces the `MATLMVM` matrix inside the preconditioner with the one provided by the user.

Input Parameters:
- `pc` - An `PCLMVM` preconditioner
- `B`  - An `MATLMVM` type matrix

Level: intermediate

-seealso: [](ch_ksp), `PCLMVMGetMatLMVM()`

# External Links
$(_doc_external("Ksp/PCLMVMSetMatLMVM"))
"""
function PCLMVMSetMatLMVM(petsclib::PetscLibType, pc::PC, B::PetscMat) end

@for_petsc function PCLMVMSetMatLMVM(petsclib::$UnionPetscLib, pc::PC, B::PetscMat )

    @chk ccall(
               (:PCLMVMSetMatLMVM, $petsc_library),
               PetscErrorCode,
               (PC, CMat),
               pc, B,
              )


	return nothing
end 

"""
	PCLMVMGetMatLMVM(petsclib::PetscLibType,pc::PC, B::PetscMat) 
Returns a pointer to the underlying `MATLMVM` matrix.

Input Parameter:
- `pc` - An `PCLMVM` preconditioner

Output Parameter:
- `B` - `MATLMVM` matrix used by the preconditioner

Level: intermediate

-seealso: [](ch_ksp), `PCLMVMSetMatLMVM()`

# External Links
$(_doc_external("Ksp/PCLMVMGetMatLMVM"))
"""
function PCLMVMGetMatLMVM(petsclib::PetscLibType, pc::PC, B::PetscMat) end

@for_petsc function PCLMVMGetMatLMVM(petsclib::$UnionPetscLib, pc::PC, B::PetscMat )
	B_ = Ref(B.ptr)

    @chk ccall(
               (:PCLMVMGetMatLMVM, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{CMat}),
               pc, B_,
              )

	B.ptr = C_NULL

	return nothing
end 

"""
	PCLMVMSetIS(petsclib::PetscLibType,pc::PC, inactive::IS) 
Sets the index sets that reduce the `PC` application.

Input Parameters:
- `pc`       - An `PCLMVM` preconditioner
- `inactive` - Index set defining the variables removed from the problem

Level: intermediate

-seealso: [](ch_ksp), `PCLMVMClearIS()`

# External Links
$(_doc_external("Ksp/PCLMVMSetIS"))
"""
function PCLMVMSetIS(petsclib::PetscLibType, pc::PC, inactive::IS) end

@for_petsc function PCLMVMSetIS(petsclib::$UnionPetscLib, pc::PC, inactive::IS )

    @chk ccall(
               (:PCLMVMSetIS, $petsc_library),
               PetscErrorCode,
               (PC, IS),
               pc, inactive,
              )


	return nothing
end 

"""
	PCLMVMClearIS(petsclib::PetscLibType,pc::PC) 
Removes the inactive variable index set from a `PCLMVM`

Input Parameter:
- `pc` - An `PCLMVM` preconditioner

Level: intermediate

-seealso: [](ch_ksp), `PCLMVMSetIS()`

# External Links
$(_doc_external("Ksp/PCLMVMClearIS"))
"""
function PCLMVMClearIS(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCLMVMClearIS(petsclib::$UnionPetscLib, pc::PC )

    @chk ccall(
               (:PCLMVMClearIS, $petsc_library),
               PetscErrorCode,
               (PC,),
               pc,
              )


	return nothing
end 

"""
	PCRedistributeGetKSP(petsclib::PetscLibType,pc::PC, innerksp::PetscKSP) 
Gets the `KSP` created by the `PCREDISTRIBUTE`

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `innerksp` - the inner `KSP`

Level: advanced

-seealso: [](ch_ksp), `KSP`, `PCREDISTRIBUTE`

# External Links
$(_doc_external("Ksp/PCRedistributeGetKSP"))
"""
function PCRedistributeGetKSP(petsclib::PetscLibType, pc::PC, innerksp::PetscKSP) end

@for_petsc function PCRedistributeGetKSP(petsclib::$UnionPetscLib, pc::PC, innerksp::PetscKSP )
	innerksp_ = Ref(innerksp.ptr)

    @chk ccall(
               (:PCRedistributeGetKSP, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{CKSP}),
               pc, innerksp_,
              )

	innerksp.ptr = C_NULL

	return nothing
end 

"""
	PCFieldSplitRestrictIS(petsclib::PetscLibType,pc::PC, isy::IS) 
Restricts the fieldsplit `IS`s to be within a given `IS`.

Input Parameters:
- `pc`  - the preconditioner context
- `isy` - the index set that defines the indices to which the fieldsplit is to be restricted

Level: advanced

-seealso: [](sec_block_matrices), `PCFIELDSPLIT`, `PCFieldSplitSetFields()`, `PCFieldSplitSetIS()`

# External Links
$(_doc_external("Ksp/PCFieldSplitRestrictIS"))
"""
function PCFieldSplitRestrictIS(petsclib::PetscLibType, pc::PC, isy::IS) end

@for_petsc function PCFieldSplitRestrictIS(petsclib::$UnionPetscLib, pc::PC, isy::IS )

    @chk ccall(
               (:PCFieldSplitRestrictIS, $petsc_library),
               PetscErrorCode,
               (PC, IS),
               pc, isy,
              )


	return nothing
end 

"""
	PCFieldSplitSetFields(petsclib::PetscLibType,pc::PC, splitname::String, n::PetscInt, fields::Vector{PetscInt}, fields_col::Vector{PetscInt}) 
Sets the fields that define one particular split in `PCFIELDSPLIT`

Logically Collective

Input Parameters:
- `pc`         - the preconditioner context
- `splitname`  - name of this split, if `NULL` the number of the split is used
- `n`          - the number of fields in this split
- `fields`     - the fields in this split
- `fields_col` - generally the same as `fields`, if it does not match `fields` then the submatrix that is solved for this set of fields comes from an off-diagonal block
of the matrix and `fields_col` provides the column indices for that block

Options Database Key:
- `-pc_fieldsplit_%d_fields <a,b,..>` - indicates the fields to be used in the `%d`'th split

Level: intermediate

-seealso: [](sec_block_matrices), `PC`, `PCFieldSplitGetSubKSP()`, `PCFIELDSPLIT`, `PCFieldSplitSetBlockSize()`, `PCFieldSplitSetIS()`, `PCFieldSplitRestrictIS()`,
`MatSetBlockSize()`, `MatCreateNest()`

# External Links
$(_doc_external("Ksp/PCFieldSplitSetFields"))
"""
function PCFieldSplitSetFields(petsclib::PetscLibType, pc::PC, splitname::String, n::PetscInt, fields::Vector{PetscInt}, fields_col::Vector{PetscInt}) end

@for_petsc function PCFieldSplitSetFields(petsclib::$UnionPetscLib, pc::PC, splitname::String, n::$PetscInt, fields::Vector{$PetscInt}, fields_col::Vector{$PetscInt} )

    @chk ccall(
               (:PCFieldSplitSetFields, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{Cchar}, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}),
               pc, splitname, n, fields, fields_col,
              )


	return nothing
end 

"""
	PCFieldSplitSetDiagUseAmat(petsclib::PetscLibType,pc::PC, flg::PetscBool) 
set flag indicating whether to extract diagonal blocks from Amat (rather than Pmat) to build
the sub-matrices associated with each split. Where `KSPSetOperators`(ksp,Amat,Pmat) was used to supply the operators.

Logically Collective

Input Parameters:
- `pc`  - the preconditioner object
- `flg` - boolean flag indicating whether or not to use Amat to extract the diagonal blocks from

Options Database Key:
- `-pc_fieldsplit_diag_use_amat` - use the Amat to provide the diagonal blocks

Level: intermediate

-seealso: [](sec_block_matrices), `PC`, `PCSetOperators()`, `KSPSetOperators()`, `PCFieldSplitGetDiagUseAmat()`, `PCFieldSplitSetOffDiagUseAmat()`, `PCFIELDSPLIT`

# External Links
$(_doc_external("Ksp/PCFieldSplitSetDiagUseAmat"))
"""
function PCFieldSplitSetDiagUseAmat(petsclib::PetscLibType, pc::PC, flg::PetscBool) end

@for_petsc function PCFieldSplitSetDiagUseAmat(petsclib::$UnionPetscLib, pc::PC, flg::PetscBool )

    @chk ccall(
               (:PCFieldSplitSetDiagUseAmat, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, flg,
              )


	return nothing
end 

"""
	flg::PetscBool = PCFieldSplitGetDiagUseAmat(petsclib::PetscLibType,pc::PC) 
get the flag indicating whether to extract diagonal blocks from Amat (rather than Pmat) to build
the sub-matrices associated with each split.  Where `KSPSetOperators`(ksp,Amat,Pmat) was used to supply the operators.

Logically Collective

Input Parameter:
- `pc` - the preconditioner object

Output Parameter:
- `flg` - boolean flag indicating whether or not to use Amat to extract the diagonal blocks from

Level: intermediate

-seealso: [](sec_block_matrices), `PC`, `PCSetOperators()`, `KSPSetOperators()`, `PCFieldSplitSetDiagUseAmat()`, `PCFieldSplitGetOffDiagUseAmat()`, `PCFIELDSPLIT`

# External Links
$(_doc_external("Ksp/PCFieldSplitGetDiagUseAmat"))
"""
function PCFieldSplitGetDiagUseAmat(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCFieldSplitGetDiagUseAmat(petsclib::$UnionPetscLib, pc::PC )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PCFieldSplitGetDiagUseAmat, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PetscBool}),
               pc, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	PCFieldSplitSetOffDiagUseAmat(petsclib::PetscLibType,pc::PC, flg::PetscBool) 
set flag indicating whether to extract off
the sub-matrices associated with each split.  Where `KSPSetOperators`(ksp,Amat,Pmat) was used to supply the operators.

Logically Collective

Input Parameters:
- `pc`  - the preconditioner object
- `flg` - boolean flag indicating whether or not to use Amat to extract the off-diagonal blocks from

Options Database Key:
- `-pc_fieldsplit_off_diag_use_amat <bool>` - use the Amat to extract the off-diagonal blocks

Level: intermediate

-seealso: [](sec_block_matrices), `PC`, `PCSetOperators()`, `KSPSetOperators()`, `PCFieldSplitGetOffDiagUseAmat()`, `PCFieldSplitSetDiagUseAmat()`, `PCFIELDSPLIT`

# External Links
$(_doc_external("Ksp/PCFieldSplitSetOffDiagUseAmat"))
"""
function PCFieldSplitSetOffDiagUseAmat(petsclib::PetscLibType, pc::PC, flg::PetscBool) end

@for_petsc function PCFieldSplitSetOffDiagUseAmat(petsclib::$UnionPetscLib, pc::PC, flg::PetscBool )

    @chk ccall(
               (:PCFieldSplitSetOffDiagUseAmat, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, flg,
              )


	return nothing
end 

"""
	flg::PetscBool = PCFieldSplitGetOffDiagUseAmat(petsclib::PetscLibType,pc::PC) 
get the flag indicating whether to extract off
the sub-matrices associated with each split.  Where `KSPSetOperators`(ksp,Amat,Pmat) was used to supply the operators.

Logically Collective

Input Parameter:
- `pc` - the preconditioner object

Output Parameter:
- `flg` - boolean flag indicating whether or not to use Amat to extract the off-diagonal blocks from

Level: intermediate

-seealso: [](sec_block_matrices), `PC`, `PCSetOperators()`, `KSPSetOperators()`, `PCFieldSplitSetOffDiagUseAmat()`, `PCFieldSplitGetDiagUseAmat()`, `PCFIELDSPLIT`

# External Links
$(_doc_external("Ksp/PCFieldSplitGetOffDiagUseAmat"))
"""
function PCFieldSplitGetOffDiagUseAmat(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCFieldSplitGetOffDiagUseAmat(petsclib::$UnionPetscLib, pc::PC )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PCFieldSplitGetOffDiagUseAmat, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PetscBool}),
               pc, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	PCFieldSplitSetIS(petsclib::PetscLibType,pc::PC, splitname::String, is::IS) 
Sets the exact elements for a split in a `PCFIELDSPLIT`

Logically Collective

Input Parameters:
- `pc`        - the preconditioner context
- `splitname` - name of this split, if `NULL` the number of the split is used
- `is`        - the index set that defines the elements in this split

Level: intermediate

-seealso: [](sec_block_matrices), `PC`, `PCFieldSplitGetSubKSP()`, `PCFIELDSPLIT`, `PCFieldSplitSetBlockSize()`, `PCFieldSplitSetFields()`

# External Links
$(_doc_external("Ksp/PCFieldSplitSetIS"))
"""
function PCFieldSplitSetIS(petsclib::PetscLibType, pc::PC, splitname::String, is::IS) end

@for_petsc function PCFieldSplitSetIS(petsclib::$UnionPetscLib, pc::PC, splitname::String, is::IS )

    @chk ccall(
               (:PCFieldSplitSetIS, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{Cchar}, IS),
               pc, splitname, is,
              )


	return nothing
end 

"""
	PCFieldSplitGetIS(petsclib::PetscLibType,pc::PC, splitname::String, is::IS) 
Retrieves the elements for a split as an `IS`

Logically Collective

Input Parameters:
- `pc`        - the preconditioner context
- `splitname` - name of this split

Output Parameter:
- `is` - the index set that defines the elements in this split, or `NULL` if the split is not found

Level: intermediate

-seealso: [](sec_block_matrices), `PC`, `PCFieldSplitGetSubKSP()`, `PCFIELDSPLIT`, `PCFieldSplitSetIS()`, `PCFieldSplitGetISByIndex()`

# External Links
$(_doc_external("Ksp/PCFieldSplitGetIS"))
"""
function PCFieldSplitGetIS(petsclib::PetscLibType, pc::PC, splitname::String, is::IS) end

@for_petsc function PCFieldSplitGetIS(petsclib::$UnionPetscLib, pc::PC, splitname::String, is::IS )

    @chk ccall(
               (:PCFieldSplitGetIS, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{Cchar}, Ptr{IS}),
               pc, splitname, is,
              )


	return nothing
end 

"""
	PCFieldSplitGetISByIndex(petsclib::PetscLibType,pc::PC, index::PetscInt, is::IS) 
Retrieves the elements for a given split as an `IS`

Logically Collective

Input Parameters:
- `pc`    - the preconditioner context
- `index` - index of this split

Output Parameter:
- `is` - the index set that defines the elements in this split

Level: intermediate

-seealso: [](sec_block_matrices), `PC`, `PCFieldSplitGetSubKSP()`, `PCFIELDSPLIT`, `PCFieldSplitGetIS()`, `PCFieldSplitSetIS()`,


# External Links
$(_doc_external("Ksp/PCFieldSplitGetISByIndex"))
"""
function PCFieldSplitGetISByIndex(petsclib::PetscLibType, pc::PC, index::PetscInt, is::IS) end

@for_petsc function PCFieldSplitGetISByIndex(petsclib::$UnionPetscLib, pc::PC, index::$PetscInt, is::IS )

    @chk ccall(
               (:PCFieldSplitGetISByIndex, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, Ptr{IS}),
               pc, index, is,
              )


	return nothing
end 

"""
	PCFieldSplitSetBlockSize(petsclib::PetscLibType,pc::PC, bs::PetscInt) 
Sets the block size for defining where fields start in the
fieldsplit preconditioner when calling `PCFieldSplitSetFields()`. If not set the matrix block size is used.

Logically Collective

Input Parameters:
- `pc` - the preconditioner context
- `bs` - the block size

Level: intermediate

-seealso: [](sec_block_matrices), `PC`, `PCFieldSplitGetSubKSP()`, `PCFIELDSPLIT`, `PCFieldSplitSetFields()`, `PCFieldSplitSetIS()`

# External Links
$(_doc_external("Ksp/PCFieldSplitSetBlockSize"))
"""
function PCFieldSplitSetBlockSize(petsclib::PetscLibType, pc::PC, bs::PetscInt) end

@for_petsc function PCFieldSplitSetBlockSize(petsclib::$UnionPetscLib, pc::PC, bs::$PetscInt )

    @chk ccall(
               (:PCFieldSplitSetBlockSize, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, bs,
              )


	return nothing
end 

"""
	n::PetscInt = PCFieldSplitGetSubKSP(petsclib::PetscLibType,pc::PC, subksp::Vector{PetscKSP}) 
Gets the `KSP` contexts for all splits

Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameters:
- `n`      - the number of splits
- `subksp` - the array of `KSP` contexts

Level: advanced

-seealso: [](sec_block_matrices), `PC`, `PCFIELDSPLIT`, `PCFieldSplitSetFields()`, `PCFieldSplitSetIS()`, `PCFieldSplitSchurGetSubKSP()`

# External Links
$(_doc_external("Ksp/PCFieldSplitGetSubKSP"))
"""
function PCFieldSplitGetSubKSP(petsclib::PetscLibType, pc::PC, subksp::Vector{PetscKSP}) end

@for_petsc function PCFieldSplitGetSubKSP(petsclib::$UnionPetscLib, pc::PC, subksp::Vector{PetscKSP} )
	n_ = Ref{$PetscInt}()
	subksp_ = Ref(pointer(subksp))

    @chk ccall(
               (:PCFieldSplitGetSubKSP, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{$PetscInt}, Ptr{Ptr{CKSP}}),
               pc, n_, subksp_,
              )

	n = n_[]

	return n
end 

"""
	n::PetscInt = PCFieldSplitSchurGetSubKSP(petsclib::PetscLibType,pc::PC, subksp::Vector{PetscKSP}) 
Gets the `KSP` contexts used inside the Schur complement based `PCFIELDSPLIT`

Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameters:
- `n`      - the number of splits
- `subksp` - the array of `KSP` contexts

Level: advanced

-seealso: [](sec_block_matrices), `PC`, `PCFIELDSPLIT`, `PCFieldSplitSetFields()`, `PCFieldSplitSetIS()`, `PCFieldSplitGetSubKSP()`

# External Links
$(_doc_external("Ksp/PCFieldSplitSchurGetSubKSP"))
"""
function PCFieldSplitSchurGetSubKSP(petsclib::PetscLibType, pc::PC, subksp::Vector{PetscKSP}) end

@for_petsc function PCFieldSplitSchurGetSubKSP(petsclib::$UnionPetscLib, pc::PC, subksp::Vector{PetscKSP} )
	n_ = Ref{$PetscInt}()
	subksp_ = Ref(pointer(subksp))

    @chk ccall(
               (:PCFieldSplitSchurGetSubKSP, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{$PetscInt}, Ptr{Ptr{CKSP}}),
               pc, n_, subksp_,
              )

	n = n_[]

	return n
end 

"""
	PCFieldSplitSetSchurPre(petsclib::PetscLibType,pc::PC, ptype::PCFieldSplitSchurPreType, pre::PetscMat) 
Indicates from what operator the preconditioner is constructed for the Schur complement.
The default is the A11 matrix.

Collective

Input Parameters:
- `pc`    - the preconditioner context
- `ptype` - which matrix to use for preconditioning the Schur complement: `PC_FIELDSPLIT_SCHUR_PRE_A11` (default),
`PC_FIELDSPLIT_SCHUR_PRE_SELF`, `PC_FIELDSPLIT_SCHUR_PRE_USER`,
`PC_FIELDSPLIT_SCHUR_PRE_SELFP`, and `PC_FIELDSPLIT_SCHUR_PRE_FULL`
- `pre`   - matrix to use for preconditioning, or `NULL`

Options Database Keys:
- `-pc_fieldsplit_schur_precondition <self,selfp,user,a11,full>` - default is `a11`. See notes for meaning of various arguments
- `-fieldsplit_1_pc_type <pctype>`                               - the preconditioner algorithm that is used to construct the preconditioner from the operator

Level: intermediate

-seealso: [](sec_block_matrices), `PC`, `PCFieldSplitGetSchurPre()`, `PCFieldSplitGetSubKSP()`, `PCFIELDSPLIT`, `PCFieldSplitSetFields()`, `PCFieldSplitSchurPreType`,
`MatSchurComplementSetAinvType()`, `PCLSC`, `PCFieldSplitSetSchurFactType()`

# External Links
$(_doc_external("Ksp/PCFieldSplitSetSchurPre"))
"""
function PCFieldSplitSetSchurPre(petsclib::PetscLibType, pc::PC, ptype::PCFieldSplitSchurPreType, pre::PetscMat) end

@for_petsc function PCFieldSplitSetSchurPre(petsclib::$UnionPetscLib, pc::PC, ptype::PCFieldSplitSchurPreType, pre::PetscMat )

    @chk ccall(
               (:PCFieldSplitSetSchurPre, $petsc_library),
               PetscErrorCode,
               (PC, PCFieldSplitSchurPreType, CMat),
               pc, ptype, pre,
              )


	return nothing
end 

"""
	PCFieldSplitGetSchurPre(petsclib::PetscLibType,pc::PC, ptype::PCFieldSplitSchurPreType, pre::PetscMat) 
For Schur complement fieldsplit, determine how the Schur complement will be
preconditioned.  See `PCFieldSplitSetSchurPre()` for details.

Logically Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameters:
- `ptype` - which matrix to use for preconditioning the Schur complement: `PC_FIELDSPLIT_SCHUR_PRE_A11`, `PC_FIELDSPLIT_SCHUR_PRE_SELF`, `PC_FIELDSPLIT_SCHUR_PRE_USER`
- `pre`   - matrix to use for preconditioning (with `PC_FIELDSPLIT_SCHUR_PRE_USER`), or `NULL`

Level: intermediate

-seealso: [](sec_block_matrices), `PC`, `PCFieldSplitSetSchurPre()`, `PCFieldSplitGetSubKSP()`, `PCFIELDSPLIT`, `PCFieldSplitSetFields()`, `PCFieldSplitSchurPreType`, `PCLSC`

# External Links
$(_doc_external("Ksp/PCFieldSplitGetSchurPre"))
"""
function PCFieldSplitGetSchurPre(petsclib::PetscLibType, pc::PC, ptype::PCFieldSplitSchurPreType, pre::PetscMat) end

@for_petsc function PCFieldSplitGetSchurPre(petsclib::$UnionPetscLib, pc::PC, ptype::PCFieldSplitSchurPreType, pre::PetscMat )
	pre_ = Ref(pre.ptr)

    @chk ccall(
               (:PCFieldSplitGetSchurPre, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PCFieldSplitSchurPreType}, Ptr{CMat}),
               pc, ptype, pre_,
              )

	pre.ptr = C_NULL

	return nothing
end 

"""
	PCFieldSplitSchurGetS(petsclib::PetscLibType,pc::PC, S::PetscMat) 
extract the `MATSCHURCOMPLEMENT` object used by this `PCFIELDSPLIT` in case it needs to be configured separately

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `S` - the Schur complement matrix

Level: advanced

-seealso: [](sec_block_matrices), `PC`, `PCFieldSplitGetSubKSP()`, `PCFIELDSPLIT`, `PCFieldSplitSchurPreType`, `PCFieldSplitSetSchurPre()`, `MATSCHURCOMPLEMENT`, `PCFieldSplitSchurRestoreS()`,
`MatCreateSchurComplement()`, `MatSchurComplementGetKSP()`, `MatSchurComplementComputeExplicitOperator()`, `MatGetSchurComplement()`

# External Links
$(_doc_external("Ksp/PCFieldSplitSchurGetS"))
"""
function PCFieldSplitSchurGetS(petsclib::PetscLibType, pc::PC, S::PetscMat) end

@for_petsc function PCFieldSplitSchurGetS(petsclib::$UnionPetscLib, pc::PC, S::PetscMat )
	S_ = Ref(S.ptr)

    @chk ccall(
               (:PCFieldSplitSchurGetS, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{CMat}),
               pc, S_,
              )

	S.ptr = C_NULL

	return nothing
end 

"""
	PCFieldSplitSchurRestoreS(petsclib::PetscLibType,pc::PC, S::PetscMat) 
returns the `MATSCHURCOMPLEMENT` matrix used by this `PC`

Not Collective

Input Parameters:
- `pc` - the preconditioner context
- `S`  - the Schur complement matrix

Level: advanced

-seealso: [](sec_block_matrices), `PC`, `PCFieldSplitGetSubKSP()`, `PCFIELDSPLIT`, `PCFieldSplitSchurPreType`, `PCFieldSplitSetSchurPre()`, `MatSchurComplement`, `PCFieldSplitSchurGetS()`

# External Links
$(_doc_external("Ksp/PCFieldSplitSchurRestoreS"))
"""
function PCFieldSplitSchurRestoreS(petsclib::PetscLibType, pc::PC, S::PetscMat) end

@for_petsc function PCFieldSplitSchurRestoreS(petsclib::$UnionPetscLib, pc::PC, S::PetscMat )
	S_ = Ref(S.ptr)

    @chk ccall(
               (:PCFieldSplitSchurRestoreS, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{CMat}),
               pc, S_,
              )

	S.ptr = C_NULL

	return nothing
end 

"""
	PCFieldSplitSetSchurFactType(petsclib::PetscLibType,pc::PC, ftype::PCFieldSplitSchurFactType) 
sets which blocks of the approximate block factorization to retain in the preconditioner {cite}`murphy2000note` and {cite}`ipsen2001note`

Collective

Input Parameters:
- `pc`    - the preconditioner context
- `ftype` - which blocks of factorization to retain, `PC_FIELDSPLIT_SCHUR_FACT_FULL` is default

Options Database Key:
- `-pc_fieldsplit_schur_fact_type <diag,lower,upper,full>` - default is `full`

Level: intermediate

-seealso: [](sec_block_matrices), `PC`, `PCFieldSplitGetSubKSP()`, `PCFIELDSPLIT`, `PCFieldSplitSetFields()`, `PCFieldSplitSchurPreType`, `PCFieldSplitSetSchurScale()`,
[](sec_flexibleksp), `PCFieldSplitSetSchurPre()`

# External Links
$(_doc_external("Ksp/PCFieldSplitSetSchurFactType"))
"""
function PCFieldSplitSetSchurFactType(petsclib::PetscLibType, pc::PC, ftype::PCFieldSplitSchurFactType) end

@for_petsc function PCFieldSplitSetSchurFactType(petsclib::$UnionPetscLib, pc::PC, ftype::PCFieldSplitSchurFactType )

    @chk ccall(
               (:PCFieldSplitSetSchurFactType, $petsc_library),
               PetscErrorCode,
               (PC, PCFieldSplitSchurFactType),
               pc, ftype,
              )


	return nothing
end 

"""
	PCFieldSplitSetSchurScale(petsclib::PetscLibType,pc::PC, scale::PetscScalar) 
Controls the sign flip of S for `PC_FIELDSPLIT_SCHUR_FACT_DIAG`.

Collective

Input Parameters:
- `pc`    - the preconditioner context
- `scale` - scaling factor for the Schur complement

Options Database Key:
- `-pc_fieldsplit_schur_scale <scale>` - default is -1.0

Level: intermediate

-seealso: [](sec_block_matrices), `PC`, `PCFIELDSPLIT`, `PCFieldSplitSetFields()`, `PCFieldSplitSchurFactType`, `PCFieldSplitSetSchurFactType()`

# External Links
$(_doc_external("Ksp/PCFieldSplitSetSchurScale"))
"""
function PCFieldSplitSetSchurScale(petsclib::PetscLibType, pc::PC, scale::PetscScalar) end

@for_petsc function PCFieldSplitSetSchurScale(petsclib::$UnionPetscLib, pc::PC, scale::$PetscScalar )

    @chk ccall(
               (:PCFieldSplitSetSchurScale, $petsc_library),
               PetscErrorCode,
               (PC, $PetscScalar),
               pc, scale,
              )


	return nothing
end 

"""
	PCFieldSplitGetSchurBlocks(petsclib::PetscLibType,pc::PC, A00::PetscMat, A01::PetscMat, A10::PetscMat, A11::PetscMat) 
Gets all matrix blocks for the Schur complement

Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameters:
- `A00` - the (0,0) block
- `A01` - the (0,1) block
- `A10` - the (1,0) block
- `A11` - the (1,1) block

Level: advanced

-seealso: [](sec_block_matrices), `PC`, `PCFIELDSPLIT`, `MatSchurComplementGetSubMatrices()`, `MatSchurComplementSetSubMatrices()`

# External Links
$(_doc_external("Ksp/PCFieldSplitGetSchurBlocks"))
"""
function PCFieldSplitGetSchurBlocks(petsclib::PetscLibType, pc::PC, A00::PetscMat, A01::PetscMat, A10::PetscMat, A11::PetscMat) end

@for_petsc function PCFieldSplitGetSchurBlocks(petsclib::$UnionPetscLib, pc::PC, A00::PetscMat, A01::PetscMat, A10::PetscMat, A11::PetscMat )
	A00_ = Ref(A00.ptr)
	A01_ = Ref(A01.ptr)
	A10_ = Ref(A10.ptr)
	A11_ = Ref(A11.ptr)

    @chk ccall(
               (:PCFieldSplitGetSchurBlocks, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{CMat}, Ptr{CMat}, Ptr{CMat}, Ptr{CMat}),
               pc, A00_, A01_, A10_, A11_,
              )

	A00.ptr = C_NULL
	A01.ptr = C_NULL
	A10.ptr = C_NULL
	A11.ptr = C_NULL

	return nothing
end 

"""
	PCFieldSplitSetGKBTol(petsclib::PetscLibType,pc::PC, tolerance::PetscReal) 
Sets the solver tolerance for the generalized Golub

Collective

Input Parameters:
- `pc`        - the preconditioner context
- `tolerance` - the solver tolerance

Options Database Key:
- `-pc_fieldsplit_gkb_tol <tolerance>` - default is 1e-5

Level: intermediate

-seealso: [](sec_block_matrices), `PC`, `PCFIELDSPLIT`, `PCFieldSplitSetGKBDelay()`, `PCFieldSplitSetGKBNu()`, `PCFieldSplitSetGKBMaxit()`

# External Links
$(_doc_external("Ksp/PCFieldSplitSetGKBTol"))
"""
function PCFieldSplitSetGKBTol(petsclib::PetscLibType, pc::PC, tolerance::PetscReal) end

@for_petsc function PCFieldSplitSetGKBTol(petsclib::$UnionPetscLib, pc::PC, tolerance::$PetscReal )

    @chk ccall(
               (:PCFieldSplitSetGKBTol, $petsc_library),
               PetscErrorCode,
               (PC, $PetscReal),
               pc, tolerance,
              )


	return nothing
end 

"""
	PCFieldSplitSetGKBMaxit(petsclib::PetscLibType,pc::PC, maxit::PetscInt) 
Sets the maximum number of iterations for the generalized Golub

Collective

Input Parameters:
- `pc`    - the preconditioner context
- `maxit` - the maximum number of iterations

Options Database Key:
- `-pc_fieldsplit_gkb_maxit <maxit>` - default is 100

Level: intermediate

-seealso: [](sec_block_matrices), `PC`, `PCFIELDSPLIT`, `PCFieldSplitSetGKBDelay()`, `PCFieldSplitSetGKBTol()`, `PCFieldSplitSetGKBNu()`

# External Links
$(_doc_external("Ksp/PCFieldSplitSetGKBMaxit"))
"""
function PCFieldSplitSetGKBMaxit(petsclib::PetscLibType, pc::PC, maxit::PetscInt) end

@for_petsc function PCFieldSplitSetGKBMaxit(petsclib::$UnionPetscLib, pc::PC, maxit::$PetscInt )

    @chk ccall(
               (:PCFieldSplitSetGKBMaxit, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, maxit,
              )


	return nothing
end 

"""
	PCFieldSplitSetGKBDelay(petsclib::PetscLibType,pc::PC, delay::PetscInt) 
Sets the delay in the lower bound error estimate in the generalized Golub
preconditioner.

Collective

Input Parameters:
- `pc`    - the preconditioner context
- `delay` - the delay window in the lower bound estimate

Options Database Key:
- `-pc_fieldsplit_gkb_delay <delay>` - default is 5

Level: intermediate

-seealso: [](sec_block_matrices), `PC`, `PCFIELDSPLIT`, `PCFieldSplitSetGKBNu()`, `PCFieldSplitSetGKBTol()`, `PCFieldSplitSetGKBMaxit()`

# External Links
$(_doc_external("Ksp/PCFieldSplitSetGKBDelay"))
"""
function PCFieldSplitSetGKBDelay(petsclib::PetscLibType, pc::PC, delay::PetscInt) end

@for_petsc function PCFieldSplitSetGKBDelay(petsclib::$UnionPetscLib, pc::PC, delay::$PetscInt )

    @chk ccall(
               (:PCFieldSplitSetGKBDelay, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, delay,
              )


	return nothing
end 

"""
	PCFieldSplitSetGKBNu(petsclib::PetscLibType,pc::PC, nu::PetscReal) 
Sets the scalar value nu >= 0 in the transformation H = A00 + nu*A01*A01' of the (1,1) block in the
Golub-Kahan bidiagonalization preconditioner {cite}`arioli2013` in `PCFIELDSPLIT`

Collective

Input Parameters:
- `pc` - the preconditioner context
- `nu` - the shift parameter

Options Database Key:
- `-pc_fieldsplit_gkb_nu <nu>` - default is 1

Level: intermediate

-seealso: [](sec_block_matrices), `PC`, `PCFIELDSPLIT`, `PCFieldSplitSetGKBDelay()`, `PCFieldSplitSetGKBTol()`, `PCFieldSplitSetGKBMaxit()`

# External Links
$(_doc_external("Ksp/PCFieldSplitSetGKBNu"))
"""
function PCFieldSplitSetGKBNu(petsclib::PetscLibType, pc::PC, nu::PetscReal) end

@for_petsc function PCFieldSplitSetGKBNu(petsclib::$UnionPetscLib, pc::PC, nu::$PetscReal )

    @chk ccall(
               (:PCFieldSplitSetGKBNu, $petsc_library),
               PetscErrorCode,
               (PC, $PetscReal),
               pc, nu,
              )


	return nothing
end 

"""
	PCFieldSplitSetType(petsclib::PetscLibType,pc::PC, type::PCCompositeType) 
Sets the type, `PCCompositeType`, of a `PCFIELDSPLIT`

Collective

Input Parameters:
- `pc`   - the preconditioner context
- `type` - `PC_COMPOSITE_ADDITIVE`, `PC_COMPOSITE_MULTIPLICATIVE` (default), `PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE`, `PC_COMPOSITE_SPECIAL`, `PC_COMPOSITE_SCHUR`,
`PC_COMPOSITE_GKB`

Options Database Key:
- `-pc_fieldsplit_type <one of multiplicative, additive, symmetric_multiplicative, special, schur>` - Sets fieldsplit preconditioner type

Level: intermediate

-seealso: [](sec_block_matrices), `PC`, `PCFIELDSPLIT`, `PCCompositeType`, `PCCompositeGetType()`, `PC_COMPOSITE_ADDITIVE`, `PC_COMPOSITE_MULTIPLICATIVE`,
`PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE`, `PC_COMPOSITE_SPECIAL`, `PC_COMPOSITE_SCHUR`, `PCFieldSplitSetSchurFactType()`

# External Links
$(_doc_external("Ksp/PCFieldSplitSetType"))
"""
function PCFieldSplitSetType(petsclib::PetscLibType, pc::PC, type::PCCompositeType) end

@for_petsc function PCFieldSplitSetType(petsclib::$UnionPetscLib, pc::PC, type::PCCompositeType )

    @chk ccall(
               (:PCFieldSplitSetType, $petsc_library),
               PetscErrorCode,
               (PC, PCCompositeType),
               pc, type,
              )


	return nothing
end 

"""
	type::PCCompositeType = PCFieldSplitGetType(petsclib::PetscLibType,pc::PC) 
Gets the type, `PCCompositeType`, of a `PCFIELDSPLIT`

Not collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `type` - `PC_COMPOSITE_ADDITIVE`, `PC_COMPOSITE_MULTIPLICATIVE` (default), `PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE`, `PC_COMPOSITE_SPECIAL`, `PC_COMPOSITE_SCHUR`

Level: intermediate

-seealso: [](sec_block_matrices), `PC`, `PCCompositeSetType()`, `PCFIELDSPLIT`, `PCCompositeType`, `PC_COMPOSITE_ADDITIVE`, `PC_COMPOSITE_MULTIPLICATIVE`,
`PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE`, `PC_COMPOSITE_SPECIAL`, `PC_COMPOSITE_SCHUR`

# External Links
$(_doc_external("Ksp/PCFieldSplitGetType"))
"""
function PCFieldSplitGetType(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCFieldSplitGetType(petsclib::$UnionPetscLib, pc::PC )
	type_ = Ref{PCCompositeType}()

    @chk ccall(
               (:PCFieldSplitGetType, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PCCompositeType}),
               pc, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	PCFieldSplitSetDMSplits(petsclib::PetscLibType,pc::PC, flg::PetscBool) 
Flags whether `DMCreateFieldDecomposition()` should be used to define the splits in a `PCFIELDSPLIT`, whenever possible.

Logically Collective

Input Parameters:
- `pc`  - the preconditioner context
- `flg` - boolean indicating whether to use field splits defined by the `DM`

Options Database Key:
- `-pc_fieldsplit_dm_splits <bool>` - use the field splits defined by the `DM`

Level: intermediate

-seealso: [](sec_block_matrices), `PC`, `PCFIELDSPLIT`, `PCFieldSplitGetDMSplits()`, `DMCreateFieldDecomposition()`, `PCFieldSplitSetFields()`, `PCFieldSplitSetIS()`

# External Links
$(_doc_external("Ksp/PCFieldSplitSetDMSplits"))
"""
function PCFieldSplitSetDMSplits(petsclib::PetscLibType, pc::PC, flg::PetscBool) end

@for_petsc function PCFieldSplitSetDMSplits(petsclib::$UnionPetscLib, pc::PC, flg::PetscBool )

    @chk ccall(
               (:PCFieldSplitSetDMSplits, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, flg,
              )


	return nothing
end 

"""
	flg::PetscBool = PCFieldSplitGetDMSplits(petsclib::PetscLibType,pc::PC) 
Returns flag indicating whether `DMCreateFieldDecomposition()` should be used to define the splits in a `PCFIELDSPLIT`, whenever possible.

Logically Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `flg` - boolean indicating whether to use field splits defined by the `DM`

Level: intermediate

-seealso: [](sec_block_matrices), `PC`, `PCFIELDSPLIT`, `PCFieldSplitSetDMSplits()`, `DMCreateFieldDecomposition()`, `PCFieldSplitSetFields()`, `PCFieldSplitSetIS()`

# External Links
$(_doc_external("Ksp/PCFieldSplitGetDMSplits"))
"""
function PCFieldSplitGetDMSplits(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCFieldSplitGetDMSplits(petsclib::$UnionPetscLib, pc::PC )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PCFieldSplitGetDMSplits, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PetscBool}),
               pc, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	flg::PetscBool = PCFieldSplitGetDetectSaddlePoint(petsclib::PetscLibType,pc::PC) 
Returns flag indicating whether `PCFIELDSPLIT` will attempt to automatically determine fields based on zero diagonal entries.

Logically Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `flg` - boolean indicating whether to detect fields or not

Level: intermediate

-seealso: [](sec_block_matrices), `PC`, `PCFIELDSPLIT`, `PCFieldSplitSetDetectSaddlePoint()`

# External Links
$(_doc_external("Ksp/PCFieldSplitGetDetectSaddlePoint"))
"""
function PCFieldSplitGetDetectSaddlePoint(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCFieldSplitGetDetectSaddlePoint(petsclib::$UnionPetscLib, pc::PC )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PCFieldSplitGetDetectSaddlePoint, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PetscBool}),
               pc, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	flg::PetscBool = PCFieldSplitSetDetectSaddlePoint(petsclib::PetscLibType,pc::PC) 
Sets flag indicating whether `PCFIELDSPLIT` will attempt to automatically determine fields based on zero diagonal entries.

Logically Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `flg` - boolean indicating whether to detect fields or not

Options Database Key:
- `-pc_fieldsplit_detect_saddle_point <bool>` - detect and use the saddle point

Level: intermediate

-seealso: [](sec_block_matrices), `PC`, `PCFIELDSPLIT`, `PCFieldSplitGetDetectSaddlePoint()`, `PCFieldSplitSetType()`, `PCFieldSplitSetSchurPre()`, `PC_FIELDSPLIT_SCHUR_PRE_SELF`

# External Links
$(_doc_external("Ksp/PCFieldSplitSetDetectSaddlePoint"))
"""
function PCFieldSplitSetDetectSaddlePoint(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCFieldSplitSetDetectSaddlePoint(petsclib::$UnionPetscLib, pc::PC )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PCFieldSplitSetDetectSaddlePoint, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	PCGASMSetTotalSubdomains(petsclib::PetscLibType,pc::PC, N::PetscInt) 
sets the total number of subdomains to use across the communicator for `PCGASM`

Logically Collective

Input Parameters:
- `pc` - the preconditioner
- `N`  - total number of subdomains

Level: beginner

-seealso: [](ch_ksp), `PCGASM`, `PCGASMSetSubdomains()`, `PCGASMSetOverlap()`
`PCGASMCreateSubdomains2D()`

# External Links
$(_doc_external("Ksp/PCGASMSetTotalSubdomains"))
"""
function PCGASMSetTotalSubdomains(petsclib::PetscLibType, pc::PC, N::PetscInt) end

@for_petsc function PCGASMSetTotalSubdomains(petsclib::$UnionPetscLib, pc::PC, N::$PetscInt )

    @chk ccall(
               (:PCGASMSetTotalSubdomains, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, N,
              )


	return nothing
end 

"""
	PCGASMSetSubdomains(petsclib::PetscLibType,pc::PC, n::PetscInt, iis::Vector{IS}, ois::Vector{IS}) 
Sets the subdomains for this MPI process
for the additive Schwarz preconditioner with multiple MPI processes per subdomain, `PCGASM`

Collective

Input Parameters:
- `pc`  - the preconditioner object
- `n`   - the number of subdomains for this MPI process
- `iis` - the index sets that define the inner subdomains (or `NULL` for PETSc to determine subdomains), the `iis` array is
copied so may be freed after this call.
- `ois` - the index sets that define the outer subdomains (or `NULL` to use the same as `iis`, or to construct by expanding `iis` by
the requested overlap), the `ois` array is copied so may be freed after this call.

Level: advanced

-seealso: [](ch_ksp), `PCGASM`, `PCGASMSetOverlap()`, `PCGASMGetSubKSP()`, `PCGASMDestroySubdomains()`,
`PCGASMCreateSubdomains2D()`, `PCGASMGetSubdomains()`

# External Links
$(_doc_external("Ksp/PCGASMSetSubdomains"))
"""
function PCGASMSetSubdomains(petsclib::PetscLibType, pc::PC, n::PetscInt, iis::Vector{IS}, ois::Vector{IS}) end

@for_petsc function PCGASMSetSubdomains(petsclib::$UnionPetscLib, pc::PC, n::$PetscInt, iis::Vector{IS}, ois::Vector{IS} )

    @chk ccall(
               (:PCGASMSetSubdomains, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, Ptr{IS}, Ptr{IS}),
               pc, n, iis, ois,
              )


	return nothing
end 

"""
	PCGASMSetOverlap(petsclib::PetscLibType,pc::PC, ovl::PetscInt) 
Sets the overlap between a pair of subdomains for the
additive Schwarz preconditioner `PCGASM`.  Either all or no MPI processes in the
pc communicator must call this routine.

Logically Collective

Input Parameters:
- `pc`  - the preconditioner context
- `ovl` - the amount of overlap between subdomains (ovl >= 0, default value = 0)

Options Database Key:
- `-pc_gasm_overlap <overlap>` - Sets overlap

Level: intermediate

-seealso: [](ch_ksp), `PCGASM`, `PCGASMSetSubdomains()`, `PCGASMGetSubKSP()`,
`PCGASMCreateSubdomains2D()`, `PCGASMGetSubdomains()`

# External Links
$(_doc_external("Ksp/PCGASMSetOverlap"))
"""
function PCGASMSetOverlap(petsclib::PetscLibType, pc::PC, ovl::PetscInt) end

@for_petsc function PCGASMSetOverlap(petsclib::$UnionPetscLib, pc::PC, ovl::$PetscInt )

    @chk ccall(
               (:PCGASMSetOverlap, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, ovl,
              )


	return nothing
end 

"""
	PCGASMSetType(petsclib::PetscLibType,pc::PC, type::PCGASMType) 
Sets the type of restriction and interpolation used
for local problems in the `PCGASM` additive Schwarz method.

Logically Collective

Input Parameters:
- `pc`   - the preconditioner context
- `type` - variant of `PCGASM`, one of
-seealso: [](ch_ksp), `PCGASM`, `PCGASMSetSubdomains()`, `PCGASMGetSubKSP()`,
`PCGASMCreateSubdomains2D()`, `PCASM`, `PCASMSetType()`

# External Links
$(_doc_external("Ksp/PCGASMSetType"))
"""
function PCGASMSetType(petsclib::PetscLibType, pc::PC, type::PCGASMType) end

@for_petsc function PCGASMSetType(petsclib::$UnionPetscLib, pc::PC, type::PCGASMType )

    @chk ccall(
               (:PCGASMSetType, $petsc_library),
               PetscErrorCode,
               (PC, PCGASMType),
               pc, type,
              )


	return nothing
end 

"""
	PCGASMSetSortIndices(petsclib::PetscLibType,pc::PC, doSort::PetscBool) 
Determines whether subdomain indices are sorted.

Logically Collective

Input Parameters:
- `pc`     - the preconditioner context
- `doSort` - sort the subdomain indices

Level: intermediate

-seealso: [](ch_ksp), `PCGASM`, `PCGASMSetSubdomains()`, `PCGASMGetSubKSP()`,
`PCGASMCreateSubdomains2D()`

# External Links
$(_doc_external("Ksp/PCGASMSetSortIndices"))
"""
function PCGASMSetSortIndices(petsclib::PetscLibType, pc::PC, doSort::PetscBool) end

@for_petsc function PCGASMSetSortIndices(petsclib::$UnionPetscLib, pc::PC, doSort::PetscBool )

    @chk ccall(
               (:PCGASMSetSortIndices, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, doSort,
              )


	return nothing
end 

"""
	n_loc::PetscInt,first_loc::PetscInt = PCGASMGetSubKSP(petsclib::PetscLibType,pc::PC, ksp::Vector{PetscKSP}) 
Gets the local `KSP` contexts for all subdomains on this MPI process.

Collective iff first_local is requested

Input Parameter:
- `pc` - the preconditioner context

Output Parameters:
- `n_local`     - the number of blocks on this MPI process or `NULL`
- `first_local` - the global number of the first block on this process or `NULL`, all processes must request or all must pass `NULL`
- `ksp`         - the array of `KSP` contexts

Level: advanced

-seealso: [](ch_ksp), `PCGASM`, `PCGASMSetSubdomains()`, `PCGASMSetOverlap()`,
`PCGASMCreateSubdomains2D()`,

# External Links
$(_doc_external("Ksp/PCGASMGetSubKSP"))
"""
function PCGASMGetSubKSP(petsclib::PetscLibType, pc::PC, ksp::Vector{PetscKSP}) end

@for_petsc function PCGASMGetSubKSP(petsclib::$UnionPetscLib, pc::PC, ksp::Vector{PetscKSP} )
	n_loc_ = Ref{$PetscInt}()
	first_loc_ = Ref{$PetscInt}()
	ksp_ = Ref(pointer(ksp))

    @chk ccall(
               (:PCGASMGetSubKSP, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{CKSP}}),
               pc, n_loc_, first_loc_, ksp_,
              )

	n_loc = n_loc_[]
	first_loc = first_loc_[]

	return n_loc,first_loc
end 

"""
	n::PetscInt,iis::Vector{IS} = PCGASMCreateSubdomains(petsclib::PetscLibType,A::PetscMat, N::PetscInt) 
Creates `n` index sets defining `n` nonoverlapping subdomains on this MPI process for the `PCGASM` additive
Schwarz preconditioner for a any problem based on its matrix.

Collective

Input Parameters:
- `A` - The global matrix operator
- `N` - the number of global subdomains requested

Output Parameters:
- `n`   - the number of subdomains created on this MPI process
- `iis` - the array of index sets defining the local inner subdomains (on which the correction is applied)

Level: advanced

-seealso: [](ch_ksp), `PCGASM`, `PCGASMSetSubdomains()`, `PCGASMDestroySubdomains()`

# External Links
$(_doc_external("Ksp/PCGASMCreateSubdomains"))
"""
function PCGASMCreateSubdomains(petsclib::PetscLibType, A::PetscMat, N::PetscInt) end

@for_petsc function PCGASMCreateSubdomains(petsclib::$UnionPetscLib, A::PetscMat, N::$PetscInt )
	n_ = Ref{$PetscInt}()
	iis_ = Ref{Ptr{IS}}()

    @chk ccall(
               (:PCGASMCreateSubdomains, $petsc_library),
               PetscErrorCode,
               (CMat, $PetscInt, Ptr{$PetscInt}, Ptr{Ptr{IS}}),
               A, N, n_, iis_,
              )

	n = n_[]
	iis = unsafe_wrap(Array, iis_[], VecGetLocalSize(petsclib, x); own = false)

	return n,iis
end 

"""
	PCGASMDestroySubdomains(petsclib::PetscLibType,n::PetscInt, iis::Vector{IS}, ois::Vector{IS}) 
Destroys the index sets created with
`PCGASMCreateSubdomains()` or `PCGASMCreateSubdomains2D()`. Should be
called after setting subdomains with `PCGASMSetSubdomains()`.

Collective

Input Parameters:
- `n`   - the number of index sets
- `iis` - the array of inner subdomains
- `ois` - the array of outer subdomains, can be `NULL`

Level: intermediate

-seealso: [](ch_ksp), `PCGASM`, `PCGASMCreateSubdomains()`, `PCGASMSetSubdomains()`

# External Links
$(_doc_external("Ksp/PCGASMDestroySubdomains"))
"""
function PCGASMDestroySubdomains(petsclib::PetscLibType, n::PetscInt, iis::Vector{IS}, ois::Vector{IS}) end

@for_petsc function PCGASMDestroySubdomains(petsclib::$UnionPetscLib, n::$PetscInt, iis::Vector{IS}, ois::Vector{IS} )
	iis_ = Ref(pointer(iis))
	ois_ = Ref(pointer(ois))

    @chk ccall(
               (:PCGASMDestroySubdomains, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{Ptr{IS}}, Ptr{Ptr{IS}}),
               n, iis_, ois_,
              )


	return nothing
end 

"""
	nsub::PetscInt,iis::Vector{IS},ois::Vector{IS} = PCGASMCreateSubdomains2D(petsclib::PetscLibType,pc::PC, M::PetscInt, N::PetscInt, Mdomains::PetscInt, Ndomains::PetscInt, dof::PetscInt, overlap::PetscInt) 
Creates the index sets for the `PCGASM` overlapping Schwarz
preconditioner for a two-dimensional problem on a regular grid.

Collective

Input Parameters:
- `pc`       - the preconditioner context
- `M`        - the global number of grid points in the x direction
- `N`        - the global number of grid points in the y direction
- `Mdomains` - the global number of subdomains in the x direction
- `Ndomains` - the global number of subdomains in the y direction
- `dof`      - degrees of freedom per node
- `overlap`  - overlap in mesh lines

Output Parameters:
- `nsub` - the number of local subdomains created
- `iis`  - array of index sets defining inner (nonoverlapping) subdomains
- `ois`  - array of index sets defining outer (overlapping, if overlap > 0) subdomains

Level: advanced

-seealso: [](ch_ksp), `PCGASM`, `PCGASMSetSubdomains()`, `PCGASMGetSubKSP()`, `PCGASMSetOverlap()`, `PCASMCreateSubdomains2D()`,
`PCGASMDestroySubdomains()`

# External Links
$(_doc_external("Ksp/PCGASMCreateSubdomains2D"))
"""
function PCGASMCreateSubdomains2D(petsclib::PetscLibType, pc::PC, M::PetscInt, N::PetscInt, Mdomains::PetscInt, Ndomains::PetscInt, dof::PetscInt, overlap::PetscInt) end

@for_petsc function PCGASMCreateSubdomains2D(petsclib::$UnionPetscLib, pc::PC, M::$PetscInt, N::$PetscInt, Mdomains::$PetscInt, Ndomains::$PetscInt, dof::$PetscInt, overlap::$PetscInt )
	nsub_ = Ref{$PetscInt}()
	iis_ = Ref{Ptr{IS}}()
	ois_ = Ref{Ptr{IS}}()

    @chk ccall(
               (:PCGASMCreateSubdomains2D, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{Ptr{IS}}, Ptr{Ptr{IS}}),
               pc, M, N, Mdomains, Ndomains, dof, overlap, nsub_, iis_, ois_,
              )

	nsub = nsub_[]
	iis = unsafe_wrap(Array, iis_[], VecGetLocalSize(petsclib, x); own = false)
	ois = unsafe_wrap(Array, ois_[], VecGetLocalSize(petsclib, x); own = false)

	return nsub,iis,ois
end 

"""
	n::PetscInt = PCGASMGetSubdomains(petsclib::PetscLibType,pc::PC, iis::Vector{IS}, ois::Vector{IS}) 
Gets the subdomains supported on this MPI process
for the `PCGASM` additive Schwarz preconditioner.

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameters:
- `n`   - the number of subdomains for this MPI process (default value = 1)
- `iis` - the index sets that define the inner subdomains (without overlap) supported on this process (can be `NULL`)
- `ois` - the index sets that define the outer subdomains (with overlap) supported on this process (can be `NULL`)

Level: advanced

-seealso: [](ch_ksp), `PCGASM`, `PCGASMSetOverlap()`, `PCGASMGetSubKSP()`, `PCGASMCreateSubdomains2D()`,
`PCGASMSetSubdomains()`, `PCGASMGetSubmatrices()`, `PCGASMDestroySubdomains()`

# External Links
$(_doc_external("Ksp/PCGASMGetSubdomains"))
"""
function PCGASMGetSubdomains(petsclib::PetscLibType, pc::PC, iis::Vector{IS}, ois::Vector{IS}) end

@for_petsc function PCGASMGetSubdomains(petsclib::$UnionPetscLib, pc::PC, iis::Vector{IS}, ois::Vector{IS} )
	n_ = Ref{$PetscInt}()
	iis_ = Ref(pointer(iis))
	ois_ = Ref(pointer(ois))

    @chk ccall(
               (:PCGASMGetSubdomains, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{$PetscInt}, Ptr{Ptr{IS}}, Ptr{Ptr{IS}}),
               pc, n_, iis_, ois_,
              )

	n = n_[]

	return n
end 

"""
	n::PetscInt = PCGASMGetSubmatrices(petsclib::PetscLibType,pc::PC, mat::Vector{PetscMat}) 
Gets the local submatrices (for this MPI process
only) for the `PCGASM` additive Schwarz preconditioner.

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameters:
- `n`   - the number of matrices for this MPI process (default value = 1)
- `mat` - the matrices

Level: advanced

-seealso: [](ch_ksp), `PCGASM`, `PCGASMSetOverlap()`, `PCGASMGetSubKSP()`,
`PCGASMCreateSubdomains2D()`, `PCGASMSetSubdomains()`, `PCGASMGetSubdomains()`

# External Links
$(_doc_external("Ksp/PCGASMGetSubmatrices"))
"""
function PCGASMGetSubmatrices(petsclib::PetscLibType, pc::PC, mat::Vector{PetscMat}) end

@for_petsc function PCGASMGetSubmatrices(petsclib::$UnionPetscLib, pc::PC, mat::Vector{PetscMat} )
	n_ = Ref{$PetscInt}()
	mat_ = Ref(pointer(mat))

    @chk ccall(
               (:PCGASMGetSubmatrices, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{$PetscInt}, Ptr{Ptr{CMat}}),
               pc, n_, mat_,
              )

	n = n_[]

	return n
end 

"""
	PCGASMSetUseDMSubdomains(petsclib::PetscLibType,pc::PC, flg::PetscBool) 
Indicates whether to use `DMCreateDomainDecomposition()` to define the subdomains, whenever possible for `PCGASM`

Logically Collective

Input Parameters:
- `pc`  - the preconditioner
- `flg` - boolean indicating whether to use subdomains defined by the `DM`

Options Database Key:
- `-pc_gasm_dm_subdomains`    - configure subdomains
- `-pc_gasm_overlap`          - set overlap
- `-pc_gasm_total_subdomains` - set number of subdomains

Level: intermediate

-seealso: [](ch_ksp), `PCGASM`, `PCGASMGetUseDMSubdomains()`, `PCGASMSetSubdomains()`, `PCGASMSetOverlap()`
`PCGASMCreateSubdomains2D()`

# External Links
$(_doc_external("Ksp/PCGASMSetUseDMSubdomains"))
"""
function PCGASMSetUseDMSubdomains(petsclib::PetscLibType, pc::PC, flg::PetscBool) end

@for_petsc function PCGASMSetUseDMSubdomains(petsclib::$UnionPetscLib, pc::PC, flg::PetscBool )

    @chk ccall(
               (:PCGASMSetUseDMSubdomains, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, flg,
              )


	return nothing
end 

"""
	flg::PetscBool = PCGASMGetUseDMSubdomains(petsclib::PetscLibType,pc::PC) 
Returns flag indicating whether to use `DMCreateDomainDecomposition()` to define the subdomains, whenever possible with `PCGASM`

Not Collective

Input Parameter:
- `pc` - the preconditioner

Output Parameter:
- `flg` - boolean indicating whether to use subdomains defined by the `DM`

Level: intermediate

-seealso: [](ch_ksp), `PCGASM`, `PCGASMSetUseDMSubdomains()`, `PCGASMSetOverlap()`
`PCGASMCreateSubdomains2D()`

# External Links
$(_doc_external("Ksp/PCGASMGetUseDMSubdomains"))
"""
function PCGASMGetUseDMSubdomains(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCGASMGetUseDMSubdomains(petsclib::$UnionPetscLib, pc::PC )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PCGASMGetUseDMSubdomains, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PetscBool}),
               pc, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	PCHYPRESetDiscreteGradient(petsclib::PetscLibType,pc::PC, G::PetscMat) 
Set the discrete gradient matrix for `PCHYPRE` type of AMS or ADS

Collective

Input Parameters:
- `pc` - the preconditioning context
- `G`  - the discrete gradient

Level: intermediate

-seealso: [](ch_ksp), `PCHYPRE`, `PCHYPRESetDiscreteCurl()`

# External Links
$(_doc_external("Ksp/PCHYPRESetDiscreteGradient"))
"""
function PCHYPRESetDiscreteGradient(petsclib::PetscLibType, pc::PC, G::PetscMat) end

@for_petsc function PCHYPRESetDiscreteGradient(petsclib::$UnionPetscLib, pc::PC, G::PetscMat )

    @chk ccall(
               (:PCHYPRESetDiscreteGradient, $petsc_library),
               PetscErrorCode,
               (PC, CMat),
               pc, G,
              )


	return nothing
end 

"""
	PCHYPRESetDiscreteCurl(petsclib::PetscLibType,pc::PC, C::PetscMat) 
Set the discrete curl matrix for `PCHYPRE` type of ADS

Collective

Input Parameters:
- `pc` - the preconditioning context
- `C`  - the discrete curl

Level: intermediate

-seealso: [](ch_ksp), `PCHYPRE`, `PCHYPRESetDiscreteGradient()`

# External Links
$(_doc_external("Ksp/PCHYPRESetDiscreteCurl"))
"""
function PCHYPRESetDiscreteCurl(petsclib::PetscLibType, pc::PC, C::PetscMat) end

@for_petsc function PCHYPRESetDiscreteCurl(petsclib::$UnionPetscLib, pc::PC, C::PetscMat )

    @chk ccall(
               (:PCHYPRESetDiscreteCurl, $petsc_library),
               PetscErrorCode,
               (PC, CMat),
               pc, C,
              )


	return nothing
end 

"""
	PCHYPRESetInterpolations(petsclib::PetscLibType,pc::PC, dim::PetscInt, RT_PiFull::PetscMat, RT_Pi::Vector{PetscMat}, ND_PiFull::PetscMat, ND_Pi::Vector{PetscMat}) 
Set the interpolation matrices for `PCHYPRE` type of AMS or ADS

Collective

Input Parameters:
- `pc`        - the preconditioning context
- `dim`       - the dimension of the problem, only used in AMS
- `RT_PiFull` - Raviart-Thomas interpolation matrix
- `RT_Pi`     - x/y/z component of Raviart-Thomas interpolation matrix
- `ND_PiFull` - Nedelec interpolation matrix
- `ND_Pi`     - x/y/z component of Nedelec interpolation matrix

Level: intermediate

-seealso: [](ch_ksp), `PCHYPRE`

# External Links
$(_doc_external("Ksp/PCHYPRESetInterpolations"))
"""
function PCHYPRESetInterpolations(petsclib::PetscLibType, pc::PC, dim::PetscInt, RT_PiFull::PetscMat, RT_Pi::Vector{PetscMat}, ND_PiFull::PetscMat, ND_Pi::Vector{PetscMat}) end

@for_petsc function PCHYPRESetInterpolations(petsclib::$UnionPetscLib, pc::PC, dim::$PetscInt, RT_PiFull::PetscMat, RT_Pi::Vector{PetscMat}, ND_PiFull::PetscMat, ND_Pi::Vector{PetscMat} )

    @chk ccall(
               (:PCHYPRESetInterpolations, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, CMat, Ptr{CMat}, CMat, Ptr{CMat}),
               pc, dim, RT_PiFull, RT_Pi, ND_PiFull, ND_Pi,
              )


	return nothing
end 

"""
	PCHYPRESetAlphaPoissonMatrix(petsclib::PetscLibType,pc::PC, A::PetscMat) 
Set the vector Poisson matrix for `PCHYPRE` of type AMS

Collective

Input Parameters:
- `pc` - the preconditioning context
- `A`  - the matrix

Level: intermediate

-seealso: [](ch_ksp), `PCHYPRE`, `PCHYPRESetDiscreteGradient()`, `PCHYPRESetDiscreteCurl()`, `PCHYPRESetBetaPoissonMatrix()`

# External Links
$(_doc_external("Ksp/PCHYPRESetAlphaPoissonMatrix"))
"""
function PCHYPRESetAlphaPoissonMatrix(petsclib::PetscLibType, pc::PC, A::PetscMat) end

@for_petsc function PCHYPRESetAlphaPoissonMatrix(petsclib::$UnionPetscLib, pc::PC, A::PetscMat )

    @chk ccall(
               (:PCHYPRESetAlphaPoissonMatrix, $petsc_library),
               PetscErrorCode,
               (PC, CMat),
               pc, A,
              )


	return nothing
end 

"""
	PCHYPRESetBetaPoissonMatrix(petsclib::PetscLibType,pc::PC, A::PetscMat) 
Set the Poisson matrix for `PCHYPRE` of type AMS

Collective

Input Parameters:
- `pc` - the preconditioning context
- `A`  - the matrix, or `NULL` to turn it off

Level: intermediate

-seealso: [](ch_ksp), `PCHYPRE`, `PCHYPRESetDiscreteGradient()`, `PCHYPRESetDiscreteCurl()`, `PCHYPRESetAlphaPoissonMatrix()`

# External Links
$(_doc_external("Ksp/PCHYPRESetBetaPoissonMatrix"))
"""
function PCHYPRESetBetaPoissonMatrix(petsclib::PetscLibType, pc::PC, A::PetscMat) end

@for_petsc function PCHYPRESetBetaPoissonMatrix(petsclib::$UnionPetscLib, pc::PC, A::PetscMat )

    @chk ccall(
               (:PCHYPRESetBetaPoissonMatrix, $petsc_library),
               PetscErrorCode,
               (PC, CMat),
               pc, A,
              )


	return nothing
end 

"""
	PCHYPRESetEdgeConstantVectors(petsclib::PetscLibType,pc::PC, ozz::PetscVec, zoz::PetscVec, zzo::PetscVec) 
Set the representation of the constant vector fields in the edge element basis for `PCHYPRE` of type AMS

Collective

Input Parameters:
- `pc`  - the preconditioning context
- `ozz` - vector representing (1,0,0) (or (1,0) in 2D)
- `zoz` - vector representing (0,1,0) (or (0,1) in 2D)
- `zzo` - vector representing (0,0,1) (use NULL in 2D)

Level: intermediate

-seealso: [](ch_ksp), `PCHYPRE`, `PCHYPRESetDiscreteGradient()`, `PCHYPRESetDiscreteCurl()`, `PCHYPRESetAlphaPoissonMatrix()`

# External Links
$(_doc_external("Ksp/PCHYPRESetEdgeConstantVectors"))
"""
function PCHYPRESetEdgeConstantVectors(petsclib::PetscLibType, pc::PC, ozz::PetscVec, zoz::PetscVec, zzo::PetscVec) end

@for_petsc function PCHYPRESetEdgeConstantVectors(petsclib::$UnionPetscLib, pc::PC, ozz::PetscVec, zoz::PetscVec, zzo::PetscVec )

    @chk ccall(
               (:PCHYPRESetEdgeConstantVectors, $petsc_library),
               PetscErrorCode,
               (PC, CVec, CVec, CVec),
               pc, ozz, zoz, zzo,
              )


	return nothing
end 

"""
	PCHYPREAMSSetInteriorNodes(petsclib::PetscLibType,pc::PC, interior::PetscVec) 
Set the list of interior nodes to a zero

Collective

Input Parameters:
- `pc`       - the preconditioning context
- `interior` - vector. node is interior if its entry in the array is 1.0.

Level: intermediate

-seealso: [](ch_ksp), `PCHYPRE`, `PCHYPRESetDiscreteGradient()`, `PCHYPRESetDiscreteCurl()`, `PCHYPRESetAlphaPoissonMatrix()`

# External Links
$(_doc_external("Ksp/PCHYPREAMSSetInteriorNodes"))
"""
function PCHYPREAMSSetInteriorNodes(petsclib::PetscLibType, pc::PC, interior::PetscVec) end

@for_petsc function PCHYPREAMSSetInteriorNodes(petsclib::$UnionPetscLib, pc::PC, interior::PetscVec )

    @chk ccall(
               (:PCHYPREAMSSetInteriorNodes, $petsc_library),
               PetscErrorCode,
               (PC, CVec),
               pc, interior,
              )


	return nothing
end 

"""
	PCHYPRESetType(petsclib::PetscLibType,pc::PC, name::String) 
Sets which hypre preconditioner you wish to use

Input Parameters:
- `pc`   - the preconditioner context
- `name` - either euclid, ilu, pilut, parasails, boomeramg, ams, or ads

Options Database Key:
- `pc_hypre_type` - One of euclid, ilu, pilut, parasails, boomeramg, ams, or ads

Level: intermediate

-seealso: [](ch_ksp), `PCCreate()`, `PCSetType()`, `PCType`, `PC`, `PCHYPRE`

# External Links
$(_doc_external("Ksp/PCHYPRESetType"))
"""
function PCHYPRESetType(petsclib::PetscLibType, pc::PC, name::String) end

@for_petsc function PCHYPRESetType(petsclib::$UnionPetscLib, pc::PC, name::String )

    @chk ccall(
               (:PCHYPRESetType, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{Cchar}),
               pc, name,
              )


	return nothing
end 

"""
	n_per_level::Vector{PetscInt} = PCHYPREGetCFMarkers(petsclib::PetscLibType,pc::PC, CFMarkers::Vector{PetscBT}) 
Gets CF marker arrays for all levels (except the finest level)

Logically Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameters:
- `n_per_level` - the number of nodes per level (size of `num_levels`)
- `CFMarkers`   - the Coarse/Fine Boolean arrays (size of `num_levels` - 1)

Level: advanced

-seealso: [](ch_ksp), `PC`, `PCMG`, `PCMGGetRestriction()`, `PCMGSetInterpolation()`, `PCMGGetRScale()`, `PCMGGetInterpolation()`, `PCGetInterpolations()`

# External Links
$(_doc_external("Ksp/PCHYPREGetCFMarkers"))
"""
function PCHYPREGetCFMarkers(petsclib::PetscLibType, pc::PC, CFMarkers::Vector{PetscBT}) end

@for_petsc function PCHYPREGetCFMarkers(petsclib::$UnionPetscLib, pc::PC, CFMarkers::Vector{PetscBT} )
	n_per_level_ = Ref{Ptr{$PetscInt}}()
	CFMarkers_ = Ref(pointer(CFMarkers))

    @chk ccall(
               (:PCHYPREGetCFMarkers, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{Ptr{$PetscInt}}, Ptr{Ptr{PetscBT}}),
               pc, n_per_level_, CFMarkers_,
              )

	n_per_level = unsafe_wrap(Array, n_per_level_[], VecGetLocalSize(petsclib, x); own = false)

	return n_per_level
end 

"""
	name::String = PCHYPREGetType(petsclib::PetscLibType,pc::PC) 
Gets which hypre preconditioner you are using

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `name` - either euclid, ilu, pilut, parasails, boomeramg, ams, or ads

Level: intermediate

-seealso: [](ch_ksp), `PCCreate()`, `PCHYPRESetType()`, `PCType`, `PC`, `PCHYPRE`

# External Links
$(_doc_external("Ksp/PCHYPREGetType"))
"""
function PCHYPREGetType(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCHYPREGetType(petsclib::$UnionPetscLib, pc::PC )
	name_ = Ref{Ptr{Cchar}}()

    @chk ccall(
               (:PCHYPREGetType, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{Ptr{Cchar}}),
               pc, name_,
              )

	name = unsafe_wrap(Array, name_[], VecGetLocalSize(petsclib, x); own = false)

	return name
end 

"""
	PCMGGalerkinSetMatProductAlgorithm(petsclib::PetscLibType,pc::PC, name::String) 
Set type of sparse matrix

Logically Collective

Input Parameters:
- `pc`   - the hypre context
- `name` - one of 'cusparse', 'hypre'

Options Database Key:
- `-pc_mg_galerkin_mat_product_algorithm <cusparse,hypre>` - Type of sparse matrix-matrix product to use in hypre

Level: intermediate

-seealso: [](ch_ksp), `PCHYPRE`, `PCMGGalerkinGetMatProductAlgorithm()`

# External Links
$(_doc_external("Ksp/PCMGGalerkinSetMatProductAlgorithm"))
"""
function PCMGGalerkinSetMatProductAlgorithm(petsclib::PetscLibType, pc::PC, name::String) end

@for_petsc function PCMGGalerkinSetMatProductAlgorithm(petsclib::$UnionPetscLib, pc::PC, name::String )

    @chk ccall(
               (:PCMGGalerkinSetMatProductAlgorithm, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{Cchar}),
               pc, name,
              )


	return nothing
end 

"""
	PCMGGalerkinGetMatProductAlgorithm(petsclib::PetscLibType,pc::PC, name::String) 
Get type of sparse matrix

Not Collective

Input Parameter:
- `pc` - the multigrid context

Output Parameter:
- `name` - one of 'cusparse', 'hypre'

Level: intermediate

-seealso: [](ch_ksp), `PCHYPRE`, `PCMGGalerkinSetMatProductAlgorithm()`

# External Links
$(_doc_external("Ksp/PCMGGalerkinGetMatProductAlgorithm"))
"""
function PCMGGalerkinGetMatProductAlgorithm(petsclib::PetscLibType, pc::PC, name::String) end

@for_petsc function PCMGGalerkinGetMatProductAlgorithm(petsclib::$UnionPetscLib, pc::PC, name::String )
	name_ = Ref(pointer(name))

    @chk ccall(
               (:PCMGGalerkinGetMatProductAlgorithm, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{Ptr{Cchar}}),
               pc, name_,
              )


	return nothing
end 

"""
	PCSORGetSymmetric(petsclib::PetscLibType,pc::PC, flag::MatSORType) 
Gets the form the SOR preconditioner is using;   backward, or forward relaxation.  The local variants perform SOR on
each processor.  By default forward relaxation is used.

Logically Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `flag` - one of the following
-seealso: [](ch_ksp), `PCSOR`, `PCEisenstatSetOmega()`, `PCSORSetIterations()`, `PCSORSetOmega()`, `PCSORSetSymmetric()`

# External Links
$(_doc_external("Ksp/PCSORGetSymmetric"))
"""
function PCSORGetSymmetric(petsclib::PetscLibType, pc::PC, flag::MatSORType) end

@for_petsc function PCSORGetSymmetric(petsclib::$UnionPetscLib, pc::PC, flag::MatSORType )

    @chk ccall(
               (:PCSORGetSymmetric, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{MatSORType}),
               pc, flag,
              )


	return nothing
end 

"""
	omega::PetscReal = PCSORGetOmega(petsclib::PetscLibType,pc::PC) 
Gets the SOR relaxation coefficient, omega
(where omega = 1.0 by default).

Logically Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `omega` - relaxation coefficient (0 < omega < 2).

Options Database Key:
- `-pc_sor_omega <omega>` - Sets omega

Level: intermediate

-seealso: [](ch_ksp), `PCSOR`, `PCSORSetSymmetric()`, `PCSORSetIterations()`, `PCEisenstatSetOmega()`, `PCSORSetOmega()`

# External Links
$(_doc_external("Ksp/PCSORGetOmega"))
"""
function PCSORGetOmega(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCSORGetOmega(petsclib::$UnionPetscLib, pc::PC )
	omega_ = Ref{$PetscReal}()

    @chk ccall(
               (:PCSORGetOmega, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{$PetscReal}),
               pc, omega_,
              )

	omega = omega_[]

	return omega
end 

"""
	its::PetscInt,lits::PetscInt = PCSORGetIterations(petsclib::PetscLibType,pc::PC) 
Gets the number of inner iterations to
be used by the SOR preconditioner. The default is 1.

Logically Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameters:
- `lits` - number of local iterations, smoothings over just variables on processor
- `its`  - number of parallel iterations to use; each parallel iteration has lits local iterations

Options Database Keys:
- `-pc_sor_its <its>`   - Sets number of iterations
- `-pc_sor_lits <lits>` - Sets number of local iterations

Level: intermediate

-seealso: [](ch_ksp), `PCSOR`, `PCSORSetOmega()`, `PCSORSetSymmetric()`, `PCSORSetIterations()`

# External Links
$(_doc_external("Ksp/PCSORGetIterations"))
"""
function PCSORGetIterations(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCSORGetIterations(petsclib::$UnionPetscLib, pc::PC )
	its_ = Ref{$PetscInt}()
	lits_ = Ref{$PetscInt}()

    @chk ccall(
               (:PCSORGetIterations, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{$PetscInt}, Ptr{$PetscInt}),
               pc, its_, lits_,
              )

	its = its_[]
	lits = lits_[]

	return its,lits
end 

"""
	PCSORSetSymmetric(petsclib::PetscLibType,pc::PC, flag::MatSORType) 
Sets the SOR preconditioner to use symmetric (SSOR),
backward, or forward relaxation.  The local variants perform SOR on
each processor.  By default forward relaxation is used.

Logically Collective

Input Parameters:
- `pc`   - the preconditioner context
- `flag` - one of the following
-seealso: [](ch_ksp), `PCSOR`, `PCEisenstatSetOmega()`, `PCSORSetIterations()`, `PCSORSetOmega()`

# External Links
$(_doc_external("Ksp/PCSORSetSymmetric"))
"""
function PCSORSetSymmetric(petsclib::PetscLibType, pc::PC, flag::MatSORType) end

@for_petsc function PCSORSetSymmetric(petsclib::$UnionPetscLib, pc::PC, flag::MatSORType )

    @chk ccall(
               (:PCSORSetSymmetric, $petsc_library),
               PetscErrorCode,
               (PC, MatSORType),
               pc, flag,
              )


	return nothing
end 

"""
	PCSORSetOmega(petsclib::PetscLibType,pc::PC, omega::PetscReal) 
Sets the SOR relaxation coefficient, omega
(where omega = 1.0 by default).

Logically Collective

Input Parameters:
- `pc`    - the preconditioner context
- `omega` - relaxation coefficient (0 < omega < 2).

Options Database Key:
- `-pc_sor_omega <omega>` - Sets omega

Level: intermediate

-seealso: [](ch_ksp), `PCSOR`, `PCSORSetSymmetric()`, `PCSORSetIterations()`, `PCEisenstatSetOmega()`, `MatSetOption()`

# External Links
$(_doc_external("Ksp/PCSORSetOmega"))
"""
function PCSORSetOmega(petsclib::PetscLibType, pc::PC, omega::PetscReal) end

@for_petsc function PCSORSetOmega(petsclib::$UnionPetscLib, pc::PC, omega::$PetscReal )

    @chk ccall(
               (:PCSORSetOmega, $petsc_library),
               PetscErrorCode,
               (PC, $PetscReal),
               pc, omega,
              )


	return nothing
end 

"""
	PCSORSetIterations(petsclib::PetscLibType,pc::PC, its::PetscInt, lits::PetscInt) 
Sets the number of inner iterations to
be used by the SOR preconditioner. The default is 1.

Logically Collective

Input Parameters:
- `pc`   - the preconditioner context
- `lits` - number of local iterations, smoothings over just variables on processor
- `its`  - number of parallel iterations to use; each parallel iteration has lits local iterations

Options Database Keys:
- `-pc_sor_its <its>`   - Sets number of iterations
- `-pc_sor_lits <lits>` - Sets number of local iterations

Level: intermediate

-seealso: [](ch_ksp), `PCSOR`, `PCSORSetOmega()`, `PCSORSetSymmetric()`

# External Links
$(_doc_external("Ksp/PCSORSetIterations"))
"""
function PCSORSetIterations(petsclib::PetscLibType, pc::PC, its::PetscInt, lits::PetscInt) end

@for_petsc function PCSORSetIterations(petsclib::$UnionPetscLib, pc::PC, its::$PetscInt, lits::$PetscInt )

    @chk ccall(
               (:PCSORSetIterations, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, $PetscInt),
               pc, its, lits,
              )


	return nothing
end 

"""
	PCDeflationSetInitOnly(petsclib::PetscLibType,pc::PC, flg::PetscBool) 
Do only initialization step.
Sets initial guess to the solution on the deflation space but does not apply
the deflation preconditioner. The additional preconditioner is still applied.

Logically Collective

Input Parameters:
- `pc`  - the preconditioner context
- `flg` - default `PETSC_FALSE`

Options Database Key:
- `-pc_deflation_init_only <false>` - if true computes only the special guess

Level: intermediate

-seealso: [](ch_ksp), `PCDEFLATION`

# External Links
$(_doc_external("Ksp/PCDeflationSetInitOnly"))
"""
function PCDeflationSetInitOnly(petsclib::PetscLibType, pc::PC, flg::PetscBool) end

@for_petsc function PCDeflationSetInitOnly(petsclib::$UnionPetscLib, pc::PC, flg::PetscBool )

    @chk ccall(
               (:PCDeflationSetInitOnly, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, flg,
              )


	return nothing
end 

"""
	PCDeflationSetLevels(petsclib::PetscLibType,pc::PC, max::PetscInt) 
Set the maximum level of deflation nesting.

Logically Collective

Input Parameters:
- `pc`  - the preconditioner context
- `max` - maximum deflation level

Options Database Key:
- `-pc_deflation_max_lvl <0>` - maximum number of levels for multilevel deflation

Level: intermediate

-seealso: [](ch_ksp), `PCDeflationSetSpaceToCompute()`, `PCDeflationSetSpace()`, `PCDEFLATION`

# External Links
$(_doc_external("Ksp/PCDeflationSetLevels"))
"""
function PCDeflationSetLevels(petsclib::PetscLibType, pc::PC, max::PetscInt) end

@for_petsc function PCDeflationSetLevels(petsclib::$UnionPetscLib, pc::PC, max::$PetscInt )

    @chk ccall(
               (:PCDeflationSetLevels, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, max,
              )


	return nothing
end 

"""
	PCDeflationSetReductionFactor(petsclib::PetscLibType,pc::PC, red::PetscInt) 
Set reduction factor for the `PCDEFLATION`

Logically Collective

Input Parameters:
- `pc`  - the preconditioner context
- `red` - reduction factor (or `PETSC_DETERMINE`)

Options Database Key:
- `-pc_deflation_reduction_factor <-1>` - reduction factor on bottom level coarse problem for `PCDEFLATION`

-seealso: [](ch_ksp), `PCTELESCOPE`, `PCDEFLATION`, `PCDeflationSetLevels()`

# External Links
$(_doc_external("Ksp/PCDeflationSetReductionFactor"))
"""
function PCDeflationSetReductionFactor(petsclib::PetscLibType, pc::PC, red::PetscInt) end

@for_petsc function PCDeflationSetReductionFactor(petsclib::$UnionPetscLib, pc::PC, red::$PetscInt )

    @chk ccall(
               (:PCDeflationSetReductionFactor, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, red,
              )


	return nothing
end 

"""
	PCDeflationSetCorrectionFactor(petsclib::PetscLibType,pc::PC, fact::PetscScalar) 
Set coarse problem correction factor.
The preconditioner becomes P*M^{-1} + fact*Q.

Logically Collective

Input Parameters:
- `pc`   - the preconditioner context
- `fact` - correction factor

Options Database Keys:
- `-pc_deflation_correction        <false>` - if true apply coarse problem correction
- `-pc_deflation_correction_factor <1.0>`   - sets coarse problem correction factor

-seealso: [](ch_ksp), `PCDEFLATION`, `PCDeflationSetLevels()`, `PCDeflationSetReductionFactor()`

# External Links
$(_doc_external("Ksp/PCDeflationSetCorrectionFactor"))
"""
function PCDeflationSetCorrectionFactor(petsclib::PetscLibType, pc::PC, fact::PetscScalar) end

@for_petsc function PCDeflationSetCorrectionFactor(petsclib::$UnionPetscLib, pc::PC, fact::$PetscScalar )

    @chk ccall(
               (:PCDeflationSetCorrectionFactor, $petsc_library),
               PetscErrorCode,
               (PC, $PetscScalar),
               pc, fact,
              )


	return nothing
end 

"""
	PCDeflationSetSpaceToCompute(petsclib::PetscLibType,pc::PC, type::PCDeflationSpaceType, size::PetscInt) 
Set deflation space type and size to compute.

Logically Collective

Input Parameters:
- `pc`   - the preconditioner context
- `type` - deflation space type to compute (or `PETSC_IGNORE`)
- `size` - size of the space to compute (or `PETSC_DEFAULT`)

Options Database Keys:
- `-pc_deflation_compute_space      <haar>` - compute `PCDeflationSpaceType` deflation space
- `-pc_deflation_compute_space_size <1>`    - size of the deflation space

-seealso: [](ch_ksp), `PCDeflationSetLevels()`, `PCDEFLATION`

# External Links
$(_doc_external("Ksp/PCDeflationSetSpaceToCompute"))
"""
function PCDeflationSetSpaceToCompute(petsclib::PetscLibType, pc::PC, type::PCDeflationSpaceType, size::PetscInt) end

@for_petsc function PCDeflationSetSpaceToCompute(petsclib::$UnionPetscLib, pc::PC, type::PCDeflationSpaceType, size::$PetscInt )

    @chk ccall(
               (:PCDeflationSetSpaceToCompute, $petsc_library),
               PetscErrorCode,
               (PC, PCDeflationSpaceType, $PetscInt),
               pc, type, size,
              )


	return nothing
end 

"""
	PCDeflationSetSpace(petsclib::PetscLibType,pc::PC, W::PetscMat, transpose::PetscBool) 
Set the deflation space matrix (or its (Hermitian) transpose).

Logically Collective

Input Parameters:
- `pc`        - the preconditioner context
- `W`         - deflation matrix
- `transpose` - indicates that W is an explicit transpose of the deflation matrix

Level: intermediate

-seealso: [](ch_ksp), `PCDeflationSetLevels()`, `PCDEFLATION`, `PCDeflationSetProjectionNullSpaceMat()`

# External Links
$(_doc_external("Ksp/PCDeflationSetSpace"))
"""
function PCDeflationSetSpace(petsclib::PetscLibType, pc::PC, W::PetscMat, transpose::PetscBool) end

@for_petsc function PCDeflationSetSpace(petsclib::$UnionPetscLib, pc::PC, W::PetscMat, transpose::PetscBool )

    @chk ccall(
               (:PCDeflationSetSpace, $petsc_library),
               PetscErrorCode,
               (PC, CMat, PetscBool),
               pc, W, transpose,
              )


	return nothing
end 

"""
	PCDeflationSetProjectionNullSpaceMat(petsclib::PetscLibType,pc::PC, mat::PetscMat) 
Set the projection null space matrix (W'*A).

Collective

Input Parameters:
- `pc`  - preconditioner context
- `mat` - projection null space matrix

Level: developer

-seealso: [](ch_ksp), `PCDEFLATION`, `PCDeflationSetSpace()`

# External Links
$(_doc_external("Ksp/PCDeflationSetProjectionNullSpaceMat"))
"""
function PCDeflationSetProjectionNullSpaceMat(petsclib::PetscLibType, pc::PC, mat::PetscMat) end

@for_petsc function PCDeflationSetProjectionNullSpaceMat(petsclib::$UnionPetscLib, pc::PC, mat::PetscMat )

    @chk ccall(
               (:PCDeflationSetProjectionNullSpaceMat, $petsc_library),
               PetscErrorCode,
               (PC, CMat),
               pc, mat,
              )


	return nothing
end 

"""
	PCDeflationSetCoarseMat(petsclib::PetscLibType,pc::PC, mat::PetscMat) 
Set the coarse problem `Mat`.

Collective

Input Parameters:
- `pc`  - preconditioner context
- `mat` - coarse problem mat

Level: developer

-seealso: [](ch_ksp), `PCDEFLATION`, `PCDeflationGetCoarseKSP()`

# External Links
$(_doc_external("Ksp/PCDeflationSetCoarseMat"))
"""
function PCDeflationSetCoarseMat(petsclib::PetscLibType, pc::PC, mat::PetscMat) end

@for_petsc function PCDeflationSetCoarseMat(petsclib::$UnionPetscLib, pc::PC, mat::PetscMat )

    @chk ccall(
               (:PCDeflationSetCoarseMat, $petsc_library),
               PetscErrorCode,
               (PC, CMat),
               pc, mat,
              )


	return nothing
end 

"""
	PCDeflationGetCoarseKSP(petsclib::PetscLibType,pc::PC, ksp::PetscKSP) 
Returns the coarse problem `KSP`.

Not Collective

Input Parameter:
- `pc` - preconditioner context

Output Parameter:
- `ksp` - coarse problem `KSP` context

Level: advanced

-seealso: [](ch_ksp), `PCDEFLATION`, `PCDeflationSetCoarseMat()`

# External Links
$(_doc_external("Ksp/PCDeflationGetCoarseKSP"))
"""
function PCDeflationGetCoarseKSP(petsclib::PetscLibType, pc::PC, ksp::PetscKSP) end

@for_petsc function PCDeflationGetCoarseKSP(petsclib::$UnionPetscLib, pc::PC, ksp::PetscKSP )
	ksp_ = Ref(ksp.ptr)

    @chk ccall(
               (:PCDeflationGetCoarseKSP, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{CKSP}),
               pc, ksp_,
              )

	ksp.ptr = C_NULL

	return nothing
end 

"""
	PCDeflationGetPC(petsclib::PetscLibType,pc::PC, apc::PC) 
Returns the additional preconditioner M^{

Not Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `apc` - additional preconditioner

Level: advanced

-seealso: [](ch_ksp), `PCDEFLATION`, `PCDeflationGetCoarseKSP()`

# External Links
$(_doc_external("Ksp/PCDeflationGetPC"))
"""
function PCDeflationGetPC(petsclib::PetscLibType, pc::PC, apc::PC) end

@for_petsc function PCDeflationGetPC(petsclib::$UnionPetscLib, pc::PC, apc::PC )

    @chk ccall(
               (:PCDeflationGetPC, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PC}),
               pc, apc,
              )


	return nothing
end 

"""
	PCGAMGSetNSmooths(petsclib::PetscLibType,pc::PC, n::PetscInt) 
Set number of smoothing steps (1 is typical) used to construct the prolongation operator

Logically Collective

Input Parameters:
- `pc` - the preconditioner context
- `n`  - the number of smooths

Options Database Key:
- `-pc_gamg_agg_nsmooths <nsmooth, default=1>` - the flag

Level: intermediate

-seealso: [the Users Manual section on PCGAMG](sec_amg), [the Users Manual section on PCMG](sec_mg), [](ch_ksp), `PCMG`, `PCGAMG`

# External Links
$(_doc_external("Ksp/PCGAMGSetNSmooths"))
"""
function PCGAMGSetNSmooths(petsclib::PetscLibType, pc::PC, n::PetscInt) end

@for_petsc function PCGAMGSetNSmooths(petsclib::$UnionPetscLib, pc::PC, n::$PetscInt )

    @chk ccall(
               (:PCGAMGSetNSmooths, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, n,
              )


	return nothing
end 

"""
	PCGAMGSetAggressiveLevels(petsclib::PetscLibType,pc::PC, n::PetscInt) 
Use aggressive coarsening on first n levels

Logically Collective

Input Parameters:
- `pc` - the preconditioner context
- `n`  - 0, 1 or more

Options Database Key:
- `-pc_gamg_aggressive_coarsening <n,default = 1>` - the flag

Level: intermediate

-seealso: [the Users Manual section on PCGAMG](sec_amg), [the Users Manual section on PCMG](sec_mg), [](ch_ksp), `PCGAMG`, `PCGAMGSetThreshold()`, `PCGAMGMISkSetAggressive()`, `PCGAMGSetAggressiveSquareGraph()`, `PCGAMGMISkSetMinDegreeOrdering()`, `PCGAMGSetLowMemoryFilter()`

# External Links
$(_doc_external("Ksp/PCGAMGSetAggressiveLevels"))
"""
function PCGAMGSetAggressiveLevels(petsclib::PetscLibType, pc::PC, n::PetscInt) end

@for_petsc function PCGAMGSetAggressiveLevels(petsclib::$UnionPetscLib, pc::PC, n::$PetscInt )

    @chk ccall(
               (:PCGAMGSetAggressiveLevels, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, n,
              )


	return nothing
end 

"""
	PCGAMGMISkSetAggressive(petsclib::PetscLibType,pc::PC, n::PetscInt) 
Number (k) distance in MIS coarsening (>2 is 'aggressive')

Logically Collective

Input Parameters:
- `pc` - the preconditioner context
- `n`  - 1 or more (default = 2)

Options Database Key:
- `-pc_gamg_aggressive_mis_k <n,default=2>` - the flag

Level: intermediate

-seealso: [the Users Manual section on PCGAMG](sec_amg), [the Users Manual section on PCMG](sec_mg), [](ch_ksp), `PCGAMG`, `PCGAMGSetThreshold()`, `PCGAMGSetAggressiveLevels()`, `PCGAMGSetAggressiveSquareGraph()`, `PCGAMGMISkSetMinDegreeOrdering()`, `PCGAMGSetLowMemoryFilter()`

# External Links
$(_doc_external("Ksp/PCGAMGMISkSetAggressive"))
"""
function PCGAMGMISkSetAggressive(petsclib::PetscLibType, pc::PC, n::PetscInt) end

@for_petsc function PCGAMGMISkSetAggressive(petsclib::$UnionPetscLib, pc::PC, n::$PetscInt )

    @chk ccall(
               (:PCGAMGMISkSetAggressive, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, n,
              )


	return nothing
end 

"""
	PCGAMGSetAggressiveSquareGraph(petsclib::PetscLibType,pc::PC, b::PetscBool) 
Use graph square, A^T A, for aggressive coarsening. Coarsening is slower than the alternative (MIS

Logically Collective

Input Parameters:
- `pc` - the preconditioner context
- `b`  - default true

Options Database Key:
- `-pc_gamg_aggressive_square_graph <bool,default=true>` - the flag

Level: intermediate

-seealso: [the Users Manual section on PCGAMG](sec_amg), [the Users Manual section on PCMG](sec_mg), [](ch_ksp), `PCGAMG`, `PCGAMGSetThreshold()`, `PCGAMGSetAggressiveLevels()`, `PCGAMGMISkSetAggressive()`, `PCGAMGMISkSetMinDegreeOrdering()`, `PCGAMGSetLowMemoryFilter()`

# External Links
$(_doc_external("Ksp/PCGAMGSetAggressiveSquareGraph"))
"""
function PCGAMGSetAggressiveSquareGraph(petsclib::PetscLibType, pc::PC, b::PetscBool) end

@for_petsc function PCGAMGSetAggressiveSquareGraph(petsclib::$UnionPetscLib, pc::PC, b::PetscBool )

    @chk ccall(
               (:PCGAMGSetAggressiveSquareGraph, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, b,
              )


	return nothing
end 

"""
	PCGAMGMISkSetMinDegreeOrdering(petsclib::PetscLibType,pc::PC, b::PetscBool) 
Use minimum degree ordering in greedy MIS algorithm

Logically Collective

Input Parameters:
- `pc` - the preconditioner context
- `b`  - default false

Options Database Key:
- `-pc_gamg_mis_k_minimum_degree_ordering <bool,default=false>` - the flag

Level: intermediate

-seealso: [the Users Manual section on PCGAMG](sec_amg), [the Users Manual section on PCMG](sec_mg), [](ch_ksp), `PCGAMG`, `PCGAMGSetThreshold()`, `PCGAMGSetAggressiveLevels()`, `PCGAMGMISkSetAggressive()`, `PCGAMGSetAggressiveSquareGraph()`, `PCGAMGSetLowMemoryFilter()`

# External Links
$(_doc_external("Ksp/PCGAMGMISkSetMinDegreeOrdering"))
"""
function PCGAMGMISkSetMinDegreeOrdering(petsclib::PetscLibType, pc::PC, b::PetscBool) end

@for_petsc function PCGAMGMISkSetMinDegreeOrdering(petsclib::$UnionPetscLib, pc::PC, b::PetscBool )

    @chk ccall(
               (:PCGAMGMISkSetMinDegreeOrdering, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, b,
              )


	return nothing
end 

"""
	PCGAMGSetLowMemoryFilter(petsclib::PetscLibType,pc::PC, b::PetscBool) 
Use low memory graph/matrix filter

Logically Collective

Input Parameters:
- `pc` - the preconditioner context
- `b`  - default false

Options Database Key:
- `-pc_gamg_low_memory_threshold_filter <bool,default=false>` - the flag

Level: intermediate

-seealso: [the Users Manual section on PCGAMG](sec_amg), [the Users Manual section on PCMG](sec_mg), `PCGAMG`, `PCGAMGSetThreshold()`, `PCGAMGSetAggressiveLevels()`,
`PCGAMGMISkSetAggressive()`, `PCGAMGSetAggressiveSquareGraph()`, `PCGAMGMISkSetMinDegreeOrdering()`

# External Links
$(_doc_external("Ksp/PCGAMGSetLowMemoryFilter"))
"""
function PCGAMGSetLowMemoryFilter(petsclib::PetscLibType, pc::PC, b::PetscBool) end

@for_petsc function PCGAMGSetLowMemoryFilter(petsclib::$UnionPetscLib, pc::PC, b::PetscBool )

    @chk ccall(
               (:PCGAMGSetLowMemoryFilter, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, b,
              )


	return nothing
end 

"""
	PCGAMGSetGraphSymmetrize(petsclib::PetscLibType,pc::PC, b::PetscBool) 
Symmetrize graph used for coarsening. Defaults to true, but if matrix has symmetric attribute, then not needed since the graph is already known to be symmetric

Logically Collective

Input Parameters:
- `pc` - the preconditioner context
- `b`  - default true

Options Database Key:
- `-pc_gamg_graph_symmetrize <bool,default=true>` - the flag

Level: intermediate

-seealso: [the Users Manual section on PCGAMG](sec_amg), [the Users Manual section on PCMG](sec_mg), `PCGAMG`, `PCGAMGSetThreshold()`, `PCGAMGSetAggressiveLevels()`, `MatCreateGraph()`,
`PCGAMGMISkSetAggressive()`, `PCGAMGSetAggressiveSquareGraph()`, `PCGAMGMISkSetMinDegreeOrdering()`

# External Links
$(_doc_external("Ksp/PCGAMGSetGraphSymmetrize"))
"""
function PCGAMGSetGraphSymmetrize(petsclib::PetscLibType, pc::PC, b::PetscBool) end

@for_petsc function PCGAMGSetGraphSymmetrize(petsclib::$UnionPetscLib, pc::PC, b::PetscBool )

    @chk ccall(
               (:PCGAMGSetGraphSymmetrize, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, b,
              )


	return nothing
end 

"""
	PCGAMGSetProcEqLim(petsclib::PetscLibType,pc::PC, n::PetscInt) 
Set number of equations to aim for per process on the coarse grids via processor reduction in `PCGAMG`

Logically Collective

Input Parameters:
- `pc` - the preconditioner context
- `n`  - the number of equations

Options Database Key:
- `-pc_gamg_process_eq_limit <limit>` - set the limit

Level: intermediate

-seealso: [the Users Manual section on PCGAMG](sec_amg), [the Users Manual section on PCMG](sec_mg), [](ch_ksp), `PCGAMG`, `PCGAMGSetCoarseEqLim()`, `PCGAMGSetRankReductionFactors()`, `PCGAMGSetRepartition()`

# External Links
$(_doc_external("Ksp/PCGAMGSetProcEqLim"))
"""
function PCGAMGSetProcEqLim(petsclib::PetscLibType, pc::PC, n::PetscInt) end

@for_petsc function PCGAMGSetProcEqLim(petsclib::$UnionPetscLib, pc::PC, n::$PetscInt )

    @chk ccall(
               (:PCGAMGSetProcEqLim, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, n,
              )


	return nothing
end 

"""
	PCGAMGSetCoarseEqLim(petsclib::PetscLibType,pc::PC, n::PetscInt) 
Set maximum number of equations on the coarsest grid of `PCGAMG`

Collective

Input Parameters:
- `pc` - the preconditioner context
- `n`  - maximum number of equations to aim for

Options Database Key:
- `-pc_gamg_coarse_eq_limit <limit>` - set the limit

Level: intermediate

-seealso: [the Users Manual section on PCGAMG](sec_amg), [the Users Manual section on PCMG](sec_mg), [](ch_ksp), `PCGAMG`, `PCGAMGSetProcEqLim()`, `PCGAMGSetRankReductionFactors()`, `PCGAMGSetRepartition()`,
`PCGAMGSetParallelCoarseGridSolve()`

# External Links
$(_doc_external("Ksp/PCGAMGSetCoarseEqLim"))
"""
function PCGAMGSetCoarseEqLim(petsclib::PetscLibType, pc::PC, n::PetscInt) end

@for_petsc function PCGAMGSetCoarseEqLim(petsclib::$UnionPetscLib, pc::PC, n::$PetscInt )

    @chk ccall(
               (:PCGAMGSetCoarseEqLim, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, n,
              )


	return nothing
end 

"""
	PCGAMGSetRepartition(petsclib::PetscLibType,pc::PC, n::PetscBool) 
Repartition the degrees of freedom across the processors on the coarser grids when reducing the number of MPI ranks to use

Collective

Input Parameters:
- `pc` - the preconditioner context
- `n`  - `PETSC_TRUE` or `PETSC_FALSE`

Options Database Key:
- `-pc_gamg_repartition <true,false>` - turn on the repartitioning

Level: intermediate

-seealso: [the Users Manual section on PCGAMG](sec_amg), [the Users Manual section on PCMG](sec_mg), [](ch_ksp), `PCGAMG`, `PCGAMGSetProcEqLim()`, `PCGAMGSetRankReductionFactors()`

# External Links
$(_doc_external("Ksp/PCGAMGSetRepartition"))
"""
function PCGAMGSetRepartition(petsclib::PetscLibType, pc::PC, n::PetscBool) end

@for_petsc function PCGAMGSetRepartition(petsclib::$UnionPetscLib, pc::PC, n::PetscBool )

    @chk ccall(
               (:PCGAMGSetRepartition, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, n,
              )


	return nothing
end 

"""
	PCGAMGSetUseSAEstEig(petsclib::PetscLibType,pc::PC, b::PetscBool) 
Use the eigen estimate from smoothed aggregation for the Chebyshev smoother during the solution process

Collective

Input Parameters:
- `pc` - the preconditioner context
- `b`  - flag

Options Database Key:
- `-pc_gamg_use_sa_esteig <true,false>` - use the eigen estimate

Level: advanced

-seealso: [the Users Manual section on PCGAMG](sec_amg), [the Users Manual section on PCMG](sec_mg), [](ch_ksp), `PCGAMG`, `KSPChebyshevSetEigenvalues()`, `KSPChebyshevEstEigSet()`, `PCGAMGSetRecomputeEstEig()`

# External Links
$(_doc_external("Ksp/PCGAMGSetUseSAEstEig"))
"""
function PCGAMGSetUseSAEstEig(petsclib::PetscLibType, pc::PC, b::PetscBool) end

@for_petsc function PCGAMGSetUseSAEstEig(petsclib::$UnionPetscLib, pc::PC, b::PetscBool )

    @chk ccall(
               (:PCGAMGSetUseSAEstEig, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, b,
              )


	return nothing
end 

"""
	PCGAMGSetRecomputeEstEig(petsclib::PetscLibType,pc::PC, b::PetscBool) 
Set flag for Chebyshev smoothers to recompute the eigen estimates when a new matrix is used

Collective

Input Parameters:
- `pc` - the preconditioner context
- `b`  - flag, default is `PETSC_TRUE`

Options Database Key:
- `-pc_gamg_recompute_esteig <true>` - use the eigen estimate

Level: advanced

-seealso: [the Users Manual section on PCGAMG](sec_amg), [the Users Manual section on PCMG](sec_mg), [](ch_ksp), `PCGAMG`, `KSPChebyshevSetEigenvalues()`, `KSPChebyshevEstEigSet()`

# External Links
$(_doc_external("Ksp/PCGAMGSetRecomputeEstEig"))
"""
function PCGAMGSetRecomputeEstEig(petsclib::PetscLibType, pc::PC, b::PetscBool) end

@for_petsc function PCGAMGSetRecomputeEstEig(petsclib::$UnionPetscLib, pc::PC, b::PetscBool )

    @chk ccall(
               (:PCGAMGSetRecomputeEstEig, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, b,
              )


	return nothing
end 

"""
	PCGAMGSetEigenvalues(petsclib::PetscLibType,pc::PC, emax::PetscReal, emin::PetscReal) 
Set WHAT eigenvalues WHY?

Collective

Input Parameters:
- `pc`   - the preconditioner context
- `emax` - max eigenvalue
- `emin` - min eigenvalue

Options Database Key:
- `-pc_gamg_eigenvalues <emin,emax>` - estimates of the eigenvalues

Level: intermediate

-seealso: [the Users Manual section on PCGAMG](sec_amg), [the Users Manual section on PCMG](sec_mg), [](ch_ksp), `PCGAMG`, `PCGAMGSetUseSAEstEig()`

# External Links
$(_doc_external("Ksp/PCGAMGSetEigenvalues"))
"""
function PCGAMGSetEigenvalues(petsclib::PetscLibType, pc::PC, emax::PetscReal, emin::PetscReal) end

@for_petsc function PCGAMGSetEigenvalues(petsclib::$UnionPetscLib, pc::PC, emax::$PetscReal, emin::$PetscReal )

    @chk ccall(
               (:PCGAMGSetEigenvalues, $petsc_library),
               PetscErrorCode,
               (PC, $PetscReal, $PetscReal),
               pc, emax, emin,
              )


	return nothing
end 

"""
	PCGAMGSetReuseInterpolation(petsclib::PetscLibType,pc::PC, n::PetscBool) 
Reuse prolongation when rebuilding a `PCGAMG` algebraic multigrid preconditioner

Collective

Input Parameters:
- `pc` - the preconditioner context
- `n`  - `PETSC_TRUE` or `PETSC_FALSE`

Options Database Key:
- `-pc_gamg_reuse_interpolation <true,false>` - reuse the previous interpolation

Level: intermediate

-seealso: [the Users Manual section on PCGAMG](sec_amg), [the Users Manual section on PCMG](sec_mg), [](ch_ksp), `PCGAMG`

# External Links
$(_doc_external("Ksp/PCGAMGSetReuseInterpolation"))
"""
function PCGAMGSetReuseInterpolation(petsclib::PetscLibType, pc::PC, n::PetscBool) end

@for_petsc function PCGAMGSetReuseInterpolation(petsclib::$UnionPetscLib, pc::PC, n::PetscBool )

    @chk ccall(
               (:PCGAMGSetReuseInterpolation, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, n,
              )


	return nothing
end 

"""
	PCGAMGASMSetUseAggs(petsclib::PetscLibType,pc::PC, flg::PetscBool) 
Have the `PCGAMG` smoother on each level use `PCASM` where the aggregates defined by the coarsening process are
the subdomains for the additive Schwarz preconditioner used as the smoother

Collective

Input Parameters:
- `pc`  - the preconditioner context
- `flg` - `PETSC_TRUE` to use aggregates, `PETSC_FALSE` to not

Options Database Key:
- `-pc_gamg_asm_use_agg <true,false>` - use aggregates to define the additive Schwarz subdomains

Level: intermediate

-seealso: [the Users Manual section on PCGAMG](sec_amg), [the Users Manual section on PCMG](sec_mg), [](ch_ksp), `PCGAMG`, `PCASM`, `PCSetType`

# External Links
$(_doc_external("Ksp/PCGAMGASMSetUseAggs"))
"""
function PCGAMGASMSetUseAggs(petsclib::PetscLibType, pc::PC, flg::PetscBool) end

@for_petsc function PCGAMGASMSetUseAggs(petsclib::$UnionPetscLib, pc::PC, flg::PetscBool )

    @chk ccall(
               (:PCGAMGASMSetUseAggs, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, flg,
              )


	return nothing
end 

"""
	PCGAMGSetParallelCoarseGridSolve(petsclib::PetscLibType,pc::PC, flg::PetscBool) 
allow a parallel coarse grid solver

Collective

Input Parameters:
- `pc`  - the preconditioner context
- `flg` - `PETSC_TRUE` to not force coarse grid onto one processor

Options Database Key:
- `-pc_gamg_parallel_coarse_grid_solver` - use a parallel coarse grid solver

Level: intermediate

-seealso: [the Users Manual section on PCGAMG](sec_amg), [the Users Manual section on PCMG](sec_mg), [](ch_ksp), `PCGAMG`, `PCGAMGSetCoarseGridLayoutType()`, `PCGAMGSetCpuPinCoarseGrids()`, `PCGAMGSetRankReductionFactors()`

# External Links
$(_doc_external("Ksp/PCGAMGSetParallelCoarseGridSolve"))
"""
function PCGAMGSetParallelCoarseGridSolve(petsclib::PetscLibType, pc::PC, flg::PetscBool) end

@for_petsc function PCGAMGSetParallelCoarseGridSolve(petsclib::$UnionPetscLib, pc::PC, flg::PetscBool )

    @chk ccall(
               (:PCGAMGSetParallelCoarseGridSolve, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, flg,
              )


	return nothing
end 

"""
	PCGAMGSetCpuPinCoarseGrids(petsclib::PetscLibType,pc::PC, flg::PetscBool) 
pin the coarse grids created in `PCGAMG` to run only on the CPU since the problems may be too small to run efficiently on the GPUs

Collective

Input Parameters:
- `pc`  - the preconditioner context
- `flg` - `PETSC_TRUE` to pin coarse grids to the CPU

Options Database Key:
- `-pc_gamg_cpu_pin_coarse_grids` - pin the coarse grids to the CPU

Level: intermediate

-seealso: [the Users Manual section on PCGAMG](sec_amg), [the Users Manual section on PCMG](sec_mg), [](ch_ksp), `PCGAMG`, `PCGAMGSetCoarseGridLayoutType()`, `PCGAMGSetParallelCoarseGridSolve()`

# External Links
$(_doc_external("Ksp/PCGAMGSetCpuPinCoarseGrids"))
"""
function PCGAMGSetCpuPinCoarseGrids(petsclib::PetscLibType, pc::PC, flg::PetscBool) end

@for_petsc function PCGAMGSetCpuPinCoarseGrids(petsclib::$UnionPetscLib, pc::PC, flg::PetscBool )

    @chk ccall(
               (:PCGAMGSetCpuPinCoarseGrids, $petsc_library),
               PetscErrorCode,
               (PC, PetscBool),
               pc, flg,
              )


	return nothing
end 

"""
	PCGAMGSetCoarseGridLayoutType(petsclib::PetscLibType,pc::PC, flg::PCGAMGLayoutType) 
place coarse grids on processors with natural order (compact type)

Collective

Input Parameters:
- `pc`  - the preconditioner context
- `flg` - `PCGAMGLayoutType` type, either `PCGAMG_LAYOUT_COMPACT` or `PCGAMG_LAYOUT_SPREAD`

Options Database Key:
- `-pc_gamg_coarse_grid_layout_type` - place the coarse grids with natural ordering

Level: intermediate

-seealso: [the Users Manual section on PCGAMG](sec_amg), [the Users Manual section on PCMG](sec_mg), [](ch_ksp), `PCGAMG`, `PCGAMGSetParallelCoarseGridSolve()`, `PCGAMGSetCpuPinCoarseGrids()`, `PCGAMGLayoutType`, `PCGAMG_LAYOUT_COMPACT`, `PCGAMG_LAYOUT_SPREAD`

# External Links
$(_doc_external("Ksp/PCGAMGSetCoarseGridLayoutType"))
"""
function PCGAMGSetCoarseGridLayoutType(petsclib::PetscLibType, pc::PC, flg::PCGAMGLayoutType) end

@for_petsc function PCGAMGSetCoarseGridLayoutType(petsclib::$UnionPetscLib, pc::PC, flg::PCGAMGLayoutType )

    @chk ccall(
               (:PCGAMGSetCoarseGridLayoutType, $petsc_library),
               PetscErrorCode,
               (PC, PCGAMGLayoutType),
               pc, flg,
              )


	return nothing
end 

"""
	PCGAMGSetNlevels(petsclib::PetscLibType,pc::PC, n::PetscInt) 
Sets the maximum number of levels `PCGAMG` will use

Collective

Input Parameters:
- `pc` - the preconditioner
- `n`  - the maximum number of levels to use

Options Database Key:
- `-pc_mg_levels <n>` - set the maximum number of levels to allow

Level: intermediate

-seealso: [the Users Manual section on PCGAMG](sec_amg), [the Users Manual section on PCMG](sec_mg), [](ch_ksp), `PCGAMG`

# External Links
$(_doc_external("Ksp/PCGAMGSetNlevels"))
"""
function PCGAMGSetNlevels(petsclib::PetscLibType, pc::PC, n::PetscInt) end

@for_petsc function PCGAMGSetNlevels(petsclib::$UnionPetscLib, pc::PC, n::$PetscInt )

    @chk ccall(
               (:PCGAMGSetNlevels, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, n,
              )


	return nothing
end 

"""
	PCGAMGASMSetHEM(petsclib::PetscLibType,pc::PC, n::PetscInt) 
Sets the number of HEM matching passed

Collective

Input Parameters:
- `pc` - the preconditioner
- `n`  - number of HEM matching passed to construct ASM subdomains

Options Database Key:
- `-pc_gamg_asm_hem <n>` - set the number of HEM matching passed

Level: intermediate

-seealso: [the Users Manual section on PCGAMG](sec_amg), [the Users Manual section on PCMG](sec_mg), [](ch_ksp), `PCGAMG`

# External Links
$(_doc_external("Ksp/PCGAMGASMSetHEM"))
"""
function PCGAMGASMSetHEM(petsclib::PetscLibType, pc::PC, n::PetscInt) end

@for_petsc function PCGAMGASMSetHEM(petsclib::$UnionPetscLib, pc::PC, n::$PetscInt )

    @chk ccall(
               (:PCGAMGASMSetHEM, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt),
               pc, n,
              )


	return nothing
end 

"""
	PCGAMGSetThreshold(petsclib::PetscLibType,pc::PC, v::Vector{PetscReal}, n::PetscInt) 
Relative threshold to use for dropping edges in aggregation graph

Not Collective

Input Parameters:
- `pc` - the preconditioner context
- `v`  - array of threshold values for finest n levels; 0.0 means keep all nonzero entries in the graph; negative means keep even zero entries in the graph
- `n`  - number of threshold values provided in array

Options Database Key:
- `-pc_gamg_threshold <threshold>` - the threshold to drop edges

Level: intermediate

-seealso: [the Users Manual section on PCGAMG](sec_amg), [the Users Manual section on PCMG](sec_mg), [](ch_ksp), `PCGAMG`, `PCGAMGSetAggressiveLevels()`, `PCGAMGMISkSetAggressive()`, `PCGAMGSetMinDegreeOrderingMISk()`, `PCGAMGSetThresholdScale()`

# External Links
$(_doc_external("Ksp/PCGAMGSetThreshold"))
"""
function PCGAMGSetThreshold(petsclib::PetscLibType, pc::PC, v::Vector{PetscReal}, n::PetscInt) end

@for_petsc function PCGAMGSetThreshold(petsclib::$UnionPetscLib, pc::PC, v::Vector{$PetscReal}, n::$PetscInt )

    @chk ccall(
               (:PCGAMGSetThreshold, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{$PetscReal}, $PetscInt),
               pc, v, n,
              )


	return nothing
end 

"""
	PCGAMGSetRankReductionFactors(petsclib::PetscLibType,pc::PC, v::Vector{PetscInt}, n::PetscInt) 
Set a manual schedule for MPI rank reduction on coarse grids

Collective

Input Parameters:
- `pc` - the preconditioner context
- `v`  - array of reduction factors. 0 for first value forces a reduction to one process/device on first level in CUDA
- `n`  - number of values provided in array

Options Database Key:
- `-pc_gamg_rank_reduction_factors <factors>` - provide the schedule

Level: intermediate

-seealso: [the Users Manual section on PCGAMG](sec_amg), [the Users Manual section on PCMG](sec_mg), [](ch_ksp), `PCGAMG`, `PCGAMGSetProcEqLim()`, `PCGAMGSetCoarseEqLim()`, `PCGAMGSetParallelCoarseGridSolve()`

# External Links
$(_doc_external("Ksp/PCGAMGSetRankReductionFactors"))
"""
function PCGAMGSetRankReductionFactors(petsclib::PetscLibType, pc::PC, v::Vector{PetscInt}, n::PetscInt) end

@for_petsc function PCGAMGSetRankReductionFactors(petsclib::$UnionPetscLib, pc::PC, v::Vector{$PetscInt}, n::$PetscInt )

    @chk ccall(
               (:PCGAMGSetRankReductionFactors, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{$PetscInt}, $PetscInt),
               pc, v, n,
              )


	return nothing
end 

"""
	PCGAMGSetThresholdScale(petsclib::PetscLibType,pc::PC, v::PetscReal) 
Relative threshold reduction at each level

Not Collective

Input Parameters:
- `pc` - the preconditioner context
- `v`  - the threshold value reduction, usually < 1.0

Options Database Key:
- `-pc_gamg_threshold_scale <v>` - set the relative threshold reduction on each level

Level: advanced

-seealso: [the Users Manual section on PCGAMG](sec_amg), [the Users Manual section on PCMG](sec_mg), [](ch_ksp), `PCGAMG`, `PCGAMGSetThreshold()`

# External Links
$(_doc_external("Ksp/PCGAMGSetThresholdScale"))
"""
function PCGAMGSetThresholdScale(petsclib::PetscLibType, pc::PC, v::PetscReal) end

@for_petsc function PCGAMGSetThresholdScale(petsclib::$UnionPetscLib, pc::PC, v::$PetscReal )

    @chk ccall(
               (:PCGAMGSetThresholdScale, $petsc_library),
               PetscErrorCode,
               (PC, $PetscReal),
               pc, v,
              )


	return nothing
end 

"""
	PCGAMGSetType(petsclib::PetscLibType,pc::PC, type::PCGAMGType) 
Set the type of algorithm `PCGAMG` should use

Collective

Input Parameters:
- `pc`   - the preconditioner context
- `type` - `PCGAMGAGG`, `PCGAMGGEO`, or `PCGAMGCLASSICAL`

Options Database Key:
- `-pc_gamg_type <agg,geo,classical>` - type of algebraic multigrid to apply - only agg is supported

Level: intermediate

-seealso: [the Users Manual section on PCGAMG](sec_amg), [the Users Manual section on PCMG](sec_mg), [](ch_ksp), `PCGAMGGetType()`, `PCGAMG`, `PCGAMGType`

# External Links
$(_doc_external("Ksp/PCGAMGSetType"))
"""
function PCGAMGSetType(petsclib::PetscLibType, pc::PC, type::PCGAMGType) end

@for_petsc function PCGAMGSetType(petsclib::$UnionPetscLib, pc::PC, type::PCGAMGType )

    @chk ccall(
               (:PCGAMGSetType, $petsc_library),
               PetscErrorCode,
               (PC, PCGAMGType),
               pc, type,
              )


	return nothing
end 

"""
	type::PCGAMGType = PCGAMGGetType(petsclib::PetscLibType,pc::PC) 
Get the type of algorithm `PCGAMG` will use

Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `type` - the type of algorithm used

Level: intermediate

-seealso: [the Users Manual section on PCGAMG](sec_amg), [the Users Manual section on PCMG](sec_mg), [](ch_ksp), `PCGAMG`, `PCGAMGSetType()`, `PCGAMGType`

# External Links
$(_doc_external("Ksp/PCGAMGGetType"))
"""
function PCGAMGGetType(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCGAMGGetType(petsclib::$UnionPetscLib, pc::PC )
	type_ = Ref{PCGAMGType}()

    @chk ccall(
               (:PCGAMGGetType, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PCGAMGType}),
               pc, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	PCGAMGSetInjectionIndex(petsclib::PetscLibType,pc::PC, n::PetscInt, idx::Vector{PetscInt}) 
Array of subset of variables per vertex to inject into coarse grid space

Logically Collective

Input Parameters:
- `pc`  - the coarsen context
- `n`   - number of indices
- `idx` - array of indices

Options Database Key:
- `-pc_gamg_injection_index` - array of subset of variables per vertex to use for injection coarse grid space

Level: intermediate

-seealso: [the Users Manual section on PCGAMG](sec_amg), [the Users Manual section on PCMG](sec_mg), `PCGAMG`

# External Links
$(_doc_external("Ksp/PCGAMGSetInjectionIndex"))
"""
function PCGAMGSetInjectionIndex(petsclib::PetscLibType, pc::PC, n::PetscInt, idx::Vector{PetscInt}) end

@for_petsc function PCGAMGSetInjectionIndex(petsclib::$UnionPetscLib, pc::PC, n::$PetscInt, idx::Vector{$PetscInt} )

    @chk ccall(
               (:PCGAMGSetInjectionIndex, $petsc_library),
               PetscErrorCode,
               (PC, $PetscInt, Ptr{$PetscInt}),
               pc, n, idx,
              )


	return nothing
end 

"""
	PCGAMGInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `PCGAMG` package. It is called
from `PCInitializePackage()`.

Level: developer

-seealso: [](ch_ksp), `PetscInitialize()`

# External Links
$(_doc_external("Ksp/PCGAMGInitializePackage"))
"""
function PCGAMGInitializePackage(petsclib::PetscLibType) end

@for_petsc function PCGAMGInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PCGAMGInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PCGAMGFinalizePackage(petsclib::PetscLibType) 
This function frees everything from the `PCGAMG` package. It is
called from `PetscFinalize()` automatically.

Level: developer

-seealso: [](ch_ksp), `PetscFinalize()`

# External Links
$(_doc_external("Ksp/PCGAMGFinalizePackage"))
"""
function PCGAMGFinalizePackage(petsclib::PetscLibType) end

@for_petsc function PCGAMGFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PCGAMGFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PCGAMGRegister(petsclib::PetscLibType,type::PCGAMGType, create::external) 
Register a `PCGAMG` implementation.

Input Parameters:
- `type`   - string that will be used as the name of the `PCGAMG` type.
- `create` - function for creating the gamg context.

Level: developer

-seealso: [](ch_ksp), `PCGAMGType`, `PCGAMG`, `PCGAMGSetType()`

# External Links
$(_doc_external("Ksp/PCGAMGRegister"))
"""
function PCGAMGRegister(petsclib::PetscLibType, type::PCGAMGType, create::external) end

@for_petsc function PCGAMGRegister(petsclib::$UnionPetscLib, type::PCGAMGType, create::external )

    @chk ccall(
               (:PCGAMGRegister, $petsc_library),
               PetscErrorCode,
               (PCGAMGType, external),
               type, create,
              )


	return nothing
end 

"""
	G::PetscMat = PCGAMGCreateGraph(petsclib::PetscLibType,pc::PC, A::PetscMat) 
Creates a graph that is used by the `PCGAMGType` in the coarsening process

Input Parameters:
- `pc` - the `PCGAMG`
- `A`  - the matrix, for any level

Output Parameter:
- `G` - the graph

Level: advanced

-seealso: [](ch_ksp), `PCGAMGType`, `PCGAMG`, `PCGAMGSetType()`

# External Links
$(_doc_external("Ksp/PCGAMGCreateGraph"))
"""
function PCGAMGCreateGraph(petsclib::PetscLibType, pc::PC, A::PetscMat) end

@for_petsc function PCGAMGCreateGraph(petsclib::$UnionPetscLib, pc::PC, A::PetscMat )
	G_ = Ref{CMat}()

    @chk ccall(
               (:PCGAMGCreateGraph, $petsc_library),
               PetscErrorCode,
               (PC, CMat, Ptr{CMat}),
               pc, A, G_,
              )

	G = PetscMat(G_[], petsclib)

	return G
end 

"""
	PCGAMGClassicalSetType(petsclib::PetscLibType,pc::PC, type::PCGAMGClassicalType) 
Sets the type of classical interpolation to use with `PCGAMG`

Collective

Input Parameters:
- `pc`   - the preconditioner context
- `type` - the interpolation to use, see `PCGAMGClassicalType()`

Options Database Key:
- `-pc_gamg_classical_type <direct,standard>` - set type of classical AMG prolongation

Level: intermediate

-seealso: [](ch_ksp), `PCGAMG`, `PCGAMGClassicalType`, `PCGAMGClassicalGetType()`

# External Links
$(_doc_external("Ksp/PCGAMGClassicalSetType"))
"""
function PCGAMGClassicalSetType(petsclib::PetscLibType, pc::PC, type::PCGAMGClassicalType) end

@for_petsc function PCGAMGClassicalSetType(petsclib::$UnionPetscLib, pc::PC, type::PCGAMGClassicalType )

    @chk ccall(
               (:PCGAMGClassicalSetType, $petsc_library),
               PetscErrorCode,
               (PC, PCGAMGClassicalType),
               pc, type,
              )


	return nothing
end 

"""
	type::PCGAMGClassicalType = PCGAMGClassicalGetType(petsclib::PetscLibType,pc::PC) 
Gets the type of classical interpolation to use with `PCGAMG`

Collective

Input Parameter:
- `pc` - the preconditioner context

Output Parameter:
- `type` - the type used, see `PCGAMGClassicalType()`

Level: intermediate

-seealso: [](ch_ksp), `PCGAMG`, `PCGAMGClassicalType`, `PCGAMGClassicalSetType()`

# External Links
$(_doc_external("Ksp/PCGAMGClassicalGetType"))
"""
function PCGAMGClassicalGetType(petsclib::PetscLibType, pc::PC) end

@for_petsc function PCGAMGClassicalGetType(petsclib::$UnionPetscLib, pc::PC )
	type_ = Ref{PCGAMGClassicalType}()

    @chk ccall(
               (:PCGAMGClassicalGetType, $petsc_library),
               PetscErrorCode,
               (PC, Ptr{PCGAMGClassicalType}),
               pc, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	PCFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the `PC` package. It is
called from `PetscFinalize()`.

Level: developer

-seealso: [](ch_ksp), `PetscFinalize()`, `PCInitializePackage()`

# External Links
$(_doc_external("Ksp/PCFinalizePackage"))
"""
function PCFinalizePackage(petsclib::PetscLibType) end

@for_petsc function PCFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PCFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PCInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `PC` package. It is called
from `PetscDLLibraryRegister_petscksp()` when using dynamic libraries, and on the first call to `PCCreate()`
when using shared static libraries.

Level: developer

-seealso: [](ch_ksp), `PetscInitialize()`, `PCFinalizePackage()`

# External Links
$(_doc_external("Ksp/PCInitializePackage"))
"""
function PCInitializePackage(petsclib::PetscLibType) end

@for_petsc function PCInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PCInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PCMPIServerAddressesDestroy(petsclib::PetscLibType,ctx::Cvoid) 

# External Links
$(_doc_external("Sys/PCMPIServerAddressesDestroy"))
"""
function PCMPIServerAddressesDestroy(petsclib::PetscLibType, ctx::Cvoid) end

@for_petsc function PCMPIServerAddressesDestroy(petsclib::$UnionPetscLib, ctx::Cvoid )

    @chk ccall(
               (:PCMPIServerAddressesDestroy, $petsc_library),
               PetscErrorCode,
               (Cvoid,),
               ctx,
              )


	return nothing
end 

