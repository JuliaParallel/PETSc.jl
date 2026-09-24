"""
	KSPAppendOptionsPrefix(petsclib::PetscLibType, ksp::AbstractKSP, prefix::String) 
Appends to the prefix used for searching for all
`KSP` options in the database.

Logically Collective

Input Parameters:
- `ksp`    - the Krylov context
- `prefix` - the prefix string to prepend to all `KSP` option requests

Level: intermediate

See also: `KSP`, `KSPSetOptionsPrefix()`, `KSPGetOptionsPrefix()`, `KSPSetFromOptions()`

# External Links
$(_doc_external("KSP/KSPAppendOptionsPrefix"))
"""
function KSPAppendOptionsPrefix(petsclib::PetscLibType, ksp::AbstractKSP, prefix::String)
    error("KSPAppendOptionsPrefix: no generated method for these argument types")
end

@for_petsc function KSPAppendOptionsPrefix(petsclib::$UnionPetscLib, ksp::AbstractKSP, prefix::String )

    @chk ccall(
               (:KSPAppendOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Cchar}),
               ksp, prefix,
              )


	return nothing
end 

"""
	KSPBCGSLSetEll(petsclib::PetscLibType, ksp::AbstractKSP, ell::PetscInt) 
Sets the number of search directions to use in the `KSPBCGSL` Krylov solver

Logically Collective

Input Parameters:
- `ksp` - iterative context, `KSP`, of type `KSPBCGSL`
- `ell` - number of search directions to use

Options Database Key:
- `-ksp_bcgsl_ell ell` - Number of Krylov search directions

Level: intermediate

See also: `KSPBCGSLSetUsePseudoinverse()`, `KSP`, `KSPBCGSL`, `KSPBCGSLSetPol()`, `KSPBCGSLSetXRes()`

# External Links
$(_doc_external("KSP/KSPBCGSLSetEll"))
"""
function KSPBCGSLSetEll(petsclib::PetscLibType, ksp::AbstractKSP, ell::Integer)
    error("KSPBCGSLSetEll: no generated method for these argument types")
end

@for_petsc function KSPBCGSLSetEll(petsclib::$UnionPetscLib, ksp::AbstractKSP, ell::$PetscInt )

    @chk ccall(
               (:KSPBCGSLSetEll, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, ell,
              )


	return nothing
end 

"""
	KSPBCGSLSetPol(petsclib::PetscLibType, ksp::AbstractKSP, uMROR::PetscBool) 
Sets the type of polynomial part that will
be used in the `KSPCBGSL` `KSPSolve()`

Logically Collective

Input Parameters:
- `ksp`   - iterative context of type `KSPCBGSL`
- `uMROR` - set to `PETSC_TRUE` when the polynomial is a convex combination of an MR and an OR step.

Options Database Keys:
- `-ksp_bcgsl_cxpoly` - use enhanced polynomial
- `-ksp_bcgsl_mrpoly` - use standard polynomial

Level: intermediate

See also: `KSP`, `KSPBCGSL`, `KSPCreate()`, `KSPSetType()`, `KSPCBGSL`, `KSPBCGSLSetUsePseudoinverse()`, `KSPBCGSLSetEll()`, `KSPBCGSLSetXRes()`

# External Links
$(_doc_external("KSP/KSPBCGSLSetPol"))
"""
function KSPBCGSLSetPol(petsclib::PetscLibType, ksp::AbstractKSP, uMROR::PetscBool)
    error("KSPBCGSLSetPol: no generated method for these argument types")
end

@for_petsc function KSPBCGSLSetPol(petsclib::$UnionPetscLib, ksp::AbstractKSP, uMROR::PetscBool )

    @chk ccall(
               (:KSPBCGSLSetPol, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, uMROR,
              )


	return nothing
end 

"""
	KSPBCGSLSetUsePseudoinverse(petsclib::PetscLibType, ksp::AbstractKSP, use_pinv::PetscBool) 
Use pseudoinverse (via SVD) to solve polynomial part of the update in `KSPCBGSL` solver

Logically Collective

Input Parameters:
- `ksp`      - iterative context of type `KSPCBGSL`
- `use_pinv` - set to `PETSC_TRUE` when using pseudoinverse

Options Database Key:
- `-ksp_bcgsl_pinv (true|false)` - use pseudoinverse

Level: intermediate

See also: `KSPBCGSLSetEll()`, `KSP`, `KSPCBGSL`, `KSPBCGSLSetPol()`, `KSPBCGSLSetXRes()`

# External Links
$(_doc_external("KSP/KSPBCGSLSetUsePseudoinverse"))
"""
function KSPBCGSLSetUsePseudoinverse(petsclib::PetscLibType, ksp::AbstractKSP, use_pinv::PetscBool)
    error("KSPBCGSLSetUsePseudoinverse: no generated method for these argument types")
end

@for_petsc function KSPBCGSLSetUsePseudoinverse(petsclib::$UnionPetscLib, ksp::AbstractKSP, use_pinv::PetscBool )

    @chk ccall(
               (:KSPBCGSLSetUsePseudoinverse, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, use_pinv,
              )


	return nothing
end 

"""
	KSPBCGSLSetXRes(petsclib::PetscLibType, ksp::AbstractKSP, delta::PetscReal) 
Sets the parameter governing when
exact residuals will be used instead of computed residuals for `KSPCBGSL`.

Logically Collective

Input Parameters:
- `ksp`   - iterative context of type `KSPBCGSL`
- `delta` - computed residuals are used alone when delta is not positive

Options Database Key:
- `-ksp_bcgsl_xres delta` - Threshold used to decide when to refresh computed residuals

Level: intermediate

See also: `KSPBCGSLSetEll()`, `KSPBCGSLSetPol()`, `KSP`, `KSPCBGSL`, `KSPBCGSLSetUsePseudoinverse()`

# External Links
$(_doc_external("KSP/KSPBCGSLSetXRes"))
"""
function KSPBCGSLSetXRes(petsclib::PetscLibType, ksp::AbstractKSP, delta::Real)
    error("KSPBCGSLSetXRes: no generated method for these argument types")
end

@for_petsc function KSPBCGSLSetXRes(petsclib::$UnionPetscLib, ksp::AbstractKSP, delta::$PetscReal )

    @chk ccall(
               (:KSPBCGSLSetXRes, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscReal),
               ksp, delta,
              )


	return nothing
end 

"""
	KSPBuildResidual(petsclib::PetscLibType, ksp::AbstractKSP, t::AbstractPetscVec, v::AbstractPetscVec, M_V::AbstractPetscVec) 
Builds the residual in a vector provided.

Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameters:
- `t` - work vector.  If not provided then one is generated.
- `v` - optional location to stash residual.  If `v` is not provided, then a location is generated.
- `V` - the residual

Level: advanced

See also: `KSP`, `KSPBuildSolution()`

# External Links
$(_doc_external("KSP/KSPBuildResidual"))
"""
function KSPBuildResidual(petsclib::PetscLibType, ksp::AbstractKSP, t::AbstractPetscVec, v::AbstractPetscVec, M_V::AbstractPetscVec)
    error("KSPBuildResidual: no generated method for these argument types")
end

@for_petsc function KSPBuildResidual(petsclib::$UnionPetscLib, ksp::AbstractKSP, t::AbstractPetscVec, v::AbstractPetscVec, M_V::AbstractPetscVec )
	M_V_ = Ref(M_V.ptr)

    @chk ccall(
               (:KSPBuildResidual, $petsc_library),
               PetscErrorCode,
               (CKSP, CVec, CVec, Ptr{CVec}),
               ksp, t, v, M_V_,
              )

	M_V.ptr = M_V_[]

	return nothing
end 

"""
	KSPBuildResidualDefault(petsclib::PetscLibType, ksp::AbstractKSP, t::AbstractPetscVec, v::AbstractPetscVec, M_V::AbstractPetscVec) 
Default code to compute the residual.

Collecive on ksp

Input Parameters:
- `ksp` - iterative context
- `t`   - pointer to temporary vector
- `v`   - pointer to user vector

Output Parameter:
- `V` - pointer to a vector containing the residual

Level: advanced

See also: `KSP`, `KSPBuildSolutionDefault()`

# External Links
$(_doc_external("KSP/KSPBuildResidualDefault"))
"""
function KSPBuildResidualDefault(petsclib::PetscLibType, ksp::AbstractKSP, t::AbstractPetscVec, v::AbstractPetscVec, M_V::AbstractPetscVec)
    error("KSPBuildResidualDefault: no generated method for these argument types")
end

@for_petsc function KSPBuildResidualDefault(petsclib::$UnionPetscLib, ksp::AbstractKSP, t::AbstractPetscVec, v::AbstractPetscVec, M_V::AbstractPetscVec )
	M_V_ = Ref(M_V.ptr)

    @chk ccall(
               (:KSPBuildResidualDefault, $petsc_library),
               PetscErrorCode,
               (CKSP, CVec, CVec, Ptr{CVec}),
               ksp, t, v, M_V_,
              )

	M_V.ptr = M_V_[]

	return nothing
end 

"""
	KSPBuildSolution(petsclib::PetscLibType, ksp::AbstractKSP, v::AbstractPetscVec, M_V::AbstractPetscVec) 
Builds the approximate solution in a vector provided.

Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameter:
Provide exactly one of
- `v` - location to stash solution, optional, otherwise pass `NULL`
- `V` - the solution is returned in this location. This vector is created internally. This vector should NOT be destroyed by the user with `VecDestroy()`.

Level: developer

See also: `KSPGetSolution()`, `KSPBuildResidual()`, `KSP`

# External Links
$(_doc_external("KSP/KSPBuildSolution"))
"""
function KSPBuildSolution(petsclib::PetscLibType, ksp::AbstractKSP, v::AbstractPetscVec, M_V::AbstractPetscVec)
    error("KSPBuildSolution: no generated method for these argument types")
end

@for_petsc function KSPBuildSolution(petsclib::$UnionPetscLib, ksp::AbstractKSP, v::AbstractPetscVec, M_V::AbstractPetscVec )
	M_V_ = Ref(M_V.ptr)

    @chk ccall(
               (:KSPBuildSolution, $petsc_library),
               PetscErrorCode,
               (CKSP, CVec, Ptr{CVec}),
               ksp, v, M_V_,
              )

	M_V.ptr = M_V_[]

	return nothing
end 

"""
	KSPBuildSolutionDefault(petsclib::PetscLibType, ksp::AbstractKSP, v::AbstractPetscVec, M_V::AbstractPetscVec) 
Default code to build/move the solution.

Collective

Input Parameters:
- `ksp` - iterative context
- `v`   - pointer to the user's vector

Output Parameter:
- `V` - pointer to a vector containing the solution

Level: advanced

See also: `KSP`, `KSPGetSolution()`, `KSPBuildResidualDefault()`

# External Links
$(_doc_external("KSP/KSPBuildSolutionDefault"))
"""
function KSPBuildSolutionDefault(petsclib::PetscLibType, ksp::AbstractKSP, v::AbstractPetscVec, M_V::AbstractPetscVec)
    error("KSPBuildSolutionDefault: no generated method for these argument types")
end

@for_petsc function KSPBuildSolutionDefault(petsclib::$UnionPetscLib, ksp::AbstractKSP, v::AbstractPetscVec, M_V::AbstractPetscVec )
	M_V_ = Ref(M_V.ptr)

    @chk ccall(
               (:KSPBuildSolutionDefault, $petsc_library),
               PetscErrorCode,
               (CKSP, CVec, Ptr{CVec}),
               ksp, v, M_V_,
              )

	M_V.ptr = M_V_[]

	return nothing
end 

"""
	norm_d::PetscReal = KSPCGGetNormD(petsclib::PetscLibType, ksp::AbstractKSP) 
Get norm of the direction when the solver is used inside `SNESNEWTONTR`

Not collective

Input Parameters:
- `ksp`    - the iterative context
- `norm_d` - the norm of the direction

Level: advanced

See also: `KSP`, `KSPCG`, `KSPNASH`, `KSPSTCG`, `KSPGLTR`, `SNESNEWTONTR`

# External Links
$(_doc_external("KSP/KSPCGGetNormD"))
"""
function KSPCGGetNormD(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPCGGetNormD: no generated method for these argument types")
end

@for_petsc function KSPCGGetNormD(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	norm_d_ = Ref{$PetscReal}()

    @chk ccall(
               (:KSPCGGetNormD, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscReal}),
               ksp, norm_d_,
              )

	norm_d = norm_d_[]

	return norm_d
end 

"""
	o_fcn::PetscReal = KSPCGGetObjFcn(petsclib::PetscLibType, ksp::AbstractKSP) 
Get the conjugate gradient objective function value

Not collective

Input Parameters:
- `ksp`   - the iterative context
- `o_fcn` - the objective function value

Level: advanced

See also: `KSP`, `KSPCG`, `KSPNASH`, `KSPSTCG`, `KSPGLTR`, `KSPMonitorSet`

# External Links
$(_doc_external("KSP/KSPCGGetObjFcn"))
"""
function KSPCGGetObjFcn(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPCGGetObjFcn: no generated method for these argument types")
end

@for_petsc function KSPCGGetObjFcn(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	o_fcn_ = Ref{$PetscReal}()

    @chk ccall(
               (:KSPCGGetObjFcn, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscReal}),
               ksp, o_fcn_,
              )

	o_fcn = o_fcn_[]

	return o_fcn
end 

"""
	KSPCGSetObjectiveTarget(petsclib::PetscLibType, ksp::AbstractKSP, obj::PetscReal) 
Sets the target value for the CG quadratic model

Logically Collective

Input Parameters:
- `ksp` - the iterative context
- `obj` - the objective value (0 is the default)

Level: advanced

See also: `KSP`, `KSPCG`, `KSPNASH`, `KSPSTCG`, `KSPGLTR`, `SNESNEWTONTR`

# External Links
$(_doc_external("KSP/KSPCGSetObjectiveTarget"))
"""
function KSPCGSetObjectiveTarget(petsclib::PetscLibType, ksp::AbstractKSP, obj::Real)
    error("KSPCGSetObjectiveTarget: no generated method for these argument types")
end

@for_petsc function KSPCGSetObjectiveTarget(petsclib::$UnionPetscLib, ksp::AbstractKSP, obj::$PetscReal )

    @chk ccall(
               (:KSPCGSetObjectiveTarget, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscReal),
               ksp, obj,
              )


	return nothing
end 

"""
	KSPCGSetRadius(petsclib::PetscLibType, ksp::AbstractKSP, radius::PetscReal) 
Sets the radius of the trust region used by the `KSPCG` when the solver is used inside `SNESNEWTONTR`

Logically Collective

Input Parameters:
- `ksp`    - the iterative context
- `radius` - the trust region radius (0 is the default that disable the use of the radius)

Level: advanced

See also: `KSP`, `KSPCG`, `KSPNASH`, `KSPSTCG`, `KSPGLTR`, `SNESNEWTONTR`

# External Links
$(_doc_external("KSP/KSPCGSetRadius"))
"""
function KSPCGSetRadius(petsclib::PetscLibType, ksp::AbstractKSP, radius::Real)
    error("KSPCGSetRadius: no generated method for these argument types")
end

@for_petsc function KSPCGSetRadius(petsclib::$UnionPetscLib, ksp::AbstractKSP, radius::$PetscReal )

    @chk ccall(
               (:KSPCGSetRadius, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscReal),
               ksp, radius,
              )


	return nothing
end 

"""
	KSPCGSetType(petsclib::PetscLibType, ksp::AbstractKSP, type::KSPCGType) 
Sets the variant of the conjugate gradient method to
use for solving a linear system with a complex coefficient matrix.
This option is irrelevant when solving a real system.

Logically Collective

Input Parameters:
- `ksp`  - the iterative context
- `type` - the variant of CG to use, one of
``
KSP_CG_HERMITIAN - complex, Hermitian matrix (default)
KSP_CG_SYMMETRIC - complex, symmetric matrix
``

Options Database Keys:
- `-ksp_cg_type hermitian` - Indicates Hermitian matrix
- `-ksp_cg_type symmetric` - Indicates symmetric matrix

Level: intermediate

See also: `KSP`, `KSPCG`

# External Links
$(_doc_external("KSP/KSPCGSetType"))
"""
function KSPCGSetType(petsclib::PetscLibType, ksp::AbstractKSP, type::KSPCGType)
    error("KSPCGSetType: no generated method for these argument types")
end

@for_petsc function KSPCGSetType(petsclib::$UnionPetscLib, ksp::AbstractKSP, type::KSPCGType )

    @chk ccall(
               (:KSPCGSetType, $petsc_library),
               PetscErrorCode,
               (CKSP, KSPCGType),
               ksp, type,
              )


	return nothing
end 

"""
	KSPCGUseSingleReduction(petsclib::PetscLibType, ksp::AbstractKSP, flg::PetscBool) 
Merge the two inner products needed in `KSPCG` into a single `MPI_Allreduce()` call.

Logically Collective

Input Parameters:
- `ksp` - the iterative context
- `flg` - turn on or off the single reduction

Options Database Key:
- `-ksp_cg_single_reduction (true|false)` - Merge inner products into single `MPI_Allreduce()`

Level: intermediate

See also: `KSP`, `KSPCG`, `KSPGMRES`, `KSPPIPECG`, `KSPPIPECR`, `KSPGROPPCG`

# External Links
$(_doc_external("KSP/KSPCGUseSingleReduction"))
"""
function KSPCGUseSingleReduction(petsclib::PetscLibType, ksp::AbstractKSP, flg::PetscBool)
    error("KSPCGUseSingleReduction: no generated method for these argument types")
end

@for_petsc function KSPCGUseSingleReduction(petsclib::$UnionPetscLib, ksp::AbstractKSP, flg::PetscBool )

    @chk ccall(
               (:KSPCGUseSingleReduction, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flg,
              )


	return nothing
end 

"""
	kspest::KSP = KSPChebyshevEstEigGetKSP(petsclib::PetscLibType, ksp::AbstractKSP) 
Get the Krylov method context used to estimate the eigenvalues for the Chebyshev method.

Input Parameter:
- `ksp` - the Krylov space context

Output Parameter:
- `kspest` - the eigenvalue estimation Krylov space context

Level: advanced

See also: `KSPCHEBYSHEV`, `KSPChebyshevEstEigSet()`

# External Links
$(_doc_external("KSP/KSPChebyshevEstEigGetKSP"))
"""
function KSPChebyshevEstEigGetKSP(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPChebyshevEstEigGetKSP: no generated method for these argument types")
end

@for_petsc function KSPChebyshevEstEigGetKSP(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	kspest_ = Ref{CKSP}()

    @chk ccall(
               (:KSPChebyshevEstEigGetKSP, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{CKSP}),
               ksp, kspest_,
              )

	kspest = KSP(kspest_[], petsclib; own = false)

	return kspest
end 

"""
	KSPChebyshevEstEigSet(petsclib::PetscLibType, ksp::AbstractKSP, a::PetscReal, b::PetscReal, c::PetscReal, d::PetscReal) 
Automatically estimate the eigenvalues to use for Chebyshev

Logically Collective

Input Parameters:
- `ksp` - the Krylov space context
- `a`   - multiple of min eigenvalue estimate to use for min Chebyshev bound (or `PETSC_DECIDE`)
- `b`   - multiple of max eigenvalue estimate to use for min Chebyshev bound (or `PETSC_DECIDE`)
- `c`   - multiple of min eigenvalue estimate to use for max Chebyshev bound (or `PETSC_DECIDE`)
- `d`   - multiple of max eigenvalue estimate to use for max Chebyshev bound (or `PETSC_DECIDE`)

Options Database Key:
- `-ksp_chebyshev_esteig a,b,c,d` - estimate eigenvalues using a Krylov method, then use this transform for Chebyshev eigenvalue bounds

See also: `KSPCHEBYSHEV`, `KSPChebyshevEstEigSetUseNoisy()`, `KSPChebyshevEstEigGetKSP()`

# External Links
$(_doc_external("KSP/KSPChebyshevEstEigSet"))
"""
function KSPChebyshevEstEigSet(petsclib::PetscLibType, ksp::AbstractKSP, a::Real, b::Real, c::Real, d::Real)
    error("KSPChebyshevEstEigSet: no generated method for these argument types")
end

@for_petsc function KSPChebyshevEstEigSet(petsclib::$UnionPetscLib, ksp::AbstractKSP, a::$PetscReal, b::$PetscReal, c::$PetscReal, d::$PetscReal )

    @chk ccall(
               (:KSPChebyshevEstEigSet, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscReal, $PetscReal, $PetscReal, $PetscReal),
               ksp, a, b, c, d,
              )


	return nothing
end 

"""
	KSPChebyshevEstEigSetUseNoisy(petsclib::PetscLibType, ksp::AbstractKSP, use::PetscBool) 
use a noisy random number generated right-hand side to estimate the extreme eigenvalues instead of the given right-hand side

Logically Collective

Input Parameters:
- `ksp` - linear solver context
- `use` - `PETSC_TRUE` to use noisy

Options Database Key:
- `-ksp_chebyshev_esteig_noisy (true|false)` - Use noisy right-hand side for estimate

Level: intermediate

See also: `KSPCHEBYSHEV`, `KSPChebyshevEstEigSet()`, `KSPChebyshevEstEigGetKSP()`

# External Links
$(_doc_external("KSP/KSPChebyshevEstEigSetUseNoisy"))
"""
function KSPChebyshevEstEigSetUseNoisy(petsclib::PetscLibType, ksp::AbstractKSP, use::PetscBool)
    error("KSPChebyshevEstEigSetUseNoisy: no generated method for these argument types")
end

@for_petsc function KSPChebyshevEstEigSetUseNoisy(petsclib::$UnionPetscLib, ksp::AbstractKSP, use::PetscBool )

    @chk ccall(
               (:KSPChebyshevEstEigSetUseNoisy, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, use,
              )


	return nothing
end 

"""
	kind::KSPChebyshevKind = KSPChebyshevGetKind(petsclib::PetscLibType, ksp::AbstractKSP) 
get the kind of Chebyshev polynomial to use

Logically Collective

Input Parameters:
- `ksp`  - Linear solver context
- `kind` - The kind of Chebyshev polynomial used

Level: intermediate

See also: `KSPCHEBYSHEV`, `KSPChebyshevKind`, `KSPChebyshevSetKind()`, `KSP_CHEBYSHEV_FIRST`, `KSP_CHEBYSHEV_FOURTH`, `KSP_CHEBYSHEV_OPT_FOURTH`

# External Links
$(_doc_external("KSP/KSPChebyshevGetKind"))
"""
function KSPChebyshevGetKind(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPChebyshevGetKind: no generated method for these argument types")
end

@for_petsc function KSPChebyshevGetKind(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	kind_ = Ref{KSPChebyshevKind}()

    @chk ccall(
               (:KSPChebyshevGetKind, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPChebyshevKind}),
               ksp, kind_,
              )

	kind = kind_[]

	return kind
end 

"""
	KSPChebyshevSetEigenvalues(petsclib::PetscLibType, ksp::AbstractKSP, emax::PetscReal, emin::PetscReal) 
Sets estimates for the extreme eigenvalues of the preconditioned problem.

Logically Collective

Input Parameters:
- `ksp`  - the Krylov space context
- `emax` - the eigenvalue maximum estimate
- `emin` - the eigenvalue minimum estimate

Options Database Key:
- `-ksp_chebyshev_eigenvalues emin,emax` - extreme eigenvalues

Level: intermediate

See also: `KSPCHEBYSHEV`, `KSPChebyshevEstEigSet()`

# External Links
$(_doc_external("KSP/KSPChebyshevSetEigenvalues"))
"""
function KSPChebyshevSetEigenvalues(petsclib::PetscLibType, ksp::AbstractKSP, emax::Real, emin::Real)
    error("KSPChebyshevSetEigenvalues: no generated method for these argument types")
end

@for_petsc function KSPChebyshevSetEigenvalues(petsclib::$UnionPetscLib, ksp::AbstractKSP, emax::$PetscReal, emin::$PetscReal )

    @chk ccall(
               (:KSPChebyshevSetEigenvalues, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscReal, $PetscReal),
               ksp, emax, emin,
              )


	return nothing
end 

"""
	KSPChebyshevSetKind(petsclib::PetscLibType, ksp::AbstractKSP, kind::KSPChebyshevKind) 
set the kind of Chebyshev polynomial to use

Logically Collective

Input Parameters:
- `ksp`  - Linear solver context
- `kind` - The kind of Chebyshev polynomial to use, see `KSPChebyshevKind`, one of `KSP_CHEBYSHEV_FIRST`, `KSP_CHEBYSHEV_FOURTH`, or `KSP_CHEBYSHEV_OPT_FOURTH`

Options Database Key:
- `-ksp_chebyshev_kind (first|fourth|opt_fourth)` - which kind of Chebyshev polynomial to use

Level: intermediate

See also: `KSPCHEBYSHEV`, `KSPChebyshevKind`, `KSPChebyshevGetKind()`, `KSP_CHEBYSHEV_FIRST`, `KSP_CHEBYSHEV_FOURTH`, `KSP_CHEBYSHEV_OPT_FOURTH`

# External Links
$(_doc_external("KSP/KSPChebyshevSetKind"))
"""
function KSPChebyshevSetKind(petsclib::PetscLibType, ksp::AbstractKSP, kind::KSPChebyshevKind)
    error("KSPChebyshevSetKind: no generated method for these argument types")
end

@for_petsc function KSPChebyshevSetKind(petsclib::$UnionPetscLib, ksp::AbstractKSP, kind::KSPChebyshevKind )

    @chk ccall(
               (:KSPChebyshevSetKind, $petsc_library),
               PetscErrorCode,
               (CKSP, KSPChebyshevKind),
               ksp, kind,
              )


	return nothing
end 

"""
	KSPCheckSolve(petsclib::PetscLibType, ksp::AbstractKSP, pc::AbstractPC, vec::AbstractPetscVec) 
Checks if the `PCSetUp()` or `KSPSolve()` failed and set the error flag for the outer `PC`. A `KSP_DIVERGED_ITS` is
not considered a failure in this context

Collective

Input Parameters:
- `ksp` - the linear solver `KSP` context.
- `pc`  - the preconditioner context
- `vec` - a vector that will be initialized with infinity to indicate lack of convergence

Level: developer

See also: `KSP`, `KSPCreate()`, `KSPSetType()`, `KSPCheckNorm()`, `KSPCheckDot()`

# External Links
$(_doc_external("KSP/KSPCheckSolve"))
"""
function KSPCheckSolve(petsclib::PetscLibType, ksp::AbstractKSP, pc::AbstractPC, vec::AbstractPetscVec)
    error("KSPCheckSolve: no generated method for these argument types")
end

@for_petsc function KSPCheckSolve(petsclib::$UnionPetscLib, ksp::AbstractKSP, pc::AbstractPC, vec::AbstractPetscVec )

    @chk ccall(
               (:KSPCheckSolve, $petsc_library),
               PetscErrorCode,
               (CKSP, CPC, CVec),
               ksp, pc, vec,
              )


	return nothing
end 

"""
	cr::PetscReal,rRsq::PetscReal,ce::PetscReal,eRsq::PetscReal = KSPComputeConvergenceRate(petsclib::PetscLibType, ksp::AbstractKSP) 
Compute the convergence rate for the iteration <https:/en.wikipedia.org/wiki/Coefficient_of_determination>

Not Collective

Input Parameter:
- `ksp` - The `KSP`

Output Parameters:
- `cr`   - The residual contraction rate
- `rRsq` - The coefficient of determination, R^2, indicating the linearity of the data
- `ce`   - The error contraction rate
- `eRsq` - The coefficient of determination, R^2, indicating the linearity of the data

Level: advanced

See also: `KSP`, `KSPConvergedRateView()`

# External Links
$(_doc_external("KSP/KSPComputeConvergenceRate"))
"""
function KSPComputeConvergenceRate(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPComputeConvergenceRate: no generated method for these argument types")
end

@for_petsc function KSPComputeConvergenceRate(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	cr_ = Ref{$PetscReal}()
	rRsq_ = Ref{$PetscReal}()
	ce_ = Ref{$PetscReal}()
	eRsq_ = Ref{$PetscReal}()

    @chk ccall(
               (:KSPComputeConvergenceRate, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               ksp, cr_, rRsq_, ce_, eRsq_,
              )

	cr = cr_[]
	rRsq = rRsq_[]
	ce = ce_[]
	eRsq = eRsq_[]

	return cr,rRsq,ce,eRsq
end 

"""
	neig::PetscInt = KSPComputeEigenvalues(petsclib::PetscLibType, ksp::AbstractKSP, n::PetscInt, r::Vector{PetscReal}, c::Vector{PetscReal}) 
Computes the extreme eigenvalues for the
preconditioned operator. Called after or during `KSPSolve()`.

Not Collective

Input Parameters:
- `ksp` - iterative solver obtained from `KSPCreate()`
- `n`   - size of arrays `r` and `c`. The number of eigenvalues computed `neig` will, in general, be less than this.

Output Parameters:
- `r`    - real part of computed eigenvalues, provided by user with a dimension of at least `n`
- `c`    - complex part of computed eigenvalues, provided by user with a dimension of at least `n`
- `neig` - actual number of eigenvalues computed (will be less than or equal to `n`)

Options Database Key:
- `-ksp_view_eigenvalues` - Prints eigenvalues to stdout

Level: advanced

See also: `KSPSetComputeEigenvalues()`, `KSPSetComputeSingularValues()`, `KSPMonitorSingularValue()`, `KSPComputeExtremeSingularValues()`, `KSP`, `KSPComputeRitz()`

# External Links
$(_doc_external("KSP/KSPComputeEigenvalues"))
"""
function KSPComputeEigenvalues(petsclib::PetscLibType, ksp::AbstractKSP, n::Integer, r::AbstractVector{<:Number}, c::AbstractVector{<:Number})
    error("KSPComputeEigenvalues: no generated method for these argument types")
end

@for_petsc function KSPComputeEigenvalues(petsclib::$UnionPetscLib, ksp::AbstractKSP, n::$PetscInt, r::Vector{$PetscReal}, c::Vector{$PetscReal} )
	neig_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPComputeEigenvalues, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscInt}),
               ksp, n, r, c, neig_,
              )

	neig = neig_[]

	return neig
end 

"""
	KSPComputeEigenvaluesExplicitly(petsclib::PetscLibType, ksp::AbstractKSP, nmax::PetscInt, r::Vector{PetscReal}, c::Vector{PetscReal}) 
Computes all of the eigenvalues of the
preconditioned operator using LAPACK.

Collective

Input Parameters:
- `ksp`  - iterative context obtained from `KSPCreate()`
- `nmax` - size of arrays `r` and `c`

Output Parameters:
- `r` - real part of computed eigenvalues, provided by user with a dimension at least of `n`
- `c` - complex part of computed eigenvalues, provided by user with a dimension at least of `n`

Level: advanced

See also: `KSP`, `KSPComputeEigenvalues()`, `KSPMonitorSingularValue()`, `KSPComputeExtremeSingularValues()`, `KSPSetOperators()`, `KSPSolve()`

# External Links
$(_doc_external("KSP/KSPComputeEigenvaluesExplicitly"))
"""
function KSPComputeEigenvaluesExplicitly(petsclib::PetscLibType, ksp::AbstractKSP, nmax::Integer, r::AbstractVector{<:Number}, c::AbstractVector{<:Number})
    error("KSPComputeEigenvaluesExplicitly: no generated method for these argument types")
end

@for_petsc function KSPComputeEigenvaluesExplicitly(petsclib::$UnionPetscLib, ksp::AbstractKSP, nmax::$PetscInt, r::Vector{$PetscReal}, c::Vector{$PetscReal} )

    @chk ccall(
               (:KSPComputeEigenvaluesExplicitly, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}),
               ksp, nmax, r, c,
              )


	return nothing
end 

"""
	emax::PetscReal,emin::PetscReal = KSPComputeExtremeSingularValues(petsclib::PetscLibType, ksp::AbstractKSP) 
Computes the extreme singular values
for the preconditioned operator. Called after or during `KSPSolve()`.

Not Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameters:
- `emax` - maximum estimated singular value
- `emin` - minimum estimated singular value

Options Database Key:
- `-ksp_view_singularvalues` - compute extreme singular values and print when `KSPSolve()` completes.

Level: advanced

See also: `KSPSetComputeSingularValues()`, `KSPMonitorSingularValue()`, `KSPComputeEigenvalues()`, `KSP`, `KSPComputeRitz()`

# External Links
$(_doc_external("KSP/KSPComputeExtremeSingularValues"))
"""
function KSPComputeExtremeSingularValues(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPComputeExtremeSingularValues: no generated method for these argument types")
end

@for_petsc function KSPComputeExtremeSingularValues(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	emax_ = Ref{$PetscReal}()
	emin_ = Ref{$PetscReal}()

    @chk ccall(
               (:KSPComputeExtremeSingularValues, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscReal}, Ptr{$PetscReal}),
               ksp, emax_, emin_,
              )

	emax = emax_[]
	emin = emin_[]

	return emax,emin
end 

"""
	mat::PetscMat = KSPComputeOperator(petsclib::PetscLibType, ksp::AbstractKSP, mattype::String) 
Computes the explicit preconditioned operator, including diagonal scaling and null
space removal if applicable.

Collective

Input Parameters:
- `ksp`     - the Krylov subspace context
- `mattype` - the matrix type to be used

Output Parameter:
- `mat` - the explicit preconditioned operator

Level: advanced

See also: `KSP`, `KSPSetOperators()`, `KSPComputeEigenvaluesExplicitly()`, `PCComputeOperator()`, `KSPSetDiagonalScale()`, `KSPSetNullSpace()`, `MatType`

# External Links
$(_doc_external("KSP/KSPComputeOperator"))
"""
function KSPComputeOperator(petsclib::PetscLibType, ksp::AbstractKSP, mattype::String)
    error("KSPComputeOperator: no generated method for these argument types")
end

@for_petsc function KSPComputeOperator(petsclib::$UnionPetscLib, ksp::AbstractKSP, mattype::String )
	mat_ = Ref{CMat}()

    @chk ccall(
               (:KSPComputeOperator, $petsc_library),
               PetscErrorCode,
               (CKSP, MatType, Ptr{CMat}),
               ksp, mattype, mat_,
              )

	mat = PetscMat(mat_[], petsclib)

	return mat
end 

"""
	nrit::PetscInt = KSPComputeRitz(petsclib::PetscLibType, ksp::AbstractKSP, ritz::PetscBool, small::PetscBool, S::Vector{<:AbstractPetscVec}, tetar::Vector{PetscReal}, tetai::Vector{PetscReal}) 
Computes the Ritz or harmonic Ritz pairs associated with the
smallest or largest in modulus, for the preconditioned operator.

Not Collective

Input Parameters:
- `ksp`   - iterative solver obtained from `KSPCreate()`
- `ritz`  - `PETSC_TRUE` or `PETSC_FALSE` for Ritz pairs or harmonic Ritz pairs, respectively
- `small` - `PETSC_TRUE` or `PETSC_FALSE` for smallest or largest (harmonic) Ritz values, respectively

Output Parameters:
- `nrit`  - On input number of (harmonic) Ritz pairs to compute; on output, actual number of computed (harmonic) Ritz pairs
- `S`     - an array of the Ritz vectors, pass in an array of vectors of size `nrit`
- `tetar` - real part of the Ritz values, pass in an array of size `nrit`
- `tetai` - imaginary part of the Ritz values, pass in an array of size `nrit`

Level: advanced

See also: `KSPSetComputeRitz()`, `KSP`, `KSPGMRES`, `KSPComputeEigenvalues()`, `KSPSetComputeSingularValues()`, `KSPMonitorSingularValue()`

# External Links
$(_doc_external("KSP/KSPComputeRitz"))
"""
function KSPComputeRitz(petsclib::PetscLibType, ksp::AbstractKSP, ritz::PetscBool, small::PetscBool, S::Vector{<:AbstractPetscVec}, tetar::AbstractVector{<:Number}, tetai::AbstractVector{<:Number})
    error("KSPComputeRitz: no generated method for these argument types")
end

@for_petsc function KSPComputeRitz(petsclib::$UnionPetscLib, ksp::AbstractKSP, ritz::PetscBool, small::PetscBool, S::Vector{<:AbstractPetscVec}, tetar::Vector{$PetscReal}, tetai::Vector{$PetscReal} )
	nrit_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPComputeRitz, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool, PetscBool, Ptr{$PetscInt}, Ptr{CVec}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               ksp, ritz, small, nrit_, S, tetar, tetai,
              )

	nrit = nrit_[]

	return nrit
end 

"""
	reason::KSPConvergedReason = KSPConvergedDefault(petsclib::PetscLibType, ksp::AbstractKSP, n::PetscInt, rnorm::PetscReal, ctx::Ptr{Cvoid}) 
Default code to determine convergence of the linear iterative solvers

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - residual norm (may be estimated, depending on the method may be the preconditioned residual norm)
- `ctx`   - convergence context which must be created by `KSPConvergedDefaultCreate()`

Output Parameter:
- `reason` - the convergence reason; it is positive if the iteration has converged,
negative if the iteration has diverged, and `KSP_CONVERGED_ITERATING` otherwise

Options Database Keys:
- `-ksp_max_it`                                  - maximum number of linear iterations
- `-ksp_min_it`                                  - minimum number of linear iterations, defaults to 0
- `-ksp_rtol rtol`                               - relative tolerance used in default determination of convergence, i.e. if residual norm decreases by this factor than convergence is declared
- `-ksp_atol abstol`                             - absolute tolerance used in default convergence test, i.e. if residual norm is less than this then convergence is declared
- `-ksp_divtol tol`                              - if residual norm increases by this factor than divergence is declared
- `-ksp_converged_use_initial_residual_norm`     - see `KSPConvergedDefaultSetUIRNorm()`
- `-ksp_converged_use_min_initial_residual_norm` - see `KSPConvergedDefaultSetUMIRNorm()`
- `-ksp_converged_maxits`                        - see `KSPConvergedDefaultSetConvergedMaxits()`

Level: advanced

See also: `KSP`, `KSPSetConvergenceTest()`, `KSPSetTolerances()`, `KSPConvergedSkip()`, `KSPConvergedReason`, `KSPGetConvergedReason()`,
`KSPSetMinimumIterations()`, `KSPConvergenceTestFn`,
`KSPConvergedDefaultSetUIRNorm()`, `KSPConvergedDefaultSetUMIRNorm()`, `KSPConvergedDefaultSetConvergedMaxits()`, `KSPConvergedDefaultCreate()`, `KSPConvergedDefaultDestroy()`

# External Links
$(_doc_external("KSP/KSPConvergedDefault"))
"""
function KSPConvergedDefault(petsclib::PetscLibType, ksp::AbstractKSP, n::Integer, rnorm::Real, ctx::Ptr{Cvoid})
    error("KSPConvergedDefault: no generated method for these argument types")
end

@for_petsc function KSPConvergedDefault(petsclib::$UnionPetscLib, ksp::AbstractKSP, n::$PetscInt, rnorm::$PetscReal, ctx::Ptr{Cvoid} )
	reason_ = Ref{KSPConvergedReason}()

    @chk ccall(
               (:KSPConvergedDefault, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{KSPConvergedReason}, Ptr{Cvoid}),
               ksp, n, rnorm, reason_, ctx,
              )

	reason = reason_[]

	return reason
end 

"""
	ctx::Ptr{Cvoid} = KSPConvergedDefaultCreate(petsclib::PetscLibType) 
Creates and initializes the context used by the `KSPConvergedDefault()` function

Not Collective

Output Parameter:
- `ctx` - convergence context

Level: intermediate

See also: `KSP`, `KSPConvergedDefault()`, `KSPConvergedDefaultDestroy()`, `KSPSetConvergenceTest()`, `KSPSetTolerances()`,
`KSPConvergedSkip()`, `KSPConvergedReason`, `KSPGetConvergedReason()`, `KSPConvergedDefaultSetUIRNorm()`, `KSPConvergedDefaultSetUMIRNorm()`,
`KSPConvergedDefaultSetConvergedMaxits()`

# External Links
$(_doc_external("KSP/KSPConvergedDefaultCreate"))
"""
function KSPConvergedDefaultCreate(petsclib::PetscLibType)
    error("KSPConvergedDefaultCreate: no generated method for these argument types")
end

@for_petsc function KSPConvergedDefaultCreate(petsclib::$UnionPetscLib)
	ctx_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:KSPConvergedDefaultCreate, $petsc_library),
               PetscErrorCode,
               (Ptr{Ptr{Cvoid}},),
               ctx_,
              )

	ctx = ctx_[]

	return ctx
end 

"""
	KSPConvergedDefaultDestroy(petsclib::PetscLibType, ctx::Ptr{Cvoid}) 
Frees the space used by the `KSPConvergedDefault()` function context

Not Collective

Input Parameter:
- `ctx` - convergence context

Level: intermediate

See also: `KSP`, `KSPConvergedDefault()`, `KSPConvergedDefaultCreate()`, `KSPSetConvergenceTest()`, `KSPSetTolerances()`, `KSPConvergedSkip()`,
`KSPConvergedReason`, `KSPGetConvergedReason()`, `KSPConvergedDefaultSetUIRNorm()`, `KSPConvergedDefaultSetUMIRNorm()`

# External Links
$(_doc_external("KSP/KSPConvergedDefaultDestroy"))
"""
function KSPConvergedDefaultDestroy(petsclib::PetscLibType, ctx::Ptr{Cvoid})
    error("KSPConvergedDefaultDestroy: no generated method for these argument types")
end

@for_petsc function KSPConvergedDefaultDestroy(petsclib::$UnionPetscLib, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:KSPConvergedDefaultDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{Cvoid},),
               ctx,
              )


	return nothing
end 

"""
	KSPConvergedDefaultSetConvergedMaxits(petsclib::PetscLibType, ksp::AbstractKSP, flg::PetscBool) 
allows the default convergence test to declare convergence and return `KSP_CONVERGED_ITS` if the maximum number of iterations is reached

Collective

Input Parameters:
- `ksp` - iterative context
- `flg` - boolean flag

Options Database Key:
- `-ksp_converged_maxits (true|false)` - Declare convergence if the maximum number of iterations is reached

Level: intermediate

See also: `KSP`, `KSPSetConvergenceTest()`, `KSPSetTolerances()`, `KSPConvergedSkip()`, `KSPConvergedReason`, `KSPGetConvergedReason()`, `KSPConvergedDefaultSetUMIRNorm()`, `KSPConvergedDefaultSetUIRNorm()`

# External Links
$(_doc_external("KSP/KSPConvergedDefaultSetConvergedMaxits"))
"""
function KSPConvergedDefaultSetConvergedMaxits(petsclib::PetscLibType, ksp::AbstractKSP, flg::PetscBool)
    error("KSPConvergedDefaultSetConvergedMaxits: no generated method for these argument types")
end

@for_petsc function KSPConvergedDefaultSetConvergedMaxits(petsclib::$UnionPetscLib, ksp::AbstractKSP, flg::PetscBool )

    @chk ccall(
               (:KSPConvergedDefaultSetConvergedMaxits, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flg,
              )


	return nothing
end 

"""
	KSPConvergedDefaultSetUIRNorm(petsclib::PetscLibType, ksp::AbstractKSP) 
makes the default convergence test use  || B*(b - A*(initial guess))||
instead of  || B*b ||. In the case of right preconditioner or if `KSPSetNormType`(ksp,`KSP_NORM_UNPRECONDITIONED`)
is used there is no B in the above formula.

Collective

Input Parameters:
- `ksp` - iterative context

Options Database Key:
- `-ksp_converged_use_initial_residual_norm (true|false)` - Use initial residual norm for computing relative convergence

Level: intermediate

See also: `KSP`, `KSPSetConvergenceTest()`, `KSPSetTolerances()`, `KSPConvergedSkip()`, `KSPConvergedReason`, `KSPGetConvergedReason()`, `KSPConvergedDefaultSetUMIRNorm()`, `KSPConvergedDefaultSetConvergedMaxits()`

# External Links
$(_doc_external("KSP/KSPConvergedDefaultSetUIRNorm"))
"""
function KSPConvergedDefaultSetUIRNorm(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPConvergedDefaultSetUIRNorm: no generated method for these argument types")
end

@for_petsc function KSPConvergedDefaultSetUIRNorm(petsclib::$UnionPetscLib, ksp::AbstractKSP )

    @chk ccall(
               (:KSPConvergedDefaultSetUIRNorm, $petsc_library),
               PetscErrorCode,
               (CKSP,),
               ksp,
              )


	return nothing
end 

"""
	KSPConvergedDefaultSetUMIRNorm(petsclib::PetscLibType, ksp::AbstractKSP) 
makes the default convergence test use \\min(|| B*(b - A*(initial guess))||,|| B*b ||)
In the case of right preconditioner or if `KSPSetNormType`(ksp,`KSP_NORM_UNPRECONDITIONED`)
is used there is no B in the above formula.

Collective

Input Parameters:
- `ksp` - iterative context

Options Database Key:
- `-ksp_converged_use_min_initial_residual_norm (true|false)` - Use minimum of initial residual norm and b for computing relative convergence

Level: intermediate

See also: `KSP`, `KSPSetConvergenceTest()`, `KSPSetTolerances()`, `KSPConvergedSkip()`, `KSPConvergedReason`, `KSPGetConvergedReason()`, `KSPConvergedDefaultSetUIRNorm()`, `KSPConvergedDefaultSetConvergedMaxits()`

# External Links
$(_doc_external("KSP/KSPConvergedDefaultSetUMIRNorm"))
"""
function KSPConvergedDefaultSetUMIRNorm(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPConvergedDefaultSetUMIRNorm: no generated method for these argument types")
end

@for_petsc function KSPConvergedDefaultSetUMIRNorm(petsclib::$UnionPetscLib, ksp::AbstractKSP )

    @chk ccall(
               (:KSPConvergedDefaultSetUMIRNorm, $petsc_library),
               PetscErrorCode,
               (CKSP,),
               ksp,
              )


	return nothing
end 

"""
	KSPConvergedRateView(petsclib::PetscLibType, ksp::AbstractKSP, viewer::PetscViewer) 
Displays the convergence rate <https://en.wikipedia.org/wiki/Coefficient_of_determination> of `KSPSolve()` to a viewer

Collective

Input Parameters:
- `ksp`    - iterative solver obtained from `KSPCreate()`
- `viewer` - the `PetscViewer` to display the reason

Options Database Key:
- `-ksp_converged_rate` - print reason for convergence or divergence and the convergence rate (or 0.0 for divergence)

Level: intermediate

See also: `KSPConvergedReasonView()`, `KSPGetConvergedRate()`, `KSPSetTolerances()`, `KSPConvergedDefault()`

# External Links
$(_doc_external("KSP/KSPConvergedRateView"))
"""
function KSPConvergedRateView(petsclib::PetscLibType, ksp::AbstractKSP, viewer::PetscViewer)
    error("KSPConvergedRateView: no generated method for these argument types")
end

@for_petsc function KSPConvergedRateView(petsclib::$UnionPetscLib, ksp::AbstractKSP, viewer::PetscViewer )

    @chk ccall(
               (:KSPConvergedRateView, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscViewer),
               ksp, viewer,
              )


	return nothing
end 

"""
	KSPConvergedReasonView(petsclib::PetscLibType, ksp::AbstractKSP, viewer::PetscViewer) 
Displays the reason a `KSP` solve converged or diverged, `KSPConvergedReason` to a `PetscViewer`

Collective

Input Parameters:
- `ksp`    - iterative solver obtained from `KSPCreate()`
- `viewer` - the `PetscViewer` on which to display the reason

Options Database Keys:
- `-ksp_converged_reason`          - print reason for converged or diverged, also prints number of iterations
- `-ksp_converged_reason ::failed` - only print reason and number of iterations when diverged

Level: beginner

See also: `KSPConvergedReasonViewFromOptions()`, `KSPCreate()`, `KSPSetUp()`, `KSPDestroy()`, `KSPSetTolerances()`, `KSPConvergedDefault()`,
`KSPSolveTranspose()`, `KSPGetIterationNumber()`, `KSP`, `KSPGetConvergedReason()`, `PetscViewerPushFormat()`, `PetscViewerPopFormat()`

# External Links
$(_doc_external("KSP/KSPConvergedReasonView"))
"""
function KSPConvergedReasonView(petsclib::PetscLibType, ksp::AbstractKSP, viewer::PetscViewer)
    error("KSPConvergedReasonView: no generated method for these argument types")
end

@for_petsc function KSPConvergedReasonView(petsclib::$UnionPetscLib, ksp::AbstractKSP, viewer::PetscViewer )

    @chk ccall(
               (:KSPConvergedReasonView, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscViewer),
               ksp, viewer,
              )


	return nothing
end 

"""
	KSPConvergedReasonViewCancel(petsclib::PetscLibType, ksp::AbstractKSP) 
Clears all the `KSPConvergedReason` view functions for a `KSP` object set with `KSPConvergedReasonViewSet()`
as well as the default viewer.

Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Level: intermediate

See also: `KSPCreate()`, `KSPDestroy()`, `KSPReset()`, `KSPConvergedReasonViewSet()`

# External Links
$(_doc_external("KSP/KSPConvergedReasonViewCancel"))
"""
function KSPConvergedReasonViewCancel(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPConvergedReasonViewCancel: no generated method for these argument types")
end

@for_petsc function KSPConvergedReasonViewCancel(petsclib::$UnionPetscLib, ksp::AbstractKSP )

    @chk ccall(
               (:KSPConvergedReasonViewCancel, $petsc_library),
               PetscErrorCode,
               (CKSP,),
               ksp,
              )


	return nothing
end 

"""
	KSPConvergedReasonViewFromOptions(petsclib::PetscLibType, ksp::AbstractKSP) 
Processes command line options to determine if/how a `KSPConvergedReason` is to be viewed.

Collective

Input Parameter:
- `ksp` - the `KSP` object

Level: intermediate

See also: `KSPConvergedReasonView()`, `KSPConvergedReasonViewSet()`

# External Links
$(_doc_external("KSP/KSPConvergedReasonViewFromOptions"))
"""
function KSPConvergedReasonViewFromOptions(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPConvergedReasonViewFromOptions: no generated method for these argument types")
end

@for_petsc function KSPConvergedReasonViewFromOptions(petsclib::$UnionPetscLib, ksp::AbstractKSP )

    @chk ccall(
               (:KSPConvergedReasonViewFromOptions, $petsc_library),
               PetscErrorCode,
               (CKSP,),
               ksp,
              )


	return nothing
end 

"""
	KSPConvergedReasonViewSet(petsclib::PetscLibType, ksp::AbstractKSP, f::Ptr{Cvoid}, ctx::Ptr{Cvoid}, reasonviewdestroy::Ptr{Cvoid}) 
Sets an ADDITIONAL function that is to be used at the
end of the linear solver to display the convergence reason of the linear solver.

Logically Collective

Input Parameters:
- `ksp`               - the `KSP` context
- `f`                 - the `ksp` converged reason view function, see `KSPConvergedReasonViewFn`
- `ctx`               - [optional] context for private data for the `KSPConvergedReason` view routine (use `NULL` if context is not needed)
- `reasonviewdestroy` - [optional] routine that frees `ctx` (may be `NULL`), see `PetscCtxDestroyFn` for the calling sequence

Options Database Keys:
- `-ksp_converged_reason`             - sets a default `KSPConvergedReasonView()`
- `-ksp_converged_reason_view_cancel` - cancels all converged reason viewers that have been hardwired into a code by
calls to `KSPConvergedReasonViewSet()`, but does not cancel those set via the options database.

Level: intermediate

See also: `KSPConvergedReasonView()`, `KSPConvergedReasonViewFn`, `KSPConvergedReasonViewCancel()`, `PetscCtxDestroyFn`

# External Links
$(_doc_external("KSP/KSPConvergedReasonViewSet"))
"""
function KSPConvergedReasonViewSet(petsclib::PetscLibType, ksp::AbstractKSP, f::Ptr{Cvoid}, ctx::Ptr{Cvoid}, reasonviewdestroy::Ptr{Cvoid})
    error("KSPConvergedReasonViewSet: no generated method for these argument types")
end

@for_petsc function KSPConvergedReasonViewSet(petsclib::$UnionPetscLib, ksp::AbstractKSP, f::Ptr{Cvoid}, ctx::Ptr{Cvoid}, reasonviewdestroy::Ptr{Cvoid} )

    @chk ccall(
               (:KSPConvergedReasonViewSet, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
               ksp, f, ctx, reasonviewdestroy,
              )


	return nothing
end 

"""
	reason::KSPConvergedReason = KSPConvergedSkip(petsclib::PetscLibType, ksp::AbstractKSP, n::PetscInt, rnorm::PetscReal, dtx::Ptr{Cvoid}) 
Convergence test that do not return as converged
until the maximum number of iterations is reached.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm residual value (may be estimated)
- `dtx`   - unused convergence context

Output Parameter:
- `reason` - `KSP_CONVERGED_ITERATING` or `KSP_CONVERGED_ITS`

Options Database Key:
- `-ksp_convergence_test skip` - skips the test

Level: advanced

See also: `KSP`, `KSPCG`, `KSPBCGS`, `KSPConvergenceTestFn`, `KSPSetConvergenceTest()`, `KSPSetTolerances()`, `KSPSetNormType()`,
`KSPConvergedReason`

# External Links
$(_doc_external("KSP/KSPConvergedSkip"))
"""
function KSPConvergedSkip(petsclib::PetscLibType, ksp::AbstractKSP, n::Integer, rnorm::Real, dtx::Ptr{Cvoid})
    error("KSPConvergedSkip: no generated method for these argument types")
end

@for_petsc function KSPConvergedSkip(petsclib::$UnionPetscLib, ksp::AbstractKSP, n::$PetscInt, rnorm::$PetscReal, dtx::Ptr{Cvoid} )
	reason_ = Ref{KSPConvergedReason}()

    @chk ccall(
               (:KSPConvergedSkip, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{KSPConvergedReason}, Ptr{Cvoid}),
               ksp, n, rnorm, reason_, dtx,
              )

	reason = reason_[]

	return reason
end 

"""
	inksp::KSP = KSPCreate(petsclib::PetscLibType, comm::MPI_Comm) 
Creates the `KSP` context. This `KSP` context is used in PETSc to solve linear systems with `KSPSolve()`

Collective

Input Parameter:
- `comm` - MPI communicator

Output Parameter:
- `inksp` - location to put the `KSP` context

Level: beginner

See also: `KSPSetUp()`, `KSPSolve()`, `KSPDestroy()`, `KSP`, `KSPGMRES`, `KSPType`, `KSPSetType()`

# External Links
$(_doc_external("KSP/KSPCreate"))
"""
function KSPCreate(petsclib::PetscLibType, comm::MPI_Comm)
    error("KSPCreate: no generated method for these argument types")
end

@for_petsc function KSPCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	inksp_ = Ref{CKSP}()

    @chk ccall(
               (:KSPCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{CKSP}),
               comm, inksp_,
              )

	inksp = KSP(inksp_[], petsclib)

	return inksp
end 

# override for KSPCreateVecs; C signature: KSPCreateVecs(KSP ksp, PetscInt rightn, Vec* right[], PetscInt leftn, Vec* left[])
"""
	right::Vector{PetscVec},left::Vector{PetscVec} = KSPCreateVecs(petsclib::PetscLibType, ksp::AbstractKSP, rightn::PetscInt, leftn::PetscInt) 
Gets a number of work vectors suitably sized for the operator in the `KSP`

Collective

Input Parameters:
- `ksp`    - iterative context
- `rightn` - number of right work vectors to allocate
- `leftn`  - number of left work vectors to allocate

Output Parameters:
- `right` - the array of vectors created
- `left`  - the array of left vectors

Level: advanced

See also: `MatCreateVecs()`, `VecDestroyVecs()`, `KSPSetWorkVecs()`

# External Links
$(_doc_external("KSP/KSPCreateVecs"))
"""
function KSPCreateVecs(petsclib::PetscLibType, ksp::AbstractKSP, rightn::PetscInt, leftn::PetscInt) end

@for_petsc function KSPCreateVecs(petsclib::$UnionPetscLib, ksp::AbstractKSP, rightn::$PetscInt, leftn::$PetscInt )
	right_ = Ref{Ptr{CVec}}()
	left_ = Ref{Ptr{CVec}}()

    @chk ccall(
               (:KSPCreateVecs, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, Ptr{Ptr{CVec}}, $PetscInt, Ptr{Ptr{CVec}}),
               ksp, rightn, right_, leftn, left_,
              )

	a_v = unsafe_wrap(Array, right_[], rightn; own = false)
    if rightn != 0
        v = PetscVec(a_v[1], petsclib)
        right = ntuple(i -> similar(v), rightn)
    else
        right = nothing
    end

    
    a_v = unsafe_wrap(Array, left_[], leftn; own = false)
    if leftn != 0
        v = PetscVec(a_v[1], petsclib)
        left = ntuple(i -> similar(v), leftn)
    else
        left = nothing
    end

	return right,left
end

"""
	KSPDestroy(petsclib::PetscLibType, ksp::AbstractKSP) 
Destroys a `KSP` context.

Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Level: beginner

See also: `KSPCreate()`, `KSPSetUp()`, `KSPSolve()`, `KSP`

# External Links
$(_doc_external("KSP/KSPDestroy"))
"""
function KSPDestroy(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPDestroy: no generated method for these argument types")
end

@for_petsc function KSPDestroy(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	ksp_ = Ref(ksp.ptr)

    @chk ccall(
               (:KSPDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{CKSP},),
               ksp_,
              )

	ksp.ptr = C_NULL

	return nothing
end 

"""
	KSPDestroyDefault(petsclib::PetscLibType, ksp::AbstractKSP) 
Destroys an iterative context variable for methods with no separate context.  Preferred calling sequence `KSPDestroy()`.

Collective

Input Parameter:
- `ksp` - the iterative context

Level: advanced

See also: `KSP`, `KSPDestroy()`

# External Links
$(_doc_external("KSP/KSPDestroyDefault"))
"""
function KSPDestroyDefault(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPDestroyDefault: no generated method for these argument types")
end

@for_petsc function KSPDestroyDefault(petsclib::$UnionPetscLib, ksp::AbstractKSP )

    @chk ccall(
               (:KSPDestroyDefault, $petsc_library),
               PetscErrorCode,
               (CKSP,),
               ksp,
              )


	return nothing
end 

"""
	mmax::PetscInt = KSPFCGGetMmax(petsclib::PetscLibType, ksp::AbstractKSP) 
get the maximum number of previous directions `KSPFCG` will store

Not Collective

Input Parameter:
- `ksp` - the Krylov space context

Output Parameter:
- `mmax` - the maximum number of previous directions allowed for orthogonalization

Level: intermediate

See also: `KSPFCG`, `KSPFCGGetTruncationType()`, `KSPFCGGetNprealloc()`, `KSPFCGSetMmax()`

# External Links
$(_doc_external("KSP/KSPFCGGetMmax"))
"""
function KSPFCGGetMmax(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPFCGGetMmax: no generated method for these argument types")
end

@for_petsc function KSPFCGGetMmax(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	mmax_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPFCGGetMmax, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscInt}),
               ksp, mmax_,
              )

	mmax = mmax_[]

	return mmax
end 

"""
	nprealloc::PetscInt = KSPFCGGetNprealloc(petsclib::PetscLibType, ksp::AbstractKSP) 
get the number of directions preallocate by `KSPFCG`

Not Collective

Input Parameter:
- `ksp` - the Krylov space context

Output Parameter:
- `nprealloc` - the number of directions preallocated

Level: advanced

See also: `KSPFCG`, `KSPFCGGetTruncationType()`, `KSPFCGSetNprealloc()`, `KSPFCGSetMmax()`, `KSPFCGGetMmax()`

# External Links
$(_doc_external("KSP/KSPFCGGetNprealloc"))
"""
function KSPFCGGetNprealloc(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPFCGGetNprealloc: no generated method for these argument types")
end

@for_petsc function KSPFCGGetNprealloc(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	nprealloc_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPFCGGetNprealloc, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscInt}),
               ksp, nprealloc_,
              )

	nprealloc = nprealloc_[]

	return nprealloc
end 

"""
	truncstrat::KSPFCDTruncationType = KSPFCGGetTruncationType(petsclib::PetscLibType, ksp::AbstractKSP) 
get the truncation strategy employed by `KSPFCG`

Not Collective

Input Parameter:
- `ksp` - the Krylov space context

Output Parameter:
- `truncstrat` - the strategy type

Level: intermediate

See also: `KSPFCG`, `KSPFCGSetTruncationType()`, `KSPFCDTruncationType`, `KSP_FCD_TRUNC_TYPE_STANDARD`, `KSP_FCD_TRUNC_TYPE_NOTAY`

# External Links
$(_doc_external("KSP/KSPFCGGetTruncationType"))
"""
function KSPFCGGetTruncationType(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPFCGGetTruncationType: no generated method for these argument types")
end

@for_petsc function KSPFCGGetTruncationType(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	truncstrat_ = Ref{KSPFCDTruncationType}()

    @chk ccall(
               (:KSPFCGGetTruncationType, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPFCDTruncationType}),
               ksp, truncstrat_,
              )

	truncstrat = truncstrat_[]

	return truncstrat
end 

"""
	KSPFCGSetMmax(petsclib::PetscLibType, ksp::AbstractKSP, mmax::PetscInt) 
set the maximum number of previous directions `KSPFCG` will store for orthogonalization

Logically Collective

Input Parameters:
- `ksp`  - the Krylov space context
- `mmax` - the maximum number of previous directions to orthogonalize against

Options Database Key:
- `-ksp_fcg_mmax N` - maximum number of search directions

Level: intermediate

See also: `KSPFCG`, `KSPFCGGetTruncationType()`, `KSPFCDTruncationType`, `KSPFCGSetTruncationType()`, `KSPFCGGetNprealloc()`, `KSPFCGGetMmax()`

# External Links
$(_doc_external("KSP/KSPFCGSetMmax"))
"""
function KSPFCGSetMmax(petsclib::PetscLibType, ksp::AbstractKSP, mmax::Integer)
    error("KSPFCGSetMmax: no generated method for these argument types")
end

@for_petsc function KSPFCGSetMmax(petsclib::$UnionPetscLib, ksp::AbstractKSP, mmax::$PetscInt )

    @chk ccall(
               (:KSPFCGSetMmax, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, mmax,
              )


	return nothing
end 

"""
	KSPFCGSetNprealloc(petsclib::PetscLibType, ksp::AbstractKSP, nprealloc::PetscInt) 
set the number of directions to preallocate with `KSPFCG`

Logically Collective

Input Parameters:
- `ksp`       - the Krylov space context
- `nprealloc` - the number of vectors to preallocate

Options Database Key:
- `-ksp_fcg_nprealloc N` - number of directions to preallocate

Level: advanced

See also: `KSPFCG`, `KSPFCGGetTruncationType()`, `KSPFCGGetNprealloc()`, `KSPFCGSetMmax()`, `KSPFCGGetMmax()`

# External Links
$(_doc_external("KSP/KSPFCGSetNprealloc"))
"""
function KSPFCGSetNprealloc(petsclib::PetscLibType, ksp::AbstractKSP, nprealloc::Integer)
    error("KSPFCGSetNprealloc: no generated method for these argument types")
end

@for_petsc function KSPFCGSetNprealloc(petsclib::$UnionPetscLib, ksp::AbstractKSP, nprealloc::$PetscInt )

    @chk ccall(
               (:KSPFCGSetNprealloc, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, nprealloc,
              )


	return nothing
end 

"""
	KSPFCGSetTruncationType(petsclib::PetscLibType, ksp::AbstractKSP, truncstrat::KSPFCDTruncationType) 
specify how many of its stored previous directions `KSPFCG` uses during orthogonalization

Logically Collective

Input Parameters:
- `ksp`        - the Krylov space context
- `truncstrat` - the choice of strategy
``
KSP_FCD_TRUNC_TYPE_STANDARD uses all (up to `mmax`) stored directions
KSP_FCD_TRUNC_TYPE_NOTAY uses the last `max(1,mod(i,mmax))` stored directions at iteration i = 0, 1, ...
``

Options Database Key:
- `-ksp_fcg_truncation_type (standard|notay)` - specify how many of its stored previous directions `KSPFCG` uses during orthogonalization

Level: intermediate

See also: `KSPFCG`, `KSPFCDTruncationType`, `KSPFCGGetTruncationType()`, `KSPFCGSetNprealloc()`, `KSPFCGSetMmax()`, `KSPFCGGetMmax()`,
`KSP_FCD_TRUNC_TYPE_STANDARD`, `KSP_FCD_TRUNC_TYPE_NOTAY`

# External Links
$(_doc_external("KSP/KSPFCGSetTruncationType"))
"""
function KSPFCGSetTruncationType(petsclib::PetscLibType, ksp::AbstractKSP, truncstrat::KSPFCDTruncationType)
    error("KSPFCGSetTruncationType: no generated method for these argument types")
end

@for_petsc function KSPFCGSetTruncationType(petsclib::$UnionPetscLib, ksp::AbstractKSP, truncstrat::KSPFCDTruncationType )

    @chk ccall(
               (:KSPFCGSetTruncationType, $petsc_library),
               PetscErrorCode,
               (CKSP, KSPFCDTruncationType),
               ksp, truncstrat,
              )


	return nothing
end 

"""
	pc::PC = KSPFETIDPGetInnerBDDC(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets the `PCBDDC` preconditioner used to set up the `KSPFETIDP` matrix for the Lagrange multipliers

Input Parameter:
- `ksp` - the `KSPFETIDP` Krylov solver

Output Parameter:
- `pc` - the `PCBDDC` preconditioner

Level: advanced

See also: `MATIS`, `PCBDDC`, `KSPFETIDP`, `KSPFETIDPSetInnerBDDC()`, `KSPFETIDPGetInnerKSP()`

# External Links
$(_doc_external("KSP/KSPFETIDPGetInnerBDDC"))
"""
function KSPFETIDPGetInnerBDDC(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPFETIDPGetInnerBDDC: no generated method for these argument types")
end

@for_petsc function KSPFETIDPGetInnerBDDC(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	pc_ = Ref{CPC}()

    @chk ccall(
               (:KSPFETIDPGetInnerBDDC, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{CPC}),
               ksp, pc_,
              )

	pc = PC(pc_[], petsclib; own = false)

	return pc
end 

"""
	innerksp::KSP = KSPFETIDPGetInnerKSP(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets the `KSP` object for the Lagrange multipliers from inside a `KSPFETIDP`

Input Parameter:
- `ksp` - the `KSPFETIDP`

Output Parameter:
- `innerksp` - the `KSP` for the multipliers

Level: advanced

See also: `KSPFETIDP`, `MATIS`, `PCBDDC`, `KSPFETIDPSetInnerBDDC()`, `KSPFETIDPGetInnerBDDC()`

# External Links
$(_doc_external("KSP/KSPFETIDPGetInnerKSP"))
"""
function KSPFETIDPGetInnerKSP(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPFETIDPGetInnerKSP: no generated method for these argument types")
end

@for_petsc function KSPFETIDPGetInnerKSP(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	innerksp_ = Ref{CKSP}()

    @chk ccall(
               (:KSPFETIDPGetInnerKSP, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{CKSP}),
               ksp, innerksp_,
              )

	innerksp = KSP(innerksp_[], petsclib; own = false)

	return innerksp
end 

"""
	KSPFETIDPSetInnerBDDC(petsclib::PetscLibType, ksp::AbstractKSP, pc::AbstractPC) 
Provides the `PCBDDC` preconditioner used to set up the `KSPFETIDP` matrix for the Lagrange multipliers

Collective

Input Parameters:
- `ksp` - the `KSPFETIDP` Krylov solver
- `pc`  - the `PCBDDC` preconditioner

Level: advanced

See also: `MATIS`, `PCBDDC`, `KSPFETIDPGetInnerBDDC()`, `KSPFETIDPGetInnerKSP()`

# External Links
$(_doc_external("KSP/KSPFETIDPSetInnerBDDC"))
"""
function KSPFETIDPSetInnerBDDC(petsclib::PetscLibType, ksp::AbstractKSP, pc::AbstractPC)
    error("KSPFETIDPSetInnerBDDC: no generated method for these argument types")
end

@for_petsc function KSPFETIDPSetInnerBDDC(petsclib::$UnionPetscLib, ksp::AbstractKSP, pc::AbstractPC )

    @chk ccall(
               (:KSPFETIDPSetInnerBDDC, $petsc_library),
               PetscErrorCode,
               (CKSP, CPC),
               ksp, pc,
              )


	return nothing
end 

"""
	KSPFETIDPSetPressureOperator(petsclib::PetscLibType, ksp::AbstractKSP, P::AbstractPetscMat) 
Sets the operator used to set up the pressure preconditioner for the saddle point `KSPFETIDP` solver,

Collective

Input Parameters:
- `ksp` - the `KSPFETIDP` solver
- `P`   - the linear operator to be preconditioned, usually the mass matrix.

Level: advanced

See also: `KSPFETIDP`, `MATIS`, `PCBDDC`, `KSPFETIDPGetInnerBDDC()`, `KSPFETIDPGetInnerKSP()`, `KSPSetOperators()`

# External Links
$(_doc_external("KSP/KSPFETIDPSetPressureOperator"))
"""
function KSPFETIDPSetPressureOperator(petsclib::PetscLibType, ksp::AbstractKSP, P::AbstractPetscMat)
    error("KSPFETIDPSetPressureOperator: no generated method for these argument types")
end

@for_petsc function KSPFETIDPSetPressureOperator(petsclib::$UnionPetscLib, ksp::AbstractKSP, P::AbstractPetscMat )

    @chk ccall(
               (:KSPFETIDPSetPressureOperator, $petsc_library),
               PetscErrorCode,
               (CKSP, CMat),
               ksp, P,
              )


	return nothing
end 

"""
	KSPFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the `KSP` package. It is
called from `PetscFinalize()`.

Level: developer

See also: `PetscFinalize()`, `KSPInitializePackage()`

# External Links
$(_doc_external("KSP/KSPFinalizePackage"))
"""
function KSPFinalizePackage(petsclib::PetscLibType)
    error("KSPFinalizePackage: no generated method for these argument types")
end

@for_petsc function KSPFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:KSPFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	KSPFlexibleModifyPCKSP(petsclib::PetscLibType, ksp::AbstractKSP, total_its::PetscInt, loc_its::PetscInt, res_norm::PetscReal, ctx::Ptr{Cvoid}) 
modifies the attributes of the `PCKSP` preconditioner, see .

Input Parameters:
- `ksp`       - the ksp context being used.
- `total_its` - the total number of `KSP` iterations that have occurred.
- `loc_its`   - the number of `KSP` iterations since last restart.
- `res_norm`  - the current residual norm.
- `ctx`       - context, unused in this routine

Level: intermediate

See also: `KSPFGMRES`, `KSPFCG`, `KSPPIPEFCG`, `KSPGCR`, `KSPPIPEGCR`, `KSPFlexibleModifyPCFn`, `KSPFlexibleSetModifyPC()`

# External Links
$(_doc_external("KSP/KSPFlexibleModifyPCKSP"))
"""
function KSPFlexibleModifyPCKSP(petsclib::PetscLibType, ksp::AbstractKSP, total_its::Integer, loc_its::Integer, res_norm::Real, ctx::Ptr{Cvoid})
    error("KSPFlexibleModifyPCKSP: no generated method for these argument types")
end

@for_petsc function KSPFlexibleModifyPCKSP(petsclib::$UnionPetscLib, ksp::AbstractKSP, total_its::$PetscInt, loc_its::$PetscInt, res_norm::$PetscReal, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:KSPFlexibleModifyPCKSP, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscInt, $PetscReal, Ptr{Cvoid}),
               ksp, total_its, loc_its, res_norm, ctx,
              )


	return nothing
end 

"""
	KSPFlexibleModifyPCNoChange(petsclib::PetscLibType, ksp::AbstractKSP, total_its::PetscInt, loc_its::PetscInt, res_norm::PetscReal, ctx::Ptr{Cvoid}) 
this is the default used by the flexible Krylov methods - it doesn't change the preconditioner. 

Input Parameters:
- `ksp`       - the ksp context being used.
- `total_its` - the total number of `KSP` iterations that have occurred.
- `loc_its`   - the number of `KSP` iterations since last restart.
- `res_norm`  - the current residual norm.
- `ctx`       - context variable, unused in this routine

Level: intermediate

See also: `KSPFGMRES`, `KSPFCG`, `KSPPIPEFCG`, `KSPGCR`, `KSPPIPEGCR`, `KSPFlexibleModifyPCFn`, `KSPFlexibleSetModifyPC()`, `KSPFlexibleModifyPCKSP()`

# External Links
$(_doc_external("KSP/KSPFlexibleModifyPCNoChange"))
"""
function KSPFlexibleModifyPCNoChange(petsclib::PetscLibType, ksp::AbstractKSP, total_its::Integer, loc_its::Integer, res_norm::Real, ctx::Ptr{Cvoid})
    error("KSPFlexibleModifyPCNoChange: no generated method for these argument types")
end

@for_petsc function KSPFlexibleModifyPCNoChange(petsclib::$UnionPetscLib, ksp::AbstractKSP, total_its::$PetscInt, loc_its::$PetscInt, res_norm::$PetscReal, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:KSPFlexibleModifyPCNoChange, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscInt, $PetscReal, Ptr{Cvoid}),
               ksp, total_its, loc_its, res_norm, ctx,
              )


	return nothing
end 

"""
	KSPFlexibleSetModifyPC(petsclib::PetscLibType, ksp::AbstractKSP, fcn::Ptr{Cvoid}, ctx::Ptr{Cvoid}, destroy::Ptr{Cvoid}) 
Sets the routine used by flexible `KSP` methods to modify the preconditioner. 

Logically Collective

Input Parameters:
- `ksp`     - iterative context obtained from `KSPCreate()`
- `fcn`     - function to modify the `PC`, see `KSPFlexibleModifyPCFn`
- `ctx`     - optional context
- `destroy` - optional context destroy routine

Level: intermediate

See also: `KSPFGMRES`, `KSPFCG`, `KSPPIPEFCG`, `KSPGCR`, `KSPPIPEGCR`, `KSPFlexibleModifyPCFn`, `KSPFlexibleModifyPCNoChange()`, `KSPFlexibleModifyPCKSP()`

# External Links
$(_doc_external("KSP/KSPFlexibleSetModifyPC"))
"""
function KSPFlexibleSetModifyPC(petsclib::PetscLibType, ksp::AbstractKSP, fcn::Ptr{Cvoid}, ctx::Ptr{Cvoid}, destroy::Ptr{Cvoid})
    error("KSPFlexibleSetModifyPC: no generated method for these argument types")
end

@for_petsc function KSPFlexibleSetModifyPC(petsclib::$UnionPetscLib, ksp::AbstractKSP, fcn::Ptr{Cvoid}, ctx::Ptr{Cvoid}, destroy::Ptr{Cvoid} )

    @chk ccall(
               (:KSPFlexibleSetModifyPC, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
               ksp, fcn, ctx, destroy,
              )


	return nothing
end 

"""
	restart::PetscInt = KSPGCRGetRestart(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets number of iterations at which `KSPGCR` restarts.

Not Collective

Input Parameter:
- `ksp` - the Krylov space context

Output Parameter:
- `restart` - integer restart value

Level: intermediate

See also: `KSPGCR`, `KSPSetTolerances()`, `KSPGCRSetRestart()`, `KSPGMRESGetRestart()`

# External Links
$(_doc_external("KSP/KSPGCRGetRestart"))
"""
function KSPGCRGetRestart(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGCRGetRestart: no generated method for these argument types")
end

@for_petsc function KSPGCRGetRestart(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	restart_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPGCRGetRestart, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscInt}),
               ksp, restart_,
              )

	restart = restart_[]

	return restart
end 

"""
	KSPGCRSetRestart(petsclib::PetscLibType, ksp::AbstractKSP, restart::PetscInt) 
Sets number of iterations at which `KSPGCR` restarts.

Not Collective

Input Parameters:
- `ksp`     - the Krylov space context
- `restart` - integer restart value

Options Database Key:
- `-ksp_gcr_restart restart` - the number of stored vectors to orthogonalize against

Level: intermediate

See also: `KSPGCR`, `KSPSetTolerances()`, `KSPGCRGetRestart()`, `KSPGMRESSetRestart()`

# External Links
$(_doc_external("KSP/KSPGCRSetRestart"))
"""
function KSPGCRSetRestart(petsclib::PetscLibType, ksp::AbstractKSP, restart::Integer)
    error("KSPGCRSetRestart: no generated method for these argument types")
end

@for_petsc function KSPGCRSetRestart(petsclib::$UnionPetscLib, ksp::AbstractKSP, restart::$PetscInt )

    @chk ccall(
               (:KSPGCRSetRestart, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, restart,
              )


	return nothing
end 

"""
	lambda::PetscReal = KSPGLTRGetLambda(petsclib::PetscLibType, ksp::AbstractKSP) 
Get the multiplier on the trust-region constraint when using `KSPGLTR`

Not Collective

Input Parameter:
- `ksp` - the iterative context

Output Parameter:
- `lambda` - the multiplier

Level: advanced

See also: `KSP`, `KSPGLTR`, `KSPGLTRGetMinEig()`

# External Links
$(_doc_external("KSP/KSPGLTRGetLambda"))
"""
function KSPGLTRGetLambda(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGLTRGetLambda: no generated method for these argument types")
end

@for_petsc function KSPGLTRGetLambda(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	lambda_ = Ref{$PetscReal}()

    @chk ccall(
               (:KSPGLTRGetLambda, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscReal}),
               ksp, lambda_,
              )

	lambda = lambda_[]

	return lambda
end 

"""
	e_min::PetscReal = KSPGLTRGetMinEig(petsclib::PetscLibType, ksp::AbstractKSP) 
Get minimum eigenvalue computed by `KSPGLTR`

Collective

Input Parameter:
- `ksp` - the iterative context

Output Parameter:
- `e_min` - the minimum eigenvalue

Level: advanced

See also: `KSP`, `KSPGLTR`, `KSPGLTRGetLambda()`

# External Links
$(_doc_external("KSP/KSPGLTRGetMinEig"))
"""
function KSPGLTRGetMinEig(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGLTRGetMinEig: no generated method for these argument types")
end

@for_petsc function KSPGLTRGetMinEig(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	e_min_ = Ref{$PetscReal}()

    @chk ccall(
               (:KSPGLTRGetMinEig, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscReal}),
               ksp, e_min_,
              )

	e_min = e_min_[]

	return e_min
end 

"""
	KSPGMRESClassicalGramSchmidtOrthogonalization(petsclib::PetscLibType, ksp::AbstractKSP, it::PetscInt) 
This is the basic orthogonalization routine
using classical Gram-Schmidt with possible iterative refinement to improve the stability

Collective, No Fortran Support

Input Parameters:
- `ksp` - `KSP` object, must be associated with `KSPGMRES`, `KSPFGMRES`, or `KSPLGMRES` Krylov method
- `it`  - one less than the current GMRES restart iteration, i.e. the size of the Krylov space

Options Database Keys:
- `-ksp_gmres_classicalgramschmidt (true|false)`                                - Activates `KSPGMRESClassicalGramSchmidtOrthogonalization()`
- `-ksp_gmres_cgs_refinement_type (refine_never|refine_ifneeded|refine_always)` - determine if iterative refinement is
used to increase the stability of the classical Gram-Schmidt  orthogonalization.

Level: intermediate

See also: `KSPGMRESCGSRefinementType`, `KSPGMRESSetOrthogonalization()`, `KSPGMRESSetCGSRefinementType()`,
`KSPGMRESGetCGSRefinementType()`, `KSPGMRESGetOrthogonalization()`, `KSPGMRESModifiedGramSchmidtOrthogonalization()`

# External Links
$(_doc_external("KSP/KSPGMRESClassicalGramSchmidtOrthogonalization"))
"""
function KSPGMRESClassicalGramSchmidtOrthogonalization(petsclib::PetscLibType, ksp::AbstractKSP, it::Integer)
    error("KSPGMRESClassicalGramSchmidtOrthogonalization: no generated method for these argument types")
end

@for_petsc function KSPGMRESClassicalGramSchmidtOrthogonalization(petsclib::$UnionPetscLib, ksp::AbstractKSP, it::$PetscInt )

    @chk ccall(
               (:KSPGMRESClassicalGramSchmidtOrthogonalization, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, it,
              )


	return nothing
end 

"""
	type::KSPGMRESCGSRefinementType = KSPGMRESGetCGSRefinementType(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets the type of iterative refinement to use
in the classical Gram-Schmidt orthogonalization used by `KSPGMRES` and other PETSc GMRES implementations.

Not Collective

Input Parameter:
- `ksp` - the Krylov space solver context

Output Parameter:
- `type` - the type of refinement

Level: intermediate

See also: `KSPGMRES`, `KSPGMRESSetOrthogonalization()`, `KSPGMRESCGSRefinementType`, `KSPGMRESClassicalGramSchmidtOrthogonalization()`, `KSPGMRESSetCGSRefinementType()`,
`KSPGMRESGetOrthogonalization()`

# External Links
$(_doc_external("KSP/KSPGMRESGetCGSRefinementType"))
"""
function KSPGMRESGetCGSRefinementType(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGMRESGetCGSRefinementType: no generated method for these argument types")
end

@for_petsc function KSPGMRESGetCGSRefinementType(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	type_ = Ref{KSPGMRESCGSRefinementType}()

    @chk ccall(
               (:KSPGMRESGetCGSRefinementType, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPGMRESCGSRefinementType}),
               ksp, type_,
              )

	type = type_[]

	return type
end 

"""
	KSPGMRESGetOrthogonalization(petsclib::PetscLibType, ksp::AbstractKSP, noname::Ptr{Cvoid}) 
Gets the orthogonalization routine used by `KSPGMRES` and `KSPFGMRES`.

Not Collective

Input Parameter:
- `ksp` - iterative context obtained from `KSPCreate()`

Output Parameter:
- `fcn` - orthogonalization function

Calling sequence of `fcn`:
- `ksp` - the solver context
- `it`  - the current iteration

Level: intermediate

See also: `KSPGMRESSetRestart()`, `KSPGMRESSetPreAllocateVectors()`, `KSPGMRESSetCGSRefinementType()`, `KSPGMRESSetOrthogonalization()`,
`KSPGMRESModifiedGramSchmidtOrthogonalization()`, `KSPGMRESClassicalGramSchmidtOrthogonalization()`, `KSPGMRESGetCGSRefinementType()`

# External Links
$(_doc_external("KSP/KSPGMRESGetOrthogonalization"))
"""
function KSPGMRESGetOrthogonalization(petsclib::PetscLibType, ksp::AbstractKSP, noname::Ptr{Cvoid})
    error("KSPGMRESGetOrthogonalization: no generated method for these argument types")
end

@for_petsc function KSPGMRESGetOrthogonalization(petsclib::$UnionPetscLib, ksp::AbstractKSP, noname::Ptr{Cvoid} )

    @chk ccall(
               (:KSPGMRESGetOrthogonalization, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Cvoid}),
               ksp, noname,
              )


	return nothing
end 

"""
	restart::PetscInt = KSPGMRESGetRestart(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets number of iterations at which GMRES (`KSPGMRES`, `KSPFGMRES`, `KSPPGMRES`, `KSPAGMRES`, `KSPDGMRES`, `KSPPIPEFGMRES`,
and `KSPLGMRES`) restarts.

Not Collective

Input Parameter:
- `ksp` - the Krylov space solver context

Output Parameter:
- `restart` - integer restart value

Level: intermediate

See also: `KSPGMRES`, `KSPSetTolerances()`, `KSPGMRESSetOrthogonalization()`, `KSPGMRESSetPreAllocateVectors()`, `KSPGMRESSetRestart()`,
`KSPFGMRES`, `KSPLGMRES`, `KSPPGMRES`, `KSPAGMRES`, `KSPDGMRES`, `KSPPIPEFGMRES`

# External Links
$(_doc_external("KSP/KSPGMRESGetRestart"))
"""
function KSPGMRESGetRestart(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGMRESGetRestart: no generated method for these argument types")
end

@for_petsc function KSPGMRESGetRestart(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	restart_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPGMRESGetRestart, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscInt}),
               ksp, restart_,
              )

	restart = restart_[]

	return restart
end 

"""
	KSPGMRESModifiedGramSchmidtOrthogonalization(petsclib::PetscLibType, ksp::AbstractKSP, it::PetscInt) 
This is the basic orthogonalization routine
using modified Gram-Schmidt.

Collective, No Fortran Support

Input Parameters:
- `ksp` - `KSP` object, must be associated with `KSPGMRES`, `KSPFGMRES`, or `KSPLGMRES` Krylov method
- `it`  - one less than the current GMRES restart iteration, i.e. the size of the Krylov space

Options Database Key:
- `-ksp_gmres_modifiedgramschmidt` - Activates `KSPGMRESModifiedGramSchmidtOrthogonalization()`

Level: intermediate

See also: `KSPGMRESSetOrthogonalization()`, `KSPGMRESClassicalGramSchmidtOrthogonalization()`, `KSPGMRESGetOrthogonalization()`

# External Links
$(_doc_external("KSP/KSPGMRESModifiedGramSchmidtOrthogonalization"))
"""
function KSPGMRESModifiedGramSchmidtOrthogonalization(petsclib::PetscLibType, ksp::AbstractKSP, it::Integer)
    error("KSPGMRESModifiedGramSchmidtOrthogonalization: no generated method for these argument types")
end

@for_petsc function KSPGMRESModifiedGramSchmidtOrthogonalization(petsclib::$UnionPetscLib, ksp::AbstractKSP, it::$PetscInt )

    @chk ccall(
               (:KSPGMRESModifiedGramSchmidtOrthogonalization, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, it,
              )


	return nothing
end 

"""
	KSPGMRESMonitorKrylov(petsclib::PetscLibType, ksp::AbstractKSP, its::PetscInt, fgnorm::PetscReal, Viewers::Ptr{Cvoid}) 
Calls `VecView()` to monitor each new direction in the `KSPGMRES` accumulated Krylov space.

Collective

Input Parameters:
- `ksp`     - the `KSP` context
- `its`     - iteration number
- `fgnorm`  - 2-norm of residual (or gradient)
- `Viewers` - a collection of viewers created with `PetscViewersCreate()`

Options Database Key:
- `-ksp_gmres_krylov_monitor (true|false)` - Plot the Krylov directions

Level: intermediate

See also: `KSPGMRES`, `KSPMonitorSet()`, `KSPMonitorResidual()`, `VecView()`, `PetscViewersCreate()`, `PetscViewersDestroy()`

# External Links
$(_doc_external("KSP/KSPGMRESMonitorKrylov"))
"""
function KSPGMRESMonitorKrylov(petsclib::PetscLibType, ksp::AbstractKSP, its::Integer, fgnorm::Real, Viewers::Ptr{Cvoid})
    error("KSPGMRESMonitorKrylov: no generated method for these argument types")
end

@for_petsc function KSPGMRESMonitorKrylov(petsclib::$UnionPetscLib, ksp::AbstractKSP, its::$PetscInt, fgnorm::$PetscReal, Viewers::Ptr{Cvoid} )

    @chk ccall(
               (:KSPGMRESMonitorKrylov, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{Cvoid}),
               ksp, its, fgnorm, Viewers,
              )


	return nothing
end 

"""
	KSPGMRESSetBreakdownTolerance(petsclib::PetscLibType, ksp::AbstractKSP, tol::PetscReal) 
Sets the tolerance for determining divergence breakdown in `KSPGMRES` at restart.

Logically Collective

Input Parameters:
- `ksp` - the Krylov space solver context
- `tol` - the tolerance

Options Database Key:
- `-ksp_gmres_breakdown_tolerance tol` - set tolerance for determining divergence breakdown

Level: intermediate

See also: `KSPGMRES`, `KSPSetTolerances()`, `KSPGMRESSetHapTol()`, `KSPConvergedReason`

# External Links
$(_doc_external("KSP/KSPGMRESSetBreakdownTolerance"))
"""
function KSPGMRESSetBreakdownTolerance(petsclib::PetscLibType, ksp::AbstractKSP, tol::Real)
    error("KSPGMRESSetBreakdownTolerance: no generated method for these argument types")
end

@for_petsc function KSPGMRESSetBreakdownTolerance(petsclib::$UnionPetscLib, ksp::AbstractKSP, tol::$PetscReal )

    @chk ccall(
               (:KSPGMRESSetBreakdownTolerance, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscReal),
               ksp, tol,
              )


	return nothing
end 

"""
	KSPGMRESSetCGSRefinementType(petsclib::PetscLibType, ksp::AbstractKSP, type::KSPGMRESCGSRefinementType) 
Sets the type of iterative refinement to use
in the classical Gram-Schmidt orthogonalization used by `KSPGMRES` and other PETSc GMRES implementations.

Logically Collective

Input Parameters:
- `ksp`  - the Krylov space solver context
- `type` - the type of refinement
``
KSP_GMRES_CGS_REFINE_NEVER
KSP_GMRES_CGS_REFINE_IFNEEDED
KSP_GMRES_CGS_REFINE_ALWAYS
``

Options Database Key:
- `-ksp_gmres_cgs_refinement_type (refine_never|refine_ifneeded|refine_always)` - refinement type

Level: intermediate

See also: `KSPGMRES`, `KSPGMRESSetOrthogonalization()`, `KSPGMRESCGSRefinementType`, `KSPGMRESClassicalGramSchmidtOrthogonalization()`, `KSPGMRESGetCGSRefinementType()`,
`KSPGMRESGetOrthogonalization()`

# External Links
$(_doc_external("KSP/KSPGMRESSetCGSRefinementType"))
"""
function KSPGMRESSetCGSRefinementType(petsclib::PetscLibType, ksp::AbstractKSP, type::KSPGMRESCGSRefinementType)
    error("KSPGMRESSetCGSRefinementType: no generated method for these argument types")
end

@for_petsc function KSPGMRESSetCGSRefinementType(petsclib::$UnionPetscLib, ksp::AbstractKSP, type::KSPGMRESCGSRefinementType )

    @chk ccall(
               (:KSPGMRESSetCGSRefinementType, $petsc_library),
               PetscErrorCode,
               (CKSP, KSPGMRESCGSRefinementType),
               ksp, type,
              )


	return nothing
end 

"""
	KSPGMRESSetHapTol(petsclib::PetscLibType, ksp::AbstractKSP, tol::PetscReal) 
Sets the tolerance for detecting a happy breakdown in GMRES (`KSPGMRES`, `KSPFGMRES` and `KSPLGMRES` and others)

Logically Collective

Input Parameters:
- `ksp` - the Krylov space solver context
- `tol` - the tolerance for detecting a happy breakdown

Options Database Key:
- `-ksp_gmres_haptol tol` - set tolerance for determining happy breakdown

Level: intermediate

See also: `KSPGMRES`, `KSPSetTolerances()`

# External Links
$(_doc_external("KSP/KSPGMRESSetHapTol"))
"""
function KSPGMRESSetHapTol(petsclib::PetscLibType, ksp::AbstractKSP, tol::Real)
    error("KSPGMRESSetHapTol: no generated method for these argument types")
end

@for_petsc function KSPGMRESSetHapTol(petsclib::$UnionPetscLib, ksp::AbstractKSP, tol::$PetscReal )

    @chk ccall(
               (:KSPGMRESSetHapTol, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscReal),
               ksp, tol,
              )


	return nothing
end 

"""
	KSPGMRESSetOrthogonalization(petsclib::PetscLibType, ksp::AbstractKSP, fcn::external) 
Sets the orthogonalization routine used by `KSPGMRES` and `KSPFGMRES`.

Logically Collective

Input Parameters:
- `ksp` - iterative context obtained from `KSPCreate()`
- `fcn` - orthogonalization function

Calling sequence of `fcn`:
- `ksp` - the solver context
- `it`  - the current iteration

Options Database Keys:
- `-ksp_gmres_classicalgramschmidt` - Activates KSPGMRESClassicalGramSchmidtOrthogonalization() (default)
- `-ksp_gmres_modifiedgramschmidt`  - Activates KSPGMRESModifiedGramSchmidtOrthogonalization()

Level: intermediate

See also: `KSPGMRESSetRestart()`, `KSPGMRESSetPreAllocateVectors()`,
`KSPGMRESSetCGSRefinementType()`, `KSPGMRESModifiedGramSchmidtOrthogonalization()`,
`KSPGMRESClassicalGramSchmidtOrthogonalization()`, `KSPGMRESGetCGSRefinementType()`

# External Links
$(_doc_external("KSP/KSPGMRESSetOrthogonalization"))
"""
function KSPGMRESSetOrthogonalization(petsclib::PetscLibType, ksp::AbstractKSP, fcn::external)
    error("KSPGMRESSetOrthogonalization: no generated method for these argument types")
end

@for_petsc function KSPGMRESSetOrthogonalization(petsclib::$UnionPetscLib, ksp::AbstractKSP, fcn::external )

    @chk ccall(
               (:KSPGMRESSetOrthogonalization, $petsc_library),
               PetscErrorCode,
               (CKSP, external),
               ksp, fcn,
              )


	return nothing
end 

"""
	KSPGMRESSetPreAllocateVectors(petsclib::PetscLibType, ksp::AbstractKSP) 
Causes `KSPGMRES` and `KSPFGMRES` to preallocate all its
needed work vectors at initial setup rather than the default, which
is to allocate several at a time when needed.

Logically Collective

Input Parameter:
- `ksp` - iterative context obtained from `KSPCreate()`

Options Database Key:
- `-ksp_gmres_preallocate` - Activates `KSPGmresSetPreAllocateVectors()`

Level: intermediate

See also: `KSPGMRESSetRestart()`, `KSPGMRESSetOrthogonalization()`, `KSPGMRESGetOrthogonalization()`,
`VecMDot()`, `VecMAXPY()`

# External Links
$(_doc_external("KSP/KSPGMRESSetPreAllocateVectors"))
"""
function KSPGMRESSetPreAllocateVectors(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGMRESSetPreAllocateVectors: no generated method for these argument types")
end

@for_petsc function KSPGMRESSetPreAllocateVectors(petsclib::$UnionPetscLib, ksp::AbstractKSP )

    @chk ccall(
               (:KSPGMRESSetPreAllocateVectors, $petsc_library),
               PetscErrorCode,
               (CKSP,),
               ksp,
              )


	return nothing
end 

"""
	KSPGMRESSetRestart(petsclib::PetscLibType, ksp::AbstractKSP, restart::PetscInt) 
Sets number of iterations at which GMRES (`KSPGMRES`, `KSPFGMRES`, `KSPPGMRES`, `KSPAGMRES`, `KSPDGMRES`, `KSPPIPEFGMRES`,
and `KSPLGMRES`) restarts.

Logically Collective

Input Parameters:
- `ksp`     - the Krylov space solver context
- `restart` - integer restart value, this corresponds to the number of iterations of GMRES to perform before restarting

Options Database Key:
- `-ksp_gmres_restart restart` - integer restart value

Level: intermediate

See also: `KSPGMRES`, `KSPSetTolerances()`, `KSPGMRESSetOrthogonalization()`, `KSPGMRESSetPreAllocateVectors()`, `KSPGMRESGetRestart()`,
`KSPFGMRES`, `KSPLGMRES`, `KSPPGMRES`, `KSPAGMRES`, `KSPDGMRES`, `KSPPIPEFGMRES`

# External Links
$(_doc_external("KSP/KSPGMRESSetRestart"))
"""
function KSPGMRESSetRestart(petsclib::PetscLibType, ksp::AbstractKSP, restart::Integer)
    error("KSPGMRESSetRestart: no generated method for these argument types")
end

@for_petsc function KSPGMRESSetRestart(petsclib::$UnionPetscLib, ksp::AbstractKSP, restart::$PetscInt )

    @chk ccall(
               (:KSPGMRESSetRestart, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, restart,
              )


	return nothing
end 

"""
	converge::Ptr{Cvoid},ctx::Ptr{Cvoid},destroy::Ptr{Cvoid} = KSPGetAndClearConvergenceTest(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets the function to be used to determine convergence. Removes the current test without calling destroy on the test context

Logically Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameters:
- `converge` - pointer to convergence test function, see `KSPConvergenceTestFn`
- `ctx`      - context for private data for the convergence routine
- `destroy`  - a routine for destroying the context

Level: advanced

See also: `KSP`, `KSPConvergedDefault()`, `KSPGetConvergenceContext()`, `KSPSetTolerances()`, `KSPSetConvergenceTest()`, `KSPGetConvergenceTest()`

# External Links
$(_doc_external("KSP/KSPGetAndClearConvergenceTest"))
"""
function KSPGetAndClearConvergenceTest(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetAndClearConvergenceTest: no generated method for these argument types")
end

@for_petsc function KSPGetAndClearConvergenceTest(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	converge_ = Ref{Ptr{Cvoid}}()
	ctx_ = Ref{Ptr{Cvoid}}()
	destroy_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:KSPGetAndClearConvergenceTest, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Ptr{Cvoid}}, Ptr{Cvoid}, Ptr{Ptr{Cvoid}}),
               ksp, converge_, ctx_, destroy_,
              )

	converge = converge_[]
	ctx = ctx_[]
	destroy = destroy_[]

	return converge,ctx,destroy
end 

"""
	ctx::Ptr{Cvoid} = KSPGetApplicationContext(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets the user-defined context for the linear solver set with `KSPSetApplicationContext()`

Not Collective

Input Parameter:
- `ksp` - `KSP` context

Output Parameter:
- `ctx` - a pointer to the application context

Level: intermediate

See also: `KSP`, `KSPSetApplicationContext()`

# External Links
$(_doc_external("KSP/KSPGetApplicationContext"))
"""
function KSPGetApplicationContext(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetApplicationContext: no generated method for these argument types")
end

@for_petsc function KSPGetApplicationContext(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	ctx_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:KSPGetApplicationContext, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Cvoid}),
               ksp, ctx_,
              )

	ctx = ctx_[]

	return ctx
end 

"""
	flg::PetscBool = KSPGetComputeEigenvalues(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets the flag indicating that the extreme eigenvalues
values will be calculated via a Lanczos or Arnoldi process as the linear
system is solved.

Not Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameter:
- `flg` - `PETSC_TRUE` or `PETSC_FALSE`

Level: advanced

See also: `KSPComputeEigenvalues()`, `KSPComputeEigenvaluesExplicitly()`, `KSP`, `KSPSetComputeRitz()`

# External Links
$(_doc_external("KSP/KSPGetComputeEigenvalues"))
"""
function KSPGetComputeEigenvalues(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetComputeEigenvalues: no generated method for these argument types")
end

@for_petsc function KSPGetComputeEigenvalues(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:KSPGetComputeEigenvalues, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{PetscBool}),
               ksp, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	flg::PetscBool = KSPGetComputeSingularValues(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets the flag indicating whether the extreme singular
values will be calculated via a Lanczos or Arnoldi process as the linear
system is solved.

Not Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameter:
- `flg` - `PETSC_TRUE` or `PETSC_FALSE`

Options Database Key:
- `-ksp_monitor_singular_value` - Activates `KSPSetComputeSingularValues()`

Level: advanced

See also: `KSPComputeExtremeSingularValues()`, `KSPMonitorSingularValue()`, `KSP`

# External Links
$(_doc_external("KSP/KSPGetComputeSingularValues"))
"""
function KSPGetComputeSingularValues(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetComputeSingularValues: no generated method for these argument types")
end

@for_petsc function KSPGetComputeSingularValues(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:KSPGetComputeSingularValues, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{PetscBool}),
               ksp, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	flg::PetscBool = KSPGetConvergedNegativeCurvature(petsclib::PetscLibType, ksp::AbstractKSP) 
Get the flag to declare convergence if negative curvature is detected

Collective

Input Parameter:
- `ksp` - iterative context

Output Parameter:
- `flg` - the Boolean value

Level: advanced

See also: `KSP`, `KSPConvergedReason`, `KSPSetConvergedNegativeCurvature()`

# External Links
$(_doc_external("KSP/KSPGetConvergedNegativeCurvature"))
"""
function KSPGetConvergedNegativeCurvature(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetConvergedNegativeCurvature: no generated method for these argument types")
end

@for_petsc function KSPGetConvergedNegativeCurvature(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:KSPGetConvergedNegativeCurvature, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{PetscBool}),
               ksp, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	reason::KSPConvergedReason = KSPGetConvergedReason(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets the reason the `KSP` iteration was stopped.

Not Collective

Input Parameter:
- `ksp` - the `KSP` context

Output Parameter:
- `reason` - negative value indicates diverged, positive value converged, see `KSPConvergedReason` for the possible values

Options Database Key:
- `-ksp_converged_reason` - prints the reason to standard out when the solve ends

Level: intermediate

See also: `KSPConvergedReason`, `KSP`, `KSPSetConvergenceTest()`, `KSPConvergedDefault()`, `KSPSetTolerances()`,
`KSPConvergedReasonView()`, `KSPGetConvergedReasonString()`

# External Links
$(_doc_external("KSP/KSPGetConvergedReason"))
"""
function KSPGetConvergedReason(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetConvergedReason: no generated method for these argument types")
end

@for_petsc function KSPGetConvergedReason(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	reason_ = Ref{KSPConvergedReason}()

    @chk ccall(
               (:KSPGetConvergedReason, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPConvergedReason}),
               ksp, reason_,
              )

	reason = reason_[]

	return reason
end 

"""
	strreason::String = KSPGetConvergedReasonString(petsclib::PetscLibType, ksp::AbstractKSP) 
Return a human readable string for a `KSPConvergedReason`

Not Collective

Input Parameter:
- `ksp` - the `KSP` context

Output Parameter:
- `strreason` - a human readable string that describes ksp converged reason

Level: beginner

See also: `KSP`, `KSPGetConvergedReason()`

# External Links
$(_doc_external("KSP/KSPGetConvergedReasonString"))
"""
function KSPGetConvergedReasonString(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetConvergedReasonString: no generated method for these argument types")
end

@for_petsc function KSPGetConvergedReasonString(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	strreason_ = Ref{Ptr{Cchar}}()

    @chk ccall(
               (:KSPGetConvergedReasonString, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Ptr{Cchar}}),
               ksp, strreason_,
              )

	strreason = strreason_[] == C_NULL ? "" : unsafe_string(strreason_[])

	return strreason
end 

"""
	ctx::Ptr{Cvoid} = KSPGetConvergenceContext(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets the convergence context set with `KSPSetConvergenceTest()`.

Not Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameter:
- `ctx` - monitoring context

Level: advanced

See also: `KSP`, `KSPConvergedDefault()`, `KSPSetConvergenceTest()`, `KSPGetConvergenceTest()`

# External Links
$(_doc_external("KSP/KSPGetConvergenceContext"))
"""
function KSPGetConvergenceContext(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetConvergenceContext: no generated method for these argument types")
end

@for_petsc function KSPGetConvergenceContext(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	ctx_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:KSPGetConvergenceContext, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Cvoid}),
               ksp, ctx_,
              )

	ctx = ctx_[]

	return ctx
end 

"""
	converge::Ptr{Cvoid},ctx::Ptr{Cvoid},destroy::Ptr{Cvoid} = KSPGetConvergenceTest(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets the function to be used to determine convergence.

Logically Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameters:
- `converge` - pointer to convergence test function, see `KSPConvergenceTestFn`
- `ctx`      - context for private data for the convergence routine (may be `NULL`)
- `destroy`  - a routine for destroying the context (may be `NULL`)

Level: advanced

See also: `KSP`, `KSPConvergedDefault()`, `KSPGetConvergenceContext()`, `KSPSetTolerances()`, `KSPSetConvergenceTest()`, `KSPGetAndClearConvergenceTest()`

# External Links
$(_doc_external("KSP/KSPGetConvergenceTest"))
"""
function KSPGetConvergenceTest(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetConvergenceTest: no generated method for these argument types")
end

@for_petsc function KSPGetConvergenceTest(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	converge_ = Ref{Ptr{Cvoid}}()
	ctx_ = Ref{Ptr{Cvoid}}()
	destroy_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:KSPGetConvergenceTest, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Ptr{Cvoid}}, Ptr{Cvoid}, Ptr{Ptr{Cvoid}}),
               ksp, converge_, ctx_, destroy_,
              )

	converge = converge_[]
	ctx = ctx_[]
	destroy = destroy_[]

	return converge,ctx,destroy
end 

"""
	dm::PetscDM = KSPGetDM(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets the `DM` that may be used by some preconditioners and that may be used to construct the linear system

Not Collective

Input Parameter:
- `ksp` - the `KSP`

Output Parameter:
- `dm` - the `DM`

Level: intermediate

See also: `KSP`, `DM`, `KSPSetDM()`, `KSPSetDMActive()`

# External Links
$(_doc_external("KSP/KSPGetDM"))
"""
function KSPGetDM(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetDM: no generated method for these argument types")
end

@for_petsc function KSPGetDM(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	dm_ = Ref{CDM}()

    @chk ccall(
               (:KSPGetDM, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{CDM}),
               ksp, dm_,
              )

	dm = PetscDM(dm_[], petsclib; own = false)

	return dm
end 

"""
	scale::PetscBool = KSPGetDiagonalScale(petsclib::PetscLibType, ksp::AbstractKSP) 
Checks if `KSP` solver scales the matrix and right-hand side, that is if `KSPSetDiagonalScale()` has been called

Not Collective

Input Parameter:
- `ksp` - the `KSP` context

Output Parameter:
- `scale` - `PETSC_TRUE` or `PETSC_FALSE`

Level: intermediate

See also: `KSP`, `KSPSetDiagonalScale()`, `KSPSetDiagonalScaleFix()`

# External Links
$(_doc_external("KSP/KSPGetDiagonalScale"))
"""
function KSPGetDiagonalScale(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetDiagonalScale: no generated method for these argument types")
end

@for_petsc function KSPGetDiagonalScale(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	scale_ = Ref{PetscBool}()

    @chk ccall(
               (:KSPGetDiagonalScale, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{PetscBool}),
               ksp, scale_,
              )

	scale = scale_[]

	return scale
end 

"""
	fix::PetscBool = KSPGetDiagonalScaleFix(petsclib::PetscLibType, ksp::AbstractKSP) 
Determines if `KSP` diagonally scales the system back after solving. That is `KSPSetDiagonalScaleFix()` has been called

Not Collective

Input Parameter:
- `ksp` - the `KSP` context

Output Parameter:
- `fix` - `PETSC_TRUE` to scale back after the system solve, `PETSC_FALSE` to not
rescale (default)

Level: intermediate

See also: `KSPGetDiagonalScale()`, `KSPSetDiagonalScale()`, `KSPSetDiagonalScaleFix()`, `KSP`

# External Links
$(_doc_external("KSP/KSPGetDiagonalScaleFix"))
"""
function KSPGetDiagonalScaleFix(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetDiagonalScaleFix: no generated method for these argument types")
end

@for_petsc function KSPGetDiagonalScaleFix(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	fix_ = Ref{PetscBool}()

    @chk ccall(
               (:KSPGetDiagonalScaleFix, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{PetscBool}),
               ksp, fix_,
              )

	fix = fix_[]

	return fix
end 

"""
	a::Vector{PetscReal},na::PetscInt = KSPGetErrorHistory(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets the array used to hold the error history and the number of residuals it contains.

Not Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameters:
- `a`  - pointer to array to hold history (or `NULL`)
- `na` - number of used entries in a (or `NULL`)

Level: advanced

See also: `KSPSetErrorHistory()`, `KSPGetResidualHistory()`, `KSP`

# External Links
$(_doc_external("KSP/KSPGetErrorHistory"))
"""
function KSPGetErrorHistory(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetErrorHistory: no generated method for these argument types")
end

@for_petsc function KSPGetErrorHistory(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	a_ = Ref{Ptr{$PetscReal}}()
	na_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPGetErrorHistory, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Ptr{$PetscReal}}, Ptr{$PetscInt}),
               ksp, a_, na_,
              )

	na = na_[]
	a = a_[] == C_NULL ? $PetscReal[] : unsafe_wrap(Array, a_[], na; own = false)

	return a,na
end 

"""
	flag::PetscBool = KSPGetErrorIfNotConverged(petsclib::PetscLibType, ksp::AbstractKSP) 
Will `KSPSolve()` generate an error if the solver does not converge?

Not Collective

Input Parameter:
- `ksp` - iterative solver obtained from KSPCreate()

Output Parameter:
- `flag` - `PETSC_TRUE` if it will generate an error, else `PETSC_FALSE`

Level: intermediate

See also: `KSPSetErrorIfNotConverged()`, `KSP`

# External Links
$(_doc_external("KSP/KSPGetErrorIfNotConverged"))
"""
function KSPGetErrorIfNotConverged(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetErrorIfNotConverged: no generated method for these argument types")
end

@for_petsc function KSPGetErrorIfNotConverged(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	flag_ = Ref{PetscBool}()

    @chk ccall(
               (:KSPGetErrorIfNotConverged, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{PetscBool}),
               ksp, flag_,
              )

	flag = flag_[]

	return flag
end 

"""
	guess::KSPGuess = KSPGetGuess(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets the initial guess generator for the `KSP`.

Not Collective

Input Parameter:
- `ksp` - the Krylov context

Output Parameter:
- `guess` - the object

Level: developer

See also: `KSPGuess`, `KSP`, `KSPSetOptionsPrefix()`, `KSPAppendOptionsPrefix()`, `KSPSetUseFischerGuess()`, `KSPSetGuess()`

# External Links
$(_doc_external("KSP/KSPGetGuess"))
"""
function KSPGetGuess(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetGuess: no generated method for these argument types")
end

@for_petsc function KSPGetGuess(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	guess_ = Ref{KSPGuess}()

    @chk ccall(
               (:KSPGetGuess, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPGuess}),
               ksp, guess_,
              )

	guess = guess_[]

	return guess
end 

"""
	flag::PetscBool = KSPGetInitialGuessKnoll(petsclib::PetscLibType, ksp::AbstractKSP) 
Determines whether the `KSP` solver is using the Knoll trick (using PCApply(pc,b,...) to compute
the initial guess

Not Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameter:
- `flag` - `PETSC_TRUE` if using Knoll trick, else `PETSC_FALSE`

Level: advanced

See also: `KSPSetInitialGuessKnoll()`, `KSPSetInitialGuessNonzero()`, `KSPGetInitialGuessNonzero()`, `KSP`

# External Links
$(_doc_external("KSP/KSPGetInitialGuessKnoll"))
"""
function KSPGetInitialGuessKnoll(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetInitialGuessKnoll: no generated method for these argument types")
end

@for_petsc function KSPGetInitialGuessKnoll(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	flag_ = Ref{PetscBool}()

    @chk ccall(
               (:KSPGetInitialGuessKnoll, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{PetscBool}),
               ksp, flag_,
              )

	flag = flag_[]

	return flag
end 

"""
	flag::PetscBool = KSPGetInitialGuessNonzero(petsclib::PetscLibType, ksp::AbstractKSP) 
Determines whether the `KSP` solver is using
a zero initial guess.

Not Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameter:
- `flag` - `PETSC_TRUE` if guess is nonzero, else `PETSC_FALSE`

Level: intermediate

See also: `KSPSetInitialGuessNonzero()`, `KSP`

# External Links
$(_doc_external("KSP/KSPGetInitialGuessNonzero"))
"""
function KSPGetInitialGuessNonzero(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetInitialGuessNonzero: no generated method for these argument types")
end

@for_petsc function KSPGetInitialGuessNonzero(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	flag_ = Ref{PetscBool}()

    @chk ccall(
               (:KSPGetInitialGuessNonzero, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{PetscBool}),
               ksp, flag_,
              )

	flag = flag_[]

	return flag
end 

"""
	its::PetscInt = KSPGetIterationNumber(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets the current iteration number; if the `KSPSolve()` is complete, returns the number of iterations used.

Not Collective

Input Parameter:
- `ksp` - the iterative context

Output Parameter:
- `its` - number of iterations

Level: intermediate

See also: `KSP`, `KSPGetResidualNorm()`, `KSPBuildResidual()`, `KSPGetTotalIterations()`

# External Links
$(_doc_external("KSP/KSPGetIterationNumber"))
"""
function KSPGetIterationNumber(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetIterationNumber: no generated method for these argument types")
end

@for_petsc function KSPGetIterationNumber(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	its_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPGetIterationNumber, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscInt}),
               ksp, its_,
              )

	its = its_[]

	return its
end 

"""
	bs::PetscInt = KSPGetMatSolveBatchSize(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets the maximum number of columns treated simultaneously in `KSPMatSolve()`.

Input Parameter:
- `ksp` - iterative solver context

Output Parameter:
- `bs` - batch size

Level: advanced

See also: `KSPMatSolve()`, `KSPSetMatSolveBatchSize()`, `-mat_mumps_icntl_27`, `-matproduct_batch_size`

# External Links
$(_doc_external("KSP/KSPGetMatSolveBatchSize"))
"""
function KSPGetMatSolveBatchSize(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetMatSolveBatchSize: no generated method for these argument types")
end

@for_petsc function KSPGetMatSolveBatchSize(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	bs_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPGetMatSolveBatchSize, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscInt}),
               ksp, bs_,
              )

	bs = bs_[]

	return bs
end 

"""
	minit::PetscInt = KSPGetMinimumIterations(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets the minimum number of iterations to use, regardless of the tolerances, that was set with `KSPSetMinimumIterations()` or `-ksp_min_it`

Not Collective

Input Parameter:
- `ksp` - the Krylov subspace context

Output Parameter:
- `minit` - minimum number of iterations to use

Level: intermediate

See also: `KSPGetTolerances()`, `KSPConvergedDefault()`, `KSPSetConvergenceTest()`, `KSP`, `KSPSetTolerances()`, `KSPSetMinimumIterations()`

# External Links
$(_doc_external("KSP/KSPGetMinimumIterations"))
"""
function KSPGetMinimumIterations(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetMinimumIterations: no generated method for these argument types")
end

@for_petsc function KSPGetMinimumIterations(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	minit_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPGetMinimumIterations, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscInt}),
               ksp, minit_,
              )

	minit = minit_[]

	return minit
end 

"""
	ctx::Ptr{Cvoid} = KSPGetMonitorContext(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets the monitoring context, as set by `KSPMonitorSet()` for the FIRST monitor only.

Not Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameter:
- `ctx` - monitoring context

Level: intermediate

See also: `KSPMonitorResidual()`, `KSP`

# External Links
$(_doc_external("KSP/KSPGetMonitorContext"))
"""
function KSPGetMonitorContext(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetMonitorContext: no generated method for these argument types")
end

@for_petsc function KSPGetMonitorContext(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	ctx_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:KSPGetMonitorContext, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Cvoid}),
               ksp, ctx_,
              )

	ctx = ctx_[]

	return ctx
end 

"""
	level::PetscInt = KSPGetNestLevel(petsclib::PetscLibType, ksp::AbstractKSP) 
gets the amount of nesting the `KSP` has

Not Collective

Input Parameter:
- `ksp` - the `KSP`

Output Parameter:
- `level` - the nest level

Level: developer

See also: `KSPSetUp()`, `KSPSolve()`, `KSPDestroy()`, `KSP`, `KSPGMRES`, `KSPType`, `KSPSetNestLevel()`, `PCSetKSPNestLevel()`, `PCGetKSPNestLevel()`

# External Links
$(_doc_external("KSP/KSPGetNestLevel"))
"""
function KSPGetNestLevel(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetNestLevel: no generated method for these argument types")
end

@for_petsc function KSPGetNestLevel(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	level_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPGetNestLevel, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscInt}),
               ksp, level_,
              )

	level = level_[]

	return level
end 

"""
	normtype::KSPNormType = KSPGetNormType(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets the `KSPNormType` that is used for convergence testing during `KSPSolve()` for this `KSP` context

Not Collective

Input Parameter:
- `ksp` - Krylov solver context

Output Parameter:
- `normtype` - the `KSPNormType` that is used for convergence testing

Level: advanced

See also: `KSPNormType`, `KSPSetNormType()`, `KSPConvergedSkip()`

# External Links
$(_doc_external("KSP/KSPGetNormType"))
"""
function KSPGetNormType(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetNormType: no generated method for these argument types")
end

@for_petsc function KSPGetNormType(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	normtype_ = Ref{KSPNormType}()

    @chk ccall(
               (:KSPGetNormType, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPNormType}),
               ksp, normtype_,
              )

	normtype = normtype_[]

	return normtype
end 

"""
	Amat::PetscMat,Pmat::PetscMat = KSPGetOperators(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets the matrix associated with the linear system
and a (possibly) different one used to construct the preconditioner from the `KSP` context

Collective

Input Parameter:
- `ksp` - the `KSP` context

Output Parameters:
- `Amat` - the matrix that defines the linear system
- `Pmat` - the matrix to be used in constructing the preconditioner, usually the same as `Amat`.

Level: intermediate

See also: `KSP`, `KSPSolve()`, `KSPGetPC()`, `PCSetOperators()`, `KSPSetOperators()`, `KSPGetOperatorsSet()`

# External Links
$(_doc_external("KSP/KSPGetOperators"))
"""
function KSPGetOperators(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetOperators: no generated method for these argument types")
end

@for_petsc function KSPGetOperators(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	Amat_ = Ref{CMat}()
	Pmat_ = Ref{CMat}()

    @chk ccall(
               (:KSPGetOperators, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{CMat}, Ptr{CMat}),
               ksp, Amat_, Pmat_,
              )

	Amat = PetscMat(Amat_[], petsclib; own = false)
	Pmat = PetscMat(Pmat_[], petsclib; own = false)

	return Amat,Pmat
end 

"""
	mat::PetscBool,pmat::PetscBool = KSPGetOperatorsSet(petsclib::PetscLibType, ksp::AbstractKSP) 
Determines if the matrix associated with the linear system and
possibly a different one from which the preconditioner will be built have been set in the `KSP` with `KSPSetOperators()`

Not Collective, though the results on all processes will be the same

Input Parameter:
- `ksp` - the `KSP` context

Output Parameters:
- `mat`  - the matrix associated with the linear system was set
- `pmat` - matrix from which the preconditioner will be built, usually the same as `mat` was set

Level: intermediate

See also: `KSP`, `PCSetOperators()`, `KSPGetOperators()`, `KSPSetOperators()`, `PCGetOperators()`, `PCGetOperatorsSet()`

# External Links
$(_doc_external("KSP/KSPGetOperatorsSet"))
"""
function KSPGetOperatorsSet(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetOperatorsSet: no generated method for these argument types")
end

@for_petsc function KSPGetOperatorsSet(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	mat_ = Ref{PetscBool}()
	pmat_ = Ref{PetscBool}()

    @chk ccall(
               (:KSPGetOperatorsSet, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{PetscBool}, Ptr{PetscBool}),
               ksp, mat_, pmat_,
              )

	mat = mat_[]
	pmat = pmat_[]

	return mat,pmat
end 

"""
	prefix::String = KSPGetOptionsPrefix(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets the prefix used for searching for all
`KSP` options in the database.

Not Collective

Input Parameter:
- `ksp` - the Krylov context

Output Parameter:
- `prefix` - pointer to the prefix string used is returned

Level: advanced

See also: `KSP`, `KSPSetFromOptions()`, `KSPSetOptionsPrefix()`, `KSPAppendOptionsPrefix()`

# External Links
$(_doc_external("KSP/KSPGetOptionsPrefix"))
"""
function KSPGetOptionsPrefix(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetOptionsPrefix: no generated method for these argument types")
end

@for_petsc function KSPGetOptionsPrefix(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	prefix_ = Ref{Ptr{Cchar}}()

    @chk ccall(
               (:KSPGetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Ptr{Cchar}}),
               ksp, prefix_,
              )

	prefix = prefix_[] == C_NULL ? "" : unsafe_string(prefix_[])

	return prefix
end 

"""
	pc::PC = KSPGetPC(petsclib::PetscLibType, ksp::AbstractKSP) 
Returns a pointer to the preconditioner context with the `KSP`

Not Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameter:
- `pc` - preconditioner context

Level: beginner

See also: `KSPSetPC()`, `KSP`, `PC`

# External Links
$(_doc_external("KSP/KSPGetPC"))
"""
function KSPGetPC(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetPC: no generated method for these argument types")
end

@for_petsc function KSPGetPC(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	pc_ = Ref{CPC}()

    @chk ccall(
               (:KSPGetPC, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{CPC}),
               ksp, pc_,
              )

	pc = PC(pc_[], petsclib; own = false)

	return pc
end 

"""
	side::PCSide = KSPGetPCSide(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets the preconditioning side.

Not Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameter:
- `side` - the preconditioning side, where side is one of
``
PC_LEFT      - left preconditioning (default)
PC_RIGHT     - right preconditioning
PC_SYMMETRIC - symmetric preconditioning
``

Level: intermediate

See also: `KSPSetPCSide()`, `KSP`

# External Links
$(_doc_external("KSP/KSPGetPCSide"))
"""
function KSPGetPCSide(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetPCSide: no generated method for these argument types")
end

@for_petsc function KSPGetPCSide(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	side_ = Ref{PCSide}()

    @chk ccall(
               (:KSPGetPCSide, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{PCSide}),
               ksp, side_,
              )

	side = side_[]

	return side
end 

"""
	a::Vector{PetscReal},na::PetscInt = KSPGetResidualHistory(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets the array used to hold the residual history and the number of residuals it contains.

Not Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameters:
- `a`  - pointer to array to hold history (or `NULL`)
- `na` - number of used entries in a (or `NULL`). Note this has different meanings depending on the `reset` argument to `KSPSetResidualHistory()`

Level: advanced

See also: `KSPSetResidualHistory()`, `KSP`, `KSPGetIterationNumber()`, `KSPSTCG`, `KSPBCGSL`

# External Links
$(_doc_external("KSP/KSPGetResidualHistory"))
"""
function KSPGetResidualHistory(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetResidualHistory: no generated method for these argument types")
end

@for_petsc function KSPGetResidualHistory(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	a_ = Ref{Ptr{$PetscReal}}()
	na_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPGetResidualHistory, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Ptr{$PetscReal}}, Ptr{$PetscInt}),
               ksp, a_, na_,
              )

	na = na_[]
	a = a_[] == C_NULL ? $PetscReal[] : unsafe_wrap(Array, a_[], na; own = false)

	return a,na
end 

"""
	rnorm::PetscReal = KSPGetResidualNorm(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets the last (possibly approximate and/or preconditioned) residual norm that has been computed.

Not Collective

Input Parameter:
- `ksp` - the iterative context

Output Parameter:
- `rnorm` - residual norm

Level: intermediate

See also: `KSP`, `KSPSetNormType()`, `KSPBuildResidual()`, `KSPNormType`

# External Links
$(_doc_external("KSP/KSPGetResidualNorm"))
"""
function KSPGetResidualNorm(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetResidualNorm: no generated method for these argument types")
end

@for_petsc function KSPGetResidualNorm(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	rnorm_ = Ref{$PetscReal}()

    @chk ccall(
               (:KSPGetResidualNorm, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscReal}),
               ksp, rnorm_,
              )

	rnorm = rnorm_[]

	return rnorm
end 

"""
	flag::PetscBool = KSPGetReusePreconditioner(petsclib::PetscLibType, ksp::AbstractKSP) 
Determines if the `KSP` reuses the current preconditioner even if the `Mat` operator in the `KSP` has changed.

Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameter:
- `flag` - the boolean flag indicating if the current preconditioner should be reused

Level: intermediate

See also: `KSPCreate()`, `KSPSolve()`, `KSPDestroy()`, `KSPSetReusePreconditioner()`, `KSP`

# External Links
$(_doc_external("KSP/KSPGetReusePreconditioner"))
"""
function KSPGetReusePreconditioner(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetReusePreconditioner: no generated method for these argument types")
end

@for_petsc function KSPGetReusePreconditioner(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	flag_ = Ref{PetscBool}()

    @chk ccall(
               (:KSPGetReusePreconditioner, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{PetscBool}),
               ksp, flag_,
              )

	flag = flag_[]

	return flag
end 

"""
	r::PetscVec = KSPGetRhs(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets the right-hand-side vector for the linear system to
be solved.

Not Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameter:
- `r` - right-hand-side vector

Level: developer

See also: `KSPGetSolution()`, `KSPSolve()`, `KSP`

# External Links
$(_doc_external("KSP/KSPGetRhs"))
"""
function KSPGetRhs(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetRhs: no generated method for these argument types")
end

@for_petsc function KSPGetRhs(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	r_ = Ref{CVec}()

    @chk ccall(
               (:KSPGetRhs, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{CVec}),
               ksp, r_,
              )

	r = PetscVec(r_[], petsclib; own = false)

	return r
end 

"""
	v::PetscVec = KSPGetSolution(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets the location of the solution for the
linear system to be solved.

Not Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameter:
- `v` - solution vector

Level: developer

See also: `KSPGetRhs()`, `KSPBuildSolution()`, `KSPSolve()`, `KSP`

# External Links
$(_doc_external("KSP/KSPGetSolution"))
"""
function KSPGetSolution(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetSolution: no generated method for these argument types")
end

@for_petsc function KSPGetSolution(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	v_ = Ref{CVec}()

    @chk ccall(
               (:KSPGetSolution, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{CVec}),
               ksp, v_,
              )

	v = PetscVec(v_[], petsclib; own = false)

	return v
end 

"""
	rtol::PetscReal,abstol::PetscReal,dtol::PetscReal,maxits::PetscInt = KSPGetTolerances(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets the relative, absolute, divergence, and maximum
iteration tolerances used by the default `KSP` convergence tests.

Not Collective

Input Parameter:
- `ksp` - the Krylov subspace context

Output Parameters:
- `rtol`   - the relative convergence tolerance
- `abstol` - the absolute convergence tolerance
- `dtol`   - the divergence tolerance
- `maxits` - maximum number of iterations

Level: intermediate

See also: `KSPSetTolerances()`, `KSP`, `KSPSetMinimumIterations()`, `KSPGetMinimumIterations()`

# External Links
$(_doc_external("KSP/KSPGetTolerances"))
"""
function KSPGetTolerances(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetTolerances: no generated method for these argument types")
end

@for_petsc function KSPGetTolerances(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	rtol_ = Ref{$PetscReal}()
	abstol_ = Ref{$PetscReal}()
	dtol_ = Ref{$PetscReal}()
	maxits_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPGetTolerances, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscInt}),
               ksp, rtol_, abstol_, dtol_, maxits_,
              )

	rtol = rtol_[]
	abstol = abstol_[]
	dtol = dtol_[]
	maxits = maxits_[]

	return rtol,abstol,dtol,maxits
end 

"""
	its::PetscInt = KSPGetTotalIterations(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets the total number of iterations this `KSP` object has performed since was created, counted over all linear solves

Not Collective

Input Parameter:
- `ksp` - the iterative context

Output Parameter:
- `its` - total number of iterations

Level: intermediate

See also: `KSP`, `KSPBuildResidual()`, `KSPGetResidualNorm()`, `KSPGetIterationNumber()`

# External Links
$(_doc_external("KSP/KSPGetTotalIterations"))
"""
function KSPGetTotalIterations(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetTotalIterations: no generated method for these argument types")
end

@for_petsc function KSPGetTotalIterations(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	its_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPGetTotalIterations, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscInt}),
               ksp, its_,
              )

	its = its_[]

	return its
end 

"""
	type::String = KSPGetType(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets the `KSP` type as a string from the `KSP` object.

Not Collective

Input Parameter:
- `ksp` - Krylov context

Output Parameter:
- `type` - name of the `KSP` method

Level: intermediate

See also: `KSPType`, `KSP`, `KSPSetType()`, `PetscObjectTypeCompare()`, `PetscObjectTypeCompareAny()`

# External Links
$(_doc_external("KSP/KSPGetType"))
"""
function KSPGetType(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPGetType: no generated method for these argument types")
end

@for_petsc function KSPGetType(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	type_ = Ref{KSPType}()

    @chk ccall(
               (:KSPGetType, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPType}),
               ksp, type_,
              )

	type = type_[] == C_NULL ? "" : unsafe_string(type_[])

	return type
end 

"""
	KSPHPDDMGetDeflationMat(petsclib::PetscLibType, ksp::AbstractKSP, U::AbstractPetscMat) 

# External Links
$(_doc_external("KSP/KSPHPDDMGetDeflationMat"))
"""
function KSPHPDDMGetDeflationMat(petsclib::PetscLibType, ksp::AbstractKSP, U::AbstractPetscMat)
    error("KSPHPDDMGetDeflationMat: no generated method for these argument types")
end

@for_petsc function KSPHPDDMGetDeflationMat(petsclib::$UnionPetscLib, ksp::AbstractKSP, U::AbstractPetscMat )
	U_ = Ref(U.ptr)

    @chk ccall(
               (:KSPHPDDMGetDeflationMat, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{CMat}),
               ksp, U_,
              )

	U.ptr = U_[]

	return nothing
end 

"""
	type::KSPHPDDMType = KSPHPDDMGetType(petsclib::PetscLibType, ksp::AbstractKSP) 

# External Links
$(_doc_external("KSP/KSPHPDDMGetType"))
"""
function KSPHPDDMGetType(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPHPDDMGetType: no generated method for these argument types")
end

@for_petsc function KSPHPDDMGetType(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	type_ = Ref{KSPHPDDMType}()

    @chk ccall(
               (:KSPHPDDMGetType, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPHPDDMType}),
               ksp, type_,
              )

	type = type_[]

	return type
end 

"""
	KSPHPDDMSetDeflationMat(petsclib::PetscLibType, ksp::AbstractKSP, U::AbstractPetscMat) 

# External Links
$(_doc_external("KSP/KSPHPDDMSetDeflationMat"))
"""
function KSPHPDDMSetDeflationMat(petsclib::PetscLibType, ksp::AbstractKSP, U::AbstractPetscMat)
    error("KSPHPDDMSetDeflationMat: no generated method for these argument types")
end

@for_petsc function KSPHPDDMSetDeflationMat(petsclib::$UnionPetscLib, ksp::AbstractKSP, U::AbstractPetscMat )

    @chk ccall(
               (:KSPHPDDMSetDeflationMat, $petsc_library),
               PetscErrorCode,
               (CKSP, CMat),
               ksp, U,
              )


	return nothing
end 

"""
	KSPHPDDMSetType(petsclib::PetscLibType, ksp::AbstractKSP, type::KSPHPDDMType) 

# External Links
$(_doc_external("KSP/KSPHPDDMSetType"))
"""
function KSPHPDDMSetType(petsclib::PetscLibType, ksp::AbstractKSP, type::KSPHPDDMType)
    error("KSPHPDDMSetType: no generated method for these argument types")
end

@for_petsc function KSPHPDDMSetType(petsclib::$UnionPetscLib, ksp::AbstractKSP, type::KSPHPDDMType )

    @chk ccall(
               (:KSPHPDDMSetType, $petsc_library),
               PetscErrorCode,
               (CKSP, KSPHPDDMType),
               ksp, type,
              )


	return nothing
end 

"""
	KSPInitialResidual(petsclib::PetscLibType, ksp::AbstractKSP, vsoln::AbstractPetscVec, vt1::AbstractPetscVec, vt2::AbstractPetscVec, vres::AbstractPetscVec, vb::AbstractPetscVec) 
Computes the residual. Either b - A*C*u = b - A*x with right
preconditioning or C*(b - A*x) with left preconditioning; the latter
residual is often called the "preconditioned residual".

Collective

Input Parameters:
- `ksp`   - the `KSP` solver object
- `vsoln` - solution to use in computing residual
- `vt1`   - temporary work vector
- `vt2`   - temporary work vector
- `vb`    - right-hand-side vector

Output Parameter:
- `vres` - calculated residual

Level: developer

See also: `KSP`, `KSPSolve()`, `KSPMonitor()`

# External Links
$(_doc_external("KSP/KSPInitialResidual"))
"""
function KSPInitialResidual(petsclib::PetscLibType, ksp::AbstractKSP, vsoln::AbstractPetscVec, vt1::AbstractPetscVec, vt2::AbstractPetscVec, vres::AbstractPetscVec, vb::AbstractPetscVec)
    error("KSPInitialResidual: no generated method for these argument types")
end

@for_petsc function KSPInitialResidual(petsclib::$UnionPetscLib, ksp::AbstractKSP, vsoln::AbstractPetscVec, vt1::AbstractPetscVec, vt2::AbstractPetscVec, vres::AbstractPetscVec, vb::AbstractPetscVec )

    @chk ccall(
               (:KSPInitialResidual, $petsc_library),
               PetscErrorCode,
               (CKSP, CVec, CVec, CVec, CVec, CVec),
               ksp, vsoln, vt1, vt2, vres, vb,
              )


	return nothing
end 

"""
	KSPInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `KSP` package. It is called
from `PetscDLLibraryRegister_petscksp()` when using dynamic libraries, and on the first call to `KSPCreate()`
when using shared or static libraries.

Level: developer

See also: `PetscInitialize()`, `KSPFinalizePackage()`

# External Links
$(_doc_external("KSP/KSPInitializePackage"))
"""
function KSPInitializePackage(petsclib::PetscLibType)
    error("KSPInitializePackage: no generated method for these argument types")
end

@for_petsc function KSPInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:KSPInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	KSPLGMRESSetAugDim(petsclib::PetscLibType, ksp::AbstractKSP, dim::PetscInt) 
Set the number of error approximations to include in the approximation space (default is 2) for `KSPLGMRES`

Collective

Input Parameters:
- `ksp` - the `KSP` context
- `dim` - the number of vectors to use

Options Database Key:
- `-ksp_lgmres_augment dim` - the number of error approximations to include

Level: intermediate

See also: `KSPLGMRES`, `KSPLGMRESSetConstant()`

# External Links
$(_doc_external("KSP/KSPLGMRESSetAugDim"))
"""
function KSPLGMRESSetAugDim(petsclib::PetscLibType, ksp::AbstractKSP, dim::Integer)
    error("KSPLGMRESSetAugDim: no generated method for these argument types")
end

@for_petsc function KSPLGMRESSetAugDim(petsclib::$UnionPetscLib, ksp::AbstractKSP, dim::$PetscInt )

    @chk ccall(
               (:KSPLGMRESSetAugDim, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, dim,
              )


	return nothing
end 

"""
	KSPLGMRESSetConstant(petsclib::PetscLibType, ksp::AbstractKSP) 
keep the error approximation space a constant size for every restart cycle

Collective

Input Parameters:
- `ksp` - the `KSP` context

Options Database Key:
- `-ksp_lgmres_constant` - set the size to be constant

Level: intermediate

See also: `KSPLGMRES`, `KSPLGMRESSetAugDim()`

# External Links
$(_doc_external("KSP/KSPLGMRESSetConstant"))
"""
function KSPLGMRESSetConstant(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPLGMRESSetConstant: no generated method for these argument types")
end

@for_petsc function KSPLGMRESSetConstant(petsclib::$UnionPetscLib, ksp::AbstractKSP )

    @chk ccall(
               (:KSPLGMRESSetConstant, $petsc_library),
               PetscErrorCode,
               (CKSP,),
               ksp,
              )


	return nothing
end 

"""
	reason::KSPConvergedReason = KSPLSQRConvergedDefault(petsclib::PetscLibType, ksp::AbstractKSP, n::PetscInt, rnorm::PetscReal, ctx::Ptr{Cvoid}) 
Determines convergence of the `KSPLSQR` Krylov method, including a check on the residual norm of the normal equations.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm residual value (may be estimated)
- `ctx`   - convergence context which must have been created by `KSPConvergedDefaultCreate()`

Output Parameter:
- `reason` - the convergence reason

Level: advanced

See also: `KSPLSQR`, `KSPSetConvergenceTest()`, `KSPSetTolerances()`, `KSPConvergedSkip()`, `KSPConvergedReason`, `KSPGetConvergedReason()`,
`KSPConvergedDefaultSetUIRNorm()`, `KSPConvergedDefaultSetUMIRNorm()`, `KSPConvergedDefaultCreate()`, `KSPConvergedDefaultDestroy()`,
`KSPConvergedDefault()`, `KSPLSQRGetNorms()`, `KSPLSQRSetExactMatNorm()`

# External Links
$(_doc_external("KSP/KSPLSQRConvergedDefault"))
"""
function KSPLSQRConvergedDefault(petsclib::PetscLibType, ksp::AbstractKSP, n::Integer, rnorm::Real, ctx::Ptr{Cvoid})
    error("KSPLSQRConvergedDefault: no generated method for these argument types")
end

@for_petsc function KSPLSQRConvergedDefault(petsclib::$UnionPetscLib, ksp::AbstractKSP, n::$PetscInt, rnorm::$PetscReal, ctx::Ptr{Cvoid} )
	reason_ = Ref{KSPConvergedReason}()

    @chk ccall(
               (:KSPLSQRConvergedDefault, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{KSPConvergedReason}, Ptr{Cvoid}),
               ksp, n, rnorm, reason_, ctx,
              )

	reason = reason_[]

	return reason
end 

"""
	arnorm::PetscReal,anorm::PetscReal = KSPLSQRGetNorms(petsclib::PetscLibType, ksp::AbstractKSP) 
Get the norm estimates that `KSPLSQR` computes internally during `KSPSolve()`.

Not Collective

Input Parameter:
- `ksp` - iterative context

Output Parameters:
- `arnorm` - good estimate of \\|(A*Pmat^{-T})*r\\|, where r = A x - b, used in specific stopping criterion
- `anorm`  - poor estimate of \\|A*Pmat^{-T}\\|_{frobenius} used in specific stopping criterion

Level: intermediate

See also: `KSPSolve()`, `KSPLSQR`, `KSPLSQRSetExactMatNorm()`

# External Links
$(_doc_external("KSP/KSPLSQRGetNorms"))
"""
function KSPLSQRGetNorms(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPLSQRGetNorms: no generated method for these argument types")
end

@for_petsc function KSPLSQRGetNorms(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	arnorm_ = Ref{$PetscReal}()
	anorm_ = Ref{$PetscReal}()

    @chk ccall(
               (:KSPLSQRGetNorms, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscReal}, Ptr{$PetscReal}),
               ksp, arnorm_, anorm_,
              )

	arnorm = arnorm_[]
	anorm = anorm_[]

	return arnorm,anorm
end 

"""
	se::PetscVec = KSPLSQRGetStandardErrorVec(petsclib::PetscLibType, ksp::AbstractKSP) 
Get vector of standard error estimates.
Only available if -ksp_lsqr_set_standard_error was set to true
or `KSPLSQRSetComputeStandardErrorVec`(ksp, `PETSC_TRUE`) was called.
Otherwise returns `NULL`.

Not Collective

Input Parameter:
- `ksp` - iterative context

Output Parameter:
- `se` - vector of standard estimates

Level: intermediate

See also: `KSPSolve()`, `KSPLSQR`, `KSPLSQRSetComputeStandardErrorVec()`

# External Links
$(_doc_external("KSP/KSPLSQRGetStandardErrorVec"))
"""
function KSPLSQRGetStandardErrorVec(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPLSQRGetStandardErrorVec: no generated method for these argument types")
end

@for_petsc function KSPLSQRGetStandardErrorVec(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	se_ = Ref{CVec}()

    @chk ccall(
               (:KSPLSQRGetStandardErrorVec, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{CVec}),
               ksp, se_,
              )

	se = PetscVec(se_[], petsclib; own = false)

	return se
end 

"""
	KSPLSQRMonitorResidual(petsclib::PetscLibType, ksp::AbstractKSP, n::PetscInt, rnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 
Prints the residual norm, as well as the normal equation residual norm, at each iteration of an iterative solver for the `KSPLSQR` solver

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-ksp_lsqr_monitor` - Activates `KSPLSQRMonitorResidual()`

Level: intermediate

See also: `KSPLSQR`, `KSPMonitorSet()`, `KSPMonitorResidual()`, `KSPMonitorTrueResidualMaxNorm()`, `KSPLSQRMonitorResidualDrawLG()`

# External Links
$(_doc_external("KSP/KSPLSQRMonitorResidual"))
"""
function KSPLSQRMonitorResidual(petsclib::PetscLibType, ksp::AbstractKSP, n::Integer, rnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("KSPLSQRMonitorResidual: no generated method for these argument types")
end

@for_petsc function KSPLSQRMonitorResidual(petsclib::$UnionPetscLib, ksp::AbstractKSP, n::$PetscInt, rnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:KSPLSQRMonitorResidual, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	KSPLSQRMonitorResidualDrawLG(petsclib::PetscLibType, ksp::AbstractKSP, n::PetscInt, rnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 
Plots the true residual norm at each iteration of an iterative solver for the `KSPLSQR` solver

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-ksp_lsqr_monitor draw::draw_lg` - Activates `KSPMonitorTrueResidualDrawLG()`

Level: intermediate

See also: `KSPLSQR`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `KSPLSQRMonitorResidual()`, `KSPLSQRMonitorResidualDrawLGCreate()`

# External Links
$(_doc_external("KSP/KSPLSQRMonitorResidualDrawLG"))
"""
function KSPLSQRMonitorResidualDrawLG(petsclib::PetscLibType, ksp::AbstractKSP, n::Integer, rnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("KSPLSQRMonitorResidualDrawLG: no generated method for these argument types")
end

@for_petsc function KSPLSQRMonitorResidualDrawLG(petsclib::$UnionPetscLib, ksp::AbstractKSP, n::$PetscInt, rnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:KSPLSQRMonitorResidualDrawLG, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	vf::Ptr{PetscViewerAndFormat} = KSPLSQRMonitorResidualDrawLGCreate(petsclib::PetscLibType, viewer::PetscViewer, format::PetscViewerFormat, ctx::Ptr{Cvoid}) 
Creates the line graph object for the `KSPLSQR` residual and normal equation residual norm

Collective

Input Parameters:
- `viewer` - The `PetscViewer`
- `format` - The viewer format
- `ctx`    - An optional application context

Output Parameter:
- `vf` - The `PetscViewerAndFormat`

Level: intermediate

See also: `KSPLSQR`, `KSPMonitorSet()`, `KSPLSQRMonitorResidual()`, `KSPLSQRMonitorResidualDrawLG()`

# External Links
$(_doc_external("KSP/KSPLSQRMonitorResidualDrawLGCreate"))
"""
function KSPLSQRMonitorResidualDrawLGCreate(petsclib::PetscLibType, viewer::PetscViewer, format::PetscViewerFormat, ctx::Ptr{Cvoid})
    error("KSPLSQRMonitorResidualDrawLGCreate: no generated method for these argument types")
end

@for_petsc function KSPLSQRMonitorResidualDrawLGCreate(petsclib::$UnionPetscLib, viewer::PetscViewer, format::PetscViewerFormat, ctx::Ptr{Cvoid} )
	vf_ = Ref{Ptr{PetscViewerAndFormat}}()

    @chk ccall(
               (:KSPLSQRMonitorResidualDrawLGCreate, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscViewerFormat, Ptr{Cvoid}, Ptr{Ptr{PetscViewerAndFormat}}),
               viewer, format, ctx, vf_,
              )

	vf = vf_[]

	return vf
end 

"""
	KSPLSQRSetComputeStandardErrorVec(petsclib::PetscLibType, ksp::AbstractKSP, flg::PetscBool) 
Compute a vector of standard error estimates during `KSPSolve()` for  `KSPLSQR`.

Logically Collective

Input Parameters:
- `ksp` - iterative context
- `flg` - compute the vector of standard estimates or not

Level: intermediate

See also: `KSPSolve()`, `KSPLSQR`, `KSPLSQRGetStandardErrorVec()`

# External Links
$(_doc_external("KSP/KSPLSQRSetComputeStandardErrorVec"))
"""
function KSPLSQRSetComputeStandardErrorVec(petsclib::PetscLibType, ksp::AbstractKSP, flg::PetscBool)
    error("KSPLSQRSetComputeStandardErrorVec: no generated method for these argument types")
end

@for_petsc function KSPLSQRSetComputeStandardErrorVec(petsclib::$UnionPetscLib, ksp::AbstractKSP, flg::PetscBool )

    @chk ccall(
               (:KSPLSQRSetComputeStandardErrorVec, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flg,
              )


	return nothing
end 

"""
	KSPLSQRSetExactMatNorm(petsclib::PetscLibType, ksp::AbstractKSP, flg::PetscBool) 
Compute exact matrix norm instead of iteratively refined estimate.

Not Collective

Input Parameters:
- `ksp` - iterative context
- `flg` - compute exact matrix norm or not

Level: intermediate

See also: `KSPSolve()`, `KSPLSQR`, `KSPLSQRGetNorms()`, `KSPLSQRConvergedDefault()`

# External Links
$(_doc_external("KSP/KSPLSQRSetExactMatNorm"))
"""
function KSPLSQRSetExactMatNorm(petsclib::PetscLibType, ksp::AbstractKSP, flg::PetscBool)
    error("KSPLSQRSetExactMatNorm: no generated method for these argument types")
end

@for_petsc function KSPLSQRSetExactMatNorm(petsclib::$UnionPetscLib, ksp::AbstractKSP, flg::PetscBool )

    @chk ccall(
               (:KSPLSQRSetExactMatNorm, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flg,
              )


	return nothing
end 

"""
	KSPLoad(petsclib::PetscLibType, newdm::AbstractKSP, viewer::PetscViewer) 
Loads a `KSP` that has been stored in a `PETSCVIEWERBINARY`  with `KSPView()`.

Collective

Input Parameters:
- `newdm`  - the newly loaded `KSP`, this needs to have been created with `KSPCreate()` or
some related function before a call to `KSPLoad()`.
- `viewer` - binary file viewer, obtained from `PetscViewerBinaryOpen()`

Level: intermediate

See also: `KSP`, `PetscViewerBinaryOpen()`, `KSPView()`, `MatLoad()`, `VecLoad()`

# External Links
$(_doc_external("KSP/KSPLoad"))
"""
function KSPLoad(petsclib::PetscLibType, newdm::AbstractKSP, viewer::PetscViewer)
    error("KSPLoad: no generated method for these argument types")
end

@for_petsc function KSPLoad(petsclib::$UnionPetscLib, newdm::AbstractKSP, viewer::PetscViewer )

    @chk ccall(
               (:KSPLoad, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscViewer),
               newdm, viewer,
              )


	return nothing
end 

"""
	qlp::PetscBool = KSPMINRESGetUseQLP(petsclib::PetscLibType, ksp::AbstractKSP) 
Get the flag that indicates if the QLP variant is being used

Logically Collective

Input Parameter:
- `ksp` - the iterative context

Output Parameter:
- `qlp` - a Boolean indicating if the QLP variant is used

Level: beginner

See also: `KSP`, `KSPMINRES`, `KSPMINRESSetUseQLP()`

# External Links
$(_doc_external("KSP/KSPMINRESGetUseQLP"))
"""
function KSPMINRESGetUseQLP(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPMINRESGetUseQLP: no generated method for these argument types")
end

@for_petsc function KSPMINRESGetUseQLP(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	qlp_ = Ref{PetscBool}()

    @chk ccall(
               (:KSPMINRESGetUseQLP, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{PetscBool}),
               ksp, qlp_,
              )

	qlp = qlp_[]

	return qlp
end 

"""
	KSPMINRESSetRadius(petsclib::PetscLibType, ksp::AbstractKSP, radius::PetscReal) 
Set the maximum solution norm allowed for use with trust region methods

Logically Collective

Input Parameters:
- `ksp`    - the iterative context
- `radius` - the value

Level: beginner

Options Database Key:
- `-ksp_minres_radius radius` - maximum allowed solution norm

See also: `KSP`, `KSPMINRES`, `KSPMINRESSetUseQLP()`

# External Links
$(_doc_external("KSP/KSPMINRESSetRadius"))
"""
function KSPMINRESSetRadius(petsclib::PetscLibType, ksp::AbstractKSP, radius::Real)
    error("KSPMINRESSetRadius: no generated method for these argument types")
end

@for_petsc function KSPMINRESSetRadius(petsclib::$UnionPetscLib, ksp::AbstractKSP, radius::$PetscReal )

    @chk ccall(
               (:KSPMINRESSetRadius, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscReal),
               ksp, radius,
              )


	return nothing
end 

"""
	KSPMINRESSetUseQLP(petsclib::PetscLibType, ksp::AbstractKSP, qlp::PetscBool) 
Use the QLP variant of `KSPMINRES`

Logically Collective

Input Parameters:
- `ksp` - the iterative context
- `qlp` - a Boolean indicating if the QLP variant should be used

Level: beginner

See also: `KSP`, `KSPMINRES`, `KSPMINRESGetUseQLP()`

# External Links
$(_doc_external("KSP/KSPMINRESSetUseQLP"))
"""
function KSPMINRESSetUseQLP(petsclib::PetscLibType, ksp::AbstractKSP, qlp::PetscBool)
    error("KSPMINRESSetUseQLP: no generated method for these argument types")
end

@for_petsc function KSPMINRESSetUseQLP(petsclib::$UnionPetscLib, ksp::AbstractKSP, qlp::PetscBool )

    @chk ccall(
               (:KSPMINRESSetUseQLP, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, qlp,
              )


	return nothing
end 

"""
	KSPMatSolve(petsclib::PetscLibType, ksp::AbstractKSP, B::AbstractPetscMat, X::AbstractPetscMat) 
Solves a linear system with multiple right-hand sides stored as a `MATDENSE`.

Input Parameters:
- `ksp` - iterative solver
- `B`   - block of right-hand sides

Output Parameter:
- `X` - block of solutions

Level: intermediate

See also: `KSPSolve()`, `MatMatSolve()`, `KSPMatSolveTranspose()`, `MATDENSE`, `KSPHPDDM`, `PCBJACOBI`, `PCASM`, `KSPSetMatSolveBatchSize()`

# External Links
$(_doc_external("KSP/KSPMatSolve"))
"""
function KSPMatSolve(petsclib::PetscLibType, ksp::AbstractKSP, B::AbstractPetscMat, X::AbstractPetscMat)
    error("KSPMatSolve: no generated method for these argument types")
end

@for_petsc function KSPMatSolve(petsclib::$UnionPetscLib, ksp::AbstractKSP, B::AbstractPetscMat, X::AbstractPetscMat )

    @chk ccall(
               (:KSPMatSolve, $petsc_library),
               PetscErrorCode,
               (CKSP, CMat, CMat),
               ksp, B, X,
              )


	return nothing
end 

"""
	KSPMatSolveTranspose(petsclib::PetscLibType, ksp::AbstractKSP, B::AbstractPetscMat, X::AbstractPetscMat) 
Solves a linear system with the transposed matrix with multiple right-hand sides stored as a `MATDENSE`.

Input Parameters:
- `ksp` - iterative solver
- `B`   - block of right-hand sides

Output Parameter:
- `X` - block of solutions

Level: intermediate

See also: `KSPSolveTranspose()`, `MatMatTransposeSolve()`, `KSPMatSolve()`, `MATDENSE`, `KSPHPDDM`, `PCBJACOBI`, `PCASM`

# External Links
$(_doc_external("KSP/KSPMatSolveTranspose"))
"""
function KSPMatSolveTranspose(petsclib::PetscLibType, ksp::AbstractKSP, B::AbstractPetscMat, X::AbstractPetscMat)
    error("KSPMatSolveTranspose: no generated method for these argument types")
end

@for_petsc function KSPMatSolveTranspose(petsclib::$UnionPetscLib, ksp::AbstractKSP, B::AbstractPetscMat, X::AbstractPetscMat )

    @chk ccall(
               (:KSPMatSolveTranspose, $petsc_library),
               PetscErrorCode,
               (CKSP, CMat, CMat),
               ksp, B, X,
              )


	return nothing
end 

"""
	KSPMonitor(petsclib::PetscLibType, ksp::AbstractKSP, it::PetscInt, rnorm::PetscReal) 
runs the user provided monitor routines, if they exist

Collective

Input Parameters:
- `ksp`   - iterative solver obtained from `KSPCreate()`
- `it`    - iteration number
- `rnorm` - relative norm of the residual

Level: developer

See also: `KSPMonitorSet()`

# External Links
$(_doc_external("KSP/KSPMonitor"))
"""
function KSPMonitor(petsclib::PetscLibType, ksp::AbstractKSP, it::Integer, rnorm::Real)
    error("KSPMonitor: no generated method for these argument types")
end

@for_petsc function KSPMonitor(petsclib::$UnionPetscLib, ksp::AbstractKSP, it::$PetscInt, rnorm::$PetscReal )

    @chk ccall(
               (:KSPMonitor, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal),
               ksp, it, rnorm,
              )


	return nothing
end 

"""
	KSPMonitorCancel(petsclib::PetscLibType, ksp::AbstractKSP) 
Clears all monitors for a `KSP` object.

Logically Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Options Database Key:
- `-ksp_monitor_cancel` - Cancels all monitors that have been hardwired into a code by calls to `KSPMonitorSet()`, but does not cancel those set via the options database.

Level: intermediate

See also: `KSPMonitorResidual()`, `KSPMonitorSet()`, `KSP`

# External Links
$(_doc_external("KSP/KSPMonitorCancel"))
"""
function KSPMonitorCancel(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPMonitorCancel: no generated method for these argument types")
end

@for_petsc function KSPMonitorCancel(petsclib::$UnionPetscLib, ksp::AbstractKSP )

    @chk ccall(
               (:KSPMonitorCancel, $petsc_library),
               PetscErrorCode,
               (CKSP,),
               ksp,
              )


	return nothing
end 

"""
	KSPMonitorDynamicTolerance(petsclib::PetscLibType, ksp::AbstractKSP, its::PetscInt, fnorm::PetscReal, ctx::Ptr{Cvoid}) 
A monitor that changes the inner tolerance of nested preconditioners in every outer iteration in an adaptive way.

Collective

Input Parameters:
- `ksp`   - iterative context
- `its`   - iteration number (not used)
- `fnorm` - the current residual norm
- `ctx`   - context used by monitor

Options Database Key:
- `-sub_ksp_dynamic_tolerance coef` - coefficient of dynamic tolerance for inner solver, default is 1.0

Level: advanced

See also: `KSP`, `KSPMonitorDynamicToleranceCreate()`, `KSPMonitorDynamicToleranceDestroy()`, `KSPMonitorDynamicToleranceSetCoefficient()`

# External Links
$(_doc_external("KSP/KSPMonitorDynamicTolerance"))
"""
function KSPMonitorDynamicTolerance(petsclib::PetscLibType, ksp::AbstractKSP, its::Integer, fnorm::Real, ctx::Ptr{Cvoid})
    error("KSPMonitorDynamicTolerance: no generated method for these argument types")
end

@for_petsc function KSPMonitorDynamicTolerance(petsclib::$UnionPetscLib, ksp::AbstractKSP, its::$PetscInt, fnorm::$PetscReal, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:KSPMonitorDynamicTolerance, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{Cvoid}),
               ksp, its, fnorm, ctx,
              )


	return nothing
end 

"""
	ctx::Ptr{Cvoid} = KSPMonitorDynamicToleranceCreate(petsclib::PetscLibType) 
Creates the context used by `KSPMonitorDynamicTolerance()`

Logically Collective

Output Parameter:
- `ctx` - a void pointer

Options Database Key:
- `-sub_ksp_dynamic_tolerance coef` - coefficient of dynamic tolerance for inner solver, default is 1.0

Level: advanced

See also: `KSP`, `KSPMonitorDynamicTolerance()`, `KSPMonitorDynamicToleranceDestroy()`, `KSPMonitorDynamicToleranceSetCoefficient()`

# External Links
$(_doc_external("KSP/KSPMonitorDynamicToleranceCreate"))
"""
function KSPMonitorDynamicToleranceCreate(petsclib::PetscLibType)
    error("KSPMonitorDynamicToleranceCreate: no generated method for these argument types")
end

@for_petsc function KSPMonitorDynamicToleranceCreate(petsclib::$UnionPetscLib)
	ctx_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:KSPMonitorDynamicToleranceCreate, $petsc_library),
               PetscErrorCode,
               (Ptr{Cvoid},),
               ctx_,
              )

	ctx = ctx_[]

	return ctx
end 

"""
	KSPMonitorDynamicToleranceDestroy(petsclib::PetscLibType, ctx::Ptr{Cvoid}) 
Destroy the monitor context used in `KSPMonitorDynamicTolerance()`

Input Parameter:
- `ctx` - the monitor context

Level: advanced

See also: `KSP`, `KSPMonitorDynamicTolerance()`, `KSPMonitorSet()`, `KSPMonitorDynamicToleranceCreate()`

# External Links
$(_doc_external("KSP/KSPMonitorDynamicToleranceDestroy"))
"""
function KSPMonitorDynamicToleranceDestroy(petsclib::PetscLibType, ctx::Ptr{Cvoid})
    error("KSPMonitorDynamicToleranceDestroy: no generated method for these argument types")
end

@for_petsc function KSPMonitorDynamicToleranceDestroy(petsclib::$UnionPetscLib, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:KSPMonitorDynamicToleranceDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{Cvoid},),
               ctx,
              )


	return nothing
end 

"""
	ctx::Ptr{Cvoid} = KSPMonitorDynamicToleranceSetCoefficient(petsclib::PetscLibType, coeff::PetscReal) 
Sets the coefficient in the context used by `KSPMonitorDynamicTolerance()`

Logically Collective

Output Parameters:
- `ctx`   - the context for `KSPMonitorDynamicTolerance()`
- `coeff` - the coefficient, default is 1.0

Options Database Key:
- `-sub_ksp_dynamic_tolerance coef` - coefficient of dynamic tolerance for inner solver, default is 1.0

Level: advanced

See also: `KSP`, `KSPMonitorDynamicTolerance()`, `KSPMonitorDynamicToleranceDestroy()`, `KSPMonitorDynamicToleranceCreate()`

# External Links
$(_doc_external("KSP/KSPMonitorDynamicToleranceSetCoefficient"))
"""
function KSPMonitorDynamicToleranceSetCoefficient(petsclib::PetscLibType, coeff::Real)
    error("KSPMonitorDynamicToleranceSetCoefficient: no generated method for these argument types")
end

@for_petsc function KSPMonitorDynamicToleranceSetCoefficient(petsclib::$UnionPetscLib, coeff::$PetscReal )
	ctx_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:KSPMonitorDynamicToleranceSetCoefficient, $petsc_library),
               PetscErrorCode,
               (Ptr{Cvoid}, $PetscReal),
               ctx_, coeff,
              )

	ctx = ctx_[]

	return ctx
end 

"""
	KSPMonitorError(petsclib::PetscLibType, ksp::AbstractKSP, n::PetscInt, rnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 
Prints the error norm, as well as the (possibly preconditioned) residual norm, at each iteration of an iterative solver.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-ksp_monitor_error` - Activates `KSPMonitorError()`

Level: intermediate

See also: `KSP`, `KSPMonitorSet()`, `KSPMonitorResidual()`, `KSPMonitorTrueResidualMaxNorm()`

# External Links
$(_doc_external("KSP/KSPMonitorError"))
"""
function KSPMonitorError(petsclib::PetscLibType, ksp::AbstractKSP, n::Integer, rnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("KSPMonitorError: no generated method for these argument types")
end

@for_petsc function KSPMonitorError(petsclib::$UnionPetscLib, ksp::AbstractKSP, n::$PetscInt, rnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:KSPMonitorError, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	KSPMonitorErrorDraw(petsclib::PetscLibType, ksp::AbstractKSP, n::PetscInt, rnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 
Plots the error at each iteration of an iterative solver.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-ksp_monitor_error draw` - Activates `KSPMonitorErrorDraw()`

Level: intermediate

See also: `PETSCVIEWERDRAW`, `KSP`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `KSPMonitorErrorDrawLG()`

# External Links
$(_doc_external("KSP/KSPMonitorErrorDraw"))
"""
function KSPMonitorErrorDraw(petsclib::PetscLibType, ksp::AbstractKSP, n::Integer, rnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("KSPMonitorErrorDraw: no generated method for these argument types")
end

@for_petsc function KSPMonitorErrorDraw(petsclib::$UnionPetscLib, ksp::AbstractKSP, n::$PetscInt, rnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:KSPMonitorErrorDraw, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	KSPMonitorErrorDrawLG(petsclib::PetscLibType, ksp::AbstractKSP, n::PetscInt, rnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 
Plots the error and residual norm at each iteration of an iterative solver.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-ksp_monitor_error draw::draw_lg` - Activates `KSPMonitorTrueResidualDrawLG()`

Level: intermediate

See also: `PETSCVIEWERDRAW`, `KSP`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `KSPMonitorErrorDraw()`

# External Links
$(_doc_external("KSP/KSPMonitorErrorDrawLG"))
"""
function KSPMonitorErrorDrawLG(petsclib::PetscLibType, ksp::AbstractKSP, n::Integer, rnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("KSPMonitorErrorDrawLG: no generated method for these argument types")
end

@for_petsc function KSPMonitorErrorDrawLG(petsclib::$UnionPetscLib, ksp::AbstractKSP, n::$PetscInt, rnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:KSPMonitorErrorDrawLG, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	vf::Ptr{PetscViewerAndFormat} = KSPMonitorErrorDrawLGCreate(petsclib::PetscLibType, viewer::PetscViewer, format::PetscViewerFormat, ctx::Ptr{Cvoid}) 
Creates the context for the error and preconditioned residual plotter `KSPMonitorErrorDrawLG()`

Collective

Input Parameters:
- `viewer` - The `PetscViewer`
- `format` - The viewer format
- `ctx`    - An optional application context

Output Parameter:
- `vf` - The viewer context

Level: intermediate

See also: `PETSCVIEWERDRAW`, `KSP`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `KSPMonitorErrorDrawLG()`

# External Links
$(_doc_external("KSP/KSPMonitorErrorDrawLGCreate"))
"""
function KSPMonitorErrorDrawLGCreate(petsclib::PetscLibType, viewer::PetscViewer, format::PetscViewerFormat, ctx::Ptr{Cvoid})
    error("KSPMonitorErrorDrawLGCreate: no generated method for these argument types")
end

@for_petsc function KSPMonitorErrorDrawLGCreate(petsclib::$UnionPetscLib, viewer::PetscViewer, format::PetscViewerFormat, ctx::Ptr{Cvoid} )
	vf_ = Ref{Ptr{PetscViewerAndFormat}}()

    @chk ccall(
               (:KSPMonitorErrorDrawLGCreate, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscViewerFormat, Ptr{Cvoid}, Ptr{Ptr{PetscViewerAndFormat}}),
               viewer, format, ctx, vf_,
              )

	vf = vf_[]

	return vf
end 

"""
	KSPMonitorLGRange(petsclib::PetscLibType, ksp::AbstractKSP, n::PetscInt, rnorm::PetscReal, monctx::Ptr{Cvoid}) 
Prints line graphs summarizing the residual norm, the fraction of elements that dominate the residual, and the convergence factor at each iteration of the `KSP` solver

Collective

Input Parameters:
- `ksp`    - iterative context
- `n`      - iteration number
- `rnorm`  - the 2-norm of the residual (or an approximation)
- `monctx` - a `PetscViewer` (typically of type `PETSCVIEWERDRAW`) containing the line graphs to update

Level: intermediate

See also: `KSP`, `KSPMonitorSet()`, `KSPMonitorResidual()`, `PETSCVIEWERDRAW`, `PetscDrawLG`

# External Links
$(_doc_external("KSP/KSPMonitorLGRange"))
"""
function KSPMonitorLGRange(petsclib::PetscLibType, ksp::AbstractKSP, n::Integer, rnorm::Real, monctx::Ptr{Cvoid})
    error("KSPMonitorLGRange: no generated method for these argument types")
end

@for_petsc function KSPMonitorLGRange(petsclib::$UnionPetscLib, ksp::AbstractKSP, n::$PetscInt, rnorm::$PetscReal, monctx::Ptr{Cvoid} )

    @chk ccall(
               (:KSPMonitorLGRange, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{Cvoid}),
               ksp, n, rnorm, monctx,
              )


	return nothing
end 

"""
	KSPMonitorRegister(petsclib::PetscLibType, name::String, vtype::String, format::PetscViewerFormat, monitor::Ptr{Cvoid}, create::Ptr{Cvoid}, destroy::Ptr{Cvoid}) 
Registers a Krylov subspace solver monitor routine that may be accessed with `KSPMonitorSetFromOptions()`

Not Collective

Input Parameters:
- `name`    - name of a new monitor type
- `vtype`   - A `PetscViewerType` for the output
- `format`  - A `PetscViewerFormat` for the output
- `monitor` - Monitor routine, see `KSPMonitorRegisterFn`
- `create`  - Creation routine, or `NULL`
- `destroy` - Destruction routine, or `NULL`

Level: advanced

See also: `KSP`, `KSPMonitorSet()`, `KSPMonitorRegisterAll()`, `KSPMonitorSetFromOptions()`

# External Links
$(_doc_external("KSP/KSPMonitorRegister"))
"""
function KSPMonitorRegister(petsclib::PetscLibType, name::String, vtype::String, format::PetscViewerFormat, monitor::Ptr{Cvoid}, create::Ptr{Cvoid}, destroy::Ptr{Cvoid})
    error("KSPMonitorRegister: no generated method for these argument types")
end

@for_petsc function KSPMonitorRegister(petsclib::$UnionPetscLib, name::String, vtype::String, format::PetscViewerFormat, monitor::Ptr{Cvoid}, create::Ptr{Cvoid}, destroy::Ptr{Cvoid} )

    @chk ccall(
               (:KSPMonitorRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, PetscViewerType, PetscViewerFormat, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
               name, vtype, format, monitor, create, destroy,
              )


	return nothing
end 

"""
	KSPMonitorResidual(petsclib::PetscLibType, ksp::AbstractKSP, n::PetscInt, rnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 
Print the (possibly preconditioned, possibly approximate) residual norm at each iteration of an iterative solver.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - (preconditioned) residual norm value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-ksp_monitor` - Activates `KSPMonitorResidual()` to print the norm value at each iteration

Level: intermediate

See also: `KSP`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `KSPMonitorResidualView()`, `KSPMonitorResidualDrawLG()`,
`KSPMonitorResidualRange()`, `KSPMonitorTrueResidualDraw()`, `KSPMonitorTrueResidualDrawLG()`, `KSPMonitorTrueResidualMax()`,
`KSPMonitorSingularValue()`, `KSPMonitorSolutionDrawLG()`, `KSPMonitorSolutionDraw()`, `KSPMonitorSolution()`,
`KSPMonitorErrorDrawLG()`, `KSPMonitorErrorDraw()`, `KSPMonitorError()`

# External Links
$(_doc_external("KSP/KSPMonitorResidual"))
"""
function KSPMonitorResidual(petsclib::PetscLibType, ksp::AbstractKSP, n::Integer, rnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("KSPMonitorResidual: no generated method for these argument types")
end

@for_petsc function KSPMonitorResidual(petsclib::$UnionPetscLib, ksp::AbstractKSP, n::$PetscInt, rnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:KSPMonitorResidual, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	KSPMonitorResidualDrawLG(petsclib::PetscLibType, ksp::AbstractKSP, n::PetscInt, rnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 
Plots the (possibly preconditioned) residual norm at each iteration of an iterative solver.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-ksp_monitor draw::draw_lg` - Activates `KSPMonitorResidualDrawLG()`

Level: intermediate

See also: `KSP`, `PETSCVIEWERDRAW`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `KSPMonitorResidualView()`, `KSPMonitorResidual()`

# External Links
$(_doc_external("KSP/KSPMonitorResidualDrawLG"))
"""
function KSPMonitorResidualDrawLG(petsclib::PetscLibType, ksp::AbstractKSP, n::Integer, rnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("KSPMonitorResidualDrawLG: no generated method for these argument types")
end

@for_petsc function KSPMonitorResidualDrawLG(petsclib::$UnionPetscLib, ksp::AbstractKSP, n::$PetscInt, rnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:KSPMonitorResidualDrawLG, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	vf::Ptr{PetscViewerAndFormat} = KSPMonitorResidualDrawLGCreate(petsclib::PetscLibType, viewer::PetscViewer, format::PetscViewerFormat, ctx::Ptr{Cvoid}) 
Creates the context for the (possibly preconditioned) residual norm monitor `KSPMonitorResidualDrawLG()`

Collective

Input Parameters:
- `viewer` - The `PetscViewer` of type `PETSCVIEWERDRAW`
- `format` - The viewer format
- `ctx`    - An optional application context

Output Parameter:
- `vf` - The viewer context

Level: intermediate

See also: `KSP`, `PETSCVIEWERDRAW`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `KSPMonitorResidualDrawLG()`,
`PetscViewerFormat`, `PetscViewer`, `PetscViewerAndFormat`

# External Links
$(_doc_external("KSP/KSPMonitorResidualDrawLGCreate"))
"""
function KSPMonitorResidualDrawLGCreate(petsclib::PetscLibType, viewer::PetscViewer, format::PetscViewerFormat, ctx::Ptr{Cvoid})
    error("KSPMonitorResidualDrawLGCreate: no generated method for these argument types")
end

@for_petsc function KSPMonitorResidualDrawLGCreate(petsclib::$UnionPetscLib, viewer::PetscViewer, format::PetscViewerFormat, ctx::Ptr{Cvoid} )
	vf_ = Ref{Ptr{PetscViewerAndFormat}}()

    @chk ccall(
               (:KSPMonitorResidualDrawLGCreate, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscViewerFormat, Ptr{Cvoid}, Ptr{Ptr{PetscViewerAndFormat}}),
               viewer, format, ctx, vf_,
              )

	vf = vf_[]

	return vf
end 

"""
	KSPMonitorResidualRange(petsclib::PetscLibType, ksp::AbstractKSP, it::PetscInt, rnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 
Prints the percentage of residual elements that are more than 10 percent of the maximum value.

Collective

Input Parameters:
- `ksp`   - iterative context
- `it`    - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-ksp_monitor_range` - Activates `KSPMonitorResidualRange()`

Level: intermediate

See also: `KSP`, `KSPMonitorSet()`, `KSPMonitorResidual()`

# External Links
$(_doc_external("KSP/KSPMonitorResidualRange"))
"""
function KSPMonitorResidualRange(petsclib::PetscLibType, ksp::AbstractKSP, it::Integer, rnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("KSPMonitorResidualRange: no generated method for these argument types")
end

@for_petsc function KSPMonitorResidualRange(petsclib::$UnionPetscLib, ksp::AbstractKSP, it::$PetscInt, rnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:KSPMonitorResidualRange, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, it, rnorm, vf,
              )


	return nothing
end 

"""
	KSPMonitorResidualShort(petsclib::PetscLibType, ksp::AbstractKSP, its::PetscInt, fnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 

# External Links
$(_doc_external("KSP/KSPMonitorResidualShort"))
"""
function KSPMonitorResidualShort(petsclib::PetscLibType, ksp::AbstractKSP, its::Integer, fnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("KSPMonitorResidualShort: no generated method for these argument types")
end

@for_petsc function KSPMonitorResidualShort(petsclib::$UnionPetscLib, ksp::AbstractKSP, its::$PetscInt, fnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:KSPMonitorResidualShort, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, its, fnorm, vf,
              )


	return nothing
end 

"""
	KSPMonitorResidualView(petsclib::PetscLibType, ksp::AbstractKSP, n::PetscInt, rnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 
Plots the (possibly preconditioned) residual at each iteration of an iterative solver.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-ksp_monitor viewertype` - Activates `KSPMonitorResidualView()`

Level: intermediate

See also: `KSP`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `KSPMonitorResidual()`, `KSPMonitorResidualDrawLG()`

# External Links
$(_doc_external("KSP/KSPMonitorResidualView"))
"""
function KSPMonitorResidualView(petsclib::PetscLibType, ksp::AbstractKSP, n::Integer, rnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("KSPMonitorResidualView: no generated method for these argument types")
end

@for_petsc function KSPMonitorResidualView(petsclib::$UnionPetscLib, ksp::AbstractKSP, n::$PetscInt, rnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:KSPMonitorResidualView, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	KSPMonitorSAWs(petsclib::PetscLibType, ksp::AbstractKSP, n::PetscInt, rnorm::PetscReal, ctx::Ptr{Cvoid}) 
monitor `KSP` solution using SAWs

Logically Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `ctx`   - created with `KSPMonitorSAWsCreate()`

Level: advanced

See also: `KSP`, `KSPMonitorSet()`, `KSPMonitorSAWsCreate()`, `KSPMonitorSAWsDestroy()`, `KSPMonitorSingularValue()`, `KSPComputeExtremeSingularValues()`, `PetscViewerSAWsOpen()`

# External Links
$(_doc_external("KSP/KSPMonitorSAWs"))
"""
function KSPMonitorSAWs(petsclib::PetscLibType, ksp::AbstractKSP, n::Integer, rnorm::Real, ctx::Ptr{Cvoid})
    error("KSPMonitorSAWs: no generated method for these argument types")
end

@for_petsc function KSPMonitorSAWs(petsclib::$UnionPetscLib, ksp::AbstractKSP, n::$PetscInt, rnorm::$PetscReal, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:KSPMonitorSAWs, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{Cvoid}),
               ksp, n, rnorm, ctx,
              )


	return nothing
end 

"""
	ctx::Ptr{Cvoid} = KSPMonitorSAWsCreate(petsclib::PetscLibType, ksp::AbstractKSP) 
create an SAWs monitor context for `KSP`

Collective

Input Parameter:
- `ksp` - `KSP` to monitor

Output Parameter:
- `ctx` - context for monitor

Level: developer

See also: `KSP`, `KSPMonitorSet()`, `KSPMonitorSAWs()`, `KSPMonitorSAWsDestroy()`

# External Links
$(_doc_external("KSP/KSPMonitorSAWsCreate"))
"""
function KSPMonitorSAWsCreate(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPMonitorSAWsCreate: no generated method for these argument types")
end

@for_petsc function KSPMonitorSAWsCreate(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	ctx_ = Ref{Ptr{Cvoid}}()

    @chk ccall(
               (:KSPMonitorSAWsCreate, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Ptr{Cvoid}}),
               ksp, ctx_,
              )

	ctx = ctx_[]

	return ctx
end 

"""
	KSPMonitorSAWsDestroy(petsclib::PetscLibType, ctx::Ptr{Cvoid}) 
destroy a monitor context created with `KSPMonitorSAWsCreate()`

Collective

Input Parameter:
- `ctx` - monitor context

Level: developer

See also: `KSP`, `KSPMonitorSet()`, `KSPMonitorSAWsCreate()`

# External Links
$(_doc_external("KSP/KSPMonitorSAWsDestroy"))
"""
function KSPMonitorSAWsDestroy(petsclib::PetscLibType, ctx::Ptr{Cvoid})
    error("KSPMonitorSAWsDestroy: no generated method for these argument types")
end

@for_petsc function KSPMonitorSAWsDestroy(petsclib::$UnionPetscLib, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:KSPMonitorSAWsDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{Cvoid},),
               ctx,
              )


	return nothing
end 

"""
	KSPMonitorSNESResidual(petsclib::PetscLibType, ksp::AbstractKSP, n::PetscInt, rnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 
Prints the `SNES` residual norm, as well as the `KSP` residual norm, at each iteration of a `KSPSolve()` called within a `SNESSolve()`.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-snes_monitor_ksp` - Activates `KSPMonitorSNESResidual()`

Level: intermediate

See also: `SNES`, `KSPMonitorSet()`, `KSPMonitorResidual()`, `KSPMonitorTrueResidualMaxNorm()`, `KSPMonitor()`, `SNESMonitor()`, `PetscViewerAndFormat()`

# External Links
$(_doc_external("SNES/KSPMonitorSNESResidual"))
"""
function KSPMonitorSNESResidual(petsclib::PetscLibType, ksp::AbstractKSP, n::Integer, rnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("KSPMonitorSNESResidual: no generated method for these argument types")
end

@for_petsc function KSPMonitorSNESResidual(petsclib::$UnionPetscLib, ksp::AbstractKSP, n::$PetscInt, rnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:KSPMonitorSNESResidual, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	KSPMonitorSNESResidualDrawLG(petsclib::PetscLibType, ksp::AbstractKSP, n::PetscInt, rnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 
Plots the linear `KSP` residual norm and the `SNES` residual norm of a `KSPSolve()` called within a `SNESSolve()`.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context, created with `KSPMonitorSNESResidualDrawLGCreate()`

Options Database Key:
- `-snes_monitor_ksp draw::draw_lg` - Activates `KSPMonitorSNESResidualDrawLG()`

Level: intermediate

See also: `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `SNESMonitor()`, `KSPMonitor()`, `KSPMonitorSNESResidualDrawLGCreate()`

# External Links
$(_doc_external("SNES/KSPMonitorSNESResidualDrawLG"))
"""
function KSPMonitorSNESResidualDrawLG(petsclib::PetscLibType, ksp::AbstractKSP, n::Integer, rnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("KSPMonitorSNESResidualDrawLG: no generated method for these argument types")
end

@for_petsc function KSPMonitorSNESResidualDrawLG(petsclib::$UnionPetscLib, ksp::AbstractKSP, n::$PetscInt, rnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:KSPMonitorSNESResidualDrawLG, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	vf::Ptr{PetscViewerAndFormat} = KSPMonitorSNESResidualDrawLGCreate(petsclib::PetscLibType, viewer::PetscViewer, format::PetscViewerFormat, ctx::Ptr{Cvoid}) 
Creates the `PetscViewer` used by `KSPMonitorSNESResidualDrawLG()`

Collective

Input Parameters:
- `viewer` - The `PetscViewer`
- `format` - The viewer format
- `ctx`    - An optional application context

Output Parameter:
- `vf` - The viewer context

Level: intermediate

See also: `KSP`, `SNES`, `PetscViewerFormat`, `PetscViewerAndFormat`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`

# External Links
$(_doc_external("SNES/KSPMonitorSNESResidualDrawLGCreate"))
"""
function KSPMonitorSNESResidualDrawLGCreate(petsclib::PetscLibType, viewer::PetscViewer, format::PetscViewerFormat, ctx::Ptr{Cvoid})
    error("KSPMonitorSNESResidualDrawLGCreate: no generated method for these argument types")
end

@for_petsc function KSPMonitorSNESResidualDrawLGCreate(petsclib::$UnionPetscLib, viewer::PetscViewer, format::PetscViewerFormat, ctx::Ptr{Cvoid} )
	vf_ = Ref{Ptr{PetscViewerAndFormat}}()

    @chk ccall(
               (:KSPMonitorSNESResidualDrawLGCreate, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscViewerFormat, Ptr{Cvoid}, Ptr{Ptr{PetscViewerAndFormat}}),
               viewer, format, ctx, vf_,
              )

	vf = vf_[]

	return vf
end 

"""
	KSPMonitorSet(petsclib::PetscLibType, ksp::AbstractKSP, monitor::Ptr{Cvoid}, ctx::Ptr{Cvoid}, monitordestroy::Ptr{Cvoid}) 
Sets an ADDITIONAL function to be called at every iteration to monitor, i.e. display in some way, perhaps by printing in the terminal,
the residual norm computed in a `KSPSolve()`

Logically Collective

Input Parameters:
- `ksp`            - iterative solver obtained from `KSPCreate()`
- `monitor`        - pointer to function (if this is `NULL`, it turns off monitoring, see `KSPMonitorFn`
- `ctx`            - [optional] context for private data for the monitor routine (use `NULL` if no context is needed)
- `monitordestroy` - [optional] routine that frees monitor context (may be `NULL`), see `PetscCtxDestroyFn` for the calling sequence

Options Database Keys:
- `-ksp_monitor`                             - sets `KSPMonitorResidual()`
- `-ksp_monitor hdf5:filename`               - sets `KSPMonitorResidualView()` and saves residual
- `-ksp_monitor draw`                        - sets `KSPMonitorResidualView()` and plots residual
- `-ksp_monitor draw::draw_lg`               - sets `KSPMonitorResidualDrawLG()` and plots residual
- `-ksp_monitor_pause_final`                 - Pauses any graphics when the solve finishes (only works for internal monitors)
- `-ksp_monitor_true_residual`               - sets `KSPMonitorTrueResidual()`
- `-ksp_monitor_true_residual draw::draw_lg` - sets `KSPMonitorTrueResidualDrawLG()` and plots residual
- `-ksp_monitor_max`                         - sets `KSPMonitorTrueResidualMax()`
- `-ksp_monitor_singular_value`              - sets `KSPMonitorSingularValue()`
- `-ksp_monitor_cancel`                      - cancels all monitors that have been hardwired into a code by calls to `KSPMonitorSet()`, but
does not cancel those set via the options database.

Level: beginner

See also: `KSPMonitorResidual()`, `KSPMonitorRegister()`, `KSPMonitorCancel()`, `KSP`, `PetscCtxDestroyFn`

# External Links
$(_doc_external("KSP/KSPMonitorSet"))
"""
function KSPMonitorSet(petsclib::PetscLibType, ksp::AbstractKSP, monitor::Ptr{Cvoid}, ctx::Ptr{Cvoid}, monitordestroy::Ptr{Cvoid})
    error("KSPMonitorSet: no generated method for these argument types")
end

@for_petsc function KSPMonitorSet(petsclib::$UnionPetscLib, ksp::AbstractKSP, monitor::Ptr{Cvoid}, ctx::Ptr{Cvoid}, monitordestroy::Ptr{Cvoid} )

    @chk ccall(
               (:KSPMonitorSet, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
               ksp, monitor, ctx, monitordestroy,
              )


	return nothing
end 

"""
	KSPMonitorSetFromOptions(petsclib::PetscLibType, ksp::AbstractKSP, opt::String, name::String, ctx::Ptr{Cvoid}) 
Sets a monitor function and viewer appropriate for the type indicated by the user in the options database

Collective

Input Parameters:
- `ksp`  - `KSP` object you wish to monitor
- `opt`  - the command line option for this monitor
- `name` - the monitor type one is seeking
- `ctx`  - An optional application context for the monitor, or `NULL`

Level: developer

See also: `KSPMonitorRegister()`, `KSPMonitorSet()`, `PetscOptionsCreateViewer()`, `PetscOptionsGetReal()`, `PetscOptionsHasName()`, `PetscOptionsGetString()`,
`PetscOptionsGetIntArray()`, `PetscOptionsGetRealArray()`, `PetscOptionsBool()`,
`PetscOptionsInt()`, `PetscOptionsString()`, `PetscOptionsReal()`,
`PetscOptionsName()`, `PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`,
`PetscOptionsStringArray()`, `PetscOptionsRealArray()`, `PetscOptionsScalar()`,
`PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`

# External Links
$(_doc_external("KSP/KSPMonitorSetFromOptions"))
"""
function KSPMonitorSetFromOptions(petsclib::PetscLibType, ksp::AbstractKSP, opt::String, name::String, ctx::Ptr{Cvoid})
    error("KSPMonitorSetFromOptions: no generated method for these argument types")
end

@for_petsc function KSPMonitorSetFromOptions(petsclib::$UnionPetscLib, ksp::AbstractKSP, opt::String, name::String, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:KSPMonitorSetFromOptions, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cvoid}),
               ksp, opt, name, ctx,
              )


	return nothing
end 

"""
	KSPMonitorSingularValue(petsclib::PetscLibType, ksp::AbstractKSP, n::PetscInt, rnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 
Prints the two norm of the true residual and estimation of the extreme singular values of the preconditioned problem at each iteration.

Logically Collective

Input Parameters:
- `ksp`   - the iterative context
- `n`     - the iteration
- `rnorm` - the two norm of the residual
- `vf`    - The viewer context

Options Database Key:
- `-ksp_monitor_singular_value` - Activates `KSPMonitorSingularValue()`

Level: intermediate

See also: `KSP`, `KSPMonitorSet()`, `KSPComputeExtremeSingularValues()`, `KSPMonitorSingularValueCreate()`

# External Links
$(_doc_external("KSP/KSPMonitorSingularValue"))
"""
function KSPMonitorSingularValue(petsclib::PetscLibType, ksp::AbstractKSP, n::Integer, rnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("KSPMonitorSingularValue: no generated method for these argument types")
end

@for_petsc function KSPMonitorSingularValue(petsclib::$UnionPetscLib, ksp::AbstractKSP, n::$PetscInt, rnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:KSPMonitorSingularValue, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	vf::Ptr{PetscViewerAndFormat} = KSPMonitorSingularValueCreate(petsclib::PetscLibType, viewer::PetscViewer, format::PetscViewerFormat, ctx::Ptr{Cvoid}) 
Creates the singular value monitor context needed by `KSPMonitorSingularValue()`

Collective

Input Parameters:
- `viewer` - The PetscViewer
- `format` - The viewer format
- `ctx`    - An optional application context

Output Parameter:
- `vf` - The viewer context

Level: intermediate

See also: `KSP`, `KSPMonitorSet()`, `KSPMonitorSingularValue()`, `PetscViewer`

# External Links
$(_doc_external("KSP/KSPMonitorSingularValueCreate"))
"""
function KSPMonitorSingularValueCreate(petsclib::PetscLibType, viewer::PetscViewer, format::PetscViewerFormat, ctx::Ptr{Cvoid})
    error("KSPMonitorSingularValueCreate: no generated method for these argument types")
end

@for_petsc function KSPMonitorSingularValueCreate(petsclib::$UnionPetscLib, viewer::PetscViewer, format::PetscViewerFormat, ctx::Ptr{Cvoid} )
	vf_ = Ref{Ptr{PetscViewerAndFormat}}()

    @chk ccall(
               (:KSPMonitorSingularValueCreate, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscViewerFormat, Ptr{Cvoid}, Ptr{Ptr{PetscViewerAndFormat}}),
               viewer, format, ctx, vf_,
              )

	vf = vf_[]

	return vf
end 

"""
	KSPMonitorSolution(petsclib::PetscLibType, ksp::AbstractKSP, n::PetscInt, rnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 
Print the solution norm at each iteration of an iterative solver.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-ksp_monitor_solution` - Activates `KSPMonitorSolution()`

Level: intermediate

See also: `KSPMonitorSet()`, `KSPMonitorTrueResidual()`

# External Links
$(_doc_external("KSP/KSPMonitorSolution"))
"""
function KSPMonitorSolution(petsclib::PetscLibType, ksp::AbstractKSP, n::Integer, rnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("KSPMonitorSolution: no generated method for these argument types")
end

@for_petsc function KSPMonitorSolution(petsclib::$UnionPetscLib, ksp::AbstractKSP, n::$PetscInt, rnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:KSPMonitorSolution, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	KSPMonitorSolutionDraw(petsclib::PetscLibType, ksp::AbstractKSP, n::PetscInt, rnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 
Plots the solution at each iteration of an iterative solver.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-ksp_monitor_solution draw` - Activates `KSPMonitorSolutionDraw()`

Level: intermediate

See also: `KSPMonitorSet()`, `KSPMonitorTrueResidual()`

# External Links
$(_doc_external("KSP/KSPMonitorSolutionDraw"))
"""
function KSPMonitorSolutionDraw(petsclib::PetscLibType, ksp::AbstractKSP, n::Integer, rnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("KSPMonitorSolutionDraw: no generated method for these argument types")
end

@for_petsc function KSPMonitorSolutionDraw(petsclib::$UnionPetscLib, ksp::AbstractKSP, n::$PetscInt, rnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:KSPMonitorSolutionDraw, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	KSPMonitorSolutionDrawLG(petsclib::PetscLibType, ksp::AbstractKSP, n::PetscInt, rnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 
Plots the solution norm at each iteration of an iterative solver.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-ksp_monitor_solution draw::draw_lg` - Activates `KSPMonitorSolutionDrawLG()`

Level: intermediate

See also: `KSP`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `KSPMonitorSolutionDrawLGCreate()`

# External Links
$(_doc_external("KSP/KSPMonitorSolutionDrawLG"))
"""
function KSPMonitorSolutionDrawLG(petsclib::PetscLibType, ksp::AbstractKSP, n::Integer, rnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("KSPMonitorSolutionDrawLG: no generated method for these argument types")
end

@for_petsc function KSPMonitorSolutionDrawLG(petsclib::$UnionPetscLib, ksp::AbstractKSP, n::$PetscInt, rnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:KSPMonitorSolutionDrawLG, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	vf::Ptr{PetscViewerAndFormat} = KSPMonitorSolutionDrawLGCreate(petsclib::PetscLibType, viewer::PetscViewer, format::PetscViewerFormat, ctx::Ptr{Cvoid}) 
Creates the context for the `KSP` monitor `KSPMonitorSolutionDrawLG()`

Collective

Input Parameters:
- `viewer` - The `PetscViewer`
- `format` - The viewer format
- `ctx`    - An optional application context

Output Parameter:
- `vf` - The viewer context

Level: intermediate

See also: `KSPMonitorSet()`, `KSPMonitorTrueResidual()`

# External Links
$(_doc_external("KSP/KSPMonitorSolutionDrawLGCreate"))
"""
function KSPMonitorSolutionDrawLGCreate(petsclib::PetscLibType, viewer::PetscViewer, format::PetscViewerFormat, ctx::Ptr{Cvoid})
    error("KSPMonitorSolutionDrawLGCreate: no generated method for these argument types")
end

@for_petsc function KSPMonitorSolutionDrawLGCreate(petsclib::$UnionPetscLib, viewer::PetscViewer, format::PetscViewerFormat, ctx::Ptr{Cvoid} )
	vf_ = Ref{Ptr{PetscViewerAndFormat}}()

    @chk ccall(
               (:KSPMonitorSolutionDrawLGCreate, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscViewerFormat, Ptr{Cvoid}, Ptr{Ptr{PetscViewerAndFormat}}),
               viewer, format, ctx, vf_,
              )

	vf = vf_[]

	return vf
end 

"""
	KSPMonitorTrueResidual(petsclib::PetscLibType, ksp::AbstractKSP, n::PetscInt, rnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 
Prints the true residual norm, as well as the (possibly preconditioned, possibly approximate) residual norm,
at each iteration of a `KSPSolve()` iterative solver.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-ksp_monitor_true_residual` - Activates `KSPMonitorTrueResidual()` to print both norm values at each iteration

Level: intermediate

See also: `KSP`, `KSPMonitorSet()`, `KSPMonitorResidual()`, `KSPMonitorTrueResidualMaxNorm()`, `PetscViewerAndFormat`

# External Links
$(_doc_external("KSP/KSPMonitorTrueResidual"))
"""
function KSPMonitorTrueResidual(petsclib::PetscLibType, ksp::AbstractKSP, n::Integer, rnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("KSPMonitorTrueResidual: no generated method for these argument types")
end

@for_petsc function KSPMonitorTrueResidual(petsclib::$UnionPetscLib, ksp::AbstractKSP, n::$PetscInt, rnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:KSPMonitorTrueResidual, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	KSPMonitorTrueResidualDrawLG(petsclib::PetscLibType, ksp::AbstractKSP, n::PetscInt, rnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 
Plots the true residual norm at each iteration of an iterative solver.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-ksp_monitor_true_residual draw::draw_lg` - Activates `KSPMonitorTrueResidualDrawLG()`

Level: intermediate

See also: `PETSCVIEWERDRAW`, `KSP`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `KSPMonitorTrueResidualDraw()`, `KSPMonitorResidual`,
`KSPMonitorTrueResidualDrawLGCreate()`

# External Links
$(_doc_external("KSP/KSPMonitorTrueResidualDrawLG"))
"""
function KSPMonitorTrueResidualDrawLG(petsclib::PetscLibType, ksp::AbstractKSP, n::Integer, rnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("KSPMonitorTrueResidualDrawLG: no generated method for these argument types")
end

@for_petsc function KSPMonitorTrueResidualDrawLG(petsclib::$UnionPetscLib, ksp::AbstractKSP, n::$PetscInt, rnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:KSPMonitorTrueResidualDrawLG, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	vf::Ptr{PetscViewerAndFormat} = KSPMonitorTrueResidualDrawLGCreate(petsclib::PetscLibType, viewer::PetscViewer, format::PetscViewerFormat, ctx::Ptr{Cvoid}) 
Creates the context for the true residual monitor `KSPMonitorTrueResidualDrawLG()`

Collective

Input Parameters:
- `viewer` - The `PetscViewer` of type `PETSCVIEWERDRAW`
- `format` - The viewer format
- `ctx`    - An optional application context

Output Parameter:
- `vf` - The viewer context

Level: intermediate

See also: `PETSCVIEWERDRAW`, `KSP`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `PetscViewerAndFormat`

# External Links
$(_doc_external("KSP/KSPMonitorTrueResidualDrawLGCreate"))
"""
function KSPMonitorTrueResidualDrawLGCreate(petsclib::PetscLibType, viewer::PetscViewer, format::PetscViewerFormat, ctx::Ptr{Cvoid})
    error("KSPMonitorTrueResidualDrawLGCreate: no generated method for these argument types")
end

@for_petsc function KSPMonitorTrueResidualDrawLGCreate(petsclib::$UnionPetscLib, viewer::PetscViewer, format::PetscViewerFormat, ctx::Ptr{Cvoid} )
	vf_ = Ref{Ptr{PetscViewerAndFormat}}()

    @chk ccall(
               (:KSPMonitorTrueResidualDrawLGCreate, $petsc_library),
               PetscErrorCode,
               (PetscViewer, PetscViewerFormat, Ptr{Cvoid}, Ptr{Ptr{PetscViewerAndFormat}}),
               viewer, format, ctx, vf_,
              )

	vf = vf_[]

	return vf
end 

"""
	KSPMonitorTrueResidualMax(petsclib::PetscLibType, ksp::AbstractKSP, n::PetscInt, rnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 
Prints the true residual max norm at each iteration of an iterative solver.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context

Options Database Key:
- `-ksp_monitor_true_residual_max` - Activates `KSPMonitorTrueResidualMax()`

Level: intermediate

See also: `KSP`, `KSPMonitorSet()`, `KSPMonitorResidual()`, `KSPMonitorTrueResidualMaxNorm()`

# External Links
$(_doc_external("KSP/KSPMonitorTrueResidualMax"))
"""
function KSPMonitorTrueResidualMax(petsclib::PetscLibType, ksp::AbstractKSP, n::Integer, rnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("KSPMonitorTrueResidualMax: no generated method for these argument types")
end

@for_petsc function KSPMonitorTrueResidualMax(petsclib::$UnionPetscLib, ksp::AbstractKSP, n::$PetscInt, rnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:KSPMonitorTrueResidualMax, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	KSPMonitorTrueResidualView(petsclib::PetscLibType, ksp::AbstractKSP, n::PetscInt, rnorm::PetscReal, vf::Vector{PetscViewerAndFormat}) 
Plots the true residual at each iteration of an iterative solver.

Collective

Input Parameters:
- `ksp`   - iterative context
- `n`     - iteration number
- `rnorm` - 2-norm (preconditioned) residual value (may be estimated).
- `vf`    - The viewer context of type `PETSCVIEWERDRAW`

Options Database Key:
- `-ksp_monitor_true_residual viewertype` - Activates `KSPMonitorTrueResidualView()`

Level: intermediate

See also: `PETSCVIEWERDRAW`, `KSP`, `KSPMonitorSet()`, `KSPMonitorTrueResidual()`, `KSPMonitorResidual()`,
`KSPMonitorTrueResidualDrawLG()`, `PetscViewerAndFormat`

# External Links
$(_doc_external("KSP/KSPMonitorTrueResidualView"))
"""
function KSPMonitorTrueResidualView(petsclib::PetscLibType, ksp::AbstractKSP, n::Integer, rnorm::Real, vf::Vector{PetscViewerAndFormat})
    error("KSPMonitorTrueResidualView: no generated method for these argument types")
end

@for_petsc function KSPMonitorTrueResidualView(petsclib::$UnionPetscLib, ksp::AbstractKSP, n::$PetscInt, rnorm::$PetscReal, vf::Vector{PetscViewerAndFormat} )

    @chk ccall(
               (:KSPMonitorTrueResidualView, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscReal, Ptr{PetscViewerAndFormat}),
               ksp, n, rnorm, vf,
              )


	return nothing
end 

"""
	mmax::PetscInt = KSPPIPEFCGGetMmax(petsclib::PetscLibType, ksp::AbstractKSP) 
get the maximum number of previous directions `KSPPIPEFCG` will store

Not Collective

Input Parameter:
- `ksp` - the Krylov space context

Output Parameter:
- `mmax` - the maximum number of previous directions allowed for orthogonalization

Level: intermediate

See also: `KSPPIPEFCG`, `KSPPIPEFCGGetTruncationType()`, `KSPPIPEFCGGetNprealloc()`, `KSPPIPEFCGSetMmax()`, `KSPFCGGetMmax()`, `KSPFCGSetMmax()`

# External Links
$(_doc_external("KSP/KSPPIPEFCGGetMmax"))
"""
function KSPPIPEFCGGetMmax(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPPIPEFCGGetMmax: no generated method for these argument types")
end

@for_petsc function KSPPIPEFCGGetMmax(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	mmax_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPPIPEFCGGetMmax, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscInt}),
               ksp, mmax_,
              )

	mmax = mmax_[]

	return mmax
end 

"""
	nprealloc::PetscInt = KSPPIPEFCGGetNprealloc(petsclib::PetscLibType, ksp::AbstractKSP) 
get the number of directions to preallocate by `KSPPIPEFCG`

Not Collective

Input Parameter:
- `ksp` - the Krylov space context

Output Parameter:
- `nprealloc` - the number of directions preallocated

Level: advanced

See also: `KSPPIPEFCG`, `KSPPIPEFCGGetTruncationType()`, `KSPPIPEFCGSetNprealloc()`, `KSPPIPEFCGSetMmax()`, `KSPPIPEFCGGetMmax()`

# External Links
$(_doc_external("KSP/KSPPIPEFCGGetNprealloc"))
"""
function KSPPIPEFCGGetNprealloc(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPPIPEFCGGetNprealloc: no generated method for these argument types")
end

@for_petsc function KSPPIPEFCGGetNprealloc(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	nprealloc_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPPIPEFCGGetNprealloc, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscInt}),
               ksp, nprealloc_,
              )

	nprealloc = nprealloc_[]

	return nprealloc
end 

"""
	truncstrat::KSPFCDTruncationType = KSPPIPEFCGGetTruncationType(petsclib::PetscLibType, ksp::AbstractKSP) 
get the truncation strategy employed by `KSPPIPEFCG`

Not Collective

Input Parameter:
- `ksp` - the Krylov space context

Output Parameter:
- `truncstrat` - the strategy type

Level: intermediate

See also: `KSPPIPEFCG`, `KSPPIPEFCGSetTruncationType()`, `KSPFCDTruncationType`, `KSP_FCD_TRUNC_TYPE_STANDARD`, `KSP_FCD_TRUNC_TYPE_NOTAY`

# External Links
$(_doc_external("KSP/KSPPIPEFCGGetTruncationType"))
"""
function KSPPIPEFCGGetTruncationType(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPPIPEFCGGetTruncationType: no generated method for these argument types")
end

@for_petsc function KSPPIPEFCGGetTruncationType(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	truncstrat_ = Ref{KSPFCDTruncationType}()

    @chk ccall(
               (:KSPPIPEFCGGetTruncationType, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPFCDTruncationType}),
               ksp, truncstrat_,
              )

	truncstrat = truncstrat_[]

	return truncstrat
end 

"""
	KSPPIPEFCGSetMmax(petsclib::PetscLibType, ksp::AbstractKSP, mmax::PetscInt) 
set the maximum number of previous directions `KSPPIPEFCG` will store for orthogonalization

Logically Collective

Input Parameters:
- `ksp`  - the Krylov space context
- `mmax` - the maximum number of previous directions to orthogonalize against

Options Database Key:
- `-ksp_pipefcg_mmax N` - maximum number of previous directions

Level: intermediate

See also: `KSPPIPEFCG`, `KSPPIPEFCGSetTruncationType()`, `KSPPIPEFCGSetNprealloc()`, `KSPFCGSetMmax()`, `KSPFCGGetMmax()`

# External Links
$(_doc_external("KSP/KSPPIPEFCGSetMmax"))
"""
function KSPPIPEFCGSetMmax(petsclib::PetscLibType, ksp::AbstractKSP, mmax::Integer)
    error("KSPPIPEFCGSetMmax: no generated method for these argument types")
end

@for_petsc function KSPPIPEFCGSetMmax(petsclib::$UnionPetscLib, ksp::AbstractKSP, mmax::$PetscInt )

    @chk ccall(
               (:KSPPIPEFCGSetMmax, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, mmax,
              )


	return nothing
end 

"""
	KSPPIPEFCGSetNprealloc(petsclib::PetscLibType, ksp::AbstractKSP, nprealloc::PetscInt) 
set the number of directions to preallocate with `KSPPIPEFCG`

Logically Collective

Input Parameters:
- `ksp`       - the Krylov space context
- `nprealloc` - the number of vectors to preallocate

Options Database Key:
- `-ksp_pipefcg_nprealloc N` - the number of vectors to preallocate

Level: advanced

See also: `KSPPIPEFCG`, `KSPPIPEFCGSetTruncationType()`, `KSPPIPEFCGGetNprealloc()`, `KSPPIPEFCGSetMmax()`, `KSPPIPEFCGGetMmax()`

# External Links
$(_doc_external("KSP/KSPPIPEFCGSetNprealloc"))
"""
function KSPPIPEFCGSetNprealloc(petsclib::PetscLibType, ksp::AbstractKSP, nprealloc::Integer)
    error("KSPPIPEFCGSetNprealloc: no generated method for these argument types")
end

@for_petsc function KSPPIPEFCGSetNprealloc(petsclib::$UnionPetscLib, ksp::AbstractKSP, nprealloc::$PetscInt )

    @chk ccall(
               (:KSPPIPEFCGSetNprealloc, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, nprealloc,
              )


	return nothing
end 

"""
	KSPPIPEFCGSetTruncationType(petsclib::PetscLibType, ksp::AbstractKSP, truncstrat::KSPFCDTruncationType) 
specify how many of its stored previous directions `KSPPIPEFCG` uses during orthogonalization

Logically Collective

Input Parameters:
- `ksp`        - the Krylov space context
- `truncstrat` - the choice of strategy
``
KSP_FCD_TRUNC_TYPE_STANDARD uses all (up to `mmax`) stored directions
KSP_FCD_TRUNC_TYPE_NOTAY uses `max(1,mod(i,mmax))` stored directions at iteration i = 0, 1, ...
``

Options Database Key:
- `-ksp_pipefcg_truncation_type (standard|notay)` - which stored search directions to orthogonalize against

Level: intermediate

See also: `KSPPIPEFCG`, `KSPPIPEFCGGetTruncationType`, `KSPFCDTruncationType`, `KSP_FCD_TRUNC_TYPE_STANDARD`, `KSP_FCD_TRUNC_TYPE_NOTAY`

# External Links
$(_doc_external("KSP/KSPPIPEFCGSetTruncationType"))
"""
function KSPPIPEFCGSetTruncationType(petsclib::PetscLibType, ksp::AbstractKSP, truncstrat::KSPFCDTruncationType)
    error("KSPPIPEFCGSetTruncationType: no generated method for these argument types")
end

@for_petsc function KSPPIPEFCGSetTruncationType(petsclib::$UnionPetscLib, ksp::AbstractKSP, truncstrat::KSPFCDTruncationType )

    @chk ccall(
               (:KSPPIPEFCGSetTruncationType, $petsc_library),
               PetscErrorCode,
               (CKSP, KSPFCDTruncationType),
               ksp, truncstrat,
              )


	return nothing
end 

"""
	KSPPIPEFGMRESSetShift(petsclib::PetscLibType, ksp::AbstractKSP, shift::PetscScalar) 
Set the shift parameter for the flexible, pipelined `KSPPIPEFGMRES` solver.

Logically Collective

Input Parameters:
- `ksp`   - the Krylov space context
- `shift` - the shift

Options Database Key:
- `-ksp_pipefgmres_shift shift` - set the shift parameter

Level: intermediate

See also: `KSPPIPEFGMRES`, `KSPComputeEigenvalues()`

# External Links
$(_doc_external("KSP/KSPPIPEFGMRESSetShift"))
"""
function KSPPIPEFGMRESSetShift(petsclib::PetscLibType, ksp::AbstractKSP, shift::Number)
    error("KSPPIPEFGMRESSetShift: no generated method for these argument types")
end

@for_petsc function KSPPIPEFGMRESSetShift(petsclib::$UnionPetscLib, ksp::AbstractKSP, shift::$PetscScalar )

    @chk ccall(
               (:KSPPIPEFGMRESSetShift, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscScalar),
               ksp, shift,
              )


	return nothing
end 

"""
	mmax::PetscInt = KSPPIPEGCRGetMmax(petsclib::PetscLibType, ksp::AbstractKSP) 
get the maximum number of previous directions `KSPPIPEGCR` will store

Not Collective

Input Parameter:
- `ksp` - the Krylov space context

Output Parameter:
- `mmax` - the maximum number of previous directions allowed for orthogonalization

Level: intermediate

See also: `KSPPIPEGCR`, `KSPPIPEGCRGetTruncationType()`, `KSPPIPEGCRGetNprealloc()`, `KSPPIPEGCRSetMmax()`

# External Links
$(_doc_external("KSP/KSPPIPEGCRGetMmax"))
"""
function KSPPIPEGCRGetMmax(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPPIPEGCRGetMmax: no generated method for these argument types")
end

@for_petsc function KSPPIPEGCRGetMmax(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	mmax_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPPIPEGCRGetMmax, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscInt}),
               ksp, mmax_,
              )

	mmax = mmax_[]

	return mmax
end 

"""
	nprealloc::PetscInt = KSPPIPEGCRGetNprealloc(petsclib::PetscLibType, ksp::AbstractKSP) 
get the number of directions preallocate by `KSPPIPEGCR`

Not Collective

Input Parameter:
- `ksp` - the Krylov space context

Output Parameter:
- `nprealloc` - the number of directions preallocated

Level: advanced

See also: `KSPPIPEGCR`, `KSPPIPEGCRGetTruncationType()`, `KSPPIPEGCRSetNprealloc()`

# External Links
$(_doc_external("KSP/KSPPIPEGCRGetNprealloc"))
"""
function KSPPIPEGCRGetNprealloc(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPPIPEGCRGetNprealloc: no generated method for these argument types")
end

@for_petsc function KSPPIPEGCRGetNprealloc(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	nprealloc_ = Ref{$PetscInt}()

    @chk ccall(
               (:KSPPIPEGCRGetNprealloc, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscInt}),
               ksp, nprealloc_,
              )

	nprealloc = nprealloc_[]

	return nprealloc
end 

"""
	truncstrat::KSPFCDTruncationType = KSPPIPEGCRGetTruncationType(petsclib::PetscLibType, ksp::AbstractKSP) 
get the truncation strategy employed by `KSPPIPEGCR`

Not Collective

Input Parameter:
- `ksp` - the Krylov space context

Output Parameter:
- `truncstrat` - the strategy type
``
KSP_FCD_TRUNC_TYPE_STANDARD uses all (up to `mmax`) stored directions
KSP_FCD_TRUNC_TYPE_NOTAY uses the last `max(1,mod(i,mmax))` directions at iteration i =0, 1, ..
``

Level: intermediate

See also: `KSPPIPEGCR`, `KSPPIPEGCRSetTruncationType()`, `KSPFCDTruncationType`, `KSP_FCD_TRUNC_TYPE_STANDARD`, `KSP_FCD_TRUNC_TYPE_NOTAY`

# External Links
$(_doc_external("KSP/KSPPIPEGCRGetTruncationType"))
"""
function KSPPIPEGCRGetTruncationType(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPPIPEGCRGetTruncationType: no generated method for these argument types")
end

@for_petsc function KSPPIPEGCRGetTruncationType(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	truncstrat_ = Ref{KSPFCDTruncationType}()

    @chk ccall(
               (:KSPPIPEGCRGetTruncationType, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{KSPFCDTruncationType}),
               ksp, truncstrat_,
              )

	truncstrat = truncstrat_[]

	return truncstrat
end 

"""
	unroll_w::PetscBool = KSPPIPEGCRGetUnrollW(petsclib::PetscLibType, ksp::AbstractKSP) 
Get information on `KSPPIPEGCR` if it uses unrolling the w vector

Logically Collective

Input Parameter:
- `ksp` - the Krylov space context

Output Parameter:
- `unroll_w` - `KSPPIPEGCR` uses unrolling (bool)

Level: intermediate

See also: `KSPPIPEGCR`, `KSPPIPEGCRGetTruncationType()`, `KSPPIPEGCRGetNprealloc()`, `KSPPIPEGCRSetUnrollW()`

# External Links
$(_doc_external("KSP/KSPPIPEGCRGetUnrollW"))
"""
function KSPPIPEGCRGetUnrollW(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPPIPEGCRGetUnrollW: no generated method for these argument types")
end

@for_petsc function KSPPIPEGCRGetUnrollW(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	unroll_w_ = Ref{PetscBool}()

    @chk ccall(
               (:KSPPIPEGCRGetUnrollW, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{PetscBool}),
               ksp, unroll_w_,
              )

	unroll_w = unroll_w_[]

	return unroll_w
end 

"""
	KSPPIPEGCRSetMmax(petsclib::PetscLibType, ksp::AbstractKSP, mmax::PetscInt) 
set the maximum number of previous directions `KSPPIPEGCR` will store for orthogonalization

Logically Collective

Input Parameters:
- `ksp`  - the Krylov space context
- `mmax` - the maximum number of previous directions to orthogonalize against

Options Database Key:
- `-ksp_pipegcr_mmax mmax` - maximum number of previous directions

Level: intermediate

See also: `KSPPIPEGCR`, `KSPPIPEGCRSetTruncationType()`, `KSPPIPEGCRSetNprealloc()`

# External Links
$(_doc_external("KSP/KSPPIPEGCRSetMmax"))
"""
function KSPPIPEGCRSetMmax(petsclib::PetscLibType, ksp::AbstractKSP, mmax::Integer)
    error("KSPPIPEGCRSetMmax: no generated method for these argument types")
end

@for_petsc function KSPPIPEGCRSetMmax(petsclib::$UnionPetscLib, ksp::AbstractKSP, mmax::$PetscInt )

    @chk ccall(
               (:KSPPIPEGCRSetMmax, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, mmax,
              )


	return nothing
end 

"""
	KSPPIPEGCRSetNprealloc(petsclib::PetscLibType, ksp::AbstractKSP, nprealloc::PetscInt) 
set the number of directions to preallocate with `KSPPIPEGCR`

Logically Collective

Input Parameters:
- `ksp`       - the Krylov space context
- `nprealloc` - the number of vectors to preallocate

Level: advanced

Options Database Key:
- `-ksp_pipegcr_nprealloc N` - number of vectors to preallocate

See also: `KSPPIPEGCR`, `KSPPIPEGCRGetTruncationType()`, `KSPPIPEGCRGetNprealloc()`

# External Links
$(_doc_external("KSP/KSPPIPEGCRSetNprealloc"))
"""
function KSPPIPEGCRSetNprealloc(petsclib::PetscLibType, ksp::AbstractKSP, nprealloc::Integer)
    error("KSPPIPEGCRSetNprealloc: no generated method for these argument types")
end

@for_petsc function KSPPIPEGCRSetNprealloc(petsclib::$UnionPetscLib, ksp::AbstractKSP, nprealloc::$PetscInt )

    @chk ccall(
               (:KSPPIPEGCRSetNprealloc, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, nprealloc,
              )


	return nothing
end 

"""
	KSPPIPEGCRSetTruncationType(petsclib::PetscLibType, ksp::AbstractKSP, truncstrat::KSPFCDTruncationType) 
specify how many of its stored previous directions `KSPPIPEGCR` uses during orthogonalization

Logically Collective

Input Parameters:
- `ksp`        - the Krylov space context
- `truncstrat` - the choice of strategy
``
KSP_FCD_TRUNC_TYPE_STANDARD uses all (up to `mmax`) stored directions
KSP_FCD_TRUNC_TYPE_NOTAY uses the last `max(1,mod(i,mmax))` directions at iteration i = 0, 1, ..
``

Options Database Key:
- `-ksp_pipegcr_truncation_type (standard|notay)` - which stored basis vectors to orthogonalize against

Level: intermediate

See also: `KSPPIPEGCR`, `KSPFCDTruncationType`, `KSPPIPEGCRGetTruncationType()`, `KSP_FCD_TRUNC_TYPE_STANDARD`, `KSP_FCD_TRUNC_TYPE_NOTAY`

# External Links
$(_doc_external("KSP/KSPPIPEGCRSetTruncationType"))
"""
function KSPPIPEGCRSetTruncationType(petsclib::PetscLibType, ksp::AbstractKSP, truncstrat::KSPFCDTruncationType)
    error("KSPPIPEGCRSetTruncationType: no generated method for these argument types")
end

@for_petsc function KSPPIPEGCRSetTruncationType(petsclib::$UnionPetscLib, ksp::AbstractKSP, truncstrat::KSPFCDTruncationType )

    @chk ccall(
               (:KSPPIPEGCRSetTruncationType, $petsc_library),
               PetscErrorCode,
               (CKSP, KSPFCDTruncationType),
               ksp, truncstrat,
              )


	return nothing
end 

"""
	KSPPIPEGCRSetUnrollW(petsclib::PetscLibType, ksp::AbstractKSP, unroll_w::PetscBool) 
Set to `PETSC_TRUE` to use `KSPPIPEGCR` with unrolling of the w vector

Logically Collective

Input Parameters:
- `ksp`      - the Krylov space context
- `unroll_w` - use unrolling

Level: intermediate

Options Database Key:
- `-ksp_pipegcr_unroll_w (true|false)` - use unrolling

See also: `KSPPIPEGCR`, `KSPPIPEGCRSetTruncationType()`, `KSPPIPEGCRSetNprealloc()`, `KSPPIPEGCRGetUnrollW()`

# External Links
$(_doc_external("KSP/KSPPIPEGCRSetUnrollW"))
"""
function KSPPIPEGCRSetUnrollW(petsclib::PetscLibType, ksp::AbstractKSP, unroll_w::PetscBool)
    error("KSPPIPEGCRSetUnrollW: no generated method for these argument types")
end

@for_petsc function KSPPIPEGCRSetUnrollW(petsclib::$UnionPetscLib, ksp::AbstractKSP, unroll_w::PetscBool )

    @chk ccall(
               (:KSPPIPEGCRSetUnrollW, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, unroll_w,
              )


	return nothing
end 

"""
	pyname::String = KSPPythonGetType(petsclib::PetscLibType, ksp::AbstractKSP) 
Get the type of a `KSP` object implemented in Python.

Not Collective

Input Parameter:
- `ksp`  - the linear solver `KSP` context.

Output Parameter:
- `pyname`  - full dotted Python name [package].module[.{class|function}]

Level: intermediate

See also: `KSPCreate()`, `KSPSetType()`, `KSPPYTHON`, `PetscPythonInitialize()`, `KSPPythonSetType()`

# External Links
$(_doc_external("KSP/KSPPythonGetType"))
"""
function KSPPythonGetType(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPPythonGetType: no generated method for these argument types")
end

@for_petsc function KSPPythonGetType(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	pyname_ = Ref{Ptr{Cchar}}()

    @chk ccall(
               (:KSPPythonGetType, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Ptr{Cchar}}),
               ksp, pyname_,
              )

	pyname = pyname_[] == C_NULL ? "" : unsafe_string(pyname_[])

	return pyname
end 

"""
	KSPPythonSetType(petsclib::PetscLibType, ksp::AbstractKSP, pyname::String) 
Initialize a `KSP` object to a type implemented in Python.

Collective

Input Parameters:
- `ksp`  - the linear solver `KSP` context.
- `pyname`  - full dotted Python name [package].module[.{class|function}]

Options Database Key:
- `-ksp_python_type pyname`  - python class

Level: intermediate

See also: `KSPCreate()`, `KSPSetType()`, `KSPPYTHON`, `PetscPythonInitialize()`

# External Links
$(_doc_external("KSP/KSPPythonSetType"))
"""
function KSPPythonSetType(petsclib::PetscLibType, ksp::AbstractKSP, pyname::String)
    error("KSPPythonSetType: no generated method for these argument types")
end

@for_petsc function KSPPythonSetType(petsclib::$UnionPetscLib, ksp::AbstractKSP, pyname::String )

    @chk ccall(
               (:KSPPythonSetType, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Cchar}),
               ksp, pyname,
              )


	return nothing
end 

"""
	quadratic::PetscReal = KSPQCGGetQuadratic(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets the value of the quadratic function, evaluated at the new iterate

Collective

Input Parameter:
- `ksp` - the iterative context

Output Parameter:
- `quadratic` - the quadratic function evaluated at the new iterate

Level: advanced

See also: `KSPQCG`

# External Links
$(_doc_external("KSP/KSPQCGGetQuadratic"))
"""
function KSPQCGGetQuadratic(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPQCGGetQuadratic: no generated method for these argument types")
end

@for_petsc function KSPQCGGetQuadratic(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	quadratic_ = Ref{$PetscReal}()

    @chk ccall(
               (:KSPQCGGetQuadratic, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscReal}),
               ksp, quadratic_,
              )

	quadratic = quadratic_[]

	return quadratic
end 

"""
	tsnorm::PetscReal = KSPQCGGetTrialStepNorm(petsclib::PetscLibType, ksp::AbstractKSP) 
Gets the norm of a trial step vector in `KSPQCG`.  The WCG step may be
constrained, so this is not necessarily the length of the ultimate step taken in `KSPQCG`.

Not Collective

Input Parameter:
- `ksp` - the iterative context

Output Parameter:
- `tsnorm` - the norm

Level: advanced

See also: `KSPQCG`, `KSPQCGSetTrustRegionRadius()`

# External Links
$(_doc_external("KSP/KSPQCGGetTrialStepNorm"))
"""
function KSPQCGGetTrialStepNorm(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPQCGGetTrialStepNorm: no generated method for these argument types")
end

@for_petsc function KSPQCGGetTrialStepNorm(petsclib::$UnionPetscLib, ksp::AbstractKSP )
	tsnorm_ = Ref{$PetscReal}()

    @chk ccall(
               (:KSPQCGGetTrialStepNorm, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscReal}),
               ksp, tsnorm_,
              )

	tsnorm = tsnorm_[]

	return tsnorm
end 

"""
	KSPQCGSetTrustRegionRadius(petsclib::PetscLibType, ksp::AbstractKSP, delta::PetscReal) 
Sets the radius of the trust region for `KSPQCG`

Logically Collective

Input Parameters:
- `ksp`   - the iterative context
- `delta` - the trust region radius (Infinity is the default)

Options Database Key:
- `-ksp_qcg_trustregionradius delta` - trust region radius

Level: advanced

See also: `KSPQCG`, `KSPQCGGetTrialStepNorm()`

# External Links
$(_doc_external("KSP/KSPQCGSetTrustRegionRadius"))
"""
function KSPQCGSetTrustRegionRadius(petsclib::PetscLibType, ksp::AbstractKSP, delta::Real)
    error("KSPQCGSetTrustRegionRadius: no generated method for these argument types")
end

@for_petsc function KSPQCGSetTrustRegionRadius(petsclib::$UnionPetscLib, ksp::AbstractKSP, delta::$PetscReal )

    @chk ccall(
               (:KSPQCGSetTrustRegionRadius, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscReal),
               ksp, delta,
              )


	return nothing
end 

"""
	KSPRegister(petsclib::PetscLibType, sname::String, fnc::external) 
Adds a method, `KSPType`, to the Krylov subspace solver package.

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - name of a new user-defined solver
- `function` - routine to create method

Level: advanced

See also: `KSP`, `KSPType`, `KSPSetType`, `KSPRegisterAll()`

# External Links
$(_doc_external("KSP/KSPRegister"))
"""
function KSPRegister(petsclib::PetscLibType, sname::String, fnc::external)
    error("KSPRegister: no generated method for these argument types")
end

@for_petsc function KSPRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:KSPRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	KSPReset(petsclib::PetscLibType, ksp::AbstractKSP) 
Removes any allocated `Vec` and `Mat` from the `KSP` data structures.

Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Level: intermediate

See also: `KSPCreate()`, `KSPSetUp()`, `KSPSolve()`, `KSP`

# External Links
$(_doc_external("KSP/KSPReset"))
"""
function KSPReset(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPReset: no generated method for these argument types")
end

@for_petsc function KSPReset(petsclib::$UnionPetscLib, ksp::AbstractKSP )

    @chk ccall(
               (:KSPReset, $petsc_library),
               PetscErrorCode,
               (CKSP,),
               ksp,
              )


	return nothing
end 

"""
	KSPResetFromOptions(petsclib::PetscLibType, ksp::AbstractKSP) 
Sets `KSP` parameters from user options ONLY if the `KSP` was previously set from options

Collective

Input Parameter:
- `ksp` - the `KSP` context

Level: advanced

See also: `KSPSetFromOptions()`, `KSPSetOptionsPrefix()`

# External Links
$(_doc_external("KSP/KSPResetFromOptions"))
"""
function KSPResetFromOptions(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPResetFromOptions: no generated method for these argument types")
end

@for_petsc function KSPResetFromOptions(petsclib::$UnionPetscLib, ksp::AbstractKSP )

    @chk ccall(
               (:KSPResetFromOptions, $petsc_library),
               PetscErrorCode,
               (CKSP,),
               ksp,
              )


	return nothing
end 

"""
	KSPResetViewers(petsclib::PetscLibType, ksp::AbstractKSP) 
Resets all the viewers set from the options database during `KSPSetFromOptions()`

Collective

Input Parameter:
- `ksp` - the `KSP` iterative solver context obtained from `KSPCreate()`

Level: beginner

See also: `KSPCreate()`, `KSPSetUp()`, `KSPSolve()`, `KSPSetFromOptions()`, `KSP`

# External Links
$(_doc_external("KSP/KSPResetViewers"))
"""
function KSPResetViewers(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPResetViewers: no generated method for these argument types")
end

@for_petsc function KSPResetViewers(petsclib::$UnionPetscLib, ksp::AbstractKSP )

    @chk ccall(
               (:KSPResetViewers, $petsc_library),
               PetscErrorCode,
               (CKSP,),
               ksp,
              )


	return nothing
end 

"""
	KSPRichardsonSetScale(petsclib::PetscLibType, ksp::AbstractKSP, scale::PetscReal) 
Set the damping factor; if this routine is not called, the factor defaults to 1.0.

Logically Collective

Input Parameters:
- `ksp`   - the iterative context
- `scale` - the damping factor

Options Database Key:
- `-ksp_richardson_scale scale` - Set the scale factor

Level: intermediate

See also: `KSPRICHARDSON`, `KSPRichardsonSetSelfScale()`

# External Links
$(_doc_external("KSP/KSPRichardsonSetScale"))
"""
function KSPRichardsonSetScale(petsclib::PetscLibType, ksp::AbstractKSP, scale::Real)
    error("KSPRichardsonSetScale: no generated method for these argument types")
end

@for_petsc function KSPRichardsonSetScale(petsclib::$UnionPetscLib, ksp::AbstractKSP, scale::$PetscReal )

    @chk ccall(
               (:KSPRichardsonSetScale, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscReal),
               ksp, scale,
              )


	return nothing
end 

"""
	KSPRichardsonSetSelfScale(petsclib::PetscLibType, ksp::AbstractKSP, scale::PetscBool) 
Sets Richardson to automatically determine optimal scaling at each iteration to minimize the 2-norm of the
preconditioned residual

Logically Collective

Input Parameters:
- `ksp`   - the iterative context
- `scale` - `PETSC_TRUE` or the default of `PETSC_FALSE`

Options Database Key:
- `-ksp_richardson_self_scale` - Use self-scaling

Level: intermediate

See also: `KSPRICHARDSON`, `KSPRichardsonSetScale()`

# External Links
$(_doc_external("KSP/KSPRichardsonSetSelfScale"))
"""
function KSPRichardsonSetSelfScale(petsclib::PetscLibType, ksp::AbstractKSP, scale::PetscBool)
    error("KSPRichardsonSetSelfScale: no generated method for these argument types")
end

@for_petsc function KSPRichardsonSetSelfScale(petsclib::$UnionPetscLib, ksp::AbstractKSP, scale::PetscBool )

    @chk ccall(
               (:KSPRichardsonSetSelfScale, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, scale,
              )


	return nothing
end 

"""
	KSPSetApplicationContext(petsclib::PetscLibType, ksp::AbstractKSP, ctx::Ptr{Cvoid}) 
Sets the optional user-defined context for the linear solver.

Logically Collective

Input Parameters:
- `ksp` - the `KSP` context
- `ctx` - application context

Level: intermediate

See also: `KSP`, `KSPGetApplicationContext()`

# External Links
$(_doc_external("KSP/KSPSetApplicationContext"))
"""
function KSPSetApplicationContext(petsclib::PetscLibType, ksp::AbstractKSP, ctx::Ptr{Cvoid})
    error("KSPSetApplicationContext: no generated method for these argument types")
end

@for_petsc function KSPSetApplicationContext(petsclib::$UnionPetscLib, ksp::AbstractKSP, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:KSPSetApplicationContext, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Cvoid}),
               ksp, ctx,
              )


	return nothing
end 

"""
	KSPSetCheckNormIteration(petsclib::PetscLibType, ksp::AbstractKSP, it::PetscInt) 
Sets the first iteration at which the norm of the residual will be
computed and used in the convergence test of `KSPSolve()` for the given `KSP` context

Logically Collective

Input Parameters:
- `ksp` - Krylov solver context
- `it`  - use -1 to check at all iterations

Level: advanced

See also: `KSP`, `KSPSetUp()`, `KSPSolve()`, `KSPDestroy()`, `KSPConvergedSkip()`, `KSPSetNormType()`, `KSPSetLagNorm()`

# External Links
$(_doc_external("KSP/KSPSetCheckNormIteration"))
"""
function KSPSetCheckNormIteration(petsclib::PetscLibType, ksp::AbstractKSP, it::Integer)
    error("KSPSetCheckNormIteration: no generated method for these argument types")
end

@for_petsc function KSPSetCheckNormIteration(petsclib::$UnionPetscLib, ksp::AbstractKSP, it::$PetscInt )

    @chk ccall(
               (:KSPSetCheckNormIteration, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, it,
              )


	return nothing
end 

"""
	KSPSetComputeEigenvalues(petsclib::PetscLibType, ksp::AbstractKSP, flg::PetscBool) 
Sets a flag so that the extreme eigenvalues
values will be calculated via a Lanczos or Arnoldi process as the linear
system is solved.

Logically Collective

Input Parameters:
- `ksp` - iterative solver obtained from `KSPCreate()`
- `flg` - `PETSC_TRUE` or `PETSC_FALSE`

Level: advanced

See also: `KSPComputeEigenvalues()`, `KSPComputeEigenvaluesExplicitly()`, `KSP`, `KSPSetComputeRitz()`

# External Links
$(_doc_external("KSP/KSPSetComputeEigenvalues"))
"""
function KSPSetComputeEigenvalues(petsclib::PetscLibType, ksp::AbstractKSP, flg::PetscBool)
    error("KSPSetComputeEigenvalues: no generated method for these argument types")
end

@for_petsc function KSPSetComputeEigenvalues(petsclib::$UnionPetscLib, ksp::AbstractKSP, flg::PetscBool )

    @chk ccall(
               (:KSPSetComputeEigenvalues, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flg,
              )


	return nothing
end 

"""
	KSPSetComputeInitialGuess(petsclib::PetscLibType, ksp::AbstractKSP, func::Ptr{Cvoid}, ctx::Ptr{Cvoid}) 
set routine to compute the initial guess of the linear system

Logically Collective

Input Parameters:
- `ksp`  - the `KSP` context
- `func` - function to compute the initial guess, see `KSPComputeInitialGuessFn` for calling sequence
- `ctx`  - optional context

Level: beginner

See also: `KSP`, `KSPSolve()`, `KSPSetComputeRHS()`, `KSPSetComputeOperators()`, `DMKSPSetComputeInitialGuess()`, `KSPSetInitialGuessNonzero()`,
`KSPComputeInitialGuessFn`

# External Links
$(_doc_external("KSP/KSPSetComputeInitialGuess"))
"""
function KSPSetComputeInitialGuess(petsclib::PetscLibType, ksp::AbstractKSP, func::Ptr{Cvoid}, ctx::Ptr{Cvoid})
    error("KSPSetComputeInitialGuess: no generated method for these argument types")
end

@for_petsc function KSPSetComputeInitialGuess(petsclib::$UnionPetscLib, ksp::AbstractKSP, func::Ptr{Cvoid}, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:KSPSetComputeInitialGuess, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Cvoid}, Ptr{Cvoid}),
               ksp, func, ctx,
              )


	return nothing
end 

"""
	KSPSetComputeOperators(petsclib::PetscLibType, ksp::AbstractKSP, func::Ptr{Cvoid}, ctx::Ptr{Cvoid}) 
set routine to compute the linear operators

Logically Collective

Input Parameters:
- `ksp`  - the `KSP` context
- `func` - function to compute the operators, see `KSPComputeOperatorsFn` for the calling sequence
- `ctx`  - optional context

Level: beginner

See also: `KSP`, `KSPSetOperators()`, `KSPSetComputeRHS()`, `DMKSPSetComputeOperators()`, `KSPSetComputeInitialGuess()`, `KSPComputeOperatorsFn`

# External Links
$(_doc_external("KSP/KSPSetComputeOperators"))
"""
function KSPSetComputeOperators(petsclib::PetscLibType, ksp::AbstractKSP, func::Ptr{Cvoid}, ctx::Ptr{Cvoid})
    error("KSPSetComputeOperators: no generated method for these argument types")
end

@for_petsc function KSPSetComputeOperators(petsclib::$UnionPetscLib, ksp::AbstractKSP, func::Ptr{Cvoid}, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:KSPSetComputeOperators, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Cvoid}, Ptr{Cvoid}),
               ksp, func, ctx,
              )


	return nothing
end 

"""
	KSPSetComputeRHS(petsclib::PetscLibType, ksp::AbstractKSP, func::Ptr{Cvoid}, ctx::Ptr{Cvoid}) 
set routine to compute the right-hand side of the linear system

Logically Collective

Input Parameters:
- `ksp`  - the `KSP` context
- `func` - function to compute the right-hand side, see `KSPComputeRHSFn` for the calling sequence
- `ctx`  - optional context

Level: beginner

See also: `KSP`, `KSPSolve()`, `DMKSPSetComputeRHS()`, `KSPSetComputeOperators()`, `KSPSetOperators()`, `KSPComputeRHSFn`

# External Links
$(_doc_external("KSP/KSPSetComputeRHS"))
"""
function KSPSetComputeRHS(petsclib::PetscLibType, ksp::AbstractKSP, func::Ptr{Cvoid}, ctx::Ptr{Cvoid})
    error("KSPSetComputeRHS: no generated method for these argument types")
end

@for_petsc function KSPSetComputeRHS(petsclib::$UnionPetscLib, ksp::AbstractKSP, func::Ptr{Cvoid}, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:KSPSetComputeRHS, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Cvoid}, Ptr{Cvoid}),
               ksp, func, ctx,
              )


	return nothing
end 

"""
	KSPSetComputeRitz(petsclib::PetscLibType, ksp::AbstractKSP, flg::PetscBool) 
Sets a flag so that the Ritz or harmonic Ritz pairs
will be calculated via a Lanczos or Arnoldi process as the linear
system is solved.

Logically Collective

Input Parameters:
- `ksp` - iterative solver obtained from `KSPCreate()`
- `flg` - `PETSC_TRUE` or `PETSC_FALSE`

Level: advanced

See also: `KSPComputeRitz()`, `KSP`, `KSPComputeEigenvalues()`, `KSPComputeExtremeSingularValues()`

# External Links
$(_doc_external("KSP/KSPSetComputeRitz"))
"""
function KSPSetComputeRitz(petsclib::PetscLibType, ksp::AbstractKSP, flg::PetscBool)
    error("KSPSetComputeRitz: no generated method for these argument types")
end

@for_petsc function KSPSetComputeRitz(petsclib::$UnionPetscLib, ksp::AbstractKSP, flg::PetscBool )

    @chk ccall(
               (:KSPSetComputeRitz, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flg,
              )


	return nothing
end 

"""
	KSPSetComputeSingularValues(petsclib::PetscLibType, ksp::AbstractKSP, flg::PetscBool) 
Sets a flag so that the extreme singular
values will be calculated via a Lanczos or Arnoldi process as the linear
system is solved.

Logically Collective

Input Parameters:
- `ksp` - iterative solver obtained from `KSPCreate()`
- `flg` - `PETSC_TRUE` or `PETSC_FALSE`

Options Database Key:
- `-ksp_monitor_singular_value` - Activates `KSPSetComputeSingularValues()`

Level: advanced

See also: `KSPComputeExtremeSingularValues()`, `KSPMonitorSingularValue()`, `KSP`, `KSPSetComputeRitz()`

# External Links
$(_doc_external("KSP/KSPSetComputeSingularValues"))
"""
function KSPSetComputeSingularValues(petsclib::PetscLibType, ksp::AbstractKSP, flg::PetscBool)
    error("KSPSetComputeSingularValues: no generated method for these argument types")
end

@for_petsc function KSPSetComputeSingularValues(petsclib::$UnionPetscLib, ksp::AbstractKSP, flg::PetscBool )

    @chk ccall(
               (:KSPSetComputeSingularValues, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flg,
              )


	return nothing
end 

"""
	KSPSetConvergedNegativeCurvature(petsclib::PetscLibType, ksp::AbstractKSP, flg::PetscBool) 
Allows to declare convergence and return `KSP_CONVERGED_NEG_CURVE` when negative curvature is detected

Collective

Input Parameters:
- `ksp` - iterative context
- `flg` - the Boolean value

Options Database Key:
- `-ksp_converged_neg_curve (true|false)` - Declare convergence if negative curvature is detected

Level: advanced

See also: `KSP`, `KSPConvergedReason`, `KSPGetConvergedNegativeCurvature()`

# External Links
$(_doc_external("KSP/KSPSetConvergedNegativeCurvature"))
"""
function KSPSetConvergedNegativeCurvature(petsclib::PetscLibType, ksp::AbstractKSP, flg::PetscBool)
    error("KSPSetConvergedNegativeCurvature: no generated method for these argument types")
end

@for_petsc function KSPSetConvergedNegativeCurvature(petsclib::$UnionPetscLib, ksp::AbstractKSP, flg::PetscBool )

    @chk ccall(
               (:KSPSetConvergedNegativeCurvature, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flg,
              )


	return nothing
end 

"""
	KSPSetConvergenceTest(petsclib::PetscLibType, ksp::AbstractKSP, converge::Ptr{Cvoid}, ctx::Ptr{Cvoid}, destroy::Ptr{Cvoid}) 
Sets the function to be used to determine convergence of `KSPSolve()`

Logically Collective

Input Parameters:
- `ksp`      - iterative solver obtained from `KSPCreate()`
- `converge` - pointer to the function, see `KSPConvergenceTestFn`
- `ctx`      - context for private data for the convergence routine (may be `NULL`)
- `destroy`  - a routine for destroying the context (may be `NULL`)

Level: advanced

See also: `KSP`, `KSPConvergenceTestFn`, `KSPConvergedDefault()`, `KSPGetConvergenceContext()`, `KSPSetTolerances()`, `KSPGetConvergenceTest()`, `KSPGetAndClearConvergenceTest()`

# External Links
$(_doc_external("KSP/KSPSetConvergenceTest"))
"""
function KSPSetConvergenceTest(petsclib::PetscLibType, ksp::AbstractKSP, converge::Ptr{Cvoid}, ctx::Ptr{Cvoid}, destroy::Ptr{Cvoid})
    error("KSPSetConvergenceTest: no generated method for these argument types")
end

@for_petsc function KSPSetConvergenceTest(petsclib::$UnionPetscLib, ksp::AbstractKSP, converge::Ptr{Cvoid}, ctx::Ptr{Cvoid}, destroy::Ptr{Cvoid} )

    @chk ccall(
               (:KSPSetConvergenceTest, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
               ksp, converge, ctx, destroy,
              )


	return nothing
end 

"""
	KSPSetDM(petsclib::PetscLibType, ksp::AbstractKSP, dm::AbstractPetscDM) 
Sets the `DM` that may be used by some preconditioners and that may be used to construct the linear system

Logically Collective

Input Parameters:
- `ksp` - the `KSP`
- `dm`  - the `DM`, cannot be `NULL` to remove a previously set `DM`

Level: intermediate

See also: `KSP`, `DM`, `KSPGetDM()`, `KSPSetDMActive()`, `KSPSetComputeOperators()`, `KSPSetComputeRHS()`, `KSPSetComputeInitialGuess()`, `DMKSPSetComputeOperators()`, `DMKSPSetComputeRHS()`, `DMKSPSetComputeInitialGuess()`

# External Links
$(_doc_external("KSP/KSPSetDM"))
"""
function KSPSetDM(petsclib::PetscLibType, ksp::AbstractKSP, dm::AbstractPetscDM)
    error("KSPSetDM: no generated method for these argument types")
end

@for_petsc function KSPSetDM(petsclib::$UnionPetscLib, ksp::AbstractKSP, dm::AbstractPetscDM )

    @chk ccall(
               (:KSPSetDM, $petsc_library),
               PetscErrorCode,
               (CKSP, CDM),
               ksp, dm,
              )


	return nothing
end 

"""
	KSPSetDMActive(petsclib::PetscLibType, ksp::AbstractKSP, active::KSPDMActive, flg::PetscBool) 
Indicates the `DM` should be used to generate the linear system matrix, the right-hand side vector, and the initial guess

Logically Collective

Input Parameters:
- `ksp`    - the `KSP`
- `active` - one of `KSP_DMACTIVE_OPERATOR`, `KSP_DMACTIVE_RHS`, or `KSP_DMACTIVE_INITIAL_GUESS`
- `flg`    - use the `DM`

Level: intermediate

See also: `KSP`, `DM`, `KSPGetDM()`, `KSPSetDM()`, `SNESSetDM()`, `KSPSetComputeOperators()`, `KSPSetComputeRHS()`, `KSPSetComputeInitialGuess()`

# External Links
$(_doc_external("KSP/KSPSetDMActive"))
"""
function KSPSetDMActive(petsclib::PetscLibType, ksp::AbstractKSP, active::KSPDMActive, flg::PetscBool)
    error("KSPSetDMActive: no generated method for these argument types")
end

@for_petsc function KSPSetDMActive(petsclib::$UnionPetscLib, ksp::AbstractKSP, active::KSPDMActive, flg::PetscBool )

    @chk ccall(
               (:KSPSetDMActive, $petsc_library),
               PetscErrorCode,
               (CKSP, KSPDMActive, PetscBool),
               ksp, active, flg,
              )


	return nothing
end 

"""
	KSPSetDiagonalScale(petsclib::PetscLibType, ksp::AbstractKSP, scale::PetscBool) 
Tells `KSP` to symmetrically diagonally scale the system
before solving. This actually CHANGES the matrix (and right-hand side).

Logically Collective

Input Parameters:
- `ksp`   - the `KSP` context
- `scale` - `PETSC_TRUE` or `PETSC_FALSE`

Options Database Keys:
- `-ksp_diagonal_scale`     - perform a diagonal scaling before the solve
- `-ksp_diagonal_scale_fix` - scale the matrix back AFTER the solve

Level: advanced

See also: `KSPGetDiagonalScale()`, `KSPSetDiagonalScaleFix()`, `KSP`

# External Links
$(_doc_external("KSP/KSPSetDiagonalScale"))
"""
function KSPSetDiagonalScale(petsclib::PetscLibType, ksp::AbstractKSP, scale::PetscBool)
    error("KSPSetDiagonalScale: no generated method for these argument types")
end

@for_petsc function KSPSetDiagonalScale(petsclib::$UnionPetscLib, ksp::AbstractKSP, scale::PetscBool )

    @chk ccall(
               (:KSPSetDiagonalScale, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, scale,
              )


	return nothing
end 

"""
	KSPSetDiagonalScaleFix(petsclib::PetscLibType, ksp::AbstractKSP, fix::PetscBool) 
Tells `KSP` to diagonally scale the system back after solving.

Logically Collective

Input Parameters:
- `ksp` - the `KSP` context
- `fix` - `PETSC_TRUE` to scale back after the system solve, `PETSC_FALSE` to not
rescale (default)

Level: intermediate

See also: `KSPGetDiagonalScale()`, `KSPSetDiagonalScale()`, `KSPGetDiagonalScaleFix()`, `KSP`

# External Links
$(_doc_external("KSP/KSPSetDiagonalScaleFix"))
"""
function KSPSetDiagonalScaleFix(petsclib::PetscLibType, ksp::AbstractKSP, fix::PetscBool)
    error("KSPSetDiagonalScaleFix: no generated method for these argument types")
end

@for_petsc function KSPSetDiagonalScaleFix(petsclib::$UnionPetscLib, ksp::AbstractKSP, fix::PetscBool )

    @chk ccall(
               (:KSPSetDiagonalScaleFix, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, fix,
              )


	return nothing
end 

"""
	KSPSetErrorHistory(petsclib::PetscLibType, ksp::AbstractKSP, a::Vector{PetscReal}, na::PetscCount, reset::PetscBool) 
Sets the array used to hold the error history. If set, this array will contain the error norms computed at each iteration of the solver.

Not Collective

Input Parameters:
- `ksp`   - iterative solver obtained from `KSPCreate()`
- `a`     - array to hold history
- `na`    - size of `a`
- `reset` - `PETSC_TRUE` indicates the history counter is reset to zero for each new linear solve

Level: advanced

See also: `KSPGetErrorHistory()`, `KSPSetResidualHistory()`, `KSP`

# External Links
$(_doc_external("KSP/KSPSetErrorHistory"))
"""
function KSPSetErrorHistory(petsclib::PetscLibType, ksp::AbstractKSP, a::AbstractVector{<:Number}, na::PetscCount, reset::PetscBool)
    error("KSPSetErrorHistory: no generated method for these argument types")
end

@for_petsc function KSPSetErrorHistory(petsclib::$UnionPetscLib, ksp::AbstractKSP, a::Vector{$PetscReal}, na::PetscCount, reset::PetscBool )

    @chk ccall(
               (:KSPSetErrorHistory, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscReal}, PetscCount, PetscBool),
               ksp, a, na, reset,
              )


	return nothing
end 

"""
	KSPSetErrorIfNotConverged(petsclib::PetscLibType, ksp::AbstractKSP, flg::PetscBool) 
Causes `KSPSolve()` to generate an error if the solver has not converged as soon as the error is detected.

Logically Collective

Input Parameters:
- `ksp` - iterative solver obtained from `KSPCreate()`
- `flg` - `PETSC_TRUE` indicates you want the error generated

Options Database Key:
- `-ksp_error_if_not_converged (true|false)` - generate an error and stop the program

Level: intermediate

See also: `KSPGetErrorIfNotConverged()`, `KSP`

# External Links
$(_doc_external("KSP/KSPSetErrorIfNotConverged"))
"""
function KSPSetErrorIfNotConverged(petsclib::PetscLibType, ksp::AbstractKSP, flg::PetscBool)
    error("KSPSetErrorIfNotConverged: no generated method for these argument types")
end

@for_petsc function KSPSetErrorIfNotConverged(petsclib::$UnionPetscLib, ksp::AbstractKSP, flg::PetscBool )

    @chk ccall(
               (:KSPSetErrorIfNotConverged, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flg,
              )


	return nothing
end 

"""
	KSPSetFromOptions(petsclib::PetscLibType, ksp::AbstractKSP) 
Sets `KSP` options from the options database.
This routine must be called before `KSPSetUp()` if the user is to be
allowed to set the Krylov type.

Collective

Input Parameter:
- `ksp` - the Krylov space context

Options Database Keys:
- `-ksp_rtol rtol`                                                          - relative tolerance used in default determination of convergence, i.e.
if residual norm decreases by this factor than convergence is declared
- `-ksp_atol abstol`                                                        - absolute tolerance used in default convergence test, i.e. if residual
norm is less than this then convergence is declared
- `-ksp_divtol tol`                                                         - if residual norm increases by this factor than divergence is declared
- `-ksp_max_it maxits`                                                      - maximum number of linear iterations
- `-ksp_min_it minits`                                                      - minimum number of linear iterations to use, defaults to zero
- `-ksp_reuse_preconditioner (true|false)`                                  - reuse the previously computed preconditioner
- `-ksp_converged_use_initial_residual_norm`                                - see `KSPConvergedDefaultSetUIRNorm()`
- `-ksp_converged_use_min_initial_residual_norm`                            - see `KSPConvergedDefaultSetUMIRNorm()`
- `-ksp_converged_maxits`                                                   - see `KSPConvergedDefaultSetConvergedMaxits()`
- `-ksp_norm_type (none|preconditioned|unpreconditioned|natural)`           - see `KSPSetNormType()`
- `-ksp_check_norm_iteration it`                                            - do not compute residual norm until iteration number it (does compute at 0th iteration)
works only for `KSPBCGS`, `KSPIBCGS`, and `KSPCG`
- `-ksp_lag_norm`                                                           - compute the norm of the residual for the ith iteration on the i+1 iteration;
this means that one can use the norm of the residual for convergence test WITHOUT
an extra `MPI_Allreduce()` limiting global synchronizations.
This will require 1 more iteration of the solver than usual.
- `-ksp_guess_type`                                                         - Type of initial guess generator for repeated linear solves
- `-ksp_fischer_guess model,size`                                           - uses the Fischer initial guess generator for repeated linear solves
- `-ksp_constant_null_space`                                                - assume the operator (matrix) has the constant vector in its null space
- `-ksp_test_null_space`                                                    - tests the null space set with `MatSetNullSpace()` to see if it truly is a null space
- `-ksp_knoll`                                                              - compute initial guess by applying the preconditioner to the right-hand side
- `-ksp_monitor_cancel`                                                     - cancel all previous convergene monitor routines set
- `-ksp_monitor`                                                            - print residual norm at each iteration
- `-ksp_monitor draw::draw_lg`                                              - plot residual norm at each iteration, see `KSPMonitorResidual()`
- `-ksp_monitor_true_residual`                                              - print the true l2 residual norm at each iteration, see `KSPMonitorTrueResidual()`
- `-all_ksp_monitor optional_filename`                                      - print residual norm at each iteration for ALL KSP solves, regardless of their prefix. This is
useful for `PCFIELDSPLIT`, `PCMG`, etc that have inner solvers and
you wish to track the convergence of all the solvers
- `-ksp_monitor_solution [ascii binary or draw][:filename][:format option]` - plot solution at each iteration
- `-ksp_monitor_singular_value`                                             - monitor extreme singular values at each iteration
- `-ksp_converged_reason`                                                   - view the convergence state at the end of the solve
- `-ksp_use_explicittranspose`                                              - transpose the system explicitly in `KSPSolveTranspose()`
- `-ksp_error_if_not_converged`                                             - stop the program as soon as an error is detected in a `KSPSolve()`, `KSP_DIVERGED_ITS`
is not treated as an error on inner solves
- `-ksp_converged_rate`                                                     - view the convergence rate at the end of the solve

Level: beginner

See also: `KSP`, `KSPSetOptionsPrefix()`, `KSPResetFromOptions()`, `KSPSetUseFischerGuess()`

# External Links
$(_doc_external("KSP/KSPSetFromOptions"))
"""
function KSPSetFromOptions(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPSetFromOptions: no generated method for these argument types")
end

@for_petsc function KSPSetFromOptions(petsclib::$UnionPetscLib, ksp::AbstractKSP )

    @chk ccall(
               (:KSPSetFromOptions, $petsc_library),
               PetscErrorCode,
               (CKSP,),
               ksp,
              )


	return nothing
end 

"""
	KSPSetGuess(petsclib::PetscLibType, ksp::AbstractKSP, guess::KSPGuess) 
Set the initial guess object `KSPGuess` to be used by the `KSP` object to generate initial guesses

Logically Collective

Input Parameters:
- `ksp`   - the Krylov context
- `guess` - the object created with `KSPGuessCreate()`

Level: advanced

See also: `KSP`, `KSPGuess`, `KSPSetOptionsPrefix()`, `KSPAppendOptionsPrefix()`, `KSPSetUseFischerGuess()`, `KSPGetGuess()`

# External Links
$(_doc_external("KSP/KSPSetGuess"))
"""
function KSPSetGuess(petsclib::PetscLibType, ksp::AbstractKSP, guess::KSPGuess)
    error("KSPSetGuess: no generated method for these argument types")
end

@for_petsc function KSPSetGuess(petsclib::$UnionPetscLib, ksp::AbstractKSP, guess::KSPGuess )

    @chk ccall(
               (:KSPSetGuess, $petsc_library),
               PetscErrorCode,
               (CKSP, KSPGuess),
               ksp, guess,
              )


	return nothing
end 

"""
	KSPSetInitialGuessKnoll(petsclib::PetscLibType, ksp::AbstractKSP, flg::PetscBool) 
Tells the iterative solver to use `PCApply()` on the right hand side vector to compute the initial guess (The Knoll trick)

Logically Collective

Input Parameters:
- `ksp` - iterative solver obtained from `KSPCreate()`
- `flg` - `PETSC_TRUE` or `PETSC_FALSE`

Level: advanced

See also: `KSPGetInitialGuessKnoll()`, `KSPGuess`, `KSPSetInitialGuessNonzero()`, `KSPGetInitialGuessNonzero()`, `KSP`

# External Links
$(_doc_external("KSP/KSPSetInitialGuessKnoll"))
"""
function KSPSetInitialGuessKnoll(petsclib::PetscLibType, ksp::AbstractKSP, flg::PetscBool)
    error("KSPSetInitialGuessKnoll: no generated method for these argument types")
end

@for_petsc function KSPSetInitialGuessKnoll(petsclib::$UnionPetscLib, ksp::AbstractKSP, flg::PetscBool )

    @chk ccall(
               (:KSPSetInitialGuessKnoll, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flg,
              )


	return nothing
end 

"""
	KSPSetInitialGuessNonzero(petsclib::PetscLibType, ksp::AbstractKSP, flg::PetscBool) 
Tells the iterative solver that the
initial guess is nonzero; otherwise `KSP` assumes the initial guess
is to be zero (and thus zeros it out before solving).

Logically Collective

Input Parameters:
- `ksp` - iterative solver obtained from `KSPCreate()`
- `flg` - `PETSC_TRUE` indicates the guess is non-zero, `PETSC_FALSE` indicates the guess is zero

Options Database Key:
- `-ksp_initial_guess_nonzero (true|false)` - use nonzero initial guess

Level: beginner

See also: `KSPGetInitialGuessNonzero()`, `KSPGuessSetType()`, `KSPGuessType`, `KSP`

# External Links
$(_doc_external("KSP/KSPSetInitialGuessNonzero"))
"""
function KSPSetInitialGuessNonzero(petsclib::PetscLibType, ksp::AbstractKSP, flg::PetscBool)
    error("KSPSetInitialGuessNonzero: no generated method for these argument types")
end

@for_petsc function KSPSetInitialGuessNonzero(petsclib::$UnionPetscLib, ksp::AbstractKSP, flg::PetscBool )

    @chk ccall(
               (:KSPSetInitialGuessNonzero, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flg,
              )


	return nothing
end 

"""
	KSPSetLagNorm(petsclib::PetscLibType, ksp::AbstractKSP, flg::PetscBool) 
Lags the residual norm calculation so that it is computed as part of the `MPI_Allreduce()` used for
computing the inner products needed for the next iteration.

Logically Collective

Input Parameters:
- `ksp` - Krylov solver context
- `flg` - `PETSC_TRUE` or `PETSC_FALSE`

Options Database Key:
- `-ksp_lag_norm` - lag the calculated residual norm

Level: advanced

See also: `KSPSetUp()`, `KSPSolve()`, `KSPDestroy()`, `KSPConvergedSkip()`, `KSPSetNormType()`, `KSPSetCheckNormIteration()`

# External Links
$(_doc_external("KSP/KSPSetLagNorm"))
"""
function KSPSetLagNorm(petsclib::PetscLibType, ksp::AbstractKSP, flg::PetscBool)
    error("KSPSetLagNorm: no generated method for these argument types")
end

@for_petsc function KSPSetLagNorm(petsclib::$UnionPetscLib, ksp::AbstractKSP, flg::PetscBool )

    @chk ccall(
               (:KSPSetLagNorm, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flg,
              )


	return nothing
end 

"""
	KSPSetMatSolveBatchSize(petsclib::PetscLibType, ksp::AbstractKSP, bs::PetscInt) 
Sets the maximum number of columns treated simultaneously in `KSPMatSolve()`.

Logically Collective

Input Parameters:
- `ksp` - the `KSP` iterative solver
- `bs`  - batch size

Level: advanced

See also: `KSPMatSolve()`, `KSPGetMatSolveBatchSize()`, `-mat_mumps_icntl_27`, `-matproduct_batch_size`

# External Links
$(_doc_external("KSP/KSPSetMatSolveBatchSize"))
"""
function KSPSetMatSolveBatchSize(petsclib::PetscLibType, ksp::AbstractKSP, bs::Integer)
    error("KSPSetMatSolveBatchSize: no generated method for these argument types")
end

@for_petsc function KSPSetMatSolveBatchSize(petsclib::$UnionPetscLib, ksp::AbstractKSP, bs::$PetscInt )

    @chk ccall(
               (:KSPSetMatSolveBatchSize, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, bs,
              )


	return nothing
end 

"""
	KSPSetMinimumIterations(petsclib::PetscLibType, ksp::AbstractKSP, minit::PetscInt) 
Sets the minimum number of iterations to use, regardless of the tolerances

Logically Collective

Input Parameters:
- `ksp`   - the Krylov subspace context
- `minit` - minimum number of iterations to use

Options Database Key:
- `-ksp_min_it minit` - Sets `minit`

Level: intermediate

See also: `KSPGetTolerances()`, `KSPConvergedDefault()`, `KSPSetConvergenceTest()`, `KSP`, `KSPSetTolerances()`, `KSPGetMinimumIterations()`

# External Links
$(_doc_external("KSP/KSPSetMinimumIterations"))
"""
function KSPSetMinimumIterations(petsclib::PetscLibType, ksp::AbstractKSP, minit::Integer)
    error("KSPSetMinimumIterations: no generated method for these argument types")
end

@for_petsc function KSPSetMinimumIterations(petsclib::$UnionPetscLib, ksp::AbstractKSP, minit::$PetscInt )

    @chk ccall(
               (:KSPSetMinimumIterations, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, minit,
              )


	return nothing
end 

"""
	KSPSetNestLevel(petsclib::PetscLibType, ksp::AbstractKSP, level::PetscInt) 
sets the amount of nesting the `KSP` has. That is the number of levels of `KSP` above this `KSP` in a linear solve.

Collective

Input Parameters:
- `ksp`   - the `KSP`
- `level` - the nest level

Level: developer

See also: `KSPSetUp()`, `KSPSolve()`, `KSPDestroy()`, `KSP`, `KSPGMRES`, `KSPType`, `KSPGetNestLevel()`, `PCSetKSPNestLevel()`, `PCGetKSPNestLevel()`

# External Links
$(_doc_external("KSP/KSPSetNestLevel"))
"""
function KSPSetNestLevel(petsclib::PetscLibType, ksp::AbstractKSP, level::Integer)
    error("KSPSetNestLevel: no generated method for these argument types")
end

@for_petsc function KSPSetNestLevel(petsclib::$UnionPetscLib, ksp::AbstractKSP, level::$PetscInt )

    @chk ccall(
               (:KSPSetNestLevel, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, level,
              )


	return nothing
end 

"""
	KSPSetNormType(petsclib::PetscLibType, ksp::AbstractKSP, normtype::KSPNormType) 
Sets the type of residual norm that is used for convergence testing in `KSPSolve()` for the given `KSP` context

Logically Collective

Input Parameters:
- `ksp`      - Krylov solver context
- `normtype` - one of
``
KSP_NORM_NONE             - skips computing the norm, this should generally only be used if you are using
the Krylov method as a smoother with a fixed small number of iterations.
Implicitly sets `KSPConvergedSkip()` as the `KSP` convergence test.
Note that certain algorithms such as `KSPGMRES` ALWAYS require the norm calculation,
for these methods the norms are still computed, they are just not used in
the convergence test.
KSP_NORM_PRECONDITIONED   - the default for left-preconditioned solves, uses the 2-norm
of the preconditioned residual  B^{-1}(b - A x).
KSP_NORM_UNPRECONDITIONED - uses the 2-norm of the true b - Ax residual.
KSP_NORM_NATURAL          - uses the A norm of the true b - Ax residual; supported by `KSPCG`, `KSPCR`, `KSPCGNE`, `KSPCGS`
``

Options Database Key:
- `-ksp_norm_type (none|preconditioned|unpreconditioned|natural)` - set `KSP` norm type

Level: advanced

See also: `KSPSetUp()`, `KSPSolve()`, `KSPDestroy()`, `KSPConvergedSkip()`, `KSPSetCheckNormIteration()`, `KSPSetPCSide()`, `KSPGetPCSide()`, `KSPNormType`

# External Links
$(_doc_external("KSP/KSPSetNormType"))
"""
function KSPSetNormType(petsclib::PetscLibType, ksp::AbstractKSP, normtype::KSPNormType)
    error("KSPSetNormType: no generated method for these argument types")
end

@for_petsc function KSPSetNormType(petsclib::$UnionPetscLib, ksp::AbstractKSP, normtype::KSPNormType )

    @chk ccall(
               (:KSPSetNormType, $petsc_library),
               PetscErrorCode,
               (CKSP, KSPNormType),
               ksp, normtype,
              )


	return nothing
end 

"""
	KSPSetOperators(petsclib::PetscLibType, ksp::AbstractKSP, Amat::AbstractPetscMat, Pmat::AbstractPetscMat) 
Sets the matrix associated with the linear system
and a (possibly) different one from which the preconditioner will be built into the `KSP` context. The matrix will then be used during `KSPSolve()`

Collective

Input Parameters:
- `ksp`  - the `KSP` context
- `Amat` - the matrix that defines the linear system
- `Pmat` - the matrix to be used in constructing the preconditioner, usually the same as `Amat`.

Level: beginner

See also: `KSP`, `Mat`, `KSPSolve()`, `KSPGetPC()`, `PCGetOperators()`, `PCSetOperators()`, `KSPGetOperators()`, `KSPSetComputeOperators()`, `KSPSetComputeInitialGuess()`, `KSPSetComputeRHS()`

# External Links
$(_doc_external("KSP/KSPSetOperators"))
"""
function KSPSetOperators(petsclib::PetscLibType, ksp::AbstractKSP, Amat::AbstractPetscMat, Pmat::AbstractPetscMat)
    error("KSPSetOperators: no generated method for these argument types")
end

@for_petsc function KSPSetOperators(petsclib::$UnionPetscLib, ksp::AbstractKSP, Amat::AbstractPetscMat, Pmat::AbstractPetscMat )

    @chk ccall(
               (:KSPSetOperators, $petsc_library),
               PetscErrorCode,
               (CKSP, CMat, CMat),
               ksp, Amat, Pmat,
              )


	return nothing
end 

"""
	KSPSetOptionsPrefix(petsclib::PetscLibType, ksp::AbstractKSP, prefix::String) 
Sets the prefix used for searching for all
`KSP` options in the database.

Logically Collective

Input Parameters:
- `ksp`    - the Krylov context
- `prefix` - the prefix string to prepend to all `KSP` option requests

Level: intermediate

See also: `KSP`, `KSPAppendOptionsPrefix()`, `KSPGetOptionsPrefix()`, `KSPSetFromOptions()`

# External Links
$(_doc_external("KSP/KSPSetOptionsPrefix"))
"""
function KSPSetOptionsPrefix(petsclib::PetscLibType, ksp::AbstractKSP, prefix::String)
    error("KSPSetOptionsPrefix: no generated method for these argument types")
end

@for_petsc function KSPSetOptionsPrefix(petsclib::$UnionPetscLib, ksp::AbstractKSP, prefix::String )

    @chk ccall(
               (:KSPSetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Cchar}),
               ksp, prefix,
              )


	return nothing
end 

"""
	KSPSetPC(petsclib::PetscLibType, ksp::AbstractKSP, pc::AbstractPC) 
Sets the preconditioner to be used to calculate the
application of the preconditioner on a vector into a `KSP`.

Collective

Input Parameters:
- `ksp` - the `KSP` iterative solver obtained from `KSPCreate()`
- `pc`  - the preconditioner object (if `NULL` it returns the `PC` currently held by the `KSP`)

Level: developer

See also: `KSPGetPC()`, `KSP`

# External Links
$(_doc_external("KSP/KSPSetPC"))
"""
function KSPSetPC(petsclib::PetscLibType, ksp::AbstractKSP, pc::AbstractPC)
    error("KSPSetPC: no generated method for these argument types")
end

@for_petsc function KSPSetPC(petsclib::$UnionPetscLib, ksp::AbstractKSP, pc::AbstractPC )

    @chk ccall(
               (:KSPSetPC, $petsc_library),
               PetscErrorCode,
               (CKSP, CPC),
               ksp, pc,
              )


	return nothing
end 

"""
	KSPSetPCSide(petsclib::PetscLibType, ksp::AbstractKSP, side::PCSide) 
Sets the preconditioning side.

Logically Collective

Input Parameter:
- `ksp` - iterative solver obtained from `KSPCreate()`

Output Parameter:
- `side` - the preconditioning side, where side is one of
``
PC_LEFT      - left preconditioning (default)
PC_RIGHT     - right preconditioning
PC_SYMMETRIC - symmetric preconditioning
``

Options Database Key:
- `-ksp_pc_side (right|left|symmetric)` - `KSP` preconditioner side

Level: intermediate

See also: `KSPGetPCSide()`, `KSPSetNormType()`, `KSPGetNormType()`, `KSP`, `KSPSetPreSolve()`, `KSPSetPostSolve()`

# External Links
$(_doc_external("KSP/KSPSetPCSide"))
"""
function KSPSetPCSide(petsclib::PetscLibType, ksp::AbstractKSP, side::PCSide)
    error("KSPSetPCSide: no generated method for these argument types")
end

@for_petsc function KSPSetPCSide(petsclib::$UnionPetscLib, ksp::AbstractKSP, side::PCSide )

    @chk ccall(
               (:KSPSetPCSide, $petsc_library),
               PetscErrorCode,
               (CKSP, PCSide),
               ksp, side,
              )


	return nothing
end 

"""
	KSPSetPostSolve(petsclib::PetscLibType, ksp::AbstractKSP, postsolve::Ptr{Cvoid}, ctx::Ptr{Cvoid}) 
Sets a function that is called at the end of each `KSPSolve()` (whether it converges or not). Used in conjunction with `KSPSetPreSolve()`.

Logically Collective

Input Parameters:
- `ksp`       - the solver object
- `postsolve` - the function to call after the solve, see` KSPPSolveFn`
- `ctx`       - an optional context needed by the function

Level: developer

See also: `KSPPSolveFn`, `KSPSetUp()`, `KSPSolve()`, `KSPDestroy()`, `KSP`, `KSPSetPreSolve()`, `PCEISENSTAT`

# External Links
$(_doc_external("KSP/KSPSetPostSolve"))
"""
function KSPSetPostSolve(petsclib::PetscLibType, ksp::AbstractKSP, postsolve::Ptr{Cvoid}, ctx::Ptr{Cvoid})
    error("KSPSetPostSolve: no generated method for these argument types")
end

@for_petsc function KSPSetPostSolve(petsclib::$UnionPetscLib, ksp::AbstractKSP, postsolve::Ptr{Cvoid}, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:KSPSetPostSolve, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Cvoid}, Ptr{Cvoid}),
               ksp, postsolve, ctx,
              )


	return nothing
end 

"""
	KSPSetPreSolve(petsclib::PetscLibType, ksp::AbstractKSP, presolve::Ptr{Cvoid}, ctx::Ptr{Cvoid}) 
Sets a function that is called at the beginning of each `KSPSolve()`. Used in conjunction with `KSPSetPostSolve()`.

Logically Collective

Input Parameters:
- `ksp`      - the solver object
- `presolve` - the function to call before the solve, see` KSPPSolveFn`
- `ctx`      - an optional context needed by the function

Level: developer

See also: `KSPPSolveFn`, `KSPSetUp()`, `KSPSolve()`, `KSPDestroy()`, `KSP`, `KSPSetPostSolve()`, `PCEISENSTAT`, `PCPreSolve()`, `PCPostSolve()`

# External Links
$(_doc_external("KSP/KSPSetPreSolve"))
"""
function KSPSetPreSolve(petsclib::PetscLibType, ksp::AbstractKSP, presolve::Ptr{Cvoid}, ctx::Ptr{Cvoid})
    error("KSPSetPreSolve: no generated method for these argument types")
end

@for_petsc function KSPSetPreSolve(petsclib::$UnionPetscLib, ksp::AbstractKSP, presolve::Ptr{Cvoid}, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:KSPSetPreSolve, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{Cvoid}, Ptr{Cvoid}),
               ksp, presolve, ctx,
              )


	return nothing
end 

"""
	KSPSetResidualHistory(petsclib::PetscLibType, ksp::AbstractKSP, a::Vector{PetscReal}, na::PetscCount, reset::PetscBool) 
Sets the array used to hold the residual history.
If set, this array will contain the residual norms computed at each
iteration of the solver.

Not Collective

Input Parameters:
- `ksp`   - iterative solver obtained from `KSPCreate()`
- `a`     - array to hold history
- `na`    - size of `a`
- `reset` - `PETSC_TRUE` indicates the history counter is reset to zero
for each new linear solve

Level: advanced

See also: `KSPGetResidualHistory()`, `KSP`

# External Links
$(_doc_external("KSP/KSPSetResidualHistory"))
"""
function KSPSetResidualHistory(petsclib::PetscLibType, ksp::AbstractKSP, a::AbstractVector{<:Number}, na::PetscCount, reset::PetscBool)
    error("KSPSetResidualHistory: no generated method for these argument types")
end

@for_petsc function KSPSetResidualHistory(petsclib::$UnionPetscLib, ksp::AbstractKSP, a::Vector{$PetscReal}, na::PetscCount, reset::PetscBool )

    @chk ccall(
               (:KSPSetResidualHistory, $petsc_library),
               PetscErrorCode,
               (CKSP, Ptr{$PetscReal}, PetscCount, PetscBool),
               ksp, a, na, reset,
              )


	return nothing
end 

"""
	KSPSetReusePreconditioner(petsclib::PetscLibType, ksp::AbstractKSP, flag::PetscBool) 
reuse the current preconditioner for future `KSPSolve()`, do not construct a new preconditioner even if the `Mat` operator
in the `KSP` has different values

Collective

Input Parameters:
- `ksp`  - iterative solver obtained from `KSPCreate()`
- `flag` - `PETSC_TRUE` to reuse the current preconditioner, or `PETSC_FALSE` to construct a new preconditioner

Options Database Key:
- `-ksp_reuse_preconditioner (true|false)` - reuse the previously computed preconditioner

Level: intermediate

See also: `KSPCreate()`, `KSPSolve()`, `KSPDestroy()`, `KSP`, `KSPGetReusePreconditioner()`,
`SNESSetLagPreconditioner()`, `SNES`

# External Links
$(_doc_external("KSP/KSPSetReusePreconditioner"))
"""
function KSPSetReusePreconditioner(petsclib::PetscLibType, ksp::AbstractKSP, flag::PetscBool)
    error("KSPSetReusePreconditioner: no generated method for these argument types")
end

@for_petsc function KSPSetReusePreconditioner(petsclib::$UnionPetscLib, ksp::AbstractKSP, flag::PetscBool )

    @chk ccall(
               (:KSPSetReusePreconditioner, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flag,
              )


	return nothing
end 

"""
	KSPSetSkipPCSetFromOptions(petsclib::PetscLibType, ksp::AbstractKSP, flag::PetscBool) 
prevents `KSPSetFromOptions()` from calling `PCSetFromOptions()`.
This is used if the same `PC` is shared by more than one `KSP` so its options are not reset for each `KSP`

Collective

Input Parameters:
- `ksp`  - iterative solver obtained from `KSPCreate()`
- `flag` - `PETSC_TRUE` to skip calling the `PCSetFromOptions()`

Level: developer

See also: `KSPCreate()`, `KSPSolve()`, `KSPDestroy()`, `PCSetReusePreconditioner()`, `KSP`

# External Links
$(_doc_external("KSP/KSPSetSkipPCSetFromOptions"))
"""
function KSPSetSkipPCSetFromOptions(petsclib::PetscLibType, ksp::AbstractKSP, flag::PetscBool)
    error("KSPSetSkipPCSetFromOptions: no generated method for these argument types")
end

@for_petsc function KSPSetSkipPCSetFromOptions(petsclib::$UnionPetscLib, ksp::AbstractKSP, flag::PetscBool )

    @chk ccall(
               (:KSPSetSkipPCSetFromOptions, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flag,
              )


	return nothing
end 

"""
	KSPSetSupportedNorm(petsclib::PetscLibType, ksp::AbstractKSP, normtype::KSPNormType, pcside::PCSide, priority::PetscInt) 
Sets a norm and preconditioner side supported by a `KSPType`

Logically Collective

Input Parameters:
- `ksp`      - Krylov method
- `normtype` - supported norm type of the type `KSPNormType`
- `pcside`   - preconditioner side, of the type `PCSide` that can be used with this `KSPNormType`
- `priority` - positive integer preference for this combination; larger values have higher priority

Level: developer

See also: `KSP`, `KSPNormType`, `PCSide`, `KSPSetNormType()`, `KSPSetPCSide()`

# External Links
$(_doc_external("KSP/KSPSetSupportedNorm"))
"""
function KSPSetSupportedNorm(petsclib::PetscLibType, ksp::AbstractKSP, normtype::KSPNormType, pcside::PCSide, priority::Integer)
    error("KSPSetSupportedNorm: no generated method for these argument types")
end

@for_petsc function KSPSetSupportedNorm(petsclib::$UnionPetscLib, ksp::AbstractKSP, normtype::KSPNormType, pcside::PCSide, priority::$PetscInt )

    @chk ccall(
               (:KSPSetSupportedNorm, $petsc_library),
               PetscErrorCode,
               (CKSP, KSPNormType, PCSide, $PetscInt),
               ksp, normtype, pcside, priority,
              )


	return nothing
end 

"""
	KSPSetTolerances(petsclib::PetscLibType, ksp::AbstractKSP, rtol::PetscReal, abstol::PetscReal, dtol::PetscReal, maxits::PetscInt) 
Sets the relative, absolute, divergence, and maximum
iteration tolerances used by the default `KSP` convergence testers.

Logically Collective

Input Parameters:
- `ksp`    - the Krylov subspace context
- `rtol`   - the relative convergence tolerance, relative decrease in the (possibly preconditioned) residual norm
- `abstol` - the absolute convergence tolerance   absolute size of the (possibly preconditioned) residual norm
- `dtol`   - the divergence tolerance,   amount (possibly preconditioned) residual norm can increase before `KSPConvergedDefault()` concludes that the method is diverging
- `maxits` - maximum number of iterations to use

Options Database Keys:
- `-ksp_atol abstol`   - Sets `abstol`
- `-ksp_rtol rtol`     - Sets `rtol`
- `-ksp_divtol dtol`   - Sets `dtol`
- `-ksp_max_it maxits` - Sets `maxits`

Level: intermediate

See also: `KSPGetTolerances()`, `KSPConvergedDefault()`, `KSPSetConvergenceTest()`, `KSP`, `KSPSetMinimumIterations()`

# External Links
$(_doc_external("KSP/KSPSetTolerances"))
"""
function KSPSetTolerances(petsclib::PetscLibType, ksp::AbstractKSP, rtol::Real, abstol::Real, dtol::Real, maxits::Integer)
    error("KSPSetTolerances: no generated method for these argument types")
end

@for_petsc function KSPSetTolerances(petsclib::$UnionPetscLib, ksp::AbstractKSP, rtol::$PetscReal, abstol::$PetscReal, dtol::$PetscReal, maxits::$PetscInt )

    @chk ccall(
               (:KSPSetTolerances, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscReal, $PetscReal, $PetscReal, $PetscInt),
               ksp, rtol, abstol, dtol, maxits,
              )


	return nothing
end 

"""
	KSPSetType(petsclib::PetscLibType, ksp::AbstractKSP, type::String) 
Sets the algorithm/method to be used to solve the linear system with the given `KSP`

Logically Collective

Input Parameters:
- `ksp`  - the Krylov space context
- `type` - a known method

Options Database Key:
- `-ksp_type type` - Sets the method; see `KSPType`

Level: intermediate

See also: `PCSetType()`, `KSPType`, `KSPRegister()`, `KSPCreate()`, `KSP`

# External Links
$(_doc_external("KSP/KSPSetType"))
"""
function KSPSetType(petsclib::PetscLibType, ksp::AbstractKSP, type::String)
    error("KSPSetType: no generated method for these argument types")
end

@for_petsc function KSPSetType(petsclib::$UnionPetscLib, ksp::AbstractKSP, type::String )

    @chk ccall(
               (:KSPSetType, $petsc_library),
               PetscErrorCode,
               (CKSP, KSPType),
               ksp, type,
              )


	return nothing
end 

"""
	KSPSetUp(petsclib::PetscLibType, ksp::AbstractKSP) 
Sets up the internal data structures for the
later use `KSPSolve()` the `KSP` linear iterative solver.

Collective

Input Parameter:
- `ksp` - iterative solver, `KSP`, obtained from `KSPCreate()`

Level: developer

See also: `KSPCreate()`, `KSPSolve()`, `KSPDestroy()`, `KSP`, `KSPSetUpOnBlocks()`

# External Links
$(_doc_external("KSP/KSPSetUp"))
"""
function KSPSetUp(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPSetUp: no generated method for these argument types")
end

@for_petsc function KSPSetUp(petsclib::$UnionPetscLib, ksp::AbstractKSP )

    @chk ccall(
               (:KSPSetUp, $petsc_library),
               PetscErrorCode,
               (CKSP,),
               ksp,
              )


	return nothing
end 

"""
	KSPSetUpOnBlocks(petsclib::PetscLibType, ksp::AbstractKSP) 
Sets up the preconditioner for each block in
the block Jacobi `PCJACOBI`, overlapping Schwarz `PCASM`, and fieldsplit `PCFIELDSPLIT` preconditioners

Collective

Input Parameter:
- `ksp` - the `KSP` context

Level: advanced

See also: `PCSetUpOnBlocks()`, `KSPSetUp()`, `PCSetUp()`, `KSP`

# External Links
$(_doc_external("KSP/KSPSetUpOnBlocks"))
"""
function KSPSetUpOnBlocks(petsclib::PetscLibType, ksp::AbstractKSP)
    error("KSPSetUpOnBlocks: no generated method for these argument types")
end

@for_petsc function KSPSetUpOnBlocks(petsclib::$UnionPetscLib, ksp::AbstractKSP )

    @chk ccall(
               (:KSPSetUpOnBlocks, $petsc_library),
               PetscErrorCode,
               (CKSP,),
               ksp,
              )


	return nothing
end 

"""
	KSPSetUseExplicitTranspose(petsclib::PetscLibType, ksp::AbstractKSP, flg::PetscBool) 
Determines the explicit transpose of the operator is formed in `KSPSolveTranspose()`. In some configurations (like GPUs) it may
be explicitly formed since the solve is much more efficient.

Logically Collective

Input Parameter:
- `ksp` - the `KSP` context

Output Parameter:
- `flg` - `PETSC_TRUE` to transpose the system in `KSPSolveTranspose()`, `PETSC_FALSE` to not transpose (default)

Level: advanced

See also: `KSPSolveTranspose()`, `KSP`

# External Links
$(_doc_external("KSP/KSPSetUseExplicitTranspose"))
"""
function KSPSetUseExplicitTranspose(petsclib::PetscLibType, ksp::AbstractKSP, flg::PetscBool)
    error("KSPSetUseExplicitTranspose: no generated method for these argument types")
end

@for_petsc function KSPSetUseExplicitTranspose(petsclib::$UnionPetscLib, ksp::AbstractKSP, flg::PetscBool )

    @chk ccall(
               (:KSPSetUseExplicitTranspose, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscBool),
               ksp, flg,
              )


	return nothing
end 

"""
	KSPSetUseFischerGuess(petsclib::PetscLibType, ksp::AbstractKSP, model::PetscInt, size::PetscInt) 
Use the Paul Fischer algorithm or its variants to compute initial guesses for a set of solves with related right-hand sides

Logically Collective

Input Parameters:
- `ksp`   - the Krylov context
- `model` - use model 1, model 2, model 3, or any other number to turn it off
- `size`  - size of subspace used to generate initial guess

Options Database Key:
- `-ksp_fischer_guess model,size` - uses the Fischer initial guess generator for repeated linear solves

Level: advanced

See also: `KSP`, `KSPSetOptionsPrefix()`, `KSPAppendOptionsPrefix()`, `KSPSetGuess()`, `KSPGetGuess()`, `KSPGuess`

# External Links
$(_doc_external("KSP/KSPSetUseFischerGuess"))
"""
function KSPSetUseFischerGuess(petsclib::PetscLibType, ksp::AbstractKSP, model::Integer, size::Integer)
    error("KSPSetUseFischerGuess: no generated method for these argument types")
end

@for_petsc function KSPSetUseFischerGuess(petsclib::$UnionPetscLib, ksp::AbstractKSP, model::$PetscInt, size::$PetscInt )

    @chk ccall(
               (:KSPSetUseFischerGuess, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt, $PetscInt),
               ksp, model, size,
              )


	return nothing
end 

"""
	KSPSetWorkVecs(petsclib::PetscLibType, ksp::AbstractKSP, nw::PetscInt) 
Sets a number of work vectors into a `KSP` object

Collective

Input Parameters:
- `ksp` - iterative context
- `nw`  - number of work vectors to allocate

Level: developer

See also: `KSP`, `KSPCreateVecs()`

# External Links
$(_doc_external("KSP/KSPSetWorkVecs"))
"""
function KSPSetWorkVecs(petsclib::PetscLibType, ksp::AbstractKSP, nw::Integer)
    error("KSPSetWorkVecs: no generated method for these argument types")
end

@for_petsc function KSPSetWorkVecs(petsclib::$UnionPetscLib, ksp::AbstractKSP, nw::$PetscInt )

    @chk ccall(
               (:KSPSetWorkVecs, $petsc_library),
               PetscErrorCode,
               (CKSP, $PetscInt),
               ksp, nw,
              )


	return nothing
end 

"""
	KSPSolve(petsclib::PetscLibType, ksp::AbstractKSP, b::Union{Ptr, AbstractPetscVec}, x::Union{Ptr, AbstractPetscVec}) 
Solves a linear system associated with `KSP` object

Collective

Input Parameters:
- `ksp` - iterative solver obtained from `KSPCreate()`
- `b`   - the right-hand side vector
- `x`   - the solution (this may be the same vector as `b`, then `b` will be overwritten with the answer)

Options Database Keys:
- `-ksp_view_eigenvalues`                      - compute preconditioned operators eigenvalues
- `-ksp_view_eigenvalues_explicit`             - compute the eigenvalues by forming the dense operator and using LAPACK
- `-ksp_view_mat binary`                       - save matrix to the default binary viewer
- `-ksp_view_pmat binary`                      - save matrix used to build preconditioner to the default binary viewer
- `-ksp_view_rhs binary`                       - save right-hand side vector to the default binary viewer
- `-ksp_view_solution binary`                  - save computed solution vector to the default binary viewer
(can be read later with src/ksp/tutorials/ex10.c for testing solvers)
- `-ksp_view_mat_explicit`                     - for matrix-free operators, computes the matrix entries and views them
- `-ksp_view_preconditioned_operator_explicit` - computes the product of the preconditioner and matrix as an explicit matrix and views it
- `-ksp_converged_reason`                      - print reason for converged or diverged, also prints number of iterations
- `-ksp_view_final_residual`                   - print 2-norm of true linear system residual at the end of the solution process
- `-ksp_view_final_residual_vec`               - print true linear system residual vector at the end of the solution process;
`-ksp_view_final_residual` must to be called first to enable this option
- `-ksp_error_if_not_converged`                - stop the program as soon as an error is detected in a `KSPSolve()`
- `-ksp_view_pre`                              - print the ksp data structure before the system solution
- `-ksp_view`                                  - print the ksp data structure at the end of the system solution

Level: beginner

See also: `KSPCreate()`, `KSPSetUp()`, `KSPDestroy()`, `KSPSetTolerances()`, `KSPConvergedDefault()`,
`KSPSolveTranspose()`, `KSPGetIterationNumber()`, `MatNullSpaceCreate()`, `MatSetNullSpace()`, `MatSetTransposeNullSpace()`, `KSP`,
`KSPConvergedReasonView()`, `KSPCheckSolve()`, `KSPSetErrorIfNotConverged()`

# External Links
$(_doc_external("KSP/KSPSolve"))
"""
function KSPSolve(petsclib::PetscLibType, ksp::AbstractKSP, b::Union{Ptr, AbstractPetscVec}, x::Union{Ptr, AbstractPetscVec})
    error("KSPSolve: no generated method for these argument types")
end

@for_petsc function KSPSolve(petsclib::$UnionPetscLib, ksp::AbstractKSP, b::Union{Ptr, AbstractPetscVec}, x::Union{Ptr, AbstractPetscVec} )

    @chk ccall(
               (:KSPSolve, $petsc_library),
               PetscErrorCode,
               (CKSP, CVec, CVec),
               ksp, b, x,
              )


	return nothing
end 

"""
	KSPSolveTranspose(petsclib::PetscLibType, ksp::AbstractKSP, b::AbstractPetscVec, x::AbstractPetscVec) 
Solves a linear system with the transpose of the matrix associated with the `KSP` object, A^T x = b.

Collective

Input Parameters:
- `ksp` - iterative solver obtained from `KSPCreate()`
- `b`   - right-hand side vector
- `x`   - solution vector

Level: developer

See also: `KSPCreate()`, `KSPSetUp()`, `KSPDestroy()`, `KSPSetTolerances()`, `KSPConvergedDefault()`,
`KSPSolve()`, `KSP`, `KSPSetOperators()`

# External Links
$(_doc_external("KSP/KSPSolveTranspose"))
"""
function KSPSolveTranspose(petsclib::PetscLibType, ksp::AbstractKSP, b::AbstractPetscVec, x::AbstractPetscVec)
    error("KSPSolveTranspose: no generated method for these argument types")
end

@for_petsc function KSPSolveTranspose(petsclib::$UnionPetscLib, ksp::AbstractKSP, b::AbstractPetscVec, x::AbstractPetscVec )

    @chk ccall(
               (:KSPSolveTranspose, $petsc_library),
               PetscErrorCode,
               (CKSP, CVec, CVec),
               ksp, b, x,
              )


	return nothing
end 

"""
	KSPUnwindPreconditioner(petsclib::PetscLibType, ksp::AbstractKSP, vsoln::AbstractPetscVec, vt1::AbstractPetscVec) 
Unwinds the preconditioning in the solution. That is,
takes solution to the preconditioned problem and gets the solution to the
original problem from it.

Collective

Input Parameters:
- `ksp`   - iterative context
- `vsoln` - solution vector
- `vt1`   - temporary work vector

Output Parameter:
- `vsoln` - contains solution on output

Level: advanced

See also: `KSP`, `KSPSetPCSide()`

# External Links
$(_doc_external("KSP/KSPUnwindPreconditioner"))
"""
function KSPUnwindPreconditioner(petsclib::PetscLibType, ksp::AbstractKSP, vsoln::AbstractPetscVec, vt1::AbstractPetscVec)
    error("KSPUnwindPreconditioner: no generated method for these argument types")
end

@for_petsc function KSPUnwindPreconditioner(petsclib::$UnionPetscLib, ksp::AbstractKSP, vsoln::AbstractPetscVec, vt1::AbstractPetscVec )

    @chk ccall(
               (:KSPUnwindPreconditioner, $petsc_library),
               PetscErrorCode,
               (CKSP, CVec, CVec),
               ksp, vsoln, vt1,
              )


	return nothing
end 

"""
	KSPView(petsclib::PetscLibType, ksp::AbstractKSP, viewer::PetscViewer) 
Prints the various parameters currently set in the `KSP` object. For example, the convergence tolerances and `KSPType`.
Also views the `PC` and `Mat` contained by the `KSP` with `PCView()` and `MatView()`.

Collective

Input Parameters:
- `ksp`    - the Krylov space context
- `viewer` - visualization context

Options Database Key:
- `-ksp_view` - print the `KSP` data structure at the end of each `KSPSolve()` call

Level: beginner

See also: `KSP`, `PetscViewer`, `PCView()`, `PetscViewerASCIIOpen()`, `KSPViewFromOptions()`

# External Links
$(_doc_external("KSP/KSPView"))
"""
function KSPView(petsclib::PetscLibType, ksp::AbstractKSP, viewer::PetscViewer)
    error("KSPView: no generated method for these argument types")
end

@for_petsc function KSPView(petsclib::$UnionPetscLib, ksp::AbstractKSP, viewer::PetscViewer )

    @chk ccall(
               (:KSPView, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscViewer),
               ksp, viewer,
              )


	return nothing
end 

"""
	KSPViewFromOptions(petsclib::PetscLibType, A::AbstractKSP, obj, name::String) 
View (print) a `KSP` object based on values in the options database. Also views the `PC` and `Mat` contained by the `KSP`
with `PCView()` and `MatView()`.

Collective

Input Parameters:
- `A`    - Krylov solver context
- `obj`  - Optional object that provides the options prefix used to query the options database
- `name` - command line option

Options Database Key:
- `-name [viewertype][:...]` - option name and values. See `PetscObjectViewFromOptions()` for the possible arguments

Level: intermediate

See also: `KSP`, `KSPView()`, `PetscObjectViewFromOptions()`, `KSPCreate()`

# External Links
$(_doc_external("KSP/KSPViewFromOptions"))
"""
function KSPViewFromOptions(petsclib::PetscLibType, A::AbstractKSP, obj, name::String)
    error("KSPViewFromOptions: no generated method for these argument types")
end

@for_petsc function KSPViewFromOptions(petsclib::$UnionPetscLib, A::AbstractKSP, obj, name::String )

    @chk ccall(
               (:KSPViewFromOptions, $petsc_library),
               PetscErrorCode,
               (CKSP, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

