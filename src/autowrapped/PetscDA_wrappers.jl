"""
	PetscDAAppendOptionsPrefix(petsclib::PetscLibType, das::PetscDA, p::String) 
Appends to the prefix used for searching for all PetscDA options in the database.

Logically Collective

Input Parameters:
- `das` - the `PetscDA` context
- `p`   - the prefix string to prepend to all `PetscDA` option requests

Level: advanced

See also: `PetscDA`, `PetscDASetFromOptions()`, `PetscDASetOptionsPrefix()`, `PetscDAGetOptionsPrefix()`

# External Links
$(_doc_external("PetscDA/PetscDAAppendOptionsPrefix"))
"""
function PetscDAAppendOptionsPrefix(petsclib::PetscLibType, das::PetscDA, p::String)
    error("PetscDAAppendOptionsPrefix: no generated method for these argument types")
end

@for_petsc function PetscDAAppendOptionsPrefix(petsclib::$UnionPetscLib, das::PetscDA, p::String )

    @chk ccall(
               (:PetscDAAppendOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (PetscDA, Ptr{Cchar}),
               das, p,
              )


	return nothing
end 

"""
	da_out::PetscDA = PetscDACreate(petsclib::PetscLibType, comm::MPI_Comm) 
Creates a new `PetscDA` object for data assimilation.

Collective

Input Parameter:
- `comm` - MPI communicator used to create the object

Output Parameter:
- `da_out` - newly created `PetscDA` object

Level: beginner

See also: `PetscDADestroy()`, `PetscDASetType()`, `PetscDASetUp()`

# External Links
$(_doc_external("PetscDA/PetscDACreate"))
"""
function PetscDACreate(petsclib::PetscLibType, comm::MPI_Comm)
    error("PetscDACreate: no generated method for these argument types")
end

@for_petsc function PetscDACreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	da_out_ = Ref{PetscDA}()

    @chk ccall(
               (:PetscDACreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscDA}),
               comm, da_out_,
              )

	da_out = da_out_[]

	return da_out
end 

"""
	PetscDADestroy(petsclib::PetscLibType, da::Union{PetscDA, Ref{PetscDA}}) 
Destroys a `PetscDA` object and releases its resources.

Collective

Input Parameter:
- `da` - pointer to the `PetscDA` object to destroy

Level: beginner

See also: `PetscDACreate()`

# External Links
$(_doc_external("PetscDA/PetscDADestroy"))
"""
function PetscDADestroy(petsclib::PetscLibType, da::Union{PetscDA, Ref{PetscDA}})
    error("PetscDADestroy: no generated method for these argument types")
end

@for_petsc function PetscDADestroy(petsclib::$UnionPetscLib, da::Union{PetscDA, Ref{PetscDA}} )
	da_ = da isa Base.RefValue ? da : Ref{PetscDA}(da)

    @chk ccall(
               (:PetscDADestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscDA},),
               da_,
              )


	return nothing
end 

"""
	PetscDAEnsembleAnalysis(petsclib::PetscLibType, da::PetscDA, observation::AbstractPetscVec, H::AbstractPetscMat) 
Executes the analysis (update) step using sparse observation matrix H

Collective

Input Parameters:
- `da`          - the `PetscDA` context
- `observation` - observation vector y in R^P
- `H`           - observation operator matrix (P x N), sparse AIJ format

See also: `PetscDA`, `PETSCDAETKF`, `PETSCDALETKF`, `PetscDAEnsembleForecast()`, `PetscDASetObsErrorVariance()`

# External Links
$(_doc_external("PetscDA/PetscDAEnsembleAnalysis"))
"""
function PetscDAEnsembleAnalysis(petsclib::PetscLibType, da::PetscDA, observation::AbstractPetscVec, H::AbstractPetscMat)
    error("PetscDAEnsembleAnalysis: no generated method for these argument types")
end

@for_petsc function PetscDAEnsembleAnalysis(petsclib::$UnionPetscLib, da::PetscDA, observation::AbstractPetscVec, H::AbstractPetscMat )

    @chk ccall(
               (:PetscDAEnsembleAnalysis, $petsc_library),
               PetscErrorCode,
               (PetscDA, CVec, CMat),
               da, observation, H,
              )


	return nothing
end 

"""
	PetscDAEnsembleApplySqrtTInverse(petsclib::PetscLibType, da::PetscDA, U::AbstractPetscMat, Y::AbstractPetscMat) 
Apply T^{-1/2} to a matrix U [Alg 6.4 line 9]

Collective

Input Parameters:
- `da` - the `PetscDA` context
- `U`  - input matrix (usually Identity, but can be general)

Output Parameter:
- `Y` - output matrix Y = T^{-1/2} * U

See also: `PetscDA`, `PETSCDAETKF`, `PETSCDALETKF`, `PetscDAEnsembleTFactor()`, `PetscDAEnsembleApplyTInverse()`

# External Links
$(_doc_external("PetscDA/PetscDAEnsembleApplySqrtTInverse"))
"""
function PetscDAEnsembleApplySqrtTInverse(petsclib::PetscLibType, da::PetscDA, U::AbstractPetscMat, Y::AbstractPetscMat)
    error("PetscDAEnsembleApplySqrtTInverse: no generated method for these argument types")
end

@for_petsc function PetscDAEnsembleApplySqrtTInverse(petsclib::$UnionPetscLib, da::PetscDA, U::AbstractPetscMat, Y::AbstractPetscMat )

    @chk ccall(
               (:PetscDAEnsembleApplySqrtTInverse, $petsc_library),
               PetscErrorCode,
               (PetscDA, CMat, CMat),
               da, U, Y,
              )


	return nothing
end 

"""
	PetscDAEnsembleApplyTInverse(petsclib::PetscLibType, da::PetscDA, sdel::AbstractPetscVec, w::AbstractPetscVec) 
Apply T^{-1} to a vector [Alg 6.4 line 8]

Collective

Input Parameters:
- `da`   - the `PetscDA` context
- `sdel` - input vector S^T-delta

Output Parameter:
- `w` - output vector w = T^{-1} * sdel

See also: `PetscDA`, `PETSCDAETKF`, `PETSCDALETKF`, `PetscDAEnsembleTFactor()`, `PetscDAEnsembleApplySqrtTInverse()`

# External Links
$(_doc_external("PetscDA/PetscDAEnsembleApplyTInverse"))
"""
function PetscDAEnsembleApplyTInverse(petsclib::PetscLibType, da::PetscDA, sdel::AbstractPetscVec, w::AbstractPetscVec)
    error("PetscDAEnsembleApplyTInverse: no generated method for these argument types")
end

@for_petsc function PetscDAEnsembleApplyTInverse(petsclib::$UnionPetscLib, da::PetscDA, sdel::AbstractPetscVec, w::AbstractPetscVec )

    @chk ccall(
               (:PetscDAEnsembleApplyTInverse, $petsc_library),
               PetscErrorCode,
               (PetscDA, CVec, CVec),
               da, sdel, w,
              )


	return nothing
end 

"""
	anomalies_out::PetscMat = PetscDAEnsembleComputeAnomalies(petsclib::PetscLibType, da::PetscDA, mean_in::AbstractPetscVec) 
Forms the state-space anomalies matrix for a `PetscDA`.

Collective

Input Parameters:
- `da`      - the `PetscDA` context
- `mean_in` - optional mean state vector (pass `NULL` to compute internally)

Output Parameter:
- `anomalies_out` - location to store the newly created anomalies matrix

See also: `PetscDA`, `PETSCDAETKF`, `PETSCDALETKF`, `PetscDAEnsembleComputeMean()`

# External Links
$(_doc_external("PetscDA/PetscDAEnsembleComputeAnomalies"))
"""
function PetscDAEnsembleComputeAnomalies(petsclib::PetscLibType, da::PetscDA, mean_in::AbstractPetscVec)
    error("PetscDAEnsembleComputeAnomalies: no generated method for these argument types")
end

@for_petsc function PetscDAEnsembleComputeAnomalies(petsclib::$UnionPetscLib, da::PetscDA, mean_in::AbstractPetscVec )
	anomalies_out_ = Ref{CMat}()

    @chk ccall(
               (:PetscDAEnsembleComputeAnomalies, $petsc_library),
               PetscErrorCode,
               (PetscDA, CVec, Ptr{CMat}),
               da, mean_in, anomalies_out_,
              )

	anomalies_out = PetscMat(anomalies_out_[], petsclib)

	return anomalies_out
end 

"""
	PetscDAEnsembleComputeMean(petsclib::PetscLibType, da::PetscDA, mean::AbstractPetscVec) 
Computes ensemble mean for a `PetscDA`

Collective

Input Parameter:
- `da` - the `PetscDA` context

Output Parameter:
- `mean` - vector that will hold the ensemble mean

Level: intermediate

See also: `PetscDA`, `PETSCDAETKF`, `PETSCDALETKF`, `PetscDAEnsembleComputeAnomalies()`

# External Links
$(_doc_external("PetscDA/PetscDAEnsembleComputeMean"))
"""
function PetscDAEnsembleComputeMean(petsclib::PetscLibType, da::PetscDA, mean::AbstractPetscVec)
    error("PetscDAEnsembleComputeMean: no generated method for these argument types")
end

@for_petsc function PetscDAEnsembleComputeMean(petsclib::$UnionPetscLib, da::PetscDA, mean::AbstractPetscVec )

    @chk ccall(
               (:PetscDAEnsembleComputeMean, $petsc_library),
               PetscErrorCode,
               (PetscDA, CVec),
               da, mean,
              )


	return nothing
end 

"""
	PetscDAEnsembleComputeNormalizedInnovationMatrix(petsclib::PetscLibType, Z::AbstractPetscMat, y_mean::AbstractPetscVec, r_inv_sqrt::AbstractPetscVec, m::PetscInt, scale::PetscScalar, S::AbstractPetscMat) 
Computes S = R^{-1/2}(Z - y_mean * 1')/sqrt(m-1) [Alg 6.4 line 5]

Collective

Input Parameters:
- `Z`          - observation ensemble matrix
- `y_mean`     - mean of observations
- `r_inv_sqrt` - R^{-1/2}
- `m`          - ensemble size
- `scale`      - 1/sqrt(m-1)

Output Parameter:
- `S` - normalized innovation matrix

Level: developer

See also: `PetscDA`, `PETSCDAETKF`, `PETSCDALETKF`, `PetscDASetSizes()`, `PetscDAGetSizes()`

# External Links
$(_doc_external("PetscDA/PetscDAEnsembleComputeNormalizedInnovationMatrix"))
"""
function PetscDAEnsembleComputeNormalizedInnovationMatrix(petsclib::PetscLibType, Z::AbstractPetscMat, y_mean::AbstractPetscVec, r_inv_sqrt::AbstractPetscVec, m::Integer, scale::Number, S::AbstractPetscMat)
    error("PetscDAEnsembleComputeNormalizedInnovationMatrix: no generated method for these argument types")
end

@for_petsc function PetscDAEnsembleComputeNormalizedInnovationMatrix(petsclib::$UnionPetscLib, Z::AbstractPetscMat, y_mean::AbstractPetscVec, r_inv_sqrt::AbstractPetscVec, m::$PetscInt, scale::$PetscScalar, S::AbstractPetscMat )

    @chk ccall(
               (:PetscDAEnsembleComputeNormalizedInnovationMatrix, $petsc_library),
               PetscErrorCode,
               (CMat, CVec, CVec, $PetscInt, $PetscScalar, CMat),
               Z, y_mean, r_inv_sqrt, m, scale, S,
              )


	return nothing
end 

"""
	PetscDAEnsembleForecast(petsclib::PetscLibType, da::PetscDA, model::external, ctx::Ptr{Cvoid}) 
Advances every ensemble member through the user-supplied forecast model.

Collective

Input Parameters:
- `da`    - the `PetscDA` context
- `model` - routine that evaluates the model map `f(input, output; ctx)`
- `ctx`   - optional context for `model`

Level: intermediate

See also: `PetscDA`, `PETSCDAETKF`, `PETSCDALETKF`, `PetscDAEnsembleAnalysis()`

# External Links
$(_doc_external("PetscDA/PetscDAEnsembleForecast"))
"""
function PetscDAEnsembleForecast(petsclib::PetscLibType, da::PetscDA, model::external, ctx::Ptr{Cvoid})
    error("PetscDAEnsembleForecast: no generated method for these argument types")
end

@for_petsc function PetscDAEnsembleForecast(petsclib::$UnionPetscLib, da::PetscDA, model::external, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:PetscDAEnsembleForecast, $petsc_library),
               PetscErrorCode,
               (PetscDA, external, Ptr{Cvoid}),
               da, model, ctx,
              )


	return nothing
end 

"""
	inflation::PetscReal = PetscDAEnsembleGetInflation(petsclib::PetscLibType, da::PetscDA) 
Gets the inflation factor for the data assimilation method.

Not Collective

Input Parameter:
- `da` - the `PetscDA` context

Output Parameter:
- `inflation` - the inflation factor

Level: intermediate

See also: `PetscDA`, `PETSCDAETKF`, `PETSCDALETKF`, `PetscDAEnsembleSetInflation()`

# External Links
$(_doc_external("PetscDA/PetscDAEnsembleGetInflation"))
"""
function PetscDAEnsembleGetInflation(petsclib::PetscLibType, da::PetscDA)
    error("PetscDAEnsembleGetInflation: no generated method for these argument types")
end

@for_petsc function PetscDAEnsembleGetInflation(petsclib::$UnionPetscLib, da::PetscDA )
	inflation_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscDAEnsembleGetInflation, $petsc_library),
               PetscErrorCode,
               (PetscDA, Ptr{$PetscReal}),
               da, inflation_,
              )

	inflation = inflation_[]

	return inflation
end 

"""
	member::PetscVec = PetscDAEnsembleGetMember(petsclib::PetscLibType, da::PetscDA, member_idx::PetscInt) 
Returns a read-only view of an ensemble member stored in the `PetscDA`.

Collective

Input Parameters:
- `da`         - the `PetscDA` context
- `member_idx` - index of the requested member (0 <= idx < ensemble_size)

Output Parameter:
- `member` - read-only vector view; call `PetscDAEnsembleRestoreMember()` when done

Level: intermediate

See also: `PetscDA`, `PETSCDAETKF`, `PETSCDALETKF`, `PetscDAEnsembleRestoreMember()`, `PetscDAEnsembleSetMember()`

# External Links
$(_doc_external("PetscDA/PetscDAEnsembleGetMember"))
"""
function PetscDAEnsembleGetMember(petsclib::PetscLibType, da::PetscDA, member_idx::Integer)
    error("PetscDAEnsembleGetMember: no generated method for these argument types")
end

@for_petsc function PetscDAEnsembleGetMember(petsclib::$UnionPetscLib, da::PetscDA, member_idx::$PetscInt )
	member_ = Ref{CVec}()

    @chk ccall(
               (:PetscDAEnsembleGetMember, $petsc_library),
               PetscErrorCode,
               (PetscDA, $PetscInt, Ptr{CVec}),
               da, member_idx, member_,
              )

	member = PetscVec(member_[], petsclib)

	return member
end 

"""
	ensemble_size::PetscInt = PetscDAEnsembleGetSize(petsclib::PetscLibType, da::PetscDA) 
Retrieves the dimension of the ensemble in a `PetscDA`.

Not Collective

Input Parameter:
- `da` - the `PetscDA` context

Output Parameters:
- `ensemble_size` - number of ensemble members

Level: beginner

See also: `PetscDA`, `PETSCDAETKF`, `PETSCDALETKF`, `PetscDASetSizes()`, `PetscDAGetSizes()`

# External Links
$(_doc_external("PetscDA/PetscDAEnsembleGetSize"))
"""
function PetscDAEnsembleGetSize(petsclib::PetscLibType, da::PetscDA)
    error("PetscDAEnsembleGetSize: no generated method for these argument types")
end

@for_petsc function PetscDAEnsembleGetSize(petsclib::$UnionPetscLib, da::PetscDA )
	ensemble_size_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDAEnsembleGetSize, $petsc_library),
               PetscErrorCode,
               (PetscDA, Ptr{$PetscInt}),
               da, ensemble_size_,
              )

	ensemble_size = ensemble_size_[]

	return ensemble_size
end 

"""
	type::PetscDASqrtType = PetscDAEnsembleGetSqrtType(petsclib::PetscLibType, da::PetscDA) 
Retrieves the current square-root implementation configured for analysis.

Not Collective

Input Parameters:
- `da` - the `PetscDA` object

Output Parameter:
- `type` - on output, the configured `PetscDASqrtType`

Level: advanced

See also: `PetscDA`, `PETSCDAETKF`, `PETSCDALETKF`, `PetscDAEnsembleSetSqrtType()`

# External Links
$(_doc_external("PetscDA/PetscDAEnsembleGetSqrtType"))
"""
function PetscDAEnsembleGetSqrtType(petsclib::PetscLibType, da::PetscDA)
    error("PetscDAEnsembleGetSqrtType: no generated method for these argument types")
end

@for_petsc function PetscDAEnsembleGetSqrtType(petsclib::$UnionPetscLib, da::PetscDA )
	type_ = Ref{PetscDASqrtType}()

    @chk ccall(
               (:PetscDAEnsembleGetSqrtType, $petsc_library),
               PetscErrorCode,
               (PetscDA, Ptr{PetscDASqrtType}),
               da, type_,
              )

	type = type_[]

	return type
end 

"""
	PetscDAEnsembleInitialize(petsclib::PetscLibType, da::PetscDA, x0::AbstractPetscVec, obs_error_std::PetscReal, rng::PetscRandom) 
Initialize ensemble members with Gaussian perturbations

Input Parameters:
- `da`            - PetscDA context
- `x0`            - Background state
- `obs_error_std` - Standard deviation for perturbations
- `rng`           - Random number generator

Level: beginner

See also: `PETSCDAETKF`, `PETSCDALETKF`, `PetscDA`

# External Links
$(_doc_external("PetscDA/PetscDAEnsembleInitialize"))
"""
function PetscDAEnsembleInitialize(petsclib::PetscLibType, da::PetscDA, x0::AbstractPetscVec, obs_error_std::Real, rng::PetscRandom)
    error("PetscDAEnsembleInitialize: no generated method for these argument types")
end

@for_petsc function PetscDAEnsembleInitialize(petsclib::$UnionPetscLib, da::PetscDA, x0::AbstractPetscVec, obs_error_std::$PetscReal, rng::PetscRandom )

    @chk ccall(
               (:PetscDAEnsembleInitialize, $petsc_library),
               PetscErrorCode,
               (PetscDA, CVec, $PetscReal, PetscRandom),
               da, x0, obs_error_std, rng,
              )


	return nothing
end 

"""
	PetscDAEnsembleRestoreMember(petsclib::PetscLibType, da::PetscDA, member_idx::PetscInt, member::AbstractPetscVec) 
Returns a column view obtained with `PetscDAEnsembleGetMember()`.

Collective

Input Parameters:
- `da`         - the `PetscDA` context
- `member_idx` - index that was previously requested
- `member`     - location that holds the view to restore

Level: intermediate

See also: `PetscDA`, `PETSCDAETKF`, `PETSCDALETKF`, `PetscDAEnsembleGetMember()`

# External Links
$(_doc_external("PetscDA/PetscDAEnsembleRestoreMember"))
"""
function PetscDAEnsembleRestoreMember(petsclib::PetscLibType, da::PetscDA, member_idx::Integer, member::AbstractPetscVec)
    error("PetscDAEnsembleRestoreMember: no generated method for these argument types")
end

@for_petsc function PetscDAEnsembleRestoreMember(petsclib::$UnionPetscLib, da::PetscDA, member_idx::$PetscInt, member::AbstractPetscVec )
	member_ = Ref(member.ptr)

    @chk ccall(
               (:PetscDAEnsembleRestoreMember, $petsc_library),
               PetscErrorCode,
               (PetscDA, $PetscInt, Ptr{CVec}),
               da, member_idx, member_,
              )

	member.ptr = member_[]

	return nothing
end 

"""
	PetscDAEnsembleSetInflation(petsclib::PetscLibType, da::PetscDA, inflation::PetscReal) 
Sets the inflation factor for the data assimilation method.

Logically Collective

Input Parameters:
- `da`        - the `PetscDA` context
- `inflation` - the inflation factor (must be >= 1.0)

Level: intermediate

See also: `PetscDA`, `PETSCDAETKF`, `PETSCDALETKF`, `PetscDAEnsembleGetInflation()`

# External Links
$(_doc_external("PetscDA/PetscDAEnsembleSetInflation"))
"""
function PetscDAEnsembleSetInflation(petsclib::PetscLibType, da::PetscDA, inflation::Real)
    error("PetscDAEnsembleSetInflation: no generated method for these argument types")
end

@for_petsc function PetscDAEnsembleSetInflation(petsclib::$UnionPetscLib, da::PetscDA, inflation::$PetscReal )

    @chk ccall(
               (:PetscDAEnsembleSetInflation, $petsc_library),
               PetscErrorCode,
               (PetscDA, $PetscReal),
               da, inflation,
              )


	return nothing
end 

"""
	PetscDAEnsembleSetMember(petsclib::PetscLibType, da::PetscDA, member_idx::PetscInt, member::AbstractPetscVec) 
Overwrites an ensemble member with user-provided state data.

Collective

Input Parameters:
- `da`         - the `PetscDA` context
- `member_idx` - index of the entry to modify
- `member`     - vector containing the new state values

Level: intermediate

See also: `PetscDA`, `PETSCDAETKF`, `PETSCDALETKF`, `PetscDAEnsembleGetMember()`

# External Links
$(_doc_external("PetscDA/PetscDAEnsembleSetMember"))
"""
function PetscDAEnsembleSetMember(petsclib::PetscLibType, da::PetscDA, member_idx::Integer, member::AbstractPetscVec)
    error("PetscDAEnsembleSetMember: no generated method for these argument types")
end

@for_petsc function PetscDAEnsembleSetMember(petsclib::$UnionPetscLib, da::PetscDA, member_idx::$PetscInt, member::AbstractPetscVec )

    @chk ccall(
               (:PetscDAEnsembleSetMember, $petsc_library),
               PetscErrorCode,
               (PetscDA, $PetscInt, CVec),
               da, member_idx, member,
              )


	return nothing
end 

"""
	PetscDAEnsembleSetSize(petsclib::PetscLibType, da::PetscDA, ensemble_size::PetscInt) 
Sets the ensemble dimensions used by a `PetscDA`.

Collective

Input Parameters:
- `da`            - the `PetscDA` context
- `ensemble_size` - number of ensemble members

Options Database Key:
- `-petscda_ensemble_size <size>` - number of ensemble members

Level: beginner

See also: `PetscDA`, `PETSCDAETKF`, `PETSCDALETKF`, `PetscDAGetSizes()`, `PetscDASetSizes()`, `PetscDASetUp()`

# External Links
$(_doc_external("PetscDA/PetscDAEnsembleSetSize"))
"""
function PetscDAEnsembleSetSize(petsclib::PetscLibType, da::PetscDA, ensemble_size::Integer)
    error("PetscDAEnsembleSetSize: no generated method for these argument types")
end

@for_petsc function PetscDAEnsembleSetSize(petsclib::$UnionPetscLib, da::PetscDA, ensemble_size::$PetscInt )

    @chk ccall(
               (:PetscDAEnsembleSetSize, $petsc_library),
               PetscErrorCode,
               (PetscDA, $PetscInt),
               da, ensemble_size,
              )


	return nothing
end 

"""
	PetscDAEnsembleSetSqrtType(petsclib::PetscLibType, da::PetscDA, type::PetscDASqrtType) 
Selects the reduced-space square-root algorithm used during analysis.

Logically Collective

Input Parameters:
- `da`   - the `PetscDA` object
- `type` - either `PETSCDA_SQRT_CHOLESKY` or `PETSCDA_SQRT_EIGEN`

Options Database Key:
- `-petscda_ensemble_sqrt_type <cholesky or eigen>` - set the `PetscDASqrtType`

Level: advanced

See also: `PetscDA`, `PETSCDAETKF`, `PETSCDALETKF`, `PetscDASqrtType`, `PetscDAEnsembleGetSqrtType()`

# External Links
$(_doc_external("PetscDA/PetscDAEnsembleSetSqrtType"))
"""
function PetscDAEnsembleSetSqrtType(petsclib::PetscLibType, da::PetscDA, type::PetscDASqrtType)
    error("PetscDAEnsembleSetSqrtType: no generated method for these argument types")
end

@for_petsc function PetscDAEnsembleSetSqrtType(petsclib::$UnionPetscLib, da::PetscDA, type::PetscDASqrtType )

    @chk ccall(
               (:PetscDAEnsembleSetSqrtType, $petsc_library),
               PetscErrorCode,
               (PetscDA, PetscDASqrtType),
               da, type,
              )


	return nothing
end 

"""
	PetscDAEnsembleTFactor(petsclib::PetscLibType, da::PetscDA, S::AbstractPetscMat) 
Compute and store factorization of T matrix

Collective

Input Parameters:
- `da` - the `PetscDA` context
- `S`  - normalized innovation matrix (obs_size x m)

See also: `PetscDA`, `PETSCDAETKF`, `PETSCDALETKF`, `PetscDAEnsembleApplyTInverse()`, `PetscDAEnsembleApplySqrtTInverse()`

# External Links
$(_doc_external("PetscDA/PetscDAEnsembleTFactor"))
"""
function PetscDAEnsembleTFactor(petsclib::PetscLibType, da::PetscDA, S::AbstractPetscMat)
    error("PetscDAEnsembleTFactor: no generated method for these argument types")
end

@for_petsc function PetscDAEnsembleTFactor(petsclib::$UnionPetscLib, da::PetscDA, S::AbstractPetscMat )

    @chk ccall(
               (:PetscDAEnsembleTFactor, $petsc_library),
               PetscErrorCode,
               (PetscDA, CMat),
               da, S,
              )


	return nothing
end 

"""
	PetscDAFinalizePackage(petsclib::PetscLibType) 
This function finalizes everything in the `PetscDA` package. It
is called from `PetscFinalize()`.

Logically Collective

Level: developer

See also: `PetscDAInitializePackage()`, `PetscInitialize()`

# External Links
$(_doc_external("PetscDA/PetscDAFinalizePackage"))
"""
function PetscDAFinalizePackage(petsclib::PetscLibType)
    error("PetscDAFinalizePackage: no generated method for these argument types")
end

@for_petsc function PetscDAFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscDAFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	ndof::PetscInt = PetscDAGetNDOF(petsclib::PetscLibType, da::PetscDA) 
Get the number of degrees of freedom per grid point

Not Collective

Input Parameter:
- `da` - the `PetscDA` context

Output Parameter:
- `ndof` - number of degrees of freedom per grid point

Level: intermediate

See also: `PetscDA`, `PetscDASetNDOF()`

# External Links
$(_doc_external("PetscDA/PetscDAGetNDOF"))
"""
function PetscDAGetNDOF(petsclib::PetscLibType, da::PetscDA)
    error("PetscDAGetNDOF: no generated method for these argument types")
end

@for_petsc function PetscDAGetNDOF(petsclib::$UnionPetscLib, da::PetscDA )
	ndof_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDAGetNDOF, $petsc_library),
               PetscErrorCode,
               (PetscDA, Ptr{$PetscInt}),
               da, ndof_,
              )

	ndof = ndof_[]

	return ndof
end 

"""
	obs_error_var::PetscVec = PetscDAGetObsErrorVariance(petsclib::PetscLibType, da::PetscDA) 
Returns a borrowed reference to the observation-error variance vector.

Not Collective

Input Parameter:
- `da` - the `PetscDA` context

Output Parameter:
- `obs_error_var` - pointer to the variance vector managed by the `PetscDA`

Level: beginner

See also: `PetscDASetObsErrorVariance()`

# External Links
$(_doc_external("PetscDA/PetscDAGetObsErrorVariance"))
"""
function PetscDAGetObsErrorVariance(petsclib::PetscLibType, da::PetscDA)
    error("PetscDAGetObsErrorVariance: no generated method for these argument types")
end

@for_petsc function PetscDAGetObsErrorVariance(petsclib::$UnionPetscLib, da::PetscDA )
	obs_error_var_ = Ref{CVec}()

    @chk ccall(
               (:PetscDAGetObsErrorVariance, $petsc_library),
               PetscErrorCode,
               (PetscDA, Ptr{CVec}),
               da, obs_error_var_,
              )

	obs_error_var = PetscVec(obs_error_var_[], petsclib)

	return obs_error_var
end 

"""
	p::Ptr{Cchar} = PetscDAGetOptionsPrefix(petsclib::PetscLibType, das::PetscDA) 
Gets the prefix used for searching for all
PetscDA options in the database

Not Collective

Input Parameter:
- `das` - the `PetscDA` context

Output Parameter:
- `p` - pointer to the prefix string used

Level: advanced

See also: `PetscDA`, `PetscDASetFromOptions()`, `PetscDASetOptionsPrefix()`, `PetscDAAppendOptionsPrefix()`

# External Links
$(_doc_external("PetscDA/PetscDAGetOptionsPrefix"))
"""
function PetscDAGetOptionsPrefix(petsclib::PetscLibType, das::PetscDA)
    error("PetscDAGetOptionsPrefix: no generated method for these argument types")
end

@for_petsc function PetscDAGetOptionsPrefix(petsclib::$UnionPetscLib, das::PetscDA )
	p_ = Ref{Ptr{Cchar}}()

    @chk ccall(
               (:PetscDAGetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (PetscDA, Ptr{Ptr{Cchar}}),
               das, p_,
              )

	p = p_[]

	return p
end 

"""
	state_size::PetscInt,obs_size::PetscInt = PetscDAGetSizes(petsclib::PetscLibType, da::PetscDA) 
Retrieves the state size and observation size from a `PetscDA`.

Not Collective

Input Parameter:
- `da` - the `PetscDA` context

Output Parameters:
- `state_size` - number of state components (may be `NULL`)
- `obs_size`   - number of observation components (may be `NULL`)

Level: beginner

See also: `PetscDASetSizes()`

# External Links
$(_doc_external("PetscDA/PetscDAGetSizes"))
"""
function PetscDAGetSizes(petsclib::PetscLibType, da::PetscDA)
    error("PetscDAGetSizes: no generated method for these argument types")
end

@for_petsc function PetscDAGetSizes(petsclib::$UnionPetscLib, da::PetscDA )
	state_size_ = Ref{$PetscInt}()
	obs_size_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDAGetSizes, $petsc_library),
               PetscErrorCode,
               (PetscDA, Ptr{$PetscInt}, Ptr{$PetscInt}),
               da, state_size_, obs_size_,
              )

	state_size = state_size_[]
	obs_size = obs_size_[]

	return state_size,obs_size
end 

"""
	type::PetscDAType = PetscDAGetType(petsclib::PetscLibType, da::PetscDA) 
Gets the name of the implementation currently associated with a `PetscDA`.

Not Collective

Input Parameter:
- `da` - the `PetscDA` context

Output Parameter:
- `type` - pointer that will receive the type name (may be `NULL`)

Level: intermediate

See also: `PetscDASetType()`

# External Links
$(_doc_external("PetscDA/PetscDAGetType"))
"""
function PetscDAGetType(petsclib::PetscLibType, da::PetscDA)
    error("PetscDAGetType: no generated method for these argument types")
end

@for_petsc function PetscDAGetType(petsclib::$UnionPetscLib, da::PetscDA )
	type_ = Ref{PetscDAType}()

    @chk ccall(
               (:PetscDAGetType, $petsc_library),
               PetscErrorCode,
               (PetscDA, Ptr{PetscDAType}),
               da, type_,
              )

	type = type_[] == C_NULL ? "" : unsafe_string(type_[])

	return type
end 

"""
	PetscDAInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `PetscDA`
package. called on the first call to `PetscDACreate()` when using static or shared
libraries.

Logically Collective

Level: developer

See also: `PetscDAFinalizePackage()`, `PetscInitialize()`

# External Links
$(_doc_external("PetscDA/PetscDAInitializePackage"))
"""
function PetscDAInitializePackage(petsclib::PetscLibType)
    error("PetscDAInitializePackage: no generated method for these argument types")
end

@for_petsc function PetscDAInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscDAInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscDALETKFGetLocalizationMatrix(petsclib::PetscLibType, n_obs_vertex::PetscInt, n_dof::PetscInt, Vecxyz::Vector{<:AbstractPetscVec}, bd::Vector{PetscReal}, H::AbstractPetscMat, Q::AbstractPetscMat) 

# External Links
$(_doc_external("PC/PetscDALETKFGetLocalizationMatrix"))
"""
function PetscDALETKFGetLocalizationMatrix(petsclib::PetscLibType, n_obs_vertex::Integer, n_dof::Integer, Vecxyz::Vector{<:AbstractPetscVec}, bd::AbstractVector{<:Number}, H::AbstractPetscMat, Q::AbstractPetscMat)
    error("PetscDALETKFGetLocalizationMatrix: no generated method for these argument types")
end

@for_petsc function PetscDALETKFGetLocalizationMatrix(petsclib::$UnionPetscLib, n_obs_vertex::$PetscInt, n_dof::$PetscInt, Vecxyz::Vector{<:AbstractPetscVec}, bd::Vector{$PetscReal}, H::AbstractPetscMat, Q::AbstractPetscMat )
	Q_ = Ref(Q.ptr)

    @chk ccall(
               (:PetscDALETKFGetLocalizationMatrix, $petsc_library),
               PetscErrorCode,
               ($PetscInt, $PetscInt, Ptr{CVec}, Ptr{$PetscReal}, CMat, Ptr{CMat}),
               n_obs_vertex, n_dof, Vecxyz, bd, H, Q_,
              )

	Q.ptr = Q_[]

	return nothing
end 

"""
	n_obs_vertex::PetscInt = PetscDALETKFGetObsPerVertex(petsclib::PetscLibType, da::PetscDA) 
Gets the number of local observations per vertex for the LETKF algorithm.

Not Collective

Input Parameter:
- `da` - the `PetscDA` context

Output Parameter:
- `n_obs_vertex` - number of observations per vertex

Level: advanced

See also: `PETSCDALETKF`, `PetscDA`, `PetscDALETKFSetObsPerVertex()`

# External Links
$(_doc_external("PetscDA/PetscDALETKFGetObsPerVertex"))
"""
function PetscDALETKFGetObsPerVertex(petsclib::PetscLibType, da::PetscDA)
    error("PetscDALETKFGetObsPerVertex: no generated method for these argument types")
end

@for_petsc function PetscDALETKFGetObsPerVertex(petsclib::$UnionPetscLib, da::PetscDA )
	n_obs_vertex_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDALETKFGetObsPerVertex, $petsc_library),
               PetscErrorCode,
               (PetscDA, Ptr{$PetscInt}),
               da, n_obs_vertex_,
              )

	n_obs_vertex = n_obs_vertex_[]

	return n_obs_vertex
end 

"""
	PetscDALETKFSetLocalization(petsclib::PetscLibType, da::PetscDA, Q::AbstractPetscMat, H::AbstractPetscMat) 
Sets the localization matrix for the LETKF algorithm.

Collective

Input Parameters:
- `da` - the `PetscDA` context
- `Q`  - the localization matrix (N x P)
- `H`  - the observation operator matrix (P x N)

Level: advanced

See also: `PETSCDALETKF`, `PetscDA`

# External Links
$(_doc_external("PetscDA/PetscDALETKFSetLocalization"))
"""
function PetscDALETKFSetLocalization(petsclib::PetscLibType, da::PetscDA, Q::AbstractPetscMat, H::AbstractPetscMat)
    error("PetscDALETKFSetLocalization: no generated method for these argument types")
end

@for_petsc function PetscDALETKFSetLocalization(petsclib::$UnionPetscLib, da::PetscDA, Q::AbstractPetscMat, H::AbstractPetscMat )

    @chk ccall(
               (:PetscDALETKFSetLocalization, $petsc_library),
               PetscErrorCode,
               (PetscDA, CMat, CMat),
               da, Q, H,
              )


	return nothing
end 

"""
	PetscDALETKFSetObsPerVertex(petsclib::PetscLibType, da::PetscDA, n_obs_vertex::PetscInt) 
Sets the number of local observations per vertex for the LETKF algorithm.

Logically Collective

Input Parameters:
- `da`           - the `PetscDA` context
- `n_obs_vertex` - number of observations per vertex

Level: advanced

See also: `PETSCDALETKF`, `PetscDA`, `PetscDALETKFSetLocalization()`

# External Links
$(_doc_external("PetscDA/PetscDALETKFSetObsPerVertex"))
"""
function PetscDALETKFSetObsPerVertex(petsclib::PetscLibType, da::PetscDA, n_obs_vertex::Integer)
    error("PetscDALETKFSetObsPerVertex: no generated method for these argument types")
end

@for_petsc function PetscDALETKFSetObsPerVertex(petsclib::$UnionPetscLib, da::PetscDA, n_obs_vertex::$PetscInt )

    @chk ccall(
               (:PetscDALETKFSetObsPerVertex, $petsc_library),
               PetscErrorCode,
               (PetscDA, $PetscInt),
               da, n_obs_vertex,
              )


	return nothing
end 

"""
	PetscDARegister(petsclib::PetscLibType, sname::String, fnc::external) 
Registers a constructor for a `PetscDA` implementation with the
dispatcher.

Not Collective

Input Parameters:
- `sname`    - name associated with the implementation
- `function` - routine that creates the implementation and installs method table

Level: developer

See also: `PetscDARegisterAll()`, `PetscDASetType()`

# External Links
$(_doc_external("PetscDA/PetscDARegister"))
"""
function PetscDARegister(petsclib::PetscLibType, sname::String, fnc::external)
    error("PetscDARegister: no generated method for these argument types")
end

@for_petsc function PetscDARegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:PetscDARegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	PetscDARegisterAll(petsclib::PetscLibType) 
Registers all data assimilation backends that were compiled in.

Not Collective

Level: developer

See also: `PetscDARegister()`

# External Links
$(_doc_external("PetscDA/PetscDARegisterAll"))
"""
function PetscDARegisterAll(petsclib::PetscLibType)
    error("PetscDARegisterAll: no generated method for these argument types")
end

@for_petsc function PetscDARegisterAll(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscDARegisterAll, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscDASetFromOptions(petsclib::PetscLibType, da::PetscDA) 
Configures a `PetscDA` object from the options database.

Collective

Input Parameter:
- `da` - the `PetscDA` context to set up

Level: intermediate

See also: `PetscDASetType()`, `PetscObjectOptionsBegin()`

# External Links
$(_doc_external("PetscDA/PetscDASetFromOptions"))
"""
function PetscDASetFromOptions(petsclib::PetscLibType, da::PetscDA)
    error("PetscDASetFromOptions: no generated method for these argument types")
end

@for_petsc function PetscDASetFromOptions(petsclib::$UnionPetscLib, da::PetscDA )

    @chk ccall(
               (:PetscDASetFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscDA,),
               da,
              )


	return nothing
end 

"""
	PetscDASetLocalSizes(petsclib::PetscLibType, da::PetscDA, local_state_size::PetscInt, local_obs_size::PetscInt) 
Sets the local state and observation dimensions used by a `PetscDA`.

Collective

Input Parameters:
- `da`               - the `PetscDA` context
- `local_state_size` - number of local state components (or `PETSC_DECIDE`)
- `local_obs_size`   - number of local observation components (or `PETSC_DECIDE`)

Level: beginner

See also: `PetscDASetSizes()`, `PetscDASetUp()`

# External Links
$(_doc_external("PetscDA/PetscDASetLocalSizes"))
"""
function PetscDASetLocalSizes(petsclib::PetscLibType, da::PetscDA, local_state_size::Integer, local_obs_size::Integer)
    error("PetscDASetLocalSizes: no generated method for these argument types")
end

@for_petsc function PetscDASetLocalSizes(petsclib::$UnionPetscLib, da::PetscDA, local_state_size::$PetscInt, local_obs_size::$PetscInt )

    @chk ccall(
               (:PetscDASetLocalSizes, $petsc_library),
               PetscErrorCode,
               (PetscDA, $PetscInt, $PetscInt),
               da, local_state_size, local_obs_size,
              )


	return nothing
end 

"""
	PetscDASetNDOF(petsclib::PetscLibType, da::PetscDA, ndof::PetscInt) 
Set the number of degrees of freedom per grid point

Logically Collective

Input Parameters:
- `da`   - the `PetscDA` context
- `ndof` - number of degrees of freedom per grid point (e.g., 2 for shallow water with h and hu)

Level: intermediate

See also: `PetscDA`, `PetscDAGetNDOF()`, `PetscDASetUp()`, `PetscDASetSizes()`

# External Links
$(_doc_external("PetscDA/PetscDASetNDOF"))
"""
function PetscDASetNDOF(petsclib::PetscLibType, da::PetscDA, ndof::Integer)
    error("PetscDASetNDOF: no generated method for these argument types")
end

@for_petsc function PetscDASetNDOF(petsclib::$UnionPetscLib, da::PetscDA, ndof::$PetscInt )

    @chk ccall(
               (:PetscDASetNDOF, $petsc_library),
               PetscErrorCode,
               (PetscDA, $PetscInt),
               da, ndof,
              )


	return nothing
end 

"""
	PetscDASetObsErrorVariance(petsclib::PetscLibType, da::PetscDA, obs_error_var::AbstractPetscVec) 
Sets the observation-error variances associated with a `PetscDA`.

Collective

Input Parameters:
- `da`            - the `PetscDA` context
- `obs_error_var` - vector containing observation error variances (assumes R is a diagonal matrix)

See also: `PetscDAGetObsErrorVariance()`

# External Links
$(_doc_external("PetscDA/PetscDASetObsErrorVariance"))
"""
function PetscDASetObsErrorVariance(petsclib::PetscLibType, da::PetscDA, obs_error_var::AbstractPetscVec)
    error("PetscDASetObsErrorVariance: no generated method for these argument types")
end

@for_petsc function PetscDASetObsErrorVariance(petsclib::$UnionPetscLib, da::PetscDA, obs_error_var::AbstractPetscVec )

    @chk ccall(
               (:PetscDASetObsErrorVariance, $petsc_library),
               PetscErrorCode,
               (PetscDA, CVec),
               da, obs_error_var,
              )


	return nothing
end 

"""
	PetscDASetOptionsPrefix(petsclib::PetscLibType, das::PetscDA, p::String) 
Sets the prefix used for searching for all
PetscDA options in the database.

Logically Collective

Input Parameters:
- `das` - the `PetscDA` context
- `p`   - the prefix string to prepend to all PetscDA option requests

Level: advanced

See also: `PetscDA`, `PetscDASetFromOptions()`, `PetscDAAppendOptionsPrefix()`, `PetscDAGetOptionsPrefix()`

# External Links
$(_doc_external("PetscDA/PetscDASetOptionsPrefix"))
"""
function PetscDASetOptionsPrefix(petsclib::PetscLibType, das::PetscDA, p::String)
    error("PetscDASetOptionsPrefix: no generated method for these argument types")
end

@for_petsc function PetscDASetOptionsPrefix(petsclib::$UnionPetscLib, das::PetscDA, p::String )

    @chk ccall(
               (:PetscDASetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (PetscDA, Ptr{Cchar}),
               das, p,
              )


	return nothing
end 

"""
	PetscDASetSizes(petsclib::PetscLibType, da::PetscDA, state_size::PetscInt, obs_size::PetscInt) 
Sets the state and observation sizes for a `PetscDA`

Collective

Input Parameters:
- `da`         - the `PetscDA` context
- `state_size` - number of state components
- `obs_size`   - number of observation components

Level: beginner

See also: `PetscDAGetSizes()`, `PetscDASetUp()`, `PetscDAEnsembleSetSize()`

# External Links
$(_doc_external("PetscDA/PetscDASetSizes"))
"""
function PetscDASetSizes(petsclib::PetscLibType, da::PetscDA, state_size::Integer, obs_size::Integer)
    error("PetscDASetSizes: no generated method for these argument types")
end

@for_petsc function PetscDASetSizes(petsclib::$UnionPetscLib, da::PetscDA, state_size::$PetscInt, obs_size::$PetscInt )

    @chk ccall(
               (:PetscDASetSizes, $petsc_library),
               PetscErrorCode,
               (PetscDA, $PetscInt, $PetscInt),
               da, state_size, obs_size,
              )


	return nothing
end 

"""
	PetscDASetType(petsclib::PetscLibType, da::PetscDA, type::PetscDAType) 
Sets the data assimilation implementation used by a `PetscDA` object.

Collective

Input Parameters:
- `da`   - the `PetscDA` context
- `type` - name of the implementation (for example `PETSCDAETKF`)

Level: intermediate

See also: `PetscDAGetType()`, `PetscDARegister()`

# External Links
$(_doc_external("PetscDA/PetscDASetType"))
"""
function PetscDASetType(petsclib::PetscLibType, da::PetscDA, type::PetscDAType)
    error("PetscDASetType: no generated method for these argument types")
end

@for_petsc function PetscDASetType(petsclib::$UnionPetscLib, da::PetscDA, type::PetscDAType )

    @chk ccall(
               (:PetscDASetType, $petsc_library),
               PetscErrorCode,
               (PetscDA, PetscDAType),
               da, type,
              )


	return nothing
end 

"""
	PetscDASetUp(petsclib::PetscLibType, da::PetscDA) 
Allocates internal data structures for a `PetscDA` based on the previously provided sizes.

Collective

Input Parameter:
- `da` - the `PetscDA` context to assemble

Level: beginner

See also: `PetscDASetSizes()`, `PetscDASetType()`

# External Links
$(_doc_external("PetscDA/PetscDASetUp"))
"""
function PetscDASetUp(petsclib::PetscLibType, da::PetscDA)
    error("PetscDASetUp: no generated method for these argument types")
end

@for_petsc function PetscDASetUp(petsclib::$UnionPetscLib, da::PetscDA )

    @chk ccall(
               (:PetscDASetUp, $petsc_library),
               PetscErrorCode,
               (PetscDA,),
               da,
              )


	return nothing
end 

"""
	PetscDAView(petsclib::PetscLibType, da::PetscDA, viewer::PetscViewer) 
Views a `PetscDA` and its implementation-specific data structure.

Collective

Input Parameters:
- `da`     - the `PetscDA` context
- `viewer` - the `PetscViewer` to use (or `NULL` for standard output)

Level: beginner

See also: `PetscDAViewFromOptions()`

# External Links
$(_doc_external("PetscDA/PetscDAView"))
"""
function PetscDAView(petsclib::PetscLibType, da::PetscDA, viewer::PetscViewer)
    error("PetscDAView: no generated method for these argument types")
end

@for_petsc function PetscDAView(petsclib::$UnionPetscLib, da::PetscDA, viewer::PetscViewer )

    @chk ccall(
               (:PetscDAView, $petsc_library),
               PetscErrorCode,
               (PetscDA, PetscViewer),
               da, viewer,
              )


	return nothing
end 

"""
	PetscDAViewFromOptions(petsclib::PetscLibType, da::PetscDA, obj, name::String) 
Processes command-line options to determine if a `PetscDA` should be viewed.

Collective

Input Parameters:
- `da`   - the `PetscDA` context
- `obj`  - optional object that provides the prefix for options
- `name` - option name to check

Options Database Key:
- `-name [viewertype][:...]` - option name and values. See `PetscObjectViewFromOptions()` for the possible arguments

Level: beginner

See also: `PetscDAView()`, `PetscObjectViewFromOptions()`

# External Links
$(_doc_external("PetscDA/PetscDAViewFromOptions"))
"""
function PetscDAViewFromOptions(petsclib::PetscLibType, da::PetscDA, obj, name::String)
    error("PetscDAViewFromOptions: no generated method for these argument types")
end

@for_petsc function PetscDAViewFromOptions(petsclib::$UnionPetscLib, da::PetscDA, obj, name::String )

    @chk ccall(
               (:PetscDAViewFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscDA, PetscObject, Ptr{Cchar}),
               da, obj, name,
              )


	return nothing
end 

"""
	ptype::PetscDataType,found::PetscBool = PetscDataTypeFromString(petsclib::PetscLibType, name::String) 
Gets the enum value of a PETSc datatype represented as a string

Not Collective

Input Parameter:
- `name` - the PETSc datatype name (for example, "double" or "real")

Output Parameters:
- `ptype` - the enum value, only valid if found is `PETSC_TRUE`
- `found` - the string matches one of the data types

Level: advanced

See also: `PetscDataType`, `PetscDataTypeToMPIDataType()`, `PetscDataTypeGetSize()`

# External Links
$(_doc_external("Sys/PetscDataTypeFromString"))
"""
function PetscDataTypeFromString(petsclib::PetscLibType, name::String)
    error("PetscDataTypeFromString: no generated method for these argument types")
end

@for_petsc function PetscDataTypeFromString(petsclib::$UnionPetscLib, name::String )
	ptype_ = Ref{PetscDataType}()
	found_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscDataTypeFromString, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{PetscDataType}, Ptr{PetscBool}),
               name, ptype_, found_,
              )

	ptype = ptype_[]
	found = found_[]

	return ptype,found
end 

"""
	size::Csize_t = PetscDataTypeGetSize(petsclib::PetscLibType, ptype::PetscDataType) 
Gets the size (in bytes) of a PETSc datatype

Not Collective

Input Parameter:
- `ptype` - the PETSc datatype name (for example `PETSC_DOUBLE`)

Output Parameter:
- `size` - the size in bytes (for example the size of `PETSC_DOUBLE` is 8)

Level: advanced

See also: `PetscDataType`, `PetscDataTypeToMPIDataType()`

# External Links
$(_doc_external("Sys/PetscDataTypeGetSize"))
"""
function PetscDataTypeGetSize(petsclib::PetscLibType, ptype::PetscDataType)
    error("PetscDataTypeGetSize: no generated method for these argument types")
end

@for_petsc function PetscDataTypeGetSize(petsclib::$UnionPetscLib, ptype::PetscDataType )
	size_ = Ref{Csize_t}()

    @chk ccall(
               (:PetscDataTypeGetSize, $petsc_library),
               PetscErrorCode,
               (PetscDataType, Ptr{Csize_t}),
               ptype, size_,
              )

	size = size_[]

	return size
end 

"""
	htype::hid_t = PetscDataTypeToHDF5DataType(petsclib::PetscLibType, ptype::PetscDataType) 
Converts the PETSc name of a datatype to its HDF5 name.

Not Collective

Input Parameter:
- `ptype` - the PETSc datatype name (for example `PETSC_DOUBLE`)

Output Parameter:
- `htype` - the HDF5 datatype

Level: advanced

See also: [](sec_viewers), `PetscDataType`, `PetscHDF5DataTypeToPetscDataType()`

# External Links
$(_doc_external("Viewer/PetscDataTypeToHDF5DataType"))
"""
function PetscDataTypeToHDF5DataType(petsclib::PetscLibType, ptype::PetscDataType)
    error("PetscDataTypeToHDF5DataType: no generated method for these argument types")
end

@for_petsc function PetscDataTypeToHDF5DataType(petsclib::$UnionPetscLib, ptype::PetscDataType )
	htype_ = Ref{hid_t}()

    @chk ccall(
               (:PetscDataTypeToHDF5DataType, $petsc_library),
               PetscErrorCode,
               (PetscDataType, Ptr{hid_t}),
               ptype, htype_,
              )

	htype = htype_[]

	return htype
end 

