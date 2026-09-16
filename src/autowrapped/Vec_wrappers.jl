"""
	VecAXPBY(petsclib::PetscLibType, y::AbstractPetscVec, alpha::PetscScalar, beta::PetscScalar, x::AbstractPetscVec) 
Computes `y = alpha x + beta y`.

Logically Collective

Input Parameters:
- `alpha` - first scalar
- `beta`  - second scalar
- `x`     - the first scaled vector
- `y`     - the second scaled vector

Output Parameter:
- `y` - output vector

Level: intermediate

See also: `Vec`, `VecAYPX()`, `VecMAXPY()`, `VecWAXPY()`, `VecAXPY()`, `VecAXPBYPCZ()`

# External Links
$(_doc_external("Vec/VecAXPBY"))
"""
function VecAXPBY(petsclib::PetscLibType, y::AbstractPetscVec, alpha::Number, beta::Number, x::AbstractPetscVec)
    error("VecAXPBY: no generated method for these argument types")
end

@for_petsc function VecAXPBY(petsclib::$UnionPetscLib, y::AbstractPetscVec, alpha::$PetscScalar, beta::$PetscScalar, x::AbstractPetscVec )

    @chk ccall(
               (:VecAXPBY, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscScalar, $PetscScalar, CVec),
               y, alpha, beta, x,
              )


	return nothing
end 

"""
	VecAXPBYPCZ(petsclib::PetscLibType, z::AbstractPetscVec, alpha::PetscScalar, beta::PetscScalar, gamma::PetscScalar, x::AbstractPetscVec, y::AbstractPetscVec) 
Computes `z = alpha x + beta y + gamma z`

Logically Collective

Input Parameters:
- `alpha` - first scalar
- `beta`  - second scalar
- `gamma` - third scalar
- `x`     - first vector
- `y`     - second vector
- `z`     - third vector

Output Parameter:
- `z` - output vector

Level: intermediate

See also: `Vec`, `VecAYPX()`, `VecMAXPY()`, `VecWAXPY()`, `VecAXPY()`, `VecAXPBY()`

# External Links
$(_doc_external("Vec/VecAXPBYPCZ"))
"""
function VecAXPBYPCZ(petsclib::PetscLibType, z::AbstractPetscVec, alpha::Number, beta::Number, gamma::Number, x::AbstractPetscVec, y::AbstractPetscVec)
    error("VecAXPBYPCZ: no generated method for these argument types")
end

@for_petsc function VecAXPBYPCZ(petsclib::$UnionPetscLib, z::AbstractPetscVec, alpha::$PetscScalar, beta::$PetscScalar, gamma::$PetscScalar, x::AbstractPetscVec, y::AbstractPetscVec )

    @chk ccall(
               (:VecAXPBYPCZ, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscScalar, $PetscScalar, $PetscScalar, CVec, CVec),
               z, alpha, beta, gamma, x, y,
              )


	return nothing
end 

"""
	VecAXPY(petsclib::PetscLibType, y::AbstractPetscVec, alpha::PetscScalar, x::AbstractPetscVec) 
Computes `y = alpha x + y`.

Logically Collective

Input Parameters:
- `alpha` - the scalar
- `x`     - vector scale by `alpha`
- `y`     - vector accumulated into

Output Parameter:
- `y` - output vector

Level: intermediate

See also: `Vec`, `VecAYPX()`, `VecMAXPY()`, `VecWAXPY()`, `VecAXPBYPCZ()`, `VecAXPBY()`

# External Links
$(_doc_external("Vec/VecAXPY"))
"""
function VecAXPY(petsclib::PetscLibType, y::AbstractPetscVec, alpha::Number, x::AbstractPetscVec)
    error("VecAXPY: no generated method for these argument types")
end

@for_petsc function VecAXPY(petsclib::$UnionPetscLib, y::AbstractPetscVec, alpha::$PetscScalar, x::AbstractPetscVec )

    @chk ccall(
               (:VecAXPY, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscScalar, CVec),
               y, alpha, x,
              )


	return nothing
end 

"""
	VecAYPX(petsclib::PetscLibType, y::AbstractPetscVec, beta::PetscScalar, x::AbstractPetscVec) 
Computes `y = x + beta y`.

Logically Collective

Input Parameters:
- `beta` - the scalar
- `x`    - the unscaled vector
- `y`    - the vector to be scaled

Output Parameter:
- `y` - output vector

Level: intermediate

See also: `Vec`, `VecMAXPY()`, `VecWAXPY()`, `VecAXPY()`, `VecAXPBYPCZ()`, `VecAXPBY()`

# External Links
$(_doc_external("Vec/VecAYPX"))
"""
function VecAYPX(petsclib::PetscLibType, y::AbstractPetscVec, beta::Number, x::AbstractPetscVec)
    error("VecAYPX: no generated method for these argument types")
end

@for_petsc function VecAYPX(petsclib::$UnionPetscLib, y::AbstractPetscVec, beta::$PetscScalar, x::AbstractPetscVec )

    @chk ccall(
               (:VecAYPX, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscScalar, CVec),
               y, beta, x,
              )


	return nothing
end 

"""
	VecAbs(petsclib::PetscLibType, v::AbstractPetscVec) 
Replaces every element in a vector with its absolute value.

Logically Collective

Input Parameter:
- `v` - the vector

Level: intermediate

See also: `Vec`, `VecExp()`, `VecSqrtAbs()`, `VecReciprocal()`, `VecLog()`, `VecPointwiseSign()`

# External Links
$(_doc_external("Vec/VecAbs"))
"""
function VecAbs(petsclib::PetscLibType, v::AbstractPetscVec)
    error("VecAbs: no generated method for these argument types")
end

@for_petsc function VecAbs(petsclib::$UnionPetscLib, v::AbstractPetscVec )

    @chk ccall(
               (:VecAbs, $petsc_library),
               PetscErrorCode,
               (CVec,),
               v,
              )


	return nothing
end 

"""
	VecAppendOptionsPrefix(petsclib::PetscLibType, v::AbstractPetscVec, prefix::String) 
Appends to the prefix used for searching for all
`Vec` options in the database.

Logically Collective

Input Parameters:
- `v`      - the `Vec` context
- `prefix` - the prefix to prepend to all option names

Level: advanced

See also: `Vec`, `VecGetOptionsPrefix()`

# External Links
$(_doc_external("Vec/VecAppendOptionsPrefix"))
"""
function VecAppendOptionsPrefix(petsclib::PetscLibType, v::AbstractPetscVec, prefix::String)
    error("VecAppendOptionsPrefix: no generated method for these argument types")
end

@for_petsc function VecAppendOptionsPrefix(petsclib::$UnionPetscLib, v::AbstractPetscVec, prefix::String )

    @chk ccall(
               (:VecAppendOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Cchar}),
               v, prefix,
              )


	return nothing
end 

"""
	VecAssemblyBegin(petsclib::PetscLibType, vec::AbstractPetscVec) 
Begins assembling the vector; that is ensuring all the vector's entries are stored on the correct MPI process. This routine should
be called after completing all calls to `VecSetValues()`.

Collective

Input Parameter:
- `vec` - the vector

Level: beginner

See also: `Vec`, `VecAssemblyEnd()`, `VecSetValues()`

# External Links
$(_doc_external("Vec/VecAssemblyBegin"))
"""
function VecAssemblyBegin(petsclib::PetscLibType, vec::AbstractPetscVec)
    error("VecAssemblyBegin: no generated method for these argument types")
end

@for_petsc function VecAssemblyBegin(petsclib::$UnionPetscLib, vec::AbstractPetscVec )

    @chk ccall(
               (:VecAssemblyBegin, $petsc_library),
               PetscErrorCode,
               (CVec,),
               vec,
              )


	return nothing
end 

"""
	VecAssemblyEnd(petsclib::PetscLibType, vec::AbstractPetscVec) 
Completes assembling the vector.  This routine should be called after `VecAssemblyBegin()`.

Collective

Input Parameter:
- `vec` - the vector

Options Database Keys:
- `-vec_view [viewertype][:...]`      - Display the vector. See `VecViewFromOptions()`/`PetscObjectViewFromOptions()` for the possible arguments
- `-vecstash_view [viewertype][:...]` - Display the vector stash. See `VecStashViewFromOptions()`/`PetscObjectViewFromOptions()` for the possible arguments

Level: beginner

See also: `Vec`, `VecAssemblyBegin()`, `VecSetValues()`, `VecViewFromOptions()`, `VecStashViewFromOptions()`,
`PetscObjectViewFromOptions()`

# External Links
$(_doc_external("Vec/VecAssemblyEnd"))
"""
function VecAssemblyEnd(petsclib::PetscLibType, vec::AbstractPetscVec)
    error("VecAssemblyEnd: no generated method for these argument types")
end

@for_petsc function VecAssemblyEnd(petsclib::$UnionPetscLib, vec::AbstractPetscVec )

    @chk ccall(
               (:VecAssemblyEnd, $petsc_library),
               PetscErrorCode,
               (CVec,),
               vec,
              )


	return nothing
end 

"""
	VecBindToCPU(petsclib::PetscLibType, v::AbstractPetscVec, flg::PetscBool) 
marks a vector to temporarily stay on the CPU and perform computations on the CPU

Logically collective

Input Parameters:
- `v`   - the vector
- `flg` - bind to the CPU if value of `PETSC_TRUE`

Level: intermediate

See also: `Vec`, `VecBoundToCPU()`

# External Links
$(_doc_external("Vec/VecBindToCPU"))
"""
function VecBindToCPU(petsclib::PetscLibType, v::AbstractPetscVec, flg::PetscBool)
    error("VecBindToCPU: no generated method for these argument types")
end

@for_petsc function VecBindToCPU(petsclib::$UnionPetscLib, v::AbstractPetscVec, flg::PetscBool )

    @chk ccall(
               (:VecBindToCPU, $petsc_library),
               PetscErrorCode,
               (CVec, PetscBool),
               v, flg,
              )


	return nothing
end 

"""
	VecBoundGradientProjection(petsclib::PetscLibType, G::AbstractPetscVec, X::AbstractPetscVec, XL::AbstractPetscVec, XU::AbstractPetscVec, GP::AbstractPetscVec) 
Projects vector according to this definition.
If XL[i] < X[i] < XU[i], then GP[i] = G[i];
If X[i] <= XL[i], then GP[i] = min(G[i],0);
If X[i] >= XU[i], then GP[i] = max(G[i],0);

Input Parameters:
- `G`  - current gradient vector
- `X`  - current solution vector with XL[i] <= X[i] <= XU[i]
- `XL` - lower bounds
- `XU` - upper bounds

Output Parameter:
- `GP` - gradient projection vector

Level: advanced

See also: `Vec`

# External Links
$(_doc_external("Vec/VecBoundGradientProjection"))
"""
function VecBoundGradientProjection(petsclib::PetscLibType, G::AbstractPetscVec, X::AbstractPetscVec, XL::AbstractPetscVec, XU::AbstractPetscVec, GP::AbstractPetscVec)
    error("VecBoundGradientProjection: no generated method for these argument types")
end

@for_petsc function VecBoundGradientProjection(petsclib::$UnionPetscLib, G::AbstractPetscVec, X::AbstractPetscVec, XL::AbstractPetscVec, XU::AbstractPetscVec, GP::AbstractPetscVec )

    @chk ccall(
               (:VecBoundGradientProjection, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, CVec, CVec, CVec),
               G, X, XL, XU, GP,
              )


	return nothing
end 

"""
	flg::PetscBool = VecBoundToCPU(petsclib::PetscLibType, v::AbstractPetscVec) 
query if a vector is bound to the CPU

Not collective

Input Parameter:
- `v` - the vector

Output Parameter:
- `flg` - the logical flag

Level: intermediate

See also: `Vec`, `VecBindToCPU()`

# External Links
$(_doc_external("Vec/VecBoundToCPU"))
"""
function VecBoundToCPU(petsclib::PetscLibType, v::AbstractPetscVec)
    error("VecBoundToCPU: no generated method for these argument types")
end

@for_petsc function VecBoundToCPU(petsclib::$UnionPetscLib, v::AbstractPetscVec )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:VecBoundToCPU, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{PetscBool}),
               v, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	Y::PetscVec,x_is::Vector{IS} = VecConcatenate(petsclib::PetscLibType, nx::PetscInt, X::Vector{<:AbstractPetscVec}) 
Creates a new vector that is a vertical concatenation of all the given array of vectors
in the order they appear in the array. The concatenated vector resides on the same
communicator and is the same type as the source vectors.

Collective

Input Parameters:
- `nx` - number of vectors to be concatenated
- `X`  - array containing the vectors to be concatenated in the order of concatenation

Output Parameters:
- `Y`    - concatenated vector
- `x_is` - array of index sets corresponding to the concatenated components of `Y` (pass `NULL` if not needed)

Level: advanced

See also: `Vec`, `VECNEST`, `VECSCATTER`, `VecScatterCreate()`

# External Links
$(_doc_external("Vec/VecConcatenate"))
"""
function VecConcatenate(petsclib::PetscLibType, nx::Integer, X::Vector{<:AbstractPetscVec})
    error("VecConcatenate: no generated method for these argument types")
end

@for_petsc function VecConcatenate(petsclib::$UnionPetscLib, nx::$PetscInt, X::Vector{<:AbstractPetscVec} )
	Y_ = Ref{CVec}()
	x_is_ = Ref{Ptr{CIS}}()

    @chk ccall(
               (:VecConcatenate, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{CVec}, Ptr{CVec}, Ptr{Ptr{CIS}}),
               nx, X, Y_, x_is_,
              )

	Y = PetscVec(Y_[], petsclib)
	x_is = x_is_[] == C_NULL ? IS{$PetscLib}[] : [IS(p, petsclib) for p in unsafe_wrap(Array, x_is_[], nx; own = false)]

	return Y,x_is
end 

"""
	VecConjugate(petsclib::PetscLibType, x::AbstractPetscVec) 
Conjugates a vector. That is, replace every entry in a vector with its complex conjugate

Logically Collective

Input Parameter:
- `x` - the vector

Level: intermediate

See also: `Vec`, `VecSet()`

# External Links
$(_doc_external("Vec/VecConjugate"))
"""
function VecConjugate(petsclib::PetscLibType, x::AbstractPetscVec)
    error("VecConjugate: no generated method for these argument types")
end

@for_petsc function VecConjugate(petsclib::$UnionPetscLib, x::AbstractPetscVec )

    @chk ccall(
               (:VecConjugate, $petsc_library),
               PetscErrorCode,
               (CVec,),
               x,
              )


	return nothing
end 

"""
	VecCopy(petsclib::PetscLibType, x::AbstractPetscVec, y::AbstractPetscVec) 
Copies a vector `y = x`

Logically Collective

Input Parameter:
- `x` - the vector

Output Parameter:
- `y` - the copy

Level: beginner

See also: `Vec`, `VecDuplicate()`

# External Links
$(_doc_external("Vec/VecCopy"))
"""
function VecCopy(petsclib::PetscLibType, x::AbstractPetscVec, y::AbstractPetscVec)
    error("VecCopy: no generated method for these argument types")
end

@for_petsc function VecCopy(petsclib::$UnionPetscLib, x::AbstractPetscVec, y::AbstractPetscVec )

    @chk ccall(
               (:VecCopy, $petsc_library),
               PetscErrorCode,
               (CVec, CVec),
               x, y,
              )


	return nothing
end 

"""
	vec::PetscVec = VecCreate(petsclib::PetscLibType, comm::MPI_Comm) 
Creates an empty vector object. The type can then be set with `VecSetType()`,
or `VecSetFromOptions().`

Collective

Input Parameter:
- `comm` - The communicator for the vector object

Output Parameter:
- `vec` - The vector object

Level: beginner

See also: `Vec`, `VecSetType()`, `VecSetSizes()`, `VecCreateMPIWithArray()`, `VecCreateMPI()`, `VecDuplicate()`,
`VecDuplicateVecs()`, `VecCreateGhost()`, `VecCreateSeq()`, `VecPlaceArray()`

# External Links
$(_doc_external("Vec/VecCreate"))
"""
function VecCreate(petsclib::PetscLibType, comm::MPI_Comm)
    error("VecCreate: no generated method for these argument types")
end

@for_petsc function VecCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	vec_ = Ref{CVec}()

    @chk ccall(
               (:VecCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{CVec}),
               comm, vec_,
              )

	vec = PetscVec(vec_[], petsclib)

	return vec
end 

"""
	vec::PetscVec = VecCreateFromOptions(petsclib::PetscLibType, comm::MPI_Comm, prefix::String, bs::PetscInt, m::PetscInt, n::PetscInt) 
Creates a vector whose type is set from the options database

Collective

Input Parameters:
- `comm`   - The communicator for the vector object
- `prefix` - [optional] prefix for the options database
- `bs`     - the block size (commonly 1)
- `m`      - the local size (or `PETSC_DECIDE`)
- `n`      - the global size (or `PETSC_DETERMINE`)

Output Parameter:
- `vec` - The vector object

Options Database Keys:
- `-vec_type` - see `VecType`, for example `seq`, `mpi`, `cuda`, defaults to `mpi`

Level: beginner

See also: `Vec`, `VecSetType()`, `VecSetSizes()`, `VecCreateMPIWithArray()`, `VecCreateMPI()`, `VecDuplicate()`,
`VecDuplicateVecs()`, `VecCreateGhost()`, `VecCreateSeq()`, `VecPlaceArray()`, `VecCreate()`, `VecType`

# External Links
$(_doc_external("Vec/VecCreateFromOptions"))
"""
function VecCreateFromOptions(petsclib::PetscLibType, comm::MPI_Comm, prefix::String, bs::Integer, m::Integer, n::Integer)
    error("VecCreateFromOptions: no generated method for these argument types")
end

@for_petsc function VecCreateFromOptions(petsclib::$UnionPetscLib, comm::MPI_Comm, prefix::String, bs::$PetscInt, m::$PetscInt, n::$PetscInt )
	vec_ = Ref{CVec}()

    @chk ccall(
               (:VecCreateFromOptions, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, $PetscInt, $PetscInt, $PetscInt, Ptr{CVec}),
               comm, prefix, bs, m, n, vec_,
              )

	vec = PetscVec(vec_[], petsclib)

	return vec
end 

"""
	vv::PetscVec = VecCreateGhost(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt, M_N::PetscInt, nghost::PetscInt, ghosts::Vector{PetscInt}) 
Creates a parallel vector with ghost padding on each processor.

Collective

Input Parameters:
- `comm`   - the MPI communicator to use
- `n`      - local vector length
- `N`      - global vector length (or `PETSC_DETERMINE` to have calculated if `n` is given)
- `nghost` - number of local ghost points
- `ghosts` - global indices of ghost points, these do not need to be in increasing order (sorted)

Output Parameter:
- `vv` - the global vector representation (without ghost points as part of vector)

Level: advanced

See also: `Vec`, `VecType`, `VecCreateSeq()`, `VecCreate()`, `VecDuplicate()`, `VecDuplicateVecs()`, `VecCreateMPI()`,
`VecGhostGetLocalForm()`, `VecGhostRestoreLocalForm()`, `VecGhostUpdateBegin()`,
`VecCreateGhostWithArray()`, `VecCreateMPIWithArray()`, `VecGhostUpdateEnd()`,
`VecCreateGhostBlock()`, `VecCreateGhostBlockWithArray()`, `VecMPISetGhost()`

# External Links
$(_doc_external("Vec/VecCreateGhost"))
"""
function VecCreateGhost(petsclib::PetscLibType, comm::MPI_Comm, n::Integer, M_N::Integer, nghost::Integer, ghosts::AbstractVector{<:Number})
    error("VecCreateGhost: no generated method for these argument types")
end

@for_petsc function VecCreateGhost(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, M_N::$PetscInt, nghost::$PetscInt, ghosts::Vector{$PetscInt} )
	vv_ = Ref{CVec}()

    @chk ccall(
               (:VecCreateGhost, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{CVec}),
               comm, n, M_N, nghost, ghosts, vv_,
              )

	vv = PetscVec(vv_[], petsclib)

	return vv
end 

"""
	vv::PetscVec = VecCreateGhostBlock(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, n::PetscInt, M_N::PetscInt, nghost::PetscInt, ghosts::Vector{PetscInt}) 
Creates a parallel vector with ghost padding on each processor.
The indicing of the ghost points is done with blocks.

Collective

Input Parameters:
- `comm`   - the MPI communicator to use
- `bs`     - the block size
- `n`      - local vector length
- `N`      - global vector length (or `PETSC_DETERMINE` to have calculated if `n` is given)
- `nghost` - number of local ghost blocks
- `ghosts` - global indices of ghost blocks, counts are by block, not by individual index, these do not need to be in increasing order (sorted)

Output Parameter:
- `vv` - the global vector representation (without ghost points as part of vector)

Level: advanced

See also: `Vec`, `VecType`, `VecCreateSeq()`, `VecCreate()`, `VecDuplicate()`, `VecDuplicateVecs()`, `VecCreateMPI()`,
`VecGhostGetLocalForm()`, `VecGhostRestoreLocalForm()`, `VecGhostUpdateBegin()`, `VecGhostUpdateEnd()`,
`VecCreateGhostWithArray()`, `VecCreateMPIWithArray()`, `VecCreateGhostBlockWithArray()`

# External Links
$(_doc_external("Vec/VecCreateGhostBlock"))
"""
function VecCreateGhostBlock(petsclib::PetscLibType, comm::MPI_Comm, bs::Integer, n::Integer, M_N::Integer, nghost::Integer, ghosts::AbstractVector{<:Number})
    error("VecCreateGhostBlock: no generated method for these argument types")
end

@for_petsc function VecCreateGhostBlock(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, n::$PetscInt, M_N::$PetscInt, nghost::$PetscInt, ghosts::Vector{$PetscInt} )
	vv_ = Ref{CVec}()

    @chk ccall(
               (:VecCreateGhostBlock, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{CVec}),
               comm, bs, n, M_N, nghost, ghosts, vv_,
              )

	vv = PetscVec(vv_[], petsclib)

	return vv
end 

"""
	vv::PetscVec = VecCreateGhostBlockWithArray(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, n::PetscInt, M_N::PetscInt, nghost::PetscInt, ghosts::Vector{PetscInt}, array::Vector{PetscScalar}) 
Creates a parallel vector with ghost padding on each processor;
the caller allocates the array space. Indices in the ghost region are based on blocks.

Collective

Input Parameters:
- `comm`   - the MPI communicator to use
- `bs`     - block size
- `n`      - local vector length
- `N`      - global vector length (or `PETSC_DETERMINE` to have calculated if `n` is given)
- `nghost` - number of local ghost blocks
- `ghosts` - global indices of ghost blocks (or `NULL` if not needed), counts are by block not by index, these do not need to be in increasing order (sorted)
- `array`  - the space to store the vector values (as long as n + nghost*bs)

Output Parameter:
- `vv` - the global vector representation (without ghost points as part of vector)

Level: advanced

See also: `Vec`, `VecType`, `VecCreate()`, `VecGhostGetLocalForm()`, `VecGhostRestoreLocalForm()`,
`VecCreateGhost()`, `VecCreateSeqWithArray()`, `VecCreateMPIWithArray()`,
`VecCreateGhostWithArray()`, `VecCreateGhostBlock()`, `VecGhostUpdateBegin()`, `VecGhostUpdateEnd()`

# External Links
$(_doc_external("Vec/VecCreateGhostBlockWithArray"))
"""
function VecCreateGhostBlockWithArray(petsclib::PetscLibType, comm::MPI_Comm, bs::Integer, n::Integer, M_N::Integer, nghost::Integer, ghosts::AbstractVector{<:Number}, array::AbstractVector{<:Number})
    error("VecCreateGhostBlockWithArray: no generated method for these argument types")
end

@for_petsc function VecCreateGhostBlockWithArray(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, n::$PetscInt, M_N::$PetscInt, nghost::$PetscInt, ghosts::Vector{$PetscInt}, array::Vector{$PetscScalar} )
	vv_ = Ref{CVec}()

    @chk ccall(
               (:VecCreateGhostBlockWithArray, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}, Ptr{CVec}),
               comm, bs, n, M_N, nghost, ghosts, array, vv_,
              )

	vv = PetscVec(vv_[], petsclib)

	return vv
end 

"""
	vv::PetscVec = VecCreateGhostWithArray(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt, M_N::PetscInt, nghost::PetscInt, ghosts::Vector{PetscInt}, array::Vector{PetscScalar}) 
Creates a parallel vector with ghost padding on each processor;
the caller allocates the array space.

Collective

Input Parameters:
- `comm`   - the MPI communicator to use
- `n`      - local vector length
- `N`      - global vector length (or `PETSC_DETERMINE` to have calculated if `n` is given)
- `nghost` - number of local ghost points
- `ghosts` - global indices of ghost points (or `NULL` if not needed), these do not need to be in increasing order (sorted)
- `array`  - the space to store the vector values (as long as n + nghost)

Output Parameter:
- `vv` - the global vector representation (without ghost points as part of vector)

Level: advanced

See also: `Vec`, `VecType`, `VecCreate()`, `VecGhostGetLocalForm()`, `VecGhostRestoreLocalForm()`,
`VecCreateGhost()`, `VecCreateSeqWithArray()`, `VecCreateMPIWithArray()`,
`VecCreateGhostBlock()`, `VecCreateGhostBlockWithArray()`, `VecMPISetGhost()`, `VecGhostUpdateBegin()`, `VecGhostUpdateEnd()`

# External Links
$(_doc_external("Vec/VecCreateGhostWithArray"))
"""
function VecCreateGhostWithArray(petsclib::PetscLibType, comm::MPI_Comm, n::Integer, M_N::Integer, nghost::Integer, ghosts::AbstractVector{<:Number}, array::AbstractVector{<:Number})
    error("VecCreateGhostWithArray: no generated method for these argument types")
end

@for_petsc function VecCreateGhostWithArray(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, M_N::$PetscInt, nghost::$PetscInt, ghosts::Vector{$PetscInt}, array::Vector{$PetscScalar} )
	vv_ = Ref{CVec}()

    @chk ccall(
               (:VecCreateGhostWithArray, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}, Ptr{CVec}),
               comm, n, M_N, nghost, ghosts, array, vv_,
              )

	vv = PetscVec(vv_[], petsclib)

	return vv
end 

"""
	w::PetscVec = VecCreateLocalVector(petsclib::PetscLibType, v::AbstractPetscVec) 
Creates a vector object suitable for use with `VecGetLocalVector()` and friends. You must call `VecDestroy()` when the
vector is no longer needed.

Not Collective.

Input Parameter:
- `v` - The vector for which the local vector is desired.

Output Parameter:
- `w` - Upon exit this contains the local vector.

Level: beginner

See also: `Vec`, `VecGetLocalVectorRead()`, `VecRestoreLocalVectorRead()`, `VecGetLocalVector()`, `VecRestoreLocalVector()`

# External Links
$(_doc_external("Vec/VecCreateLocalVector"))
"""
function VecCreateLocalVector(petsclib::PetscLibType, v::AbstractPetscVec)
    error("VecCreateLocalVector: no generated method for these argument types")
end

@for_petsc function VecCreateLocalVector(petsclib::$UnionPetscLib, v::AbstractPetscVec )
	w_ = Ref{CVec}()

    @chk ccall(
               (:VecCreateLocalVector, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{CVec}),
               v, w_,
              )

	w = PetscVec(w_[], petsclib)

	return w
end 

"""
	v::PetscVec = VecCreateMPI(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt, M_N::PetscInt) 
Creates a parallel vector.

Collective

Input Parameters:
- `comm` - the MPI communicator to use
- `n`    - local vector length (or `PETSC_DECIDE` to have calculated if `N` is given)
- `N`    - global vector length (or `PETSC_DETERMINE` to have calculated if `n` is given)

Output Parameter:
- `v` - the vector

Level: intermediate

See also: `Vec`, `VecType`, `VecCreateSeq()`, `VecCreate()`, `VecDuplicate()`, `VecDuplicateVecs()`, `VecCreateGhost()`,
`VecCreateMPIWithArray()`, `VecCreateGhostWithArray()`, `VecMPISetGhost()`, `PetscLayout`,
`VecGetOwnershipRange()`, `VecGetOwnershipRanges()`

# External Links
$(_doc_external("Vec/VecCreateMPI"))
"""
function VecCreateMPI(petsclib::PetscLibType, comm::MPI_Comm, n::Integer, M_N::Integer)
    error("VecCreateMPI: no generated method for these argument types")
end

@for_petsc function VecCreateMPI(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, M_N::$PetscInt )
	v_ = Ref{CVec}()

    @chk ccall(
               (:VecCreateMPI, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, Ptr{CVec}),
               comm, n, M_N, v_,
              )

	v = PetscVec(v_[], petsclib)

	return v
end 

"""
	v::PetscVec = VecCreateMPIKokkosWithArray(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, n::PetscInt, M_N::PetscInt, darray::Vector{PetscScalar}) 

# External Links
$(_doc_external("Vec/VecCreateMPIKokkosWithArray"))
"""
function VecCreateMPIKokkosWithArray(petsclib::PetscLibType, comm::MPI_Comm, bs::Integer, n::Integer, M_N::Integer, darray::AbstractVector{<:Number})
    error("VecCreateMPIKokkosWithArray: no generated method for these argument types")
end

@for_petsc function VecCreateMPIKokkosWithArray(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, n::$PetscInt, M_N::$PetscInt, darray::Vector{$PetscScalar} )
	v_ = Ref{CVec}()

    @chk ccall(
               (:VecCreateMPIKokkosWithArray, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscScalar}, Ptr{CVec}),
               comm, bs, n, M_N, darray, v_,
              )

	v = PetscVec(v_[], petsclib)

	return v
end 

"""
	array::ViennaCLVector,vv::PetscVec = VecCreateMPIViennaCLWithArray(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, n::PetscInt, M_N::PetscInt) 

# External Links
$(_doc_external("Vec/VecCreateMPIViennaCLWithArray"))
"""
function VecCreateMPIViennaCLWithArray(petsclib::PetscLibType, comm::MPI_Comm, bs::Integer, n::Integer, M_N::Integer)
    error("VecCreateMPIViennaCLWithArray: no generated method for these argument types")
end

@for_petsc function VecCreateMPIViennaCLWithArray(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, n::$PetscInt, M_N::$PetscInt )
	array_ = Ref{ViennaCLVector}()
	vv_ = Ref{CVec}()

    @chk ccall(
               (:VecCreateMPIViennaCLWithArray, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{ViennaCLVector}, Ptr{CVec}),
               comm, bs, n, M_N, array_, vv_,
              )

	array = array_[]
	vv = PetscVec(vv_[], petsclib)

	return array,vv
end 

"""
	viennaclvec::ViennaCLVector,vv::PetscVec = VecCreateMPIViennaCLWithArrays(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, n::PetscInt, M_N::PetscInt, cpuarray::Vector{PetscScalar}) 

# External Links
$(_doc_external("Vec/VecCreateMPIViennaCLWithArrays"))
"""
function VecCreateMPIViennaCLWithArrays(petsclib::PetscLibType, comm::MPI_Comm, bs::Integer, n::Integer, M_N::Integer, cpuarray::AbstractVector{<:Number})
    error("VecCreateMPIViennaCLWithArrays: no generated method for these argument types")
end

@for_petsc function VecCreateMPIViennaCLWithArrays(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, n::$PetscInt, M_N::$PetscInt, cpuarray::Vector{$PetscScalar} )
	viennaclvec_ = Ref{ViennaCLVector}()
	vv_ = Ref{CVec}()

    @chk ccall(
               (:VecCreateMPIViennaCLWithArrays, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscScalar}, Ptr{ViennaCLVector}, Ptr{CVec}),
               comm, bs, n, M_N, cpuarray, viennaclvec_, vv_,
              )

	viennaclvec = viennaclvec_[]
	vv = PetscVec(vv_[], petsclib)

	return viennaclvec,vv
end 

"""
	vv::PetscVec = VecCreateMPIWithArray(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, n::PetscInt, M_N::PetscInt, array::Vector{PetscScalar}) 
Creates a parallel, array-style vector,
where the user provides the array space to store the vector values.

Collective

Input Parameters:
- `comm`  - the MPI communicator to use
- `bs`    - block size, same meaning as `VecSetBlockSize()`
- `n`     - local vector length, cannot be `PETSC_DECIDE`
- `N`     - global vector length (or `PETSC_DETERMINE` to have calculated)
- `array` - the user provided array to store the vector values

Output Parameter:
- `vv` - the vector

Level: intermediate

See also: `Vec`, `VecType`, `VecCreateSeqWithArray()`, `VecCreate()`, `VecDuplicate()`, `VecDuplicateVecs()`, `VecCreateGhost()`,
`VecCreateMPI()`, `VecCreateGhostWithArray()`, `VecPlaceArray()`

# External Links
$(_doc_external("Vec/VecCreateMPIWithArray"))
"""
function VecCreateMPIWithArray(petsclib::PetscLibType, comm::MPI_Comm, bs::Integer, n::Integer, M_N::Integer, array::AbstractVector{<:Number})
    error("VecCreateMPIWithArray: no generated method for these argument types")
end

@for_petsc function VecCreateMPIWithArray(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, n::$PetscInt, M_N::$PetscInt, array::Vector{$PetscScalar} )
	vv_ = Ref{CVec}()

    @chk ccall(
               (:VecCreateMPIWithArray, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscScalar}, Ptr{CVec}),
               comm, bs, n, M_N, array, vv_,
              )

	vv = PetscVec(vv_[], petsclib)

	return vv
end 

"""
	Y::PetscVec = VecCreateNest(petsclib::PetscLibType, comm::MPI_Comm, nb::PetscInt, is::Vector{<:AbstractIS}, x::Vector{<:AbstractPetscVec}) 
Creates a new vector containing several nested subvectors, each stored separately

Collective

Input Parameters:
- `comm` - Communicator for the new `Vec`
- `nb`   - number of nested blocks
- `is`   - array of `nb` index sets describing each nested block, or `NULL` to pack subvectors contiguously
- `x`    - array of `nb` sub-vectors

Output Parameter:
- `Y` - new vector

Level: advanced

See also: `VECNEST`, `Vec`, `VecType`, `VecCreate()`, `MatCreateNest()`, `DMSetVecType()`

# External Links
$(_doc_external("Vec/VecCreateNest"))
"""
function VecCreateNest(petsclib::PetscLibType, comm::MPI_Comm, nb::Integer, is::Vector{<:AbstractIS}, x::Vector{<:AbstractPetscVec})
    error("VecCreateNest: no generated method for these argument types")
end

@for_petsc function VecCreateNest(petsclib::$UnionPetscLib, comm::MPI_Comm, nb::$PetscInt, is::Vector{<:AbstractIS}, x::Vector{<:AbstractPetscVec} )
	Y_ = Ref{CVec}()

    @chk ccall(
               (:VecCreateNest, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, Ptr{CIS}, Ptr{CVec}, Ptr{CVec}),
               comm, nb, is, x, Y_,
              )

	Y = PetscVec(Y_[], petsclib)

	return Y
end 

"""
	v::PetscVec = VecCreateSeq(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt) 
Creates a standard, sequential array-style vector.

Collective

Input Parameters:
- `comm` - the communicator, should be `PETSC_COMM_SELF`
- `n`    - the vector length

Output Parameter:
- `v` - the vector

Level: intermediate

See also: `Vec`, `VecType`, `VecCreateMPI()`, `VecCreate()`, `VecDuplicate()`, `VecDuplicateVecs()`, `VecCreateGhost()`

# External Links
$(_doc_external("Vec/VecCreateSeq"))
"""
function VecCreateSeq(petsclib::PetscLibType, comm::MPI_Comm, n::Integer)
    error("VecCreateSeq: no generated method for these argument types")
end

@for_petsc function VecCreateSeq(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt )
	v_ = Ref{CVec}()

    @chk ccall(
               (:VecCreateSeq, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, Ptr{CVec}),
               comm, n, v_,
              )

	v = PetscVec(v_[], petsclib)

	return v
end 

"""
	v::PetscVec = VecCreateSeqKokkos(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt) 

# External Links
$(_doc_external("Vec/VecCreateSeqKokkos"))
"""
function VecCreateSeqKokkos(petsclib::PetscLibType, comm::MPI_Comm, n::Integer)
    error("VecCreateSeqKokkos: no generated method for these argument types")
end

@for_petsc function VecCreateSeqKokkos(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt )
	v_ = Ref{CVec}()

    @chk ccall(
               (:VecCreateSeqKokkos, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, Ptr{CVec}),
               comm, n, v_,
              )

	v = PetscVec(v_[], petsclib)

	return v
end 

"""
	v::PetscVec = VecCreateSeqKokkosWithArray(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, n::PetscInt, darray::Vector{PetscScalar}) 

# External Links
$(_doc_external("Vec/VecCreateSeqKokkosWithArray"))
"""
function VecCreateSeqKokkosWithArray(petsclib::PetscLibType, comm::MPI_Comm, bs::Integer, n::Integer, darray::AbstractVector{<:Number})
    error("VecCreateSeqKokkosWithArray: no generated method for these argument types")
end

@for_petsc function VecCreateSeqKokkosWithArray(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, n::$PetscInt, darray::Vector{$PetscScalar} )
	v_ = Ref{CVec}()

    @chk ccall(
               (:VecCreateSeqKokkosWithArray, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, Ptr{$PetscScalar}, Ptr{CVec}),
               comm, bs, n, darray, v_,
              )

	v = PetscVec(v_[], petsclib)

	return v
end 

"""
	v::PetscVec = VecCreateSeqViennaCL(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt) 

# External Links
$(_doc_external("Vec/VecCreateSeqViennaCL"))
"""
function VecCreateSeqViennaCL(petsclib::PetscLibType, comm::MPI_Comm, n::Integer)
    error("VecCreateSeqViennaCL: no generated method for these argument types")
end

@for_petsc function VecCreateSeqViennaCL(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt )
	v_ = Ref{CVec}()

    @chk ccall(
               (:VecCreateSeqViennaCL, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, Ptr{CVec}),
               comm, n, v_,
              )

	v = PetscVec(v_[], petsclib)

	return v
end 

"""
	viennaclvec::ViennaCLVector,V::PetscVec = VecCreateSeqViennaCLWithArrays(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, n::PetscInt, cpuarray::Vector{PetscScalar}) 

# External Links
$(_doc_external("Vec/VecCreateSeqViennaCLWithArrays"))
"""
function VecCreateSeqViennaCLWithArrays(petsclib::PetscLibType, comm::MPI_Comm, bs::Integer, n::Integer, cpuarray::AbstractVector{<:Number})
    error("VecCreateSeqViennaCLWithArrays: no generated method for these argument types")
end

@for_petsc function VecCreateSeqViennaCLWithArrays(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, n::$PetscInt, cpuarray::Vector{$PetscScalar} )
	viennaclvec_ = Ref{ViennaCLVector}()
	V_ = Ref{CVec}()

    @chk ccall(
               (:VecCreateSeqViennaCLWithArrays, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, Ptr{$PetscScalar}, Ptr{ViennaCLVector}, Ptr{CVec}),
               comm, bs, n, cpuarray, viennaclvec_, V_,
              )

	viennaclvec = viennaclvec_[]
	V = PetscVec(V_[], petsclib)

	return viennaclvec,V
end 

"""
	V::PetscVec = VecCreateSeqWithArray(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, n::PetscInt, array::Vector{PetscScalar}) 
Creates a standard,sequential array-style vector,
where the user provides the array space to store the vector values.

Collective

Input Parameters:
- `comm`  - the communicator, should be `PETSC_COMM_SELF`
- `bs`    - the block size
- `n`     - the vector length
- `array` - memory where the vector elements are to be stored.

Output Parameter:
- `V` - the vector

Level: intermediate

See also: `VecCreateMPIWithArray()`, `VecCreate()`, `VecDuplicate()`, `VecDuplicateVecs()`,
`VecCreateGhost()`, `VecCreateSeq()`, `VecPlaceArray()`

# External Links
$(_doc_external("Vec/VecCreateSeqWithArray"))
"""
function VecCreateSeqWithArray(petsclib::PetscLibType, comm::MPI_Comm, bs::Integer, n::Integer, array::AbstractVector{<:Number})
    error("VecCreateSeqWithArray: no generated method for these argument types")
end

@for_petsc function VecCreateSeqWithArray(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, n::$PetscInt, array::Vector{$PetscScalar} )
	V_ = Ref{CVec}()

    @chk ccall(
               (:VecCreateSeqWithArray, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, Ptr{$PetscScalar}, Ptr{CVec}),
               comm, bs, n, array, V_,
              )

	V = PetscVec(V_[], petsclib)

	return V
end 

"""
	v::PetscVec = VecCreateShared(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt, M_N::PetscInt) 
Creates a parallel vector that uses shared memory.

Collective

Input Parameters:
- `comm` - the MPI communicator to use
- `n`    - local vector length (or `PETSC_DECIDE` to have calculated if `N` is given)
- `N`    - global vector length (or `PETSC_DECIDE` to have calculated if `n` is given)

Output Parameter:
- `v` - the vector

Level: advanced

See also: `Vec`, `VecType`, `VecCreateSeq()`, `VecCreate()`, `VecCreateMPI()`, `VecDuplicate()`, `VecDuplicateVecs()`,
`VecCreateGhost()`, `VecCreateMPIWithArray()`, `VecCreateGhostWithArray()`

# External Links
$(_doc_external("Vec/VecCreateShared"))
"""
function VecCreateShared(petsclib::PetscLibType, comm::MPI_Comm, n::Integer, M_N::Integer)
    error("VecCreateShared: no generated method for these argument types")
end

@for_petsc function VecCreateShared(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, M_N::$PetscInt )
	v_ = Ref{CVec}()

    @chk ccall(
               (:VecCreateShared, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, Ptr{CVec}),
               comm, n, M_N, v_,
              )

	v = PetscVec(v_[], petsclib)

	return v
end 

"""
	VecDestroy(petsclib::PetscLibType, v::AbstractPetscVec) 
Destroys a vector.

Collective

Input Parameter:
- `v` - the vector

Level: beginner

See also: `Vec`, `VecCreate()`, `VecDuplicate()`, `VecDestroyVecs()`

# External Links
$(_doc_external("Vec/VecDestroy"))
"""
function VecDestroy(petsclib::PetscLibType, v::AbstractPetscVec)
    error("VecDestroy: no generated method for these argument types")
end

@for_petsc function VecDestroy(petsclib::$UnionPetscLib, v::AbstractPetscVec )
	v_ = Ref(v.ptr)

    @chk ccall(
               (:VecDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{CVec},),
               v_,
              )

	v.ptr = C_NULL

	return nothing
end 

"""
	VecDestroyVecs(petsclib::PetscLibType, m::PetscInt, vv::Union{Ptr, AbstractArray{PetscVec}}) 
Frees a block of vectors obtained with `VecDuplicateVecs()`.

Collective

Input Parameters:
- `m`  - the number of vectors previously obtained, if zero no vectors are destroyed
- `vv` - pointer to pointer to array of vector pointers, if `NULL` no vectors are destroyed

Level: intermediate

See also: `Vec`, `VecDuplicateVecs()`, `VecDestroyVecsf90()`

# External Links
$(_doc_external("Vec/VecDestroyVecs"))
"""
function VecDestroyVecs(petsclib::PetscLibType, m::Integer, vv::Union{Ptr, AbstractArray{PetscVec}})
    error("VecDestroyVecs: no generated method for these argument types")
end

@for_petsc function VecDestroyVecs(petsclib::$UnionPetscLib, m::$PetscInt, vv::Union{Ptr, AbstractArray{PetscVec}} )
	vv_ = Ref{Ptr{CVec}}(vv isa Ptr ? vv : pointer(vv))

    @chk ccall(
               (:VecDestroyVecs, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{Ptr{CVec}}),
               m, vv_,
              )


	return nothing
end 

"""
	val::PetscScalar = VecDot(petsclib::PetscLibType, x::AbstractPetscVec, y::AbstractPetscVec) 
Computes the vector dot product.

Collective

Input Parameters:
- `x` - first vector
- `y` - second vector

Output Parameter:
- `val` - the dot product

Level: intermediate

See also: `Vec`, `VecMDot()`, `VecTDot()`, `VecNorm()`, `VecDotBegin()`, `VecDotEnd()`, `VecDotRealPart()`

# External Links
$(_doc_external("Vec/VecDot"))
"""
function VecDot(petsclib::PetscLibType, x::AbstractPetscVec, y::AbstractPetscVec)
    error("VecDot: no generated method for these argument types")
end

@for_petsc function VecDot(petsclib::$UnionPetscLib, x::AbstractPetscVec, y::AbstractPetscVec )
	val_ = Ref{$PetscScalar}()

    @chk ccall(
               (:VecDot, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, Ptr{$PetscScalar}),
               x, y, val_,
              )

	val = val_[]

	return val
end 

"""
	result::PetscScalar = VecDotBegin(petsclib::PetscLibType, x::AbstractPetscVec, y::AbstractPetscVec) 
Starts a split phase dot product computation.

Input Parameters:
- `x`      - the first vector
- `y`      - the second vector
- `result` - where the result will go (can be `NULL`)

Level: advanced

See also: `VecDotEnd()`, `VecNormBegin()`, `VecNormEnd()`, `VecNorm()`, `VecDot()`, `VecMDot()`,
`VecTDotBegin()`, `VecTDotEnd()`, `PetscCommSplitReductionBegin()`

# External Links
$(_doc_external("Vec/VecDotBegin"))
"""
function VecDotBegin(petsclib::PetscLibType, x::AbstractPetscVec, y::AbstractPetscVec)
    error("VecDotBegin: no generated method for these argument types")
end

@for_petsc function VecDotBegin(petsclib::$UnionPetscLib, x::AbstractPetscVec, y::AbstractPetscVec )
	result_ = Ref{$PetscScalar}()

    @chk ccall(
               (:VecDotBegin, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, Ptr{$PetscScalar}),
               x, y, result_,
              )

	result = result_[]

	return result
end 

"""
	result::PetscScalar = VecDotEnd(petsclib::PetscLibType, x::AbstractPetscVec, y::AbstractPetscVec) 
Ends a split phase dot product computation.

Input Parameters:
- `x`      - the first vector (can be `NULL`)
- `y`      - the second vector (can be `NULL`)
- `result` - where the result will go

Level: advanced

See also: `VecDotBegin()`, `VecNormBegin()`, `VecNormEnd()`, `VecNorm()`, `VecDot()`, `VecMDot()`,
`VecTDotBegin()`, `VecTDotEnd()`, `PetscCommSplitReductionBegin()`

# External Links
$(_doc_external("Vec/VecDotEnd"))
"""
function VecDotEnd(petsclib::PetscLibType, x::AbstractPetscVec, y::AbstractPetscVec)
    error("VecDotEnd: no generated method for these argument types")
end

@for_petsc function VecDotEnd(petsclib::$UnionPetscLib, x::AbstractPetscVec, y::AbstractPetscVec )
	result_ = Ref{$PetscScalar}()

    @chk ccall(
               (:VecDotEnd, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, Ptr{$PetscScalar}),
               x, y, result_,
              )

	result = result_[]

	return result
end 

"""
	dp::PetscScalar,nm::PetscReal = VecDotNorm2(petsclib::PetscLibType, s::AbstractPetscVec, t::AbstractPetscVec) 
computes the inner product of two vectors and the 2-norm squared of the second vector

Collective

Input Parameters:
- `s` - first vector
- `t` - second vector

Output Parameters:
- `dp` - s'conj(t)
- `nm` - t'conj(t)

Level: advanced

See also: `Vec`, `VecDot()`, `VecNorm()`, `VecDotBegin()`, `VecNormBegin()`, `VecDotEnd()`, `VecNormEnd()`

# External Links
$(_doc_external("Vec/VecDotNorm2"))
"""
function VecDotNorm2(petsclib::PetscLibType, s::AbstractPetscVec, t::AbstractPetscVec)
    error("VecDotNorm2: no generated method for these argument types")
end

@for_petsc function VecDotNorm2(petsclib::$UnionPetscLib, s::AbstractPetscVec, t::AbstractPetscVec )
	dp_ = Ref{$PetscScalar}()
	nm_ = Ref{$PetscReal}()

    @chk ccall(
               (:VecDotNorm2, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, Ptr{$PetscScalar}, Ptr{$PetscReal}),
               s, t, dp_, nm_,
              )

	dp = dp_[]
	nm = nm_[]

	return dp,nm
end 

"""
	val::PetscReal = VecDotRealPart(petsclib::PetscLibType, x::AbstractPetscVec, y::AbstractPetscVec) 
Computes the real part of the vector dot product.

Collective

Input Parameters:
- `x` - first vector
- `y` - second vector

Output Parameter:
- `val` - the real part of the dot product;

Level: intermediate

Notes for Users of Complex Numbers:
See `VecDot()` for more details on the definition of the dot product for complex numbers

For real numbers this returns the same value as `VecDot()`

For complex numbers in C^n (that is a vector of n components with a complex number for each component) this is equal to the usual real dot product on the
the space R^{2n} (that is a vector of 2n components with the real or imaginary part of the complex numbers for components)

See also: `Vec`, `VecMDot()`, `VecTDot()`, `VecNorm()`, `VecDotBegin()`, `VecDotEnd()`, `VecDot()`, `VecDotNorm2()`

# External Links
$(_doc_external("Vec/VecDotRealPart"))
"""
function VecDotRealPart(petsclib::PetscLibType, x::AbstractPetscVec, y::AbstractPetscVec)
    error("VecDotRealPart: no generated method for these argument types")
end

@for_petsc function VecDotRealPart(petsclib::$UnionPetscLib, x::AbstractPetscVec, y::AbstractPetscVec )
	val_ = Ref{$PetscReal}()

    @chk ccall(
               (:VecDotRealPart, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, Ptr{$PetscReal}),
               x, y, val_,
              )

	val = val_[]

	return val
end 

"""
	newv::PetscVec = VecDuplicate(petsclib::PetscLibType, v::AbstractPetscVec) 
Creates a new vector of the same type as an existing vector.

Collective

Input Parameter:
- `v` - a vector to mimic

Output Parameter:
- `newv` - location to put new vector

Level: beginner

See also: `Vec`, `VecDestroy()`, `VecDuplicateVecs()`, `VecCreate()`, `VecCopy()`

# External Links
$(_doc_external("Vec/VecDuplicate"))
"""
function VecDuplicate(petsclib::PetscLibType, v::AbstractPetscVec)
    error("VecDuplicate: no generated method for these argument types")
end

@for_petsc function VecDuplicate(petsclib::$UnionPetscLib, v::AbstractPetscVec )
	newv_ = Ref{CVec}()

    @chk ccall(
               (:VecDuplicate, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{CVec}),
               v, newv_,
              )

	newv = PetscVec(newv_[], petsclib)

	return newv
end 

"""
	M_V::Ptr{CVec} = VecDuplicateVecs(petsclib::PetscLibType, v::AbstractPetscVec, m::PetscInt) 
Creates several vectors of the same type as an existing vector.

Collective

Input Parameters:
- `m` - the number of vectors to obtain
- `v` - a vector to mimic

Output Parameter:
- `V` - location to put pointer to array of vectors

Level: intermediate

See also: `Vec`, `VecDestroyVecs()`, `VecDuplicate()`, `VecCreate()`, `VecMDot()`, `VecMAXPY()`, `KSPGMRES`,
`KSPGMRESSetPreAllocateVectors()`

# External Links
$(_doc_external("Vec/VecDuplicateVecs"))
"""
function VecDuplicateVecs(petsclib::PetscLibType, v::AbstractPetscVec, m::Integer)
    error("VecDuplicateVecs: no generated method for these argument types")
end

@for_petsc function VecDuplicateVecs(petsclib::$UnionPetscLib, v::AbstractPetscVec, m::$PetscInt )
	M_V_ = Ref{Ptr{CVec}}()

    @chk ccall(
               (:VecDuplicateVecs, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{Ptr{CVec}}),
               v, m, M_V_,
              )

	M_V = M_V_[]

	return M_V
end 

"""
	flg::PetscBool = VecEqual(petsclib::PetscLibType, vec1::AbstractPetscVec, vec2::AbstractPetscVec) 
Compares two vectors. Returns true if the two vectors are either pointing to the same memory buffer,
or if the two vectors have the same local and global layout as well as bitwise equality of all entries.
Does NOT take round-off errors into account.

Collective

Input Parameters:
- `vec1` - the first vector
- `vec2` - the second vector

Output Parameter:
- `flg` - `PETSC_TRUE` if the vectors are equal; `PETSC_FALSE` otherwise.

Level: intermediate

See also: `Vec`

# External Links
$(_doc_external("Vec/VecEqual"))
"""
function VecEqual(petsclib::PetscLibType, vec1::AbstractPetscVec, vec2::AbstractPetscVec)
    error("VecEqual: no generated method for these argument types")
end

@for_petsc function VecEqual(petsclib::$UnionPetscLib, vec1::AbstractPetscVec, vec2::AbstractPetscVec )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:VecEqual, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, Ptr{PetscBool}),
               vec1, vec2, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	norm::PetscReal,norm_loc::PetscInt,norma::PetscReal,norma_loc::PetscInt,normr::PetscReal,normr_loc::PetscInt = VecErrorWeightedNorms(petsclib::PetscLibType, U::AbstractPetscVec, Y::AbstractPetscVec, E::AbstractPetscVec, wnormtype::NormType, atol::PetscReal, vatol::AbstractPetscVec, rtol::PetscReal, vrtol::AbstractPetscVec, ignore_max::PetscReal) 
compute a weighted norm of the difference between two vectors

Collective

Input Parameters:
- `U`          - first vector to be compared
- `Y`          - second vector to be compared
- `E`          - optional third vector representing the error (if not provided, the error is ||U-Y||)
- `wnormtype`  - norm type
- `atol`       - scalar for absolute tolerance
- `vatol`      - vector representing per-entry absolute tolerances (can be `NULL`)
- `rtol`       - scalar for relative tolerance
- `vrtol`      - vector representing per-entry relative tolerances (can be `NULL`)
- `ignore_max` - ignore values smaller than this value in absolute terms.

Output Parameters:
- `norm`      - weighted norm
- `norm_loc`  - number of vector locations used for the weighted norm
- `norma`     - weighted norm based on the absolute tolerance
- `norma_loc` - number of vector locations used for the absolute weighted norm
- `normr`     - weighted norm based on the relative tolerance
- `normr_loc` - number of vector locations used for the relative weighted norm

Level: developer

See also: `Vec`, `NormType`, `TSErrorWeightedNorm()`, `TSErrorWeightedENorm()`

# External Links
$(_doc_external("Vec/VecErrorWeightedNorms"))
"""
function VecErrorWeightedNorms(petsclib::PetscLibType, U::AbstractPetscVec, Y::AbstractPetscVec, E::AbstractPetscVec, wnormtype::NormType, atol::Real, vatol::AbstractPetscVec, rtol::Real, vrtol::AbstractPetscVec, ignore_max::Real)
    error("VecErrorWeightedNorms: no generated method for these argument types")
end

@for_petsc function VecErrorWeightedNorms(petsclib::$UnionPetscLib, U::AbstractPetscVec, Y::AbstractPetscVec, E::AbstractPetscVec, wnormtype::NormType, atol::$PetscReal, vatol::AbstractPetscVec, rtol::$PetscReal, vrtol::AbstractPetscVec, ignore_max::$PetscReal )
	norm_ = Ref{$PetscReal}()
	norm_loc_ = Ref{$PetscInt}()
	norma_ = Ref{$PetscReal}()
	norma_loc_ = Ref{$PetscInt}()
	normr_ = Ref{$PetscReal}()
	normr_loc_ = Ref{$PetscInt}()

    @chk ccall(
               (:VecErrorWeightedNorms, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, CVec, NormType, $PetscReal, CVec, $PetscReal, CVec, $PetscReal, Ptr{$PetscReal}, Ptr{$PetscInt}, Ptr{$PetscReal}, Ptr{$PetscInt}, Ptr{$PetscReal}, Ptr{$PetscInt}),
               U, Y, E, wnormtype, atol, vatol, rtol, vrtol, ignore_max, norm_, norm_loc_, norma_, norma_loc_, normr_, normr_loc_,
              )

	norm = norm_[]
	norm_loc = norm_loc_[]
	norma = norma_[]
	norma_loc = norma_loc_[]
	normr = normr_[]
	normr_loc = normr_loc_[]

	return norm,norm_loc,norma,norma_loc,normr,normr_loc
end 

"""
	VecExp(petsclib::PetscLibType, v::AbstractPetscVec) 
Replaces each component of a vector by e^x_i

Not Collective

Input Parameter:
- `v` - The vector

Output Parameter:
- `v` - The vector of exponents

Level: beginner

See also: `Vec`, `VecLog()`, `VecAbs()`, `VecSqrtAbs()`, `VecReciprocal()`

# External Links
$(_doc_external("Vec/VecExp"))
"""
function VecExp(petsclib::PetscLibType, v::AbstractPetscVec)
    error("VecExp: no generated method for these argument types")
end

@for_petsc function VecExp(petsclib::$UnionPetscLib, v::AbstractPetscVec )

    @chk ccall(
               (:VecExp, $petsc_library),
               PetscErrorCode,
               (CVec,),
               v,
              )


	return nothing
end 

"""
	VecFilter(petsclib::PetscLibType, v::AbstractPetscVec, tol::PetscReal) 
Set all values in the vector with an absolute value less than or equal to the tolerance to zero

Input Parameters:
- `v`   - The vector
- `tol` - The zero tolerance

Output Parameter:
- `v` - The filtered vector

Level: intermediate

See also: `VecCreate()`, `VecSet()`, `MatFilter()`

# External Links
$(_doc_external("Vec/VecFilter"))
"""
function VecFilter(petsclib::PetscLibType, v::AbstractPetscVec, tol::Real)
    error("VecFilter: no generated method for these argument types")
end

@for_petsc function VecFilter(petsclib::$UnionPetscLib, v::AbstractPetscVec, tol::$PetscReal )

    @chk ccall(
               (:VecFilter, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscReal),
               v, tol,
              )


	return nothing
end 

"""
	VecFinalizePackage(petsclib::PetscLibType) 
This function finalizes everything in the Vec package. It is called
from PetscFinalize().

Level: developer

See also: `PetscInitialize()`

# External Links
$(_doc_external("Vec/VecFinalizePackage"))
"""
function VecFinalizePackage(petsclib::PetscLibType)
    error("VecFinalizePackage: no generated method for these argument types")
end

@for_petsc function VecFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:VecFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	VecFischer(petsclib::PetscLibType, X::AbstractPetscVec, F::AbstractPetscVec, L::AbstractPetscVec, U::AbstractPetscVec, FB::AbstractPetscVec) 
Evaluates the Fischer-Burmeister function for complementarity
problems.

Logically Collective

Input Parameters:
- `X` - current point
- `F` - function evaluated at x
- `L` - lower bounds
- `U` - upper bounds

Output Parameter:
- `FB` - The Fischer-Burmeister function vector

Level: developer

See also: `Vec`, `VecSFischer()`, `MatDFischer()`, `MatDSFischer()`

# External Links
$(_doc_external("Tao/VecFischer"))
"""
function VecFischer(petsclib::PetscLibType, X::AbstractPetscVec, F::AbstractPetscVec, L::AbstractPetscVec, U::AbstractPetscVec, FB::AbstractPetscVec)
    error("VecFischer: no generated method for these argument types")
end

@for_petsc function VecFischer(petsclib::$UnionPetscLib, X::AbstractPetscVec, F::AbstractPetscVec, L::AbstractPetscVec, U::AbstractPetscVec, FB::AbstractPetscVec )

    @chk ccall(
               (:VecFischer, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, CVec, CVec, CVec),
               X, F, L, U, FB,
              )


	return nothing
end 

"""
	VecFlag(petsclib::PetscLibType, xin::AbstractPetscVec, flg::PetscInt) 
set infinity into the local part of the vector on any subset of MPI processes

Logically Collective

Input Parameters:
- `xin` - the vector, can be `NULL` but only if on all processes
- `flg` - indicates if this processes portion of the vector should be set to infinity

Level: developer

See also: `Vec`, `PetscLayout`, `VecGetLayout()`, `VecGetSize()`, `VecGetOwnershipRange()`, `VecGetOwnershipRanges()`

# External Links
$(_doc_external("Vec/VecFlag"))
"""
function VecFlag(petsclib::PetscLibType, xin::AbstractPetscVec, flg::Integer)
    error("VecFlag: no generated method for these argument types")
end

@for_petsc function VecFlag(petsclib::$UnionPetscLib, xin::AbstractPetscVec, flg::$PetscInt )

    @chk ccall(
               (:VecFlag, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt),
               xin, flg,
              )


	return nothing
end 

"""
	a::Vector{PetscScalar} = VecGetArray(petsclib::PetscLibType, x::AbstractPetscVec) 
Returns a pointer to a contiguous array that contains this
MPI processes's portion of the vector data

Logically Collective

Input Parameter:
- `x` - the vector

Output Parameter:
- `a` - location to put pointer to the array

Level: beginner

See also: `Vec`, `VecRestoreArray()`, `VecGetArrayRead()`, `VecGetArrays()`, `VecPlaceArray()`, `VecGetArray2d()`,
`VecGetArrayPair()`, `VecRestoreArrayPair()`, `VecGetArrayWrite()`, `VecRestoreArrayWrite()`, `VecGetArrayAndMemType()`

# External Links
$(_doc_external("Vec/VecGetArray"))
"""
function VecGetArray(petsclib::PetscLibType, x::AbstractPetscVec)
    error("VecGetArray: no generated method for these argument types")
end

@for_petsc function VecGetArray(petsclib::$UnionPetscLib, x::AbstractPetscVec )
	a_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:VecGetArray, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Ptr{$PetscScalar}}),
               x, a_,
              )

	a = a_[] == C_NULL ? $PetscScalar[] : unsafe_wrap(Array, a_[], VecGetLocalSize(petsclib, x); own = false)

	return a
end 

# override for VecGetArray1d; C signature: VecGetArray1d(Vec x, PetscInt m, PetscInt mstart, PetscScalar* a[])
"""
	a::PetscArray = VecGetArray1d(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, mstart::PetscInt) 
Returns a pointer to a 1d contiguous array that contains this
processor's portion of the vector data.  You MUST call `VecRestoreArray1d()`
when you no longer need access to the array.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of two dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)

Output Parameter:
- `a` - location to put pointer to the array

Level: developer

See also: `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrays()`, `VecPlaceArray()`,
`VecRestoreArray2d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`, `VecGetArray3d()`, `VecRestoreArray3d()`,
`VecGetArray2d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecGetArray1d"))
"""
function VecGetArray1d(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, mstart::PetscInt) end

@for_petsc function VecGetArray1d(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, mstart::$PetscInt )
	a_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:VecGetArray1d, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, Ref{Ptr{$PetscScalar}}),
               x, m, mstart, a_,
              )

	data_ptr = unsafe_load(a_[])
	mat = unsafe_wrap(Array, data_ptr, m) 
	a = PetscArray(mat,a_[]) 

	return a
end

# override for VecGetArray1dWrite; C signature: VecGetArray1dWrite(Vec x, PetscInt m, PetscInt mstart, PetscScalar* a[])
"""
	a::PetscArray = VecGetArray1dWrite(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, mstart::PetscInt) 
Returns a pointer to a 1d contiguous array that will contain this
processor's portion of the vector data.  You MUST call `VecRestoreArray1dWrite()`
when you no longer need access to the array.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of two dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)

Output Parameter:
- `a` - location to put pointer to the array

Level: developer

See also: `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrays()`, `VecPlaceArray()`,
`VecRestoreArray2d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`, `VecGetArray3d()`, `VecRestoreArray3d()`,
`VecGetArray2d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecGetArray1dWrite"))
"""
function VecGetArray1dWrite(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, mstart::PetscInt) end

@for_petsc function VecGetArray1dWrite(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, mstart::$PetscInt )
	a_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:VecGetArray1dWrite, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, Ref{Ptr{$PetscScalar}}),
               x, m, mstart, a_,
              )

	data_ptr = unsafe_load(a_[])
	mat = unsafe_wrap(Array, data_ptr, m) 
	a = PetscArray(mat,data_ptr) 

	return a
end

# override for VecGetArray2d; C signature: VecGetArray2d(Vec x, PetscInt m, PetscInt n, PetscInt mstart, PetscInt nstart, PetscScalar** a[])
"""
	a::Vector{PetscScalar} = VecGetArray2d(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, mstart::PetscInt, nstart::PetscInt) 
Returns a pointer to a 2d contiguous array that contains this
processor's portion of the vector data.  You MUST call `VecRestoreArray2d()`
when you no longer need access to the array.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of two dimensional array
- `n`      - second dimension of two dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `nstart` - first index in the second coordinate direction (often 0)

Output Parameter:
- `a` - location to put pointer to the array

Level: developer

See also: `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrays()`, `VecPlaceArray()`,
`VecRestoreArray2d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`, `VecGetArray3d()`, `VecRestoreArray3d()`,
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecGetArray2d"))
"""
function VecGetArray2d(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, mstart::PetscInt, nstart::PetscInt) end

@for_petsc function VecGetArray2d(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, n::$PetscInt, mstart::$PetscInt, nstart::$PetscInt )

    arr_ptr = Ref{Ptr{Ptr{$PetscScalar}}}()

    @chk ccall(
               (:VecGetArray2d, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{$PetscScalar}}}),
               x, m, n, mstart, nstart, arr_ptr,
              )

    # Assume contiguous storage, use first row pointer
    #mat = unsafe_wrap(Array, data_ptr, (m, n))
    
    # there is a difference in C vs julia storage of arrays
    data_ptr = unsafe_load(arr_ptr[])
    sz = (m,n)
    perm = (2,1)
    mat = unsafe_wrap(Array, data_ptr, sz)
    mat = PermutedDimsArray(mat, perm)
    #arr = PetscArray(mat, a_)



	return PetscArray(mat, arr_ptr)
end

# override for VecGetArray2dRead; C signature: VecGetArray2dRead(Vec x, PetscInt m, PetscInt n, PetscInt mstart, PetscInt nstart, PetscScalar** a[])
"""
	a::Vector{PetscScalar} = VecGetArray2dRead(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, mstart::PetscInt, nstart::PetscInt) 
Returns a pointer to a 2d contiguous array that contains this
processor's portion of the vector data.  You MUST call `VecRestoreArray2dRead()`
when you no longer need access to the array.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of two dimensional array
- `n`      - second dimension of two dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `nstart` - first index in the second coordinate direction (often 0)

Output Parameter:
- `a` - location to put pointer to the array

Level: developer

See also: `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrays()`, `VecPlaceArray()`,
`VecRestoreArray2d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`, `VecGetArray3d()`, `VecRestoreArray3d()`,
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecGetArray2dRead"))
"""
function VecGetArray2dRead(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, mstart::PetscInt, nstart::PetscInt) end

@for_petsc function VecGetArray2dRead(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, n::$PetscInt, mstart::$PetscInt, nstart::$PetscInt )

    arr_ptr = Ref{Ptr{Ptr{$PetscScalar}}}()

    @chk ccall(
               (:VecGetArray2dRead, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{$PetscScalar}}}),
               x, m, n, mstart, nstart, arr_ptr,
              )

    # Assume contiguous storage, use first row pointer
    data_ptr = unsafe_load(arr_ptr[])
    mat = unsafe_wrap(Array, data_ptr, (m, n))

	return PetscArray(mat,arr_ptr)
end

# override for VecGetArray2dWrite; C signature: VecGetArray2dWrite(Vec x, PetscInt m, PetscInt n, PetscInt mstart, PetscInt nstart, PetscScalar** a[])
"""
	a::Vector{PetscScalar} = VecGetArray2dWrite(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, mstart::PetscInt, nstart::PetscInt) 
Returns a pointer to a 2d contiguous array that will contain this
processor's portion of the vector data.  You MUST call `VecRestoreArray2dWrite()`
when you no longer need access to the array.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of two dimensional array
- `n`      - second dimension of two dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `nstart` - first index in the second coordinate direction (often 0)

Output Parameter:
- `a` - location to put pointer to the array

Level: developer

See also: `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrays()`, `VecPlaceArray()`,
`VecRestoreArray2d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`, `VecGetArray3d()`, `VecRestoreArray3d()`,
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecGetArray2dWrite"))
"""
function VecGetArray2dWrite(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, mstart::PetscInt, nstart::PetscInt) end

@for_petsc function VecGetArray2dWrite(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, n::$PetscInt, mstart::$PetscInt, nstart::$PetscInt )

    arr_ptr = Ref{Ptr{Ptr{$PetscScalar}}}()

    @chk ccall(
               (:VecGetArray2dWrite, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{$PetscScalar}}}),
               x, m, n, mstart, nstart, arr_ptr,
            )

    # Assume contiguous storage, use first row pointer
    data_ptr = unsafe_load(arr_ptr[])
    mat = unsafe_wrap(Array, data_ptr, (m, n))

	return PetscArray(mat,arr_ptr)            
end

# override for VecGetArray3d; C signature: VecGetArray3d(Vec x, PetscInt m, PetscInt n, PetscInt p, PetscInt mstart, PetscInt nstart, PetscInt pstart, PetscScalar** a[])
"""
	a::PetscArray = VecGetArray3d(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt) 
Returns a pointer to a 3d contiguous array that contains this
processor's portion of the vector data.  You MUST call `VecRestoreArray3d()`
when you no longer need access to the array.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of three dimensional array
- `n`      - second dimension of three dimensional array
- `p`      - third dimension of three dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `nstart` - first index in the second coordinate direction (often 0)
- `pstart` - first index in the third coordinate direction (often 0)

Output Parameter:
- `a` - location to put pointer to the array

Level: developer

See also: `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrays()`, `VecPlaceArray()`,
`VecRestoreArray2d()`, `DMDAVecGetarray()`, `DMDAVecRestoreArray()`, `VecRestoreArray3d()`,
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecGetArray3d"))
"""
function VecGetArray3d(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt) end

@for_petsc function VecGetArray3d(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, n::$PetscInt, p::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, pstart::$PetscInt )
	a_ = Ref{Ptr{Ptr{Ptr{$PetscScalar}}}}()

    @chk ccall(
               (:VecGetArray3d, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{Ptr{$PetscScalar}}}}),
               x, m, n, p, mstart, nstart, pstart, a_,
              )

	data_ptr = unsafe_load(a_[])
	mat = unsafe_wrap(Array, data_ptr, (m,n,p)) 
	a = PetscArray(mat,data_ptr) 

	return a
end

# override for VecGetArray3dRead; C signature: VecGetArray3dRead(Vec x, PetscInt m, PetscInt n, PetscInt p, PetscInt mstart, PetscInt nstart, PetscInt pstart, PetscScalar** a[])
"""
	a::PetscArray = VecGetArray3dRead(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt) 
Returns a pointer to a 3d contiguous array that contains this
processor's portion of the vector data.  You MUST call `VecRestoreArray3dRead()`
when you no longer need access to the array.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of three dimensional array
- `n`      - second dimension of three dimensional array
- `p`      - third dimension of three dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `nstart` - first index in the second coordinate direction (often 0)
- `pstart` - first index in the third coordinate direction (often 0)

Output Parameter:
- `a` - location to put pointer to the array

Level: developer

See also: `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrays()`, `VecPlaceArray()`,
`VecRestoreArray2d()`, `DMDAVecGetarray()`, `DMDAVecRestoreArray()`, `VecGetArray3d()`, `VecRestoreArray3d()`,
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecGetArray3dRead"))
"""
function VecGetArray3dRead(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt) end

@for_petsc function VecGetArray3dRead(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, n::$PetscInt, p::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, pstart::$PetscInt )
	a_ = Ref{Ptr{Ptr{Ptr{$PetscScalar}}}}()

    @chk ccall(
               (:VecGetArray3dRead, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{Ptr{$PetscScalar}}}}),
               x, m, n, p, mstart, nstart, pstart, a_,
              )

	data_ptr = unsafe_load(a_[])
	mat = unsafe_wrap(Array, data_ptr, (m,n,p)) 
	a = PetscArray(mat,data_ptr) 

	return a
end

# override for VecGetArray3dWrite; C signature: VecGetArray3dWrite(Vec x, PetscInt m, PetscInt n, PetscInt p, PetscInt mstart, PetscInt nstart, PetscInt pstart, PetscScalar** a[])
"""
	a::PetscArray = VecGetArray3dWrite(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt) 
Returns a pointer to a 3d contiguous array that will contain this
processor's portion of the vector data.  You MUST call `VecRestoreArray3dWrite()`
when you no longer need access to the array.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of three dimensional array
- `n`      - second dimension of three dimensional array
- `p`      - third dimension of three dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `nstart` - first index in the second coordinate direction (often 0)
- `pstart` - first index in the third coordinate direction (often 0)

Output Parameter:
- `a` - location to put pointer to the array

Level: developer

See also: `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrays()`, `VecPlaceArray()`,
`VecRestoreArray2d()`, `DMDAVecGetarray()`, `DMDAVecRestoreArray()`, `VecGetArray3d()`, `VecRestoreArray3d()`,
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecGetArray3dWrite"))
"""
function VecGetArray3dWrite(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt) end

@for_petsc function VecGetArray3dWrite(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, n::$PetscInt, p::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, pstart::$PetscInt )
	a_ = Ref{Ptr{Ptr{Ptr{$PetscScalar}}}}()

    @chk ccall(
               (:VecGetArray3dWrite, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{Ptr{$PetscScalar}}}}),
               x, m, n, p, mstart, nstart, pstart, a_,
              )

	data_ptr = unsafe_load(a_[])
	mat = unsafe_wrap(Array, data_ptr, (m,n,p)) 
	a = PetscArray(mat,data_ptr) 

	return a
end

# override for VecGetArray4d; C signature: VecGetArray4d(Vec x, PetscInt m, PetscInt n, PetscInt p, PetscInt q, PetscInt mstart, PetscInt nstart, PetscInt pstart, PetscInt qstart, PetscScalar** a[])
"""
	a::PetscArray = VecGetArray4d(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, q::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, qstart::PetscInt) 
Returns a pointer to a 4d contiguous array that contains this processor's portion of the vector data.  You MUST call `VecRestoreArray4d()` when you no longer need access to the array.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of four dimensional array
- `n`      - second dimension of four dimensional array
- `p`      - third dimension of four dimensional array
- `q`      - fourth dimension of four dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `nstart` - first index in the second coordinate direction (often 0)
- `pstart` - first index in the third coordinate direction (often 0)
- `qstart` - first index in the fourth coordinate direction (often 0)

Output Parameter:
- `a` - location to put pointer to the array

Level: developer

See also: `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrays()`, `VecPlaceArray()`,
`VecRestoreArray2d()`, `DMDAVecGetarray()`, `DMDAVecRestoreArray()`, `VecGetArray3d()`, `VecRestoreArray3d()`,
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecGetArray4d"))
"""
function VecGetArray4d(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, q::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, qstart::PetscInt) end

@for_petsc function VecGetArray4d(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, n::$PetscInt, p::$PetscInt, q::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, pstart::$PetscInt, qstart::$PetscInt )
	a_ = Ref{Ptr{Ptr{Ptr{Ptr{$PetscScalar}}}}}()

    @chk ccall(
               (:VecGetArray4d, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{Ptr{Ptr{$PetscScalar}}}}}),
               x, m, n, p, q, mstart, nstart, pstart, qstart, a_,
              )

	data_ptr = unsafe_load(a_[])
	mat = unsafe_wrap(Array, data_ptr, (m,n,p,q)) 
	a = PetscArray(mat,data_ptr) 

	return a
end

# override for VecGetArray4dRead; C signature: VecGetArray4dRead(Vec x, PetscInt m, PetscInt n, PetscInt p, PetscInt q, PetscInt mstart, PetscInt nstart, PetscInt pstart, PetscInt qstart, PetscScalar** a[])
"""
	a::PetscArray = VecGetArray4dRead(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, q::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, qstart::PetscInt) 
Returns a pointer to a 4d contiguous array that contains this
processor's portion of the vector data.  You MUST call `VecRestoreArray4dRead()`
when you no longer need access to the array.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of four dimensional array
- `n`      - second dimension of four dimensional array
- `p`      - third dimension of four dimensional array
- `q`      - fourth dimension of four dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `nstart` - first index in the second coordinate direction (often 0)
- `pstart` - first index in the third coordinate direction (often 0)
- `qstart` - first index in the fourth coordinate direction (often 0)

Output Parameter:
- `a` - location to put pointer to the array

Level: beginner

See also: `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrays()`, `VecPlaceArray()`,
`VecRestoreArray2d()`, `DMDAVecGetarray()`, `DMDAVecRestoreArray()`, `VecGetArray3d()`, `VecRestoreArray3d()`,
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecGetArray4dRead"))
"""
function VecGetArray4dRead(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, q::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, qstart::PetscInt) end

@for_petsc function VecGetArray4dRead(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, n::$PetscInt, p::$PetscInt, q::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, pstart::$PetscInt, qstart::$PetscInt )
	a_ = Ref{Ptr{Ptr{Ptr{Ptr{$PetscScalar}}}}}()

    @chk ccall(
               (:VecGetArray4dRead, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{Ptr{Ptr{$PetscScalar}}}}}),
               x, m, n, p, q, mstart, nstart, pstart, qstart, a_,
              )

	data_ptr = unsafe_load(a_[])
	mat = unsafe_wrap(Array, data_ptr, (m,n,p,q)) 
	a = PetscArray(mat,data_ptr) 

	return a
end

# override for VecGetArray4dWrite; C signature: VecGetArray4dWrite(Vec x, PetscInt m, PetscInt n, PetscInt p, PetscInt q, PetscInt mstart, PetscInt nstart, PetscInt pstart, PetscInt qstart, PetscScalar** a[])
"""
	a::PetscArray = VecGetArray4dWrite(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, q::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, qstart::PetscInt) 
Returns a pointer to a 4d contiguous array that will contain this
processor's portion of the vector data.  You MUST call `VecRestoreArray4dWrite()`
when you no longer need access to the array.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of four dimensional array
- `n`      - second dimension of four dimensional array
- `p`      - third dimension of four dimensional array
- `q`      - fourth dimension of four dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `nstart` - first index in the second coordinate direction (often 0)
- `pstart` - first index in the third coordinate direction (often 0)
- `qstart` - first index in the fourth coordinate direction (often 0)

Output Parameter:
- `a` - location to put pointer to the array

Level: developer

See also: `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrays()`, `VecPlaceArray()`,
`VecRestoreArray2d()`, `DMDAVecGetarray()`, `DMDAVecRestoreArray()`, `VecGetArray3d()`, `VecRestoreArray3d()`,
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecGetArray4dWrite"))
"""
function VecGetArray4dWrite(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, q::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, qstart::PetscInt) end

@for_petsc function VecGetArray4dWrite(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, n::$PetscInt, p::$PetscInt, q::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, pstart::$PetscInt, qstart::$PetscInt )
	a_ = Ref{Ptr{Ptr{Ptr{Ptr{$PetscScalar}}}}}()

    @chk ccall(
               (:VecGetArray4dWrite, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{Ptr{Ptr{$PetscScalar}}}}}),
               x, m, n, p, q, mstart, nstart, pstart, qstart, a_,
              )

	data_ptr = unsafe_load(a_[])
	mat = unsafe_wrap(Array, data_ptr, (m,n,p,q)) 
	a = PetscArray(mat,data_ptr) 

	return a
end

"""
	a::Vector{PetscScalar},mtype::PetscMemType = VecGetArrayAndMemType(petsclib::PetscLibType, x::AbstractPetscVec) 
Like `VecGetArray()`, but if this is a standard device vector (e.g.,
`VECCUDA`), the returned pointer will be a device pointer to the device memory that contains
this MPI processes's portion of the vector data.

Logically Collective; No Fortran Support

Input Parameter:
- `x` - the vector

Output Parameters:
- `a`     - location to put pointer to the array
- `mtype` - memory type of the array

Level: beginner

See also: `Vec`, `VecRestoreArrayAndMemType()`, `VecGetArrayReadAndMemType()`, `VecGetArrayWriteAndMemType()`, `VecRestoreArray()`, `VecGetArrayRead()`, `VecGetArrays()`,
`VecPlaceArray()`, `VecGetArray2d()`, `VecGetArrayPair()`, `VecRestoreArrayPair()`, `VecGetArrayWrite()`, `VecRestoreArrayWrite()`

# External Links
$(_doc_external("Vec/VecGetArrayAndMemType"))
"""
function VecGetArrayAndMemType(petsclib::PetscLibType, x::AbstractPetscVec)
    error("VecGetArrayAndMemType: no generated method for these argument types")
end

@for_petsc function VecGetArrayAndMemType(petsclib::$UnionPetscLib, x::AbstractPetscVec )
	a_ = Ref{Ptr{$PetscScalar}}()
	mtype_ = Ref{PetscMemType}()

    @chk ccall(
               (:VecGetArrayAndMemType, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Ptr{$PetscScalar}}, Ptr{PetscMemType}),
               x, a_, mtype_,
              )

	mtype = mtype_[]
	a = a_[] == C_NULL ? $PetscScalar[] : unsafe_wrap(Array, a_[], VecGetLocalSize(petsclib, x); own = false)

	return a,mtype
end 

"""
	xv::Vector{PetscScalar},yv::Vector{PetscScalar} = VecGetArrayPair(petsclib::PetscLibType, x::AbstractPetscVec, y::AbstractPetscVec) 

# External Links
$(_doc_external("Vec/VecGetArrayPair"))
"""
function VecGetArrayPair(petsclib::PetscLibType, x::AbstractPetscVec, y::AbstractPetscVec)
    error("VecGetArrayPair: no generated method for these argument types")
end

@for_petsc function VecGetArrayPair(petsclib::$UnionPetscLib, x::AbstractPetscVec, y::AbstractPetscVec )
	xv_ = Ref{Ptr{$PetscScalar}}()
	yv_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:VecGetArrayPair, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, Ptr{Ptr{$PetscScalar}}, Ptr{Ptr{$PetscScalar}}),
               x, y, xv_, yv_,
              )

	xv = xv_[] == C_NULL ? $PetscScalar[] : unsafe_wrap(Array, xv_[], VecGetLocalSize(petsclib, x); own = false)
	yv = yv_[] == C_NULL ? $PetscScalar[] : unsafe_wrap(Array, yv_[], VecGetLocalSize(petsclib, x); own = false)

	return xv,yv
end 

"""
	a::Vector{PetscScalar} = VecGetArrayRead(petsclib::PetscLibType, x::AbstractPetscVec) 
Get read-only pointer to contiguous array containing this processor's portion of the vector data.

Not Collective

Input Parameter:
- `x` - the vector

Output Parameter:
- `a` - the array

Level: beginner

See also: `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrayPair()`, `VecRestoreArrayPair()`,
`VecGetArrayAndMemType()`

# External Links
$(_doc_external("Vec/VecGetArrayRead"))
"""
function VecGetArrayRead(petsclib::PetscLibType, x::AbstractPetscVec)
    error("VecGetArrayRead: no generated method for these argument types")
end

@for_petsc function VecGetArrayRead(petsclib::$UnionPetscLib, x::AbstractPetscVec )
	a_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:VecGetArrayRead, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Ptr{$PetscScalar}}),
               x, a_,
              )

	a = a_[] == C_NULL ? $PetscScalar[] : unsafe_wrap(Array, a_[], VecGetLocalSize(petsclib, x); own = false)

	return a
end 

"""
	a::Vector{PetscScalar},mtype::PetscMemType = VecGetArrayReadAndMemType(petsclib::PetscLibType, x::AbstractPetscVec) 
Like `VecGetArrayRead()`, but if the input vector is a device vector, it will return a read-only device pointer.
The returned pointer is guaranteed to point to up-to-date data. For host vectors, it functions as `VecGetArrayRead()`.

Not Collective; No Fortran Support

Input Parameter:
- `x` - the vector

Output Parameters:
- `a`     - the array
- `mtype` - memory type of the array

Level: beginner

See also: `Vec`, `VecRestoreArrayReadAndMemType()`, `VecGetArrayAndMemType()`, `VecGetArrayWriteAndMemType()`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrayPair()`, `VecRestoreArrayPair()`

# External Links
$(_doc_external("Vec/VecGetArrayReadAndMemType"))
"""
function VecGetArrayReadAndMemType(petsclib::PetscLibType, x::AbstractPetscVec)
    error("VecGetArrayReadAndMemType: no generated method for these argument types")
end

@for_petsc function VecGetArrayReadAndMemType(petsclib::$UnionPetscLib, x::AbstractPetscVec )
	a_ = Ref{Ptr{$PetscScalar}}()
	mtype_ = Ref{PetscMemType}()

    @chk ccall(
               (:VecGetArrayReadAndMemType, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Ptr{$PetscScalar}}, Ptr{PetscMemType}),
               x, a_, mtype_,
              )

	mtype = mtype_[]
	a = a_[] == C_NULL ? $PetscScalar[] : unsafe_wrap(Array, a_[], VecGetLocalSize(petsclib, x); own = false)

	return a,mtype
end 

"""
	a::Vector{PetscScalar} = VecGetArrayWrite(petsclib::PetscLibType, x::AbstractPetscVec) 
Returns a pointer to a contiguous array that WILL contain this
MPI processes's portion of the vector data.

Logically Collective

Input Parameter:
- `x` - the vector

Output Parameter:
- `a` - location to put pointer to the array

Level: intermediate

See also: `Vec`, `VecRestoreArray()`, `VecGetArrayRead()`, `VecGetArrays()`, `VecPlaceArray()`, `VecGetArray2d()`,
`VecGetArrayPair()`, `VecRestoreArrayPair()`, `VecGetArray()`, `VecRestoreArrayWrite()`, `VecGetArrayAndMemType()`

# External Links
$(_doc_external("Vec/VecGetArrayWrite"))
"""
function VecGetArrayWrite(petsclib::PetscLibType, x::AbstractPetscVec)
    error("VecGetArrayWrite: no generated method for these argument types")
end

@for_petsc function VecGetArrayWrite(petsclib::$UnionPetscLib, x::AbstractPetscVec )
	a_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:VecGetArrayWrite, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Ptr{$PetscScalar}}),
               x, a_,
              )

	a = a_[] == C_NULL ? $PetscScalar[] : unsafe_wrap(Array, a_[], VecGetLocalSize(petsclib, x); own = false)

	return a
end 

"""
	a::Vector{PetscScalar},mtype::PetscMemType = VecGetArrayWriteAndMemType(petsclib::PetscLibType, x::AbstractPetscVec) 
Like `VecGetArrayWrite()`, but if this is a device vector it will always return
a device pointer to the device memory that contains this processor's portion of the vector data.

Logically Collective; No Fortran Support

Input Parameter:
- `x` - the vector

Output Parameters:
- `a`     - the array
- `mtype` - memory type of the array

Level: beginner

See also: `Vec`, `VecRestoreArrayWriteAndMemType()`, `VecGetArrayReadAndMemType()`, `VecGetArrayAndMemType()`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrayPair()`, `VecRestoreArrayPair()`

# External Links
$(_doc_external("Vec/VecGetArrayWriteAndMemType"))
"""
function VecGetArrayWriteAndMemType(petsclib::PetscLibType, x::AbstractPetscVec)
    error("VecGetArrayWriteAndMemType: no generated method for these argument types")
end

@for_petsc function VecGetArrayWriteAndMemType(petsclib::$UnionPetscLib, x::AbstractPetscVec )
	a_ = Ref{Ptr{$PetscScalar}}()
	mtype_ = Ref{PetscMemType}()

    @chk ccall(
               (:VecGetArrayWriteAndMemType, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Ptr{$PetscScalar}}, Ptr{PetscMemType}),
               x, a_, mtype_,
              )

	mtype = mtype_[]
	a = a_[] == C_NULL ? $PetscScalar[] : unsafe_wrap(Array, a_[], VecGetLocalSize(petsclib, x); own = false)

	return a,mtype
end 

"""
	a::Vector{Ptr{PetscScalar}} = VecGetArrays(petsclib::PetscLibType, x::Vector{<:AbstractPetscVec}, n::PetscInt) 
Returns a pointer to the arrays in a set of vectors
that were created by a call to `VecDuplicateVecs()`.

Logically Collective; No Fortran Support

Input Parameters:
- `x` - the vectors
- `n` - the number of vectors

Output Parameter:
- `a` - location to put pointer to the array

Level: intermediate

See also: `Vec`, `VecGetArray()`, `VecRestoreArrays()`

# External Links
$(_doc_external("Vec/VecGetArrays"))
"""
function VecGetArrays(petsclib::PetscLibType, x::Vector{<:AbstractPetscVec}, n::Integer)
    error("VecGetArrays: no generated method for these argument types")
end

@for_petsc function VecGetArrays(petsclib::$UnionPetscLib, x::Vector{<:AbstractPetscVec}, n::$PetscInt )
	a_ = Ref{Ptr{Ptr{$PetscScalar}}}()

    @chk ccall(
               (:VecGetArrays, $petsc_library),
               PetscErrorCode,
               (Ptr{CVec}, $PetscInt, Ptr{Ptr{Ptr{$PetscScalar}}}),
               x, n, a_,
              )

	a = a_[] == C_NULL ? Ptr{$PetscScalar}[] : unsafe_wrap(Array, a_[], n; own = false)

	return a
end 

"""
	flg::PetscBool = VecGetBindingPropagates(petsclib::PetscLibType, v::AbstractPetscVec) 
Gets whether the state of being bound to the CPU for a GPU vector type propagates to child and some other associated objects

Input Parameter:
- `v` - the vector

Output Parameter:
- `flg` - flag indicating whether the boundtocpu flag will be propagated

Level: developer

See also: `Vec`, `VecSetBindingPropagates()`

# External Links
$(_doc_external("Vec/VecGetBindingPropagates"))
"""
function VecGetBindingPropagates(petsclib::PetscLibType, v::AbstractPetscVec)
    error("VecGetBindingPropagates: no generated method for these argument types")
end

@for_petsc function VecGetBindingPropagates(petsclib::$UnionPetscLib, v::AbstractPetscVec )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:VecGetBindingPropagates, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{PetscBool}),
               v, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	bs::PetscInt = VecGetBlockSize(petsclib::PetscLibType, v::AbstractPetscVec) 
Gets the blocksize for the vector, i.e. what is used for `VecSetValuesBlocked()`
and `VecSetValuesBlockedLocal()`.

Not Collective

Input Parameter:
- `v` - the vector

Output Parameter:
- `bs` - the blocksize

Level: advanced

See also: `Vec`, `VecSetValuesBlocked()`, `VecSetLocalToGlobalMapping()`, `VecSetBlockSize()`

# External Links
$(_doc_external("Vec/VecGetBlockSize"))
"""
function VecGetBlockSize(petsclib::PetscLibType, v::AbstractPetscVec)
    error("VecGetBlockSize: no generated method for these argument types")
end

@for_petsc function VecGetBlockSize(petsclib::$UnionPetscLib, v::AbstractPetscVec )
	bs_ = Ref{$PetscInt}()

    @chk ccall(
               (:VecGetBlockSize, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{$PetscInt}),
               v, bs_,
              )

	bs = bs_[]

	return bs
end 

"""
	dm::PetscDM = VecGetDM(petsclib::PetscLibType, v::AbstractPetscVec) 
Gets the `DM` defining the data layout of the vector

Not Collective

Input Parameter:
- `v` - The `Vec`

Output Parameter:
- `dm` - The `DM`

Level: intermediate

See also: `DM`, `VecSetDM()`, `DMGetLocalVector()`, `DMGetGlobalVector()`, `DMSetVecType()`

# External Links
$(_doc_external("DM/VecGetDM"))
"""
function VecGetDM(petsclib::PetscLibType, v::AbstractPetscVec)
    error("VecGetDM: no generated method for these argument types")
end

@for_petsc function VecGetDM(petsclib::$UnionPetscLib, v::AbstractPetscVec )
	dm_ = Ref{CDM}()

    @chk ccall(
               (:VecGetDM, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{CDM}),
               v, dm_,
              )

	dm = PetscDM(dm_[], petsclib)

	return dm
end 

"""
	map::PetscLayout = VecGetLayout(petsclib::PetscLibType, x::AbstractPetscVec) 
get `PetscLayout` describing a vector layout

Not Collective

Input Parameter:
- `x` - the vector

Output Parameter:
- `map` - the layout

Level: developer

See also: `PetscLayout`, `Vec`, `VecGetSize()`, `VecGetOwnershipRange()`, `VecGetOwnershipRanges()`

# External Links
$(_doc_external("Vec/VecGetLayout"))
"""
function VecGetLayout(petsclib::PetscLibType, x::AbstractPetscVec)
    error("VecGetLayout: no generated method for these argument types")
end

@for_petsc function VecGetLayout(petsclib::$UnionPetscLib, x::AbstractPetscVec )
	map_ = Ref{PetscLayout}()

    @chk ccall(
               (:VecGetLayout, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{PetscLayout}),
               x, map_,
              )

	map = map_[]

	return map
end 

"""
	size::PetscInt = VecGetLocalSize(petsclib::PetscLibType, x::AbstractPetscVec) 
Returns the number of elements of the vector stored
in local memory (that is on this MPI process)

Not Collective

Input Parameter:
- `x` - the vector

Output Parameter:
- `size` - the length of the local piece of the vector

Level: beginner

See also: `Vec`, `VecGetSize()`

# External Links
$(_doc_external("Vec/VecGetLocalSize"))
"""
function VecGetLocalSize(petsclib::PetscLibType, x::AbstractPetscVec)
    error("VecGetLocalSize: no generated method for these argument types")
end

@for_petsc function VecGetLocalSize(petsclib::$UnionPetscLib, x::AbstractPetscVec )
	size_ = Ref{$PetscInt}()

    @chk ccall(
               (:VecGetLocalSize, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{$PetscInt}),
               x, size_,
              )

	size = size_[]

	return size
end 

"""
	mapping::ISLocalToGlobalMapping = VecGetLocalToGlobalMapping(petsclib::PetscLibType, X::AbstractPetscVec) 
Gets the local-to-global numbering set by `VecSetLocalToGlobalMapping()`

Not Collective

Input Parameter:
- `X` - the vector

Output Parameter:
- `mapping` - the mapping

Level: advanced

See also: `Vec`, `VecSetValuesLocal()`, `VecSetLocalToGlobalMapping()`

# External Links
$(_doc_external("Vec/VecGetLocalToGlobalMapping"))
"""
function VecGetLocalToGlobalMapping(petsclib::PetscLibType, X::AbstractPetscVec)
    error("VecGetLocalToGlobalMapping: no generated method for these argument types")
end

@for_petsc function VecGetLocalToGlobalMapping(petsclib::$UnionPetscLib, X::AbstractPetscVec )
	mapping_ = Ref{ISLocalToGlobalMapping}()

    @chk ccall(
               (:VecGetLocalToGlobalMapping, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{ISLocalToGlobalMapping}),
               X, mapping_,
              )

	mapping = mapping_[]

	return mapping
end 

"""
	VecGetLocalVector(petsclib::PetscLibType, v::AbstractPetscVec, w::AbstractPetscVec) 
Maps the local portion of a vector into a
vector.

Collective

Input Parameter:
- `v` - The vector for which the local vector is desired.

Output Parameter:
- `w` - Upon exit this contains the local vector.

Level: beginner

See also: `Vec`, `VecCreateLocalVector()`, `VecRestoreLocalVector()`, `VecGetLocalVectorRead()`, `VecGetArrayRead()`, `VecGetArray()`

# External Links
$(_doc_external("Vec/VecGetLocalVector"))
"""
function VecGetLocalVector(petsclib::PetscLibType, v::AbstractPetscVec, w::AbstractPetscVec)
    error("VecGetLocalVector: no generated method for these argument types")
end

@for_petsc function VecGetLocalVector(petsclib::$UnionPetscLib, v::AbstractPetscVec, w::AbstractPetscVec )

    @chk ccall(
               (:VecGetLocalVector, $petsc_library),
               PetscErrorCode,
               (CVec, CVec),
               v, w,
              )


	return nothing
end 

"""
	VecGetLocalVectorRead(petsclib::PetscLibType, v::AbstractPetscVec, w::AbstractPetscVec) 
Maps the local portion of a vector into a
vector.

Not Collective.

Input Parameter:
- `v` - The vector for which the local vector is desired.

Output Parameter:
- `w` - Upon exit this contains the local vector.

Level: beginner

See also: `Vec`, `VecCreateLocalVector()`, `VecRestoreLocalVectorRead()`, `VecGetLocalVector()`, `VecGetArrayRead()`, `VecGetArray()`

# External Links
$(_doc_external("Vec/VecGetLocalVectorRead"))
"""
function VecGetLocalVectorRead(petsclib::PetscLibType, v::AbstractPetscVec, w::AbstractPetscVec)
    error("VecGetLocalVectorRead: no generated method for these argument types")
end

@for_petsc function VecGetLocalVectorRead(petsclib::$UnionPetscLib, v::AbstractPetscVec, w::AbstractPetscVec )

    @chk ccall(
               (:VecGetLocalVectorRead, $petsc_library),
               PetscErrorCode,
               (CVec, CVec),
               v, w,
              )


	return nothing
end 

"""
	mask::PetscOffloadMask = VecGetOffloadMask(petsclib::PetscLibType, v::AbstractPetscVec) 
Get the offload mask of a `Vec`

Not Collective

Input Parameter:
- `v` - the vector

Output Parameter:
- `mask` - corresponding `PetscOffloadMask` enum value.

Level: intermediate

See also: `Vec`, `VecCreateSeqCUDA()`, `VecCreateSeqViennaCL()`, `VecGetArray()`, `VecGetType()`

# External Links
$(_doc_external("Vec/VecGetOffloadMask"))
"""
function VecGetOffloadMask(petsclib::PetscLibType, v::AbstractPetscVec)
    error("VecGetOffloadMask: no generated method for these argument types")
end

@for_petsc function VecGetOffloadMask(petsclib::$UnionPetscLib, v::AbstractPetscVec )
	mask_ = Ref{PetscOffloadMask}()

    @chk ccall(
               (:VecGetOffloadMask, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{PetscOffloadMask}),
               v, mask_,
              )

	mask = mask_[]

	return mask
end 

"""
	prefix::Ptr{Cchar} = VecGetOptionsPrefix(petsclib::PetscLibType, v::AbstractPetscVec) 
Sets the prefix used for searching for all
Vec options in the database.

Not Collective

Input Parameter:
- `v` - the `Vec` context

Output Parameter:
- `prefix` - pointer to the prefix string used

Level: advanced

See also: `Vec`, `VecAppendOptionsPrefix()`

# External Links
$(_doc_external("Vec/VecGetOptionsPrefix"))
"""
function VecGetOptionsPrefix(petsclib::PetscLibType, v::AbstractPetscVec)
    error("VecGetOptionsPrefix: no generated method for these argument types")
end

@for_petsc function VecGetOptionsPrefix(petsclib::$UnionPetscLib, v::AbstractPetscVec )
	prefix_ = Ref{Ptr{Cchar}}()

    @chk ccall(
               (:VecGetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Ptr{Cchar}}),
               v, prefix_,
              )

	prefix = prefix_[]

	return prefix
end 

"""
	low::PetscInt,high::PetscInt = VecGetOwnershipRange(petsclib::PetscLibType, x::AbstractPetscVec) 
Returns the range of indices owned by
this process. The vector is laid out with the
first `n1` elements on the first processor, next `n2` elements on the
second, etc.  For certain parallel layouts this range may not be
well defined.

Not Collective

Input Parameter:
- `x` - the vector

Output Parameters:
- `low`  - the first local element, pass in `NULL` if not interested
- `high` - one more than the last local element, pass in `NULL` if not interested

Level: beginner

See also: `Vec`, `MatGetOwnershipRange()`, `MatGetOwnershipRanges()`, `VecGetOwnershipRanges()`, `PetscSplitOwnership()`,
`VecSetSizes()`, `VecCreateMPI()`, `PetscLayout`, `DMDAGetGhostCorners()`, `DM`

# External Links
$(_doc_external("Vec/VecGetOwnershipRange"))
"""
function VecGetOwnershipRange(petsclib::PetscLibType, x::AbstractPetscVec)
    error("VecGetOwnershipRange: no generated method for these argument types")
end

@for_petsc function VecGetOwnershipRange(petsclib::$UnionPetscLib, x::AbstractPetscVec )
	low_ = Ref{$PetscInt}()
	high_ = Ref{$PetscInt}()

    @chk ccall(
               (:VecGetOwnershipRange, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{$PetscInt}, Ptr{$PetscInt}),
               x, low_, high_,
              )

	low = low_[]
	high = high_[]

	return low,high
end 

"""
	ranges::Vector{PetscInt} = VecGetOwnershipRanges(petsclib::PetscLibType, x::AbstractPetscVec) 
Returns the range of indices owned by EACH processor,
The vector is laid out with the
first `n1` elements on the first processor, next `n2` elements on the
second, etc.  For certain parallel layouts this range may not be
well defined.

Not Collective

Input Parameter:
- `x` - the vector

Output Parameter:
- `ranges` - array of length `size` + 1 with the start and end+1 for each process

Level: beginner

See also: `Vec`, `MatGetOwnershipRange()`, `MatGetOwnershipRanges()`, `VecGetOwnershipRange()`, `PetscSplitOwnership()`,
`VecSetSizes()`, `VecCreateMPI()`, `PetscLayout`, `DMDAGetGhostCorners()`, `DM`

# External Links
$(_doc_external("Vec/VecGetOwnershipRanges"))
"""
function VecGetOwnershipRanges(petsclib::PetscLibType, x::AbstractPetscVec)
    error("VecGetOwnershipRanges: no generated method for these argument types")
end

@for_petsc function VecGetOwnershipRanges(petsclib::$UnionPetscLib, x::AbstractPetscVec )
	ranges_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:VecGetOwnershipRanges, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Ptr{$PetscInt}}),
               x, ranges_,
              )

	nproc = MPI.Comm_size(PetscObjectGetComm(petsclib, x))
	ranges = ranges_[] == C_NULL ? $PetscInt[] : unsafe_wrap(Array, ranges_[], nproc + 1; own = false)

	return ranges
end 

"""
	mbytes::Csize_t = VecGetPinnedMemoryMin(petsclib::PetscLibType, v::AbstractPetscVec) 
Get the minimum data size for which pinned memory will be used for host (CPU) allocations.

Logically Collective

Input Parameter:
- `v` - the vector

Output Parameter:
- `mbytes` - minimum data size in bytes

Level: developer

See also: `Vec`, `VecSetPinnedMemoryMin()`

# External Links
$(_doc_external("Vec/VecGetPinnedMemoryMin"))
"""
function VecGetPinnedMemoryMin(petsclib::PetscLibType, v::AbstractPetscVec)
    error("VecGetPinnedMemoryMin: no generated method for these argument types")
end

@for_petsc function VecGetPinnedMemoryMin(petsclib::$UnionPetscLib, v::AbstractPetscVec )
	mbytes_ = Ref{Csize_t}()

    @chk ccall(
               (:VecGetPinnedMemoryMin, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Csize_t}),
               v, mbytes_,
              )

	mbytes = mbytes_[]

	return mbytes
end 

"""
	size::PetscInt = VecGetSize(petsclib::PetscLibType, x::AbstractPetscVec) 
Returns the global number of elements of the vector.

Not Collective

Input Parameter:
- `x` - the vector

Output Parameter:
- `size` - the global length of the vector

Level: beginner

See also: `Vec`, `VecGetLocalSize()`

# External Links
$(_doc_external("Vec/VecGetSize"))
"""
function VecGetSize(petsclib::PetscLibType, x::AbstractPetscVec)
    error("VecGetSize: no generated method for these argument types")
end

@for_petsc function VecGetSize(petsclib::$UnionPetscLib, x::AbstractPetscVec )
	size_ = Ref{$PetscInt}()

    @chk ccall(
               (:VecGetSize, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{$PetscInt}),
               x, size_,
              )

	size = size_[]

	return size
end 

"""
	state::PetscObjectState = VecGetState(petsclib::PetscLibType, v::AbstractPetscVec) 
Gets the state of a `Vec`.

Not Collective

Input Parameter:
- `v` - the `Vec` context

Output Parameter:
- `state` - the object state

Level: advanced

See also: `Vec`, `VecCreate()`, `PetscObjectStateGet()`

# External Links
$(_doc_external("Vec/VecGetState"))
"""
function VecGetState(petsclib::PetscLibType, v::AbstractPetscVec)
    error("VecGetState: no generated method for these argument types")
end

@for_petsc function VecGetState(petsclib::$UnionPetscLib, v::AbstractPetscVec )
	state_ = Ref{PetscObjectState}()

    @chk ccall(
               (:VecGetState, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{PetscObjectState}),
               v, state_,
              )

	state = state_[]

	return state
end 

"""
	Y::PetscVec = VecGetSubVector(petsclib::PetscLibType, X::AbstractPetscVec, is::AbstractIS) 
Gets a vector representing part of another vector

Collective

Input Parameters:
- `X`  - vector from which to extract a subvector
- `is` - index set representing portion of `X` to extract

Output Parameter:
- `Y` - subvector corresponding to `is`

Level: advanced

See also: `Vec`, `IS`, `VECNEST`, `MatCreateSubMatrix()`

# External Links
$(_doc_external("Vec/VecGetSubVector"))
"""
function VecGetSubVector(petsclib::PetscLibType, X::AbstractPetscVec, is::AbstractIS)
    error("VecGetSubVector: no generated method for these argument types")
end

@for_petsc function VecGetSubVector(petsclib::$UnionPetscLib, X::AbstractPetscVec, is::AbstractIS )
	Y_ = Ref{CVec}()

    @chk ccall(
               (:VecGetSubVector, $petsc_library),
               PetscErrorCode,
               (CVec, CIS, Ptr{CVec}),
               X, is, Y_,
              )

	Y = PetscVec(Y_[], petsclib)

	return Y
end 

"""
	type::VecType = VecGetType(petsclib::PetscLibType, vec::AbstractPetscVec) 
Gets the vector type name (as a string) from a `Vec`.

Not Collective

Input Parameter:
- `vec` - The vector

Output Parameter:
- `type` - The `VecType` of the vector

Level: intermediate

See also: `Vec`, `VecType`, `VecCreate()`, `VecDuplicate()`, `VecDuplicateVecs()`, `PetscObjectTypeCompare()`, `PetscObjectTypeCompareAny()`

# External Links
$(_doc_external("Vec/VecGetType"))
"""
function VecGetType(petsclib::PetscLibType, vec::AbstractPetscVec)
    error("VecGetType: no generated method for these argument types")
end

@for_petsc function VecGetType(petsclib::$UnionPetscLib, vec::AbstractPetscVec )
	type_ = Ref{VecType}()

    @chk ccall(
               (:VecGetType, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{VecType}),
               vec, type_,
              )

	type = type_[] == C_NULL ? "" : unsafe_string(type_[])

	return type
end 

"""
	y::Vector{PetscScalar} = VecGetValues(petsclib::PetscLibType, x::AbstractPetscVec, ni::PetscInt, ix::Vector{PetscInt}) 
Gets values from certain locations of a vector. Currently
can only get values on the same processor on which they are owned

Not Collective

Input Parameters:
- `x`  - vector to get values from
- `ni` - number of elements to get
- `ix` - indices where to get them from (in global 1d numbering)

Output Parameter:
- `y` - array of values, must be passed in with a length of `ni`

Level: beginner

See also: `Vec`, `VecAssemblyBegin()`, `VecAssemblyEnd()`, `VecSetValues()`

# External Links
$(_doc_external("Vec/VecGetValues"))
"""
function VecGetValues(petsclib::PetscLibType, x::AbstractPetscVec, ni::Integer, ix::AbstractVector{<:Number})
    error("VecGetValues: no generated method for these argument types")
end

@for_petsc function VecGetValues(petsclib::$UnionPetscLib, x::AbstractPetscVec, ni::$PetscInt, ix::Vector{$PetscInt} )
	y = Vector{$PetscScalar}(undef, ni)

    @chk ccall(
               (:VecGetValues, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}),
               x, ni, ix, y,
              )


	return y
end 

"""
	values::Ptr{PetscScalar} = VecGetValuesSection(petsclib::PetscLibType, v::AbstractPetscVec, s::PetscSection, point::PetscInt) 
Gets all the values associated with a given point, according to the section, in the given `Vec`

Not Collective

Input Parameters:
- `v`     - the `Vec`
- `s`     - the organizing `PetscSection`
- `point` - the point

Output Parameter:
- `values` - the array of output values

Level: developer

See also: `PetscSection`, `PetscSectionCreate()`, `VecSetValuesSection()`

# External Links
$(_doc_external("Vec/VecGetValuesSection"))
"""
function VecGetValuesSection(petsclib::PetscLibType, v::AbstractPetscVec, s::PetscSection, point::Integer)
    error("VecGetValuesSection: no generated method for these argument types")
end

@for_petsc function VecGetValuesSection(petsclib::$UnionPetscLib, v::AbstractPetscVec, s::PetscSection, point::$PetscInt )
	values_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:VecGetValuesSection, $petsc_library),
               PetscErrorCode,
               (CVec, PetscSection, $PetscInt, Ptr{Ptr{$PetscScalar}}),
               v, s, point, values_,
              )

	values = values_[]

	return values
end 

"""
	ghost::IS = VecGhostGetGhostIS(petsclib::PetscLibType, X::AbstractPetscVec) 
Return ghosting indices of a ghost vector

Input Parameters:
- `X` - ghost vector

Output Parameter:
- `ghost` - ghosting indices

Level: beginner

See also: `VecCreateGhostWithArray()`, `VecCreateMPIWithArray()`

# External Links
$(_doc_external("Vec/VecGhostGetGhostIS"))
"""
function VecGhostGetGhostIS(petsclib::PetscLibType, X::AbstractPetscVec)
    error("VecGhostGetGhostIS: no generated method for these argument types")
end

@for_petsc function VecGhostGetGhostIS(petsclib::$UnionPetscLib, X::AbstractPetscVec )
	ghost_ = Ref{CIS}()

    @chk ccall(
               (:VecGhostGetGhostIS, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{CIS}),
               X, ghost_,
              )

	ghost = IS(ghost_[], petsclib)

	return ghost
end 

"""
	l::PetscVec = VecGhostGetLocalForm(petsclib::PetscLibType, g::AbstractPetscVec) 
Obtains the local ghosted representation of
a parallel vector (obtained with `VecCreateGhost()`, `VecCreateGhostWithArray()` or `VecCreateSeq()`).

Logically Collective

Input Parameter:
- `g` - the global vector

Output Parameter:
- `l` - the local (ghosted) representation,`NULL` if `g` is not ghosted

Level: advanced

See also: `VecGhostUpdateBegin()`, `VecGhostUpdateEnd()`, `Vec`, `VecType`, `VecCreateGhost()`, `VecGhostRestoreLocalForm()`, `VecCreateGhostWithArray()`

# External Links
$(_doc_external("Vec/VecGhostGetLocalForm"))
"""
function VecGhostGetLocalForm(petsclib::PetscLibType, g::AbstractPetscVec)
    error("VecGhostGetLocalForm: no generated method for these argument types")
end

@for_petsc function VecGhostGetLocalForm(petsclib::$UnionPetscLib, g::AbstractPetscVec )
	l_ = Ref{CVec}()

    @chk ccall(
               (:VecGhostGetLocalForm, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{CVec}),
               g, l_,
              )

	l = PetscVec(l_[], petsclib)

	return l
end 

"""
	flg::PetscBool = VecGhostIsLocalForm(petsclib::PetscLibType, g::AbstractPetscVec, l::AbstractPetscVec) 
Checks if a given vector is the local form of a global vector

Not Collective

Input Parameters:
- `g` - the global vector
- `l` - the local vector

Output Parameter:
- `flg` - `PETSC_TRUE` if `l` is the local form

Level: advanced

See also: `Vec`, `VecType`, `VecCreateGhost()`, `VecGhostRestoreLocalForm()`, `VecCreateGhostWithArray()`, `VecGhostGetLocalForm()`

# External Links
$(_doc_external("Vec/VecGhostIsLocalForm"))
"""
function VecGhostIsLocalForm(petsclib::PetscLibType, g::AbstractPetscVec, l::AbstractPetscVec)
    error("VecGhostIsLocalForm: no generated method for these argument types")
end

@for_petsc function VecGhostIsLocalForm(petsclib::$UnionPetscLib, g::AbstractPetscVec, l::AbstractPetscVec )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:VecGhostIsLocalForm, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, Ptr{PetscBool}),
               g, l, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	VecGhostRestoreLocalForm(petsclib::PetscLibType, g::AbstractPetscVec, l::AbstractPetscVec) 
Restores the local ghosted representation of
a parallel vector obtained with `VecGhostGetLocalForm()`.

Logically Collective

Input Parameters:
- `g` - the global vector
- `l` - the local (ghosted) representation

Level: advanced

See also: `VecGhostUpdateBegin()`, `VecGhostUpdateEnd()`, `Vec`, `VecType`, `VecCreateGhost()`, `VecGhostGetLocalForm()`, `VecCreateGhostWithArray()`

# External Links
$(_doc_external("Vec/VecGhostRestoreLocalForm"))
"""
function VecGhostRestoreLocalForm(petsclib::PetscLibType, g::AbstractPetscVec, l::AbstractPetscVec)
    error("VecGhostRestoreLocalForm: no generated method for these argument types")
end

@for_petsc function VecGhostRestoreLocalForm(petsclib::$UnionPetscLib, g::AbstractPetscVec, l::AbstractPetscVec )
	l_ = Ref(l.ptr)

    @chk ccall(
               (:VecGhostRestoreLocalForm, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{CVec}),
               g, l_,
              )

	l.ptr = l_[]

	return nothing
end 

"""
	VecGhostUpdateBegin(petsclib::PetscLibType, g::AbstractPetscVec, insertmode::InsertMode, scattermode::ScatterMode) 
Begins the vector scatter to update the vector from
local representation to global or global representation to local.

Neighbor-wise Collective

Input Parameters:
- `g`           - the vector (obtained with `VecCreateGhost()` or `VecDuplicate()`)
- `insertmode`  - one of `ADD_VALUES`, `MAX_VALUES`, `MIN_VALUES` or `INSERT_VALUES`
- `scattermode` - one of `SCATTER_FORWARD` (update ghosts) or `SCATTER_REVERSE` (update local values from ghosts)

Level: advanced

See also: `Vec`, `VecType`, `VecCreateGhost()`, `VecGhostUpdateEnd()`, `VecGhostGetLocalForm()`,
`VecGhostRestoreLocalForm()`, `VecCreateGhostWithArray()`

# External Links
$(_doc_external("Vec/VecGhostUpdateBegin"))
"""
function VecGhostUpdateBegin(petsclib::PetscLibType, g::AbstractPetscVec, insertmode::InsertMode, scattermode::ScatterMode)
    error("VecGhostUpdateBegin: no generated method for these argument types")
end

@for_petsc function VecGhostUpdateBegin(petsclib::$UnionPetscLib, g::AbstractPetscVec, insertmode::InsertMode, scattermode::ScatterMode )

    @chk ccall(
               (:VecGhostUpdateBegin, $petsc_library),
               PetscErrorCode,
               (CVec, InsertMode, ScatterMode),
               g, insertmode, scattermode,
              )


	return nothing
end 

"""
	VecGhostUpdateEnd(petsclib::PetscLibType, g::AbstractPetscVec, insertmode::InsertMode, scattermode::ScatterMode) 
End the vector scatter to update the vector from
local representation to global or global representation to local.

Neighbor-wise Collective

Input Parameters:
- `g`           - the vector (obtained with `VecCreateGhost()` or `VecDuplicate()`)
- `insertmode`  - one of `ADD_VALUES`, `MAX_VALUES`, `MIN_VALUES` or `INSERT_VALUES`
- `scattermode` - one of `SCATTER_FORWARD` (update ghosts) or `SCATTER_REVERSE` (update local values from ghosts)

Level: advanced

See also: `Vec`, `VecType`, `VecCreateGhost()`, `VecGhostUpdateBegin()`, `VecGhostGetLocalForm()`,
`VecGhostRestoreLocalForm()`, `VecCreateGhostWithArray()`

# External Links
$(_doc_external("Vec/VecGhostUpdateEnd"))
"""
function VecGhostUpdateEnd(petsclib::PetscLibType, g::AbstractPetscVec, insertmode::InsertMode, scattermode::ScatterMode)
    error("VecGhostUpdateEnd: no generated method for these argument types")
end

@for_petsc function VecGhostUpdateEnd(petsclib::$UnionPetscLib, g::AbstractPetscVec, insertmode::InsertMode, scattermode::ScatterMode )

    @chk ccall(
               (:VecGhostUpdateEnd, $petsc_library),
               PetscErrorCode,
               (CVec, InsertMode, ScatterMode),
               g, insertmode, scattermode,
              )


	return nothing
end 

"""
	VecISAXPY(petsclib::PetscLibType, vfull::AbstractPetscVec, is::AbstractIS, alpha::PetscScalar, vreduced::AbstractPetscVec) 
Adds a reduced vector to the appropriate elements of a full-space vector.
vfull[is[i]] += alpha*vreduced[i]

Logically Collective

Input Parameters:
- `vfull`    - the full-space vector
- `is`       - the index set for the reduced space
- `alpha`    - the scalar coefficient
- `vreduced` - the reduced-space vector

Output Parameter:
- `vfull` - the sum of the full-space vector and reduced-space vector

Level: advanced

See also: `VecISCopy()`, `VecISSet()`, `VecAXPY()`

# External Links
$(_doc_external("Vec/VecISAXPY"))
"""
function VecISAXPY(petsclib::PetscLibType, vfull::AbstractPetscVec, is::AbstractIS, alpha::Number, vreduced::AbstractPetscVec)
    error("VecISAXPY: no generated method for these argument types")
end

@for_petsc function VecISAXPY(petsclib::$UnionPetscLib, vfull::AbstractPetscVec, is::AbstractIS, alpha::$PetscScalar, vreduced::AbstractPetscVec )

    @chk ccall(
               (:VecISAXPY, $petsc_library),
               PetscErrorCode,
               (CVec, CIS, $PetscScalar, CVec),
               vfull, is, alpha, vreduced,
              )


	return nothing
end 

"""
	VecISCopy(petsclib::PetscLibType, vfull::AbstractPetscVec, is::AbstractIS, mode::ScatterMode, vreduced::AbstractPetscVec) 
Copies between a reduced vector and the appropriate elements of a full-space vector.

Logically Collective

Input Parameters:
- `vfull`    - the full-space vector
- `is`       - the index set for the reduced space
- `mode`     - the direction of copying, `SCATTER_FORWARD` or `SCATTER_REVERSE`
- `vreduced` - the reduced-space vector

Output Parameter:
- `vfull` - the sum of the full-space vector and reduced-space vector

Level: advanced

See also: `VecISSet()`, `VecISAXPY()`, `VecCopy()`

# External Links
$(_doc_external("Vec/VecISCopy"))
"""
function VecISCopy(petsclib::PetscLibType, vfull::AbstractPetscVec, is::AbstractIS, mode::ScatterMode, vreduced::AbstractPetscVec)
    error("VecISCopy: no generated method for these argument types")
end

@for_petsc function VecISCopy(petsclib::$UnionPetscLib, vfull::AbstractPetscVec, is::AbstractIS, mode::ScatterMode, vreduced::AbstractPetscVec )

    @chk ccall(
               (:VecISCopy, $petsc_library),
               PetscErrorCode,
               (CVec, CIS, ScatterMode, CVec),
               vfull, is, mode, vreduced,
              )


	return nothing
end 

"""
	VecISSet(petsclib::PetscLibType, V::AbstractPetscVec, S::AbstractIS, c::PetscScalar) 
Sets the elements of a vector, specified by an index set, to a constant

Logically Collective

Input Parameters:
- `V` - the vector
- `S` - index set for the locations in the vector
- `c` - the constant

Level: advanced

See also: `VecISCopy()`, `VecISAXPY()`, `VecISShift()`, `VecSet()`

# External Links
$(_doc_external("Vec/VecISSet"))
"""
function VecISSet(petsclib::PetscLibType, V::AbstractPetscVec, S::AbstractIS, c::Number)
    error("VecISSet: no generated method for these argument types")
end

@for_petsc function VecISSet(petsclib::$UnionPetscLib, V::AbstractPetscVec, S::AbstractIS, c::$PetscScalar )

    @chk ccall(
               (:VecISSet, $petsc_library),
               PetscErrorCode,
               (CVec, CIS, $PetscScalar),
               V, S, c,
              )


	return nothing
end 

"""
	VecISShift(petsclib::PetscLibType, V::AbstractPetscVec, S::AbstractIS, c::PetscScalar) 
Shifts the elements of a vector, specified by an index set, by a constant

Logically Collective

Input Parameters:
- `V` - the vector
- `S` - index set for the locations in the vector
- `c` - the constant

Level: advanced

See also: `VecISCopy()`, `VecISAXPY()`, `VecISSet()`, `VecShift()`

# External Links
$(_doc_external("Vec/VecISShift"))
"""
function VecISShift(petsclib::PetscLibType, V::AbstractPetscVec, S::AbstractIS, c::Number)
    error("VecISShift: no generated method for these argument types")
end

@for_petsc function VecISShift(petsclib::$UnionPetscLib, V::AbstractPetscVec, S::AbstractIS, c::$PetscScalar )

    @chk ccall(
               (:VecISShift, $petsc_library),
               PetscErrorCode,
               (CVec, CIS, $PetscScalar),
               V, S, c,
              )


	return nothing
end 

"""
	VecImaginaryPart(petsclib::PetscLibType, v::AbstractPetscVec) 
Replaces a complex vector with its imaginary part

Collective

Input Parameter:
- `v` - the vector

Level: beginner

See also: `Vec`, `VecNorm()`, `VecRealPart()`

# External Links
$(_doc_external("Vec/VecImaginaryPart"))
"""
function VecImaginaryPart(petsclib::PetscLibType, v::AbstractPetscVec)
    error("VecImaginaryPart: no generated method for these argument types")
end

@for_petsc function VecImaginaryPart(petsclib::$UnionPetscLib, v::AbstractPetscVec )

    @chk ccall(
               (:VecImaginaryPart, $petsc_library),
               PetscErrorCode,
               (CVec,),
               v,
              )


	return nothing
end 

"""
	VecInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `Vec` package. It is called
from PetscDLLibraryRegister_petscvec() when using dynamic libraries, and on the first call to `VecCreate()`
when using shared or static libraries.

Level: developer

See also: `PetscInitialize()`

# External Links
$(_doc_external("Sys/VecInitializePackage"))
"""
function VecInitializePackage(petsclib::PetscLibType)
    error("VecInitializePackage: no generated method for these argument types")
end

@for_petsc function VecInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:VecInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	a::PetscScalar = VecKokkosPlaceArray(petsclib::PetscLibType, v::AbstractPetscVec) 

# External Links
$(_doc_external("Vec/VecKokkosPlaceArray"))
"""
function VecKokkosPlaceArray(petsclib::PetscLibType, v::AbstractPetscVec)
    error("VecKokkosPlaceArray: no generated method for these argument types")
end

@for_petsc function VecKokkosPlaceArray(petsclib::$UnionPetscLib, v::AbstractPetscVec )
	a_ = Ref{$PetscScalar}()

    @chk ccall(
               (:VecKokkosPlaceArray, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{$PetscScalar}),
               v, a_,
              )

	a = a_[]

	return a
end 

"""
	VecKokkosResetArray(petsclib::PetscLibType, v::AbstractPetscVec) 

# External Links
$(_doc_external("Vec/VecKokkosResetArray"))
"""
function VecKokkosResetArray(petsclib::PetscLibType, v::AbstractPetscVec)
    error("VecKokkosResetArray: no generated method for these argument types")
end

@for_petsc function VecKokkosResetArray(petsclib::$UnionPetscLib, v::AbstractPetscVec )

    @chk ccall(
               (:VecKokkosResetArray, $petsc_library),
               PetscErrorCode,
               (CVec,),
               v,
              )


	return nothing
end 

"""
	VecLoad(petsclib::PetscLibType, vec::AbstractPetscVec, viewer::PetscViewer) 
Loads a vector that has been stored in binary or HDF5 format
with `VecView()`.

Collective

Input Parameters:
- `vec`    - the newly loaded vector, this needs to have been created with `VecCreate()` or
some related function before the call to `VecLoad()`.
- `viewer` - binary file viewer, obtained from `PetscViewerBinaryOpen()` or
HDF5 file viewer, obtained from `PetscViewerHDF5Open()`

Level: intermediate

See also: `Vec`, `PetscViewerBinaryOpen()`, `VecView()`, `MatLoad()`

# External Links
$(_doc_external("Vec/VecLoad"))
"""
function VecLoad(petsclib::PetscLibType, vec::AbstractPetscVec, viewer::PetscViewer)
    error("VecLoad: no generated method for these argument types")
end

@for_petsc function VecLoad(petsclib::$UnionPetscLib, vec::AbstractPetscVec, viewer::PetscViewer )

    @chk ccall(
               (:VecLoad, $petsc_library),
               PetscErrorCode,
               (CVec, PetscViewer),
               vec, viewer,
              )


	return nothing
end 

"""
	state::PetscInt = VecLockGet(petsclib::PetscLibType, x::AbstractPetscVec) 
Get the current lock status of a vector

Logically Collective

Input Parameter:
- `x` - the vector

Output Parameter:
- `state` - greater than zero indicates the vector is locked for read; less than zero indicates the vector is
locked for write; equal to zero means the vector is unlocked, that is, it is free to read or write.

Level: advanced

See also: `Vec`, `VecRestoreArray()`, `VecGetArrayRead()`, `VecLockReadPush()`, `VecLockReadPop()`

# External Links
$(_doc_external("Vec/VecLockGet"))
"""
function VecLockGet(petsclib::PetscLibType, x::AbstractPetscVec)
    error("VecLockGet: no generated method for these argument types")
end

@for_petsc function VecLockGet(petsclib::$UnionPetscLib, x::AbstractPetscVec )
	state_ = Ref{$PetscInt}()

    @chk ccall(
               (:VecLockGet, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{$PetscInt}),
               x, state_,
              )

	state = state_[]

	return state
end 

"""
	file::Ptr{Cchar},func::Ptr{Cchar},line::Cint = VecLockGetLocation(petsclib::PetscLibType, x::AbstractPetscVec) 
Return the source code location where a `Vec` was most recently read-locked

Not Collective

Input Parameter:
- `x` - the vector

Output Parameters:
- `file` - the source file name of the most recent `VecLockReadPush()`, or `NULL` if none is active
- `func` - the function name of the most recent `VecLockReadPush()`, or `NULL` if none is active
- `line` - the source line number of the most recent `VecLockReadPush()`, or 0 if none is active

Level: developer

See also: `Vec`, `VecLockGet()`, `VecLockReadPush()`, `VecLockReadPop()`, `VecGetArray()`

# External Links
$(_doc_external("Vec/VecLockGetLocation"))
"""
function VecLockGetLocation(petsclib::PetscLibType, x::AbstractPetscVec)
    error("VecLockGetLocation: no generated method for these argument types")
end

@for_petsc function VecLockGetLocation(petsclib::$UnionPetscLib, x::AbstractPetscVec )
	file_ = Ref{Ptr{Cchar}}()
	func_ = Ref{Ptr{Cchar}}()
	line_ = Ref{Cint}()

    @chk ccall(
               (:VecLockGetLocation, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Ptr{Cchar}}, Ptr{Ptr{Cchar}}, Ptr{Cint}),
               x, file_, func_, line_,
              )

	file = file_[]
	func = func_[]
	line = line_[]

	return file,func,line
end 

"""
	VecLockReadPop(petsclib::PetscLibType, x::AbstractPetscVec) 
Pop a read-only lock from a vector

Logically Collective

Input Parameter:
- `x` - the vector

Level: intermediate

See also: `Vec`, `VecRestoreArray()`, `VecGetArrayRead()`, `VecLockReadPush()`, `VecLockGet()`

# External Links
$(_doc_external("Vec/VecLockReadPop"))
"""
function VecLockReadPop(petsclib::PetscLibType, x::AbstractPetscVec)
    error("VecLockReadPop: no generated method for these argument types")
end

@for_petsc function VecLockReadPop(petsclib::$UnionPetscLib, x::AbstractPetscVec )

    @chk ccall(
               (:VecLockReadPop, $petsc_library),
               PetscErrorCode,
               (CVec,),
               x,
              )


	return nothing
end 

"""
	VecLockReadPush(petsclib::PetscLibType, x::AbstractPetscVec) 
Push a read-only lock on a vector to prevent it from being written to

Logically Collective

Input Parameter:
- `x` - the vector

Level: intermediate

See also: `Vec`, `VecRestoreArray()`, `VecGetArrayRead()`, `VecLockReadPop()`, `VecLockGet()`

# External Links
$(_doc_external("Vec/VecLockReadPush"))
"""
function VecLockReadPush(petsclib::PetscLibType, x::AbstractPetscVec)
    error("VecLockReadPush: no generated method for these argument types")
end

@for_petsc function VecLockReadPush(petsclib::$UnionPetscLib, x::AbstractPetscVec )

    @chk ccall(
               (:VecLockReadPush, $petsc_library),
               PetscErrorCode,
               (CVec,),
               x,
              )


	return nothing
end 

"""
	VecLockWriteSet(petsclib::PetscLibType, x::AbstractPetscVec, flg::PetscBool) 
Lock or unlock a vector for exclusive read/write access

Logically Collective

Input Parameters:
- `x`   - the vector
- `flg` - `PETSC_TRUE` to lock the vector for exclusive read/write access; `PETSC_FALSE` to unlock it.

Level: intermediate

See also: `Vec`, `VecRestoreArray()`, `VecGetArrayRead()`, `VecLockReadPush()`, `VecLockReadPop()`, `VecLockGet()`

# External Links
$(_doc_external("Vec/VecLockWriteSet"))
"""
function VecLockWriteSet(petsclib::PetscLibType, x::AbstractPetscVec, flg::PetscBool)
    error("VecLockWriteSet: no generated method for these argument types")
end

@for_petsc function VecLockWriteSet(petsclib::$UnionPetscLib, x::AbstractPetscVec, flg::PetscBool )

    @chk ccall(
               (:VecLockWriteSet, $petsc_library),
               PetscErrorCode,
               (CVec, PetscBool),
               x, flg,
              )


	return nothing
end 

"""
	VecLog(petsclib::PetscLibType, v::AbstractPetscVec) 
Replaces each component of a vector by log(x_i), the natural logarithm

Not Collective

Input Parameter:
- `v` - The vector

Output Parameter:
- `v` - The vector of logs

Level: beginner

See also: `Vec`, `VecExp()`, `VecAbs()`, `VecSqrtAbs()`, `VecReciprocal()`

# External Links
$(_doc_external("Vec/VecLog"))
"""
function VecLog(petsclib::PetscLibType, v::AbstractPetscVec)
    error("VecLog: no generated method for these argument types")
end

@for_petsc function VecLog(petsclib::$UnionPetscLib, v::AbstractPetscVec )

    @chk ccall(
               (:VecLog, $petsc_library),
               PetscErrorCode,
               (CVec,),
               v,
              )


	return nothing
end 

"""
	VecMAXPBY(petsclib::PetscLibType, y::AbstractPetscVec, nv::PetscInt, alpha::Vector{PetscScalar}, beta::PetscScalar, x::Vector{<:AbstractPetscVec}) 
Computes `y = beta y + sum alpha[i] x[i]`

Logically Collective

Input Parameters:
- `nv`    - number of scalars and `x` vectors
- `alpha` - array of scalars
- `beta`  - scalar
- `y`     - one vector
- `x`     - array of vectors

Level: intermediate

See also: `Vec`, `VecMAXPY()`, `VecAYPX()`, `VecWAXPY()`, `VecAXPY()`, `VecAXPBYPCZ()`, `VecAXPBY()`

# External Links
$(_doc_external("Vec/VecMAXPBY"))
"""
function VecMAXPBY(petsclib::PetscLibType, y::AbstractPetscVec, nv::Integer, alpha::AbstractVector{<:Number}, beta::Number, x::Vector{<:AbstractPetscVec})
    error("VecMAXPBY: no generated method for these argument types")
end

@for_petsc function VecMAXPBY(petsclib::$UnionPetscLib, y::AbstractPetscVec, nv::$PetscInt, alpha::Vector{$PetscScalar}, beta::$PetscScalar, x::Vector{<:AbstractPetscVec} )

    @chk ccall(
               (:VecMAXPBY, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{$PetscScalar}, $PetscScalar, Ptr{CVec}),
               y, nv, alpha, beta, x,
              )


	return nothing
end 

"""
	VecMAXPY(petsclib::PetscLibType, y::AbstractPetscVec, nv::PetscInt, alpha::Vector{PetscScalar}, x::Vector{<:AbstractPetscVec}) 
Computes `y = y + sum alpha[i] x[i]`

Logically Collective

Input Parameters:
- `nv`    - number of scalars and `x` vectors
- `alpha` - array of scalars
- `y`     - one vector
- `x`     - array of vectors

Level: intermediate

See also: `Vec`, `VecMAXPBY()`, `VecAYPX()`, `VecWAXPY()`, `VecAXPY()`, `VecAXPBYPCZ()`, `VecAXPBY()`, `VecDuplicateVecs()`

# External Links
$(_doc_external("Vec/VecMAXPY"))
"""
function VecMAXPY(petsclib::PetscLibType, y::AbstractPetscVec, nv::Integer, alpha::AbstractVector{<:Number}, x::Vector{<:AbstractPetscVec})
    error("VecMAXPY: no generated method for these argument types")
end

@for_petsc function VecMAXPY(petsclib::$UnionPetscLib, y::AbstractPetscVec, nv::$PetscInt, alpha::Vector{$PetscScalar}, x::Vector{<:AbstractPetscVec} )

    @chk ccall(
               (:VecMAXPY, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{$PetscScalar}, Ptr{CVec}),
               y, nv, alpha, x,
              )


	return nothing
end 

"""
	VecMDot(petsclib::PetscLibType, x::AbstractPetscVec, nv::PetscInt, y::Vector{<:AbstractPetscVec}, val::Vector{PetscScalar}) 
Computes multiple vector dot products.

Collective

Input Parameters:
- `x`  - one vector
- `nv` - number of vectors
- `y`  - array of vectors.

Output Parameter:
- `val` - array of the dot products (does not allocate the array)

Level: intermediate

Notes for Users of Complex Numbers:
For complex vectors, `VecMDot()` computes
``
val = (x,y) = y^H x,
``
where y^H denotes the conjugate transpose of y.

Use `VecMTDot()` for the indefinite form
``
val = (x,y) = y^T x,
``
where y^T denotes the transpose of y.

See also: `Vec`, `VecMTDot()`, `VecDot()`, `VecDuplicateVecs()`

# External Links
$(_doc_external("Vec/VecMDot"))
"""
function VecMDot(petsclib::PetscLibType, x::AbstractPetscVec, nv::Integer, y::Vector{<:AbstractPetscVec}, val::AbstractVector{<:Number})
    error("VecMDot: no generated method for these argument types")
end

@for_petsc function VecMDot(petsclib::$UnionPetscLib, x::AbstractPetscVec, nv::$PetscInt, y::Vector{<:AbstractPetscVec}, val::Vector{$PetscScalar} )

    @chk ccall(
               (:VecMDot, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{CVec}, Ptr{$PetscScalar}),
               x, nv, y, val,
              )


	return nothing
end 

"""
	VecMDotBegin(petsclib::PetscLibType, x::AbstractPetscVec, nv::PetscInt, y::Vector{<:AbstractPetscVec}, result::Vector{PetscScalar}) 
Starts a split phase multiple dot product computation.

Input Parameters:
- `x`      - the first vector
- `nv`     - number of vectors
- `y`      - array of vectors
- `result` - where the result will go (can be `NULL`)

Level: advanced

See also: `VecMDotEnd()`, `VecNormBegin()`, `VecNormEnd()`, `VecNorm()`, `VecDot()`, `VecMDot()`,
`VecTDotBegin()`, `VecTDotEnd()`, `VecMTDotBegin()`, `VecMTDotEnd()`, `PetscCommSplitReductionBegin()`

# External Links
$(_doc_external("Vec/VecMDotBegin"))
"""
function VecMDotBegin(petsclib::PetscLibType, x::AbstractPetscVec, nv::Integer, y::Vector{<:AbstractPetscVec}, result::AbstractVector{<:Number})
    error("VecMDotBegin: no generated method for these argument types")
end

@for_petsc function VecMDotBegin(petsclib::$UnionPetscLib, x::AbstractPetscVec, nv::$PetscInt, y::Vector{<:AbstractPetscVec}, result::Vector{$PetscScalar} )

    @chk ccall(
               (:VecMDotBegin, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{CVec}, Ptr{$PetscScalar}),
               x, nv, y, result,
              )


	return nothing
end 

"""
	VecMDotEnd(petsclib::PetscLibType, x::AbstractPetscVec, nv::PetscInt, y::Vector{<:AbstractPetscVec}, result::Vector{PetscScalar}) 
Ends a split phase multiple dot product computation.

Input Parameters:
- `x`  - the first vector (can be `NULL`)
- `nv` - number of vectors
- `y`  - array of vectors (can be `NULL`)

Output Parameter:
- `result` - where the result will go

Level: advanced

See also: `VecMDotBegin()`, `VecNormBegin()`, `VecNormEnd()`, `VecNorm()`, `VecDot()`, `VecMDot()`,
`VecTDotBegin()`, `VecTDotEnd()`, `VecMTDotBegin()`, `VecMTDotEnd()`, `PetscCommSplitReductionBegin()`

# External Links
$(_doc_external("Vec/VecMDotEnd"))
"""
function VecMDotEnd(petsclib::PetscLibType, x::AbstractPetscVec, nv::Integer, y::Vector{<:AbstractPetscVec}, result::AbstractVector{<:Number})
    error("VecMDotEnd: no generated method for these argument types")
end

@for_petsc function VecMDotEnd(petsclib::$UnionPetscLib, x::AbstractPetscVec, nv::$PetscInt, y::Vector{<:AbstractPetscVec}, result::Vector{$PetscScalar} )

    @chk ccall(
               (:VecMDotEnd, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{CVec}, Ptr{$PetscScalar}),
               x, nv, y, result,
              )


	return nothing
end 

"""
	VecMPISetGhost(petsclib::PetscLibType, vv::AbstractPetscVec, nghost::PetscInt, ghosts::Vector{PetscInt}) 
Sets the ghost points for an MPI ghost vector

Collective

Input Parameters:
- `vv`     - the MPI vector
- `nghost` - number of local ghost points
- `ghosts` - global indices of ghost points, these do not need to be in increasing order (sorted)

Level: advanced

See also: `Vec`, `VecType`, `VecCreateSeq()`, `VecCreate()`, `VecDuplicate()`, `VecDuplicateVecs()`, `VecCreateMPI()`,
`VecGhostGetLocalForm()`, `VecGhostRestoreLocalForm()`, `VecGhostUpdateBegin()`,
`VecCreateGhostWithArray()`, `VecCreateMPIWithArray()`, `VecGhostUpdateEnd()`,
`VecCreateGhostBlock()`, `VecCreateGhostBlockWithArray()`

# External Links
$(_doc_external("Vec/VecMPISetGhost"))
"""
function VecMPISetGhost(petsclib::PetscLibType, vv::AbstractPetscVec, nghost::Integer, ghosts::AbstractVector{<:Number})
    error("VecMPISetGhost: no generated method for these argument types")
end

@for_petsc function VecMPISetGhost(petsclib::$UnionPetscLib, vv::AbstractPetscVec, nghost::$PetscInt, ghosts::Vector{$PetscInt} )

    @chk ccall(
               (:VecMPISetGhost, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{$PetscInt}),
               vv, nghost, ghosts,
              )


	return nothing
end 

"""
	VecMTDot(petsclib::PetscLibType, x::AbstractPetscVec, nv::PetscInt, y::Vector{<:AbstractPetscVec}, val::Vector{PetscScalar}) 
Computes indefinite vector multiple dot products.
That is, it does NOT use the complex conjugate.

Collective

Input Parameters:
- `x`  - one vector
- `nv` - number of vectors
- `y`  - array of vectors.  Note that vectors are pointers

Output Parameter:
- `val` - array of the dot products

Level: intermediate

Notes for Users of Complex Numbers:
For complex vectors, `VecMTDot()` computes the indefinite form
``
val = (x,y) = y^T x,
``
where y^T denotes the transpose of y.

Use `VecMDot()` for the inner product
``
val = (x,y) = y^H x,
``
where y^H denotes the conjugate transpose of y.

See also: `Vec`, `VecMDot()`, `VecTDot()`

# External Links
$(_doc_external("Vec/VecMTDot"))
"""
function VecMTDot(petsclib::PetscLibType, x::AbstractPetscVec, nv::Integer, y::Vector{<:AbstractPetscVec}, val::AbstractVector{<:Number})
    error("VecMTDot: no generated method for these argument types")
end

@for_petsc function VecMTDot(petsclib::$UnionPetscLib, x::AbstractPetscVec, nv::$PetscInt, y::Vector{<:AbstractPetscVec}, val::Vector{$PetscScalar} )

    @chk ccall(
               (:VecMTDot, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{CVec}, Ptr{$PetscScalar}),
               x, nv, y, val,
              )


	return nothing
end 

"""
	VecMTDotBegin(petsclib::PetscLibType, x::AbstractPetscVec, nv::PetscInt, y::Vector{<:AbstractPetscVec}, result::Vector{PetscScalar}) 
Starts a split phase transpose multiple dot product computation.

Input Parameters:
- `x`      - the first vector
- `nv`     - number of vectors
- `y`      - array of  vectors
- `result` - where the result will go (can be `NULL`)

Level: advanced

See also: `VecMTDotEnd()`, `VecNormBegin()`, `VecNormEnd()`, `VecNorm()`, `VecDot()`, `VecMDot()`,
`VecDotBegin()`, `VecDotEnd()`, `VecMDotBegin()`, `VecMDotEnd()`, `PetscCommSplitReductionBegin()`

# External Links
$(_doc_external("Vec/VecMTDotBegin"))
"""
function VecMTDotBegin(petsclib::PetscLibType, x::AbstractPetscVec, nv::Integer, y::Vector{<:AbstractPetscVec}, result::AbstractVector{<:Number})
    error("VecMTDotBegin: no generated method for these argument types")
end

@for_petsc function VecMTDotBegin(petsclib::$UnionPetscLib, x::AbstractPetscVec, nv::$PetscInt, y::Vector{<:AbstractPetscVec}, result::Vector{$PetscScalar} )

    @chk ccall(
               (:VecMTDotBegin, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{CVec}, Ptr{$PetscScalar}),
               x, nv, y, result,
              )


	return nothing
end 

"""
	result::Vector{PetscScalar} = VecMTDotEnd(petsclib::PetscLibType, x::AbstractPetscVec, nv::PetscInt, y::Vector{<:AbstractPetscVec}) 
Ends a split phase transpose multiple dot product computation.

Input Parameters:
- `x`  - the first vector (can be `NULL`)
- `nv` - number of vectors
- `y`  - array of  vectors (can be `NULL`)

Output Parameter:
- `result` - where the result will go

Level: advanced

See also: `VecMTDotBegin()`, `VecNormBegin()`, `VecNormEnd()`, `VecNorm()`, `VecDot()`, `VecMDot()`,
`VecDotBegin()`, `VecDotEnd()`, `VecMDotBegin()`, `VecMDotEnd()`, `PetscCommSplitReductionBegin()`

# External Links
$(_doc_external("Vec/VecMTDotEnd"))
"""
function VecMTDotEnd(petsclib::PetscLibType, x::AbstractPetscVec, nv::Integer, y::Vector{<:AbstractPetscVec})
    error("VecMTDotEnd: no generated method for these argument types")
end

@for_petsc function VecMTDotEnd(petsclib::$UnionPetscLib, x::AbstractPetscVec, nv::$PetscInt, y::Vector{<:AbstractPetscVec} )
	result = Vector{$PetscScalar}(undef, nv)

    @chk ccall(
               (:VecMTDotEnd, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{CVec}, Ptr{$PetscScalar}),
               x, nv, y, result,
              )


	return result
end 

"""
	p::PetscInt,val::PetscReal = VecMax(petsclib::PetscLibType, x::AbstractPetscVec) 
Determines the vector component with maximum real part and its location.

Collective

Input Parameter:
- `x` - the vector

Output Parameters:
- `p`   - the index of `val` (pass `NULL` if you don't want this) in the vector
- `val` - the maximum component

Level: intermediate

See also: `Vec`, `VecNorm()`, `VecMin()`

# External Links
$(_doc_external("Vec/VecMax"))
"""
function VecMax(petsclib::PetscLibType, x::AbstractPetscVec)
    error("VecMax: no generated method for these argument types")
end

@for_petsc function VecMax(petsclib::$UnionPetscLib, x::AbstractPetscVec )
	p_ = Ref{$PetscInt}()
	val_ = Ref{$PetscReal}()

    @chk ccall(
               (:VecMax, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{$PetscInt}, Ptr{$PetscReal}),
               x, p_, val_,
              )

	p = p_[]
	val = val_[]

	return p,val
end 

"""
	max::PetscReal = VecMaxPointwiseDivide(petsclib::PetscLibType, x::AbstractPetscVec, y::AbstractPetscVec) 
Computes the maximum of the componentwise division `max = max_i abs(x[i]/y[i])`.

Logically Collective

Input Parameters:
- `x` - the numerators
- `y` - the denominators

Output Parameter:
- `max` - the result

Level: advanced

See also: `Vec`, `VecPointwiseDivide()`, `VecPointwiseMult()`, `VecPointwiseMax()`, `VecPointwiseMin()`, `VecPointwiseMaxAbs()`

# External Links
$(_doc_external("Vec/VecMaxPointwiseDivide"))
"""
function VecMaxPointwiseDivide(petsclib::PetscLibType, x::AbstractPetscVec, y::AbstractPetscVec)
    error("VecMaxPointwiseDivide: no generated method for these argument types")
end

@for_petsc function VecMaxPointwiseDivide(petsclib::$UnionPetscLib, x::AbstractPetscVec, y::AbstractPetscVec )
	max_ = Ref{$PetscReal}()

    @chk ccall(
               (:VecMaxPointwiseDivide, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, Ptr{$PetscReal}),
               x, y, max_,
              )

	max = max_[]

	return max
end 

"""
	mean::PetscScalar = VecMean(petsclib::PetscLibType, v::AbstractPetscVec) 
Computes the arithmetic mean of all the components of a vector.

Collective

Input Parameter:
- `v` - the vector

Output Parameter:
- `mean` - the result

Level: beginner

See also: `Vec`, `VecSum()`, `VecNorm()`

# External Links
$(_doc_external("Vec/VecMean"))
"""
function VecMean(petsclib::PetscLibType, v::AbstractPetscVec)
    error("VecMean: no generated method for these argument types")
end

@for_petsc function VecMean(petsclib::$UnionPetscLib, v::AbstractPetscVec )
	mean_ = Ref{$PetscScalar}()

    @chk ccall(
               (:VecMean, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{$PetscScalar}),
               v, mean_,
              )

	mean = mean_[]

	return mean
end 

"""
	VecMedian(petsclib::PetscLibType, Vec1::AbstractPetscVec, Vec2::AbstractPetscVec, Vec3::AbstractPetscVec, VMedian::AbstractPetscVec) 
Computes the componentwise median of three vectors
and stores the result in this vector.  Used primarily for projecting
a vector within upper and lower bounds.

Logically Collective

Input Parameters:
- `Vec1` - The first vector
- `Vec2` - The second vector
- `Vec3` - The third vector

Output Parameter:
- `VMedian` - The median vector (this can be any one of the input vectors)

Level: advanced

See also: `Vec`

# External Links
$(_doc_external("Vec/VecMedian"))
"""
function VecMedian(petsclib::PetscLibType, Vec1::AbstractPetscVec, Vec2::AbstractPetscVec, Vec3::AbstractPetscVec, VMedian::AbstractPetscVec)
    error("VecMedian: no generated method for these argument types")
end

@for_petsc function VecMedian(petsclib::$UnionPetscLib, Vec1::AbstractPetscVec, Vec2::AbstractPetscVec, Vec3::AbstractPetscVec, VMedian::AbstractPetscVec )

    @chk ccall(
               (:VecMedian, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, CVec, CVec),
               Vec1, Vec2, Vec3, VMedian,
              )


	return nothing
end 

"""
	p::PetscInt,val::PetscReal = VecMin(petsclib::PetscLibType, x::AbstractPetscVec) 
Determines the vector component with minimum real part and its location.

Collective

Input Parameter:
- `x` - the vector

Output Parameters:
- `p`   - the index of `val` (pass `NULL` if you don't want this location) in the vector
- `val` - the minimum component

Level: intermediate

See also: `Vec`, `VecMax()`

# External Links
$(_doc_external("Vec/VecMin"))
"""
function VecMin(petsclib::PetscLibType, x::AbstractPetscVec)
    error("VecMin: no generated method for these argument types")
end

@for_petsc function VecMin(petsclib::$UnionPetscLib, x::AbstractPetscVec )
	p_ = Ref{$PetscInt}()
	val_ = Ref{$PetscReal}()

    @chk ccall(
               (:VecMin, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{$PetscInt}, Ptr{$PetscReal}),
               x, p_, val_,
              )

	p = p_[]
	val = val_[]

	return p,val
end 

"""
	N::PetscInt = VecNestGetSize(petsclib::PetscLibType, X::AbstractPetscVec) 
Returns the size of the nest vector.

Not Collective

Input Parameter:
- `X` - nest vector

Output Parameter:
- `N` - number of nested vecs

Level: developer

See also: `VECNEST`, `Vec`, `VecType`, `VecNestGetSubVec()`, `VecNestGetSubVecs()`

# External Links
$(_doc_external("Vec/VecNestGetSize"))
"""
function VecNestGetSize(petsclib::PetscLibType, X::AbstractPetscVec)
    error("VecNestGetSize: no generated method for these argument types")
end

@for_petsc function VecNestGetSize(petsclib::$UnionPetscLib, X::AbstractPetscVec )
	N_ = Ref{$PetscInt}()

    @chk ccall(
               (:VecNestGetSize, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{$PetscInt}),
               X, N_,
              )

	N = N_[]

	return N
end 

"""
	sx::PetscVec = VecNestGetSubVec(petsclib::PetscLibType, X::AbstractPetscVec, idxm::PetscInt) 
Returns a single, sub-vector from a nest vector.

Not Collective

Input Parameters:
- `X`    - nest vector
- `idxm` - index of the vector within the nest

Output Parameter:
- `sx` - vector at index `idxm` within the nest

Level: developer

See also: `VECNEST`, `Vec`, `VecType`, `VecNestGetSize()`, `VecNestGetSubVecs()`

# External Links
$(_doc_external("Vec/VecNestGetSubVec"))
"""
function VecNestGetSubVec(petsclib::PetscLibType, X::AbstractPetscVec, idxm::Integer)
    error("VecNestGetSubVec: no generated method for these argument types")
end

@for_petsc function VecNestGetSubVec(petsclib::$UnionPetscLib, X::AbstractPetscVec, idxm::$PetscInt )
	sx_ = Ref{CVec}()

    @chk ccall(
               (:VecNestGetSubVec, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{CVec}),
               X, idxm, sx_,
              )

	sx = PetscVec(sx_[], petsclib)

	return sx
end 

"""
	N::PetscInt,sx::Vector{PetscVec} = VecNestGetSubVecs(petsclib::PetscLibType, X::AbstractPetscVec) 
Returns the entire array of vectors defining a nest vector.

Not Collective

Input Parameter:
- `X` - nest vector

Output Parameters:
- `N`  - number of nested vecs
- `sx` - array of vectors, can pass in `NULL`

Level: developer

See also: `VECNEST`, `Vec`, `VecType`, `VecNestGetSize()`, `VecNestGetSubVec()`, `VecNestGetSubVecsRead()`

# External Links
$(_doc_external("Vec/VecNestGetSubVecs"))
"""
function VecNestGetSubVecs(petsclib::PetscLibType, X::AbstractPetscVec)
    error("VecNestGetSubVecs: no generated method for these argument types")
end

@for_petsc function VecNestGetSubVecs(petsclib::$UnionPetscLib, X::AbstractPetscVec )
	N_ = Ref{$PetscInt}()
	sx_ = Ref{Ptr{CVec}}()

    @chk ccall(
               (:VecNestGetSubVecs, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{$PetscInt}, Ptr{Ptr{CVec}}),
               X, N_, sx_,
              )

	N = N_[]
	sx = sx_[] == C_NULL ? PetscVec{$PetscLib}[] : [PetscVec(p, petsclib) for p in unsafe_wrap(Array, sx_[], N; own = false)]

	return N,sx
end 

"""
	N::PetscInt,sx::Vector{PetscVec} = VecNestGetSubVecsRead(petsclib::PetscLibType, X::AbstractPetscVec) 
Access the subvecs of a `VECNEST` vector for read-only access

Logically collective

Input Parameter:
- `X` - nest vector

Output Parameters:
- `N`  - number of nested vecs
- `sx` - array of read-locked vectors

Level: advanced

See also: `VECNEST`, `Vec`, `VecType`, `VecNestGetSize()`, `VecNestGetSubVec()`, `VecNestRestoreSubVecsRead()`

# External Links
$(_doc_external("Vec/VecNestGetSubVecsRead"))
"""
function VecNestGetSubVecsRead(petsclib::PetscLibType, X::AbstractPetscVec)
    error("VecNestGetSubVecsRead: no generated method for these argument types")
end

@for_petsc function VecNestGetSubVecsRead(petsclib::$UnionPetscLib, X::AbstractPetscVec )
	N_ = Ref{$PetscInt}()
	sx_ = Ref{Ptr{CVec}}()

    @chk ccall(
               (:VecNestGetSubVecsRead, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{$PetscInt}, Ptr{Ptr{CVec}}),
               X, N_, sx_,
              )

	N = N_[]
	sx = sx_[] == C_NULL ? PetscVec{$PetscLib}[] : [PetscVec(p, petsclib) for p in unsafe_wrap(Array, sx_[], N; own = false)]

	return N,sx
end 

"""
	subparams::PetscVec = VecNestGetTaoTermSumParameters(petsclib::PetscLibType, params::AbstractPetscVec, index::PetscInt) 
A wrapper around `VecNestGetSubVec()` for `TAOTERMSUM`.

Not collective

Input Parameters:
- `params` - a `VECNEST` that has one nested vector for each term of a `TAOTERMSUM`
- `index`  - the index of a term

Output Parameter:
- `subparams` - the parameters of the internal terms of `TAOTERMSUM`. (may be `NULL`)

Level: intermediate

See also: [](sec_tao_term),
`TaoTerm`,
`TAOTERMSUM`,
`TaoTermSumParametersPack()`,
`TaoTermSumParametersUnpack()`,
`VECNEST`,
`VecNestGetSubVec()`

# External Links
$(_doc_external("TaoTerm/VecNestGetTaoTermSumParameters"))
"""
function VecNestGetTaoTermSumParameters(petsclib::PetscLibType, params::AbstractPetscVec, index::Integer)
    error("VecNestGetTaoTermSumParameters: no generated method for these argument types")
end

@for_petsc function VecNestGetTaoTermSumParameters(petsclib::$UnionPetscLib, params::AbstractPetscVec, index::$PetscInt )
	subparams_ = Ref{CVec}()

    @chk ccall(
               (:VecNestGetTaoTermSumParameters, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{CVec}),
               params, index, subparams_,
              )

	subparams = PetscVec(subparams_[], petsclib)

	return subparams
end 

"""
	VecNestRestoreSubVecsRead(petsclib::PetscLibType, X::AbstractPetscVec, N::PetscInt, sx::Union{Ptr, AbstractArray{PetscVec}}) 
Restore access the subvecs of a `VECNEST` vector obtained with `VecNestGetSubVecsRead()`

Logically collective

Input Parameters:
- `X`  - nest vector
- `N`  - number of nested vecs
- `sx` - array of read-locked vectors

Level: advanced

See also: `VECNEST`, `Vec`, `VecType`, `VecNestGetSize()`, `VecNestGetSubVec()`, `VecNestGetSubVecsRead()`

# External Links
$(_doc_external("Vec/VecNestRestoreSubVecsRead"))
"""
function VecNestRestoreSubVecsRead(petsclib::PetscLibType, X::AbstractPetscVec, N::Integer, sx::Union{Ptr, AbstractArray{PetscVec}})
    error("VecNestRestoreSubVecsRead: no generated method for these argument types")
end

@for_petsc function VecNestRestoreSubVecsRead(petsclib::$UnionPetscLib, X::AbstractPetscVec, N::$PetscInt, sx::Union{Ptr, AbstractArray{PetscVec}} )
	N_ = Ref{$PetscInt}(N)
	sx_ = Ref{Ptr{CVec}}(sx isa Ptr ? sx : pointer(sx))

    @chk ccall(
               (:VecNestRestoreSubVecsRead, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{$PetscInt}, Ptr{Ptr{CVec}}),
               X, N_, sx_,
              )


	return nothing
end 

"""
	VecNestSetSubVec(petsclib::PetscLibType, X::AbstractPetscVec, idxm::PetscInt, sx::AbstractPetscVec) 
Set a single component vector in a nest vector at specified index.

Not Collective

Input Parameters:
- `X`    - nest vector
- `idxm` - index of the vector within the nest vector
- `sx`   - vector at index `idxm` within the nest vector

Level: developer

See also: `VECNEST`, `Vec`, `VecType`, `VecNestSetSubVecs()`, `VecNestGetSubVec()`

# External Links
$(_doc_external("Vec/VecNestSetSubVec"))
"""
function VecNestSetSubVec(petsclib::PetscLibType, X::AbstractPetscVec, idxm::Integer, sx::AbstractPetscVec)
    error("VecNestSetSubVec: no generated method for these argument types")
end

@for_petsc function VecNestSetSubVec(petsclib::$UnionPetscLib, X::AbstractPetscVec, idxm::$PetscInt, sx::AbstractPetscVec )

    @chk ccall(
               (:VecNestSetSubVec, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, CVec),
               X, idxm, sx,
              )


	return nothing
end 

"""
	VecNestSetSubVecs(petsclib::PetscLibType, X::AbstractPetscVec, N::PetscInt, idxm::Vector{PetscInt}, sx::Vector{<:AbstractPetscVec}) 
Sets the component vectors at the specified indices in a nest vector.

Not Collective

Input Parameters:
- `X`    - nest vector
- `N`    - number of component vecs in `sx`
- `idxm` - indices of component vectors that are to be replaced
- `sx`   - array of vectors

Level: developer

See also: `VECNEST`, `Vec`, `VecType`, `VecNestGetSize()`, `VecNestGetSubVec()`

# External Links
$(_doc_external("Vec/VecNestSetSubVecs"))
"""
function VecNestSetSubVecs(petsclib::PetscLibType, X::AbstractPetscVec, N::Integer, idxm::AbstractVector{<:Number}, sx::Vector{<:AbstractPetscVec})
    error("VecNestSetSubVecs: no generated method for these argument types")
end

@for_petsc function VecNestSetSubVecs(petsclib::$UnionPetscLib, X::AbstractPetscVec, N::$PetscInt, idxm::Vector{$PetscInt}, sx::Vector{<:AbstractPetscVec} )

    @chk ccall(
               (:VecNestSetSubVecs, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{$PetscInt}, Ptr{CVec}),
               X, N, idxm, sx,
              )


	return nothing
end 

"""
	val::PetscReal = VecNorm(petsclib::PetscLibType, x::AbstractPetscVec, type::NormType) 
Computes the vector norm.

Collective

Input Parameters:
- `x`    - the vector
- `type` - the type of the norm requested

Output Parameter:
- `val` - the norm

Level: intermediate

See also: `Vec`, `NormType`, `VecDot()`, `VecTDot()`, `VecDotBegin()`, `VecDotEnd()`, `VecNormAvailable()`,
`VecNormBegin()`, `VecNormEnd()`, `NormType()`

# External Links
$(_doc_external("Vec/VecNorm"))
"""
function VecNorm(petsclib::PetscLibType, x::AbstractPetscVec, type::NormType)
    error("VecNorm: no generated method for these argument types")
end

@for_petsc function VecNorm(petsclib::$UnionPetscLib, x::AbstractPetscVec, type::NormType )
	val_ = Ref{$PetscReal}()

    @chk ccall(
               (:VecNorm, $petsc_library),
               PetscErrorCode,
               (CVec, NormType, Ptr{$PetscReal}),
               x, type, val_,
              )

	val = val_[]

	return val
end 

"""
	available::PetscBool,val::PetscReal = VecNormAvailable(petsclib::PetscLibType, x::AbstractPetscVec, type::NormType) 
Returns the vector norm if it is already known. That is, it has been previously computed and cached in the vector

Not Collective

Input Parameters:
- `x`    - the vector
- `type` - one of `NORM_1` (sum_i |x[i]|), `NORM_2` sqrt(sum_i (x[i])^2), `NORM_INFINITY` max_i |x[i]|.  Also available
`NORM_1_AND_2`, which computes both norms and stores them
in a two element array.

Output Parameters:
- `available` - `PETSC_TRUE` if the val returned is valid
- `val`       - the norm

Level: intermediate

See also: `Vec`, `VecDot()`, `VecTDot()`, `VecNorm()`, `VecDotBegin()`, `VecDotEnd()`,
`VecNormBegin()`, `VecNormEnd()`

# External Links
$(_doc_external("Vec/VecNormAvailable"))
"""
function VecNormAvailable(petsclib::PetscLibType, x::AbstractPetscVec, type::NormType)
    error("VecNormAvailable: no generated method for these argument types")
end

@for_petsc function VecNormAvailable(petsclib::$UnionPetscLib, x::AbstractPetscVec, type::NormType )
	available_ = Ref{PetscBool}()
	val_ = Ref{$PetscReal}()

    @chk ccall(
               (:VecNormAvailable, $petsc_library),
               PetscErrorCode,
               (CVec, NormType, Ptr{PetscBool}, Ptr{$PetscReal}),
               x, type, available_, val_,
              )

	available = available_[]
	val = val_[]

	return available,val
end 

"""
	result::PetscReal = VecNormBegin(petsclib::PetscLibType, x::AbstractPetscVec, ntype::NormType) 
Starts a split phase norm computation.

Input Parameters:
- `x`      - the first vector
- `ntype`  - norm type, one of `NORM_1`, `NORM_2`, `NORM_MAX`, `NORM_1_AND_2`
- `result` - where the result will go (can be `NULL`)

Level: advanced

See also: `VecNormEnd()`, `VecNorm()`, `VecDot()`, `VecMDot()`, `VecDotBegin()`, `VecDotEnd()`, `PetscCommSplitReductionBegin()`

# External Links
$(_doc_external("Vec/VecNormBegin"))
"""
function VecNormBegin(petsclib::PetscLibType, x::AbstractPetscVec, ntype::NormType)
    error("VecNormBegin: no generated method for these argument types")
end

@for_petsc function VecNormBegin(petsclib::$UnionPetscLib, x::AbstractPetscVec, ntype::NormType )
	result_ = Ref{$PetscReal}()

    @chk ccall(
               (:VecNormBegin, $petsc_library),
               PetscErrorCode,
               (CVec, NormType, Ptr{$PetscReal}),
               x, ntype, result_,
              )

	result = result_[]

	return result
end 

"""
	result::PetscReal = VecNormEnd(petsclib::PetscLibType, x::AbstractPetscVec, ntype::NormType) 
Ends a split phase norm computation.

Input Parameters:
- `x`      - the first vector
- `ntype`  - norm type, one of `NORM_1`, `NORM_2`, `NORM_MAX`, `NORM_1_AND_2`
- `result` - where the result will go

Level: advanced

See also: `VecNormBegin()`, `VecNorm()`, `VecDot()`, `VecMDot()`, `VecDotBegin()`, `VecDotEnd()`, `PetscCommSplitReductionBegin()`

# External Links
$(_doc_external("Vec/VecNormEnd"))
"""
function VecNormEnd(petsclib::PetscLibType, x::AbstractPetscVec, ntype::NormType)
    error("VecNormEnd: no generated method for these argument types")
end

@for_petsc function VecNormEnd(petsclib::$UnionPetscLib, x::AbstractPetscVec, ntype::NormType )
	result_ = Ref{$PetscReal}()

    @chk ccall(
               (:VecNormEnd, $petsc_library),
               PetscErrorCode,
               (CVec, NormType, Ptr{$PetscReal}),
               x, ntype, result_,
              )

	result = result_[]

	return result
end 

"""
	val::PetscReal = VecNormalize(petsclib::PetscLibType, x::AbstractPetscVec) 
Normalizes a vector by its 2-norm.

Collective

Input Parameter:
- `x` - the vector

Output Parameter:
- `val` - the vector norm before normalization. May be `NULL` if the value is not needed.

Level: intermediate

See also: `Vec`, `VecNorm()`, `NORM_2`, `NormType`

# External Links
$(_doc_external("Vec/VecNormalize"))
"""
function VecNormalize(petsclib::PetscLibType, x::AbstractPetscVec)
    error("VecNormalize: no generated method for these argument types")
end

@for_petsc function VecNormalize(petsclib::$UnionPetscLib, x::AbstractPetscVec )
	val_ = Ref{$PetscReal}()

    @chk ccall(
               (:VecNormalize, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{$PetscReal}),
               x, val_,
              )

	val = val_[]

	return val
end 

"""
	VecPermute(petsclib::PetscLibType, x::AbstractPetscVec, row::AbstractIS, inv::PetscBool) 
Permutes a vector in place using the given ordering.

Input Parameters:
- `x`   - The vector
- `row` - The ordering
- `inv` - The flag for inverting the permutation

Level: beginner

See also: `Vec`, `MatPermute()`

# External Links
$(_doc_external("Vec/VecPermute"))
"""
function VecPermute(petsclib::PetscLibType, x::AbstractPetscVec, row::AbstractIS, inv::PetscBool)
    error("VecPermute: no generated method for these argument types")
end

@for_petsc function VecPermute(petsclib::$UnionPetscLib, x::AbstractPetscVec, row::AbstractIS, inv::PetscBool )

    @chk ccall(
               (:VecPermute, $petsc_library),
               PetscErrorCode,
               (CVec, CIS, PetscBool),
               x, row, inv,
              )


	return nothing
end 

"""
	VecPlaceArray(petsclib::PetscLibType, vec::AbstractPetscVec, array::Vector{PetscScalar}) 
Allows one to replace the array in a vector with an
array provided by the user. This is useful to avoid copying an array
into a vector.

Logically Collective

Input Parameters:
- `vec`   - the vector
- `array` - the array

Level: developer

See also: `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecReplaceArray()`, `VecResetArray()`

# External Links
$(_doc_external("Vec/VecPlaceArray"))
"""
function VecPlaceArray(petsclib::PetscLibType, vec::AbstractPetscVec, array::AbstractVector{<:Number})
    error("VecPlaceArray: no generated method for these argument types")
end

@for_petsc function VecPlaceArray(petsclib::$UnionPetscLib, vec::AbstractPetscVec, array::Vector{$PetscScalar} )

    @chk ccall(
               (:VecPlaceArray, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{$PetscScalar}),
               vec, array,
              )


	return nothing
end 

"""
	VecPointwiseDivide(petsclib::PetscLibType, w::AbstractPetscVec, x::AbstractPetscVec, y::AbstractPetscVec) 
Computes the component-wise division `w[i] = x[i] / y[i]`.

Logically Collective

Input Parameters:
- `x` - the numerator vector
- `y` - the denominator vector

Output Parameter:
- `w` - the result

Level: advanced

See also: `Vec`, `VecPointwiseMult()`, `VecPointwiseMax()`, `VecPointwiseMin()`, `VecPointwiseMaxAbs()`, `VecMaxPointwiseDivide()`

# External Links
$(_doc_external("Vec/VecPointwiseDivide"))
"""
function VecPointwiseDivide(petsclib::PetscLibType, w::AbstractPetscVec, x::AbstractPetscVec, y::AbstractPetscVec)
    error("VecPointwiseDivide: no generated method for these argument types")
end

@for_petsc function VecPointwiseDivide(petsclib::$UnionPetscLib, w::AbstractPetscVec, x::AbstractPetscVec, y::AbstractPetscVec )

    @chk ccall(
               (:VecPointwiseDivide, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, CVec),
               w, x, y,
              )


	return nothing
end 

"""
	VecPointwiseMax(petsclib::PetscLibType, w::AbstractPetscVec, x::AbstractPetscVec, y::AbstractPetscVec) 
Computes the component-wise maximum `w[i] = max(x[i], y[i])`.

Logically Collective

Input Parameters:
- `x` - the first input vector
- `y` - the second input vector

Output Parameter:
- `w` - the result

Level: advanced

See also: `Vec`, `VecPointwiseDivide()`, `VecPointwiseMult()`, `VecPointwiseMin()`, `VecPointwiseMaxAbs()`, `VecMaxPointwiseDivide()`

# External Links
$(_doc_external("Vec/VecPointwiseMax"))
"""
function VecPointwiseMax(petsclib::PetscLibType, w::AbstractPetscVec, x::AbstractPetscVec, y::AbstractPetscVec)
    error("VecPointwiseMax: no generated method for these argument types")
end

@for_petsc function VecPointwiseMax(petsclib::$UnionPetscLib, w::AbstractPetscVec, x::AbstractPetscVec, y::AbstractPetscVec )

    @chk ccall(
               (:VecPointwiseMax, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, CVec),
               w, x, y,
              )


	return nothing
end 

"""
	VecPointwiseMaxAbs(petsclib::PetscLibType, w::AbstractPetscVec, x::AbstractPetscVec, y::AbstractPetscVec) 
Computes the component-wise maximum of the absolute values `w[i] = max(abs(x[i]), abs(y[i]))`.

Logically Collective

Input Parameters:
- `x` - the first input vector
- `y` - the second input vector

Output Parameter:
- `w` - the result

Level: advanced

See also: `Vec`, `VecPointwiseDivide()`, `VecPointwiseMult()`, `VecPointwiseMin()`, `VecPointwiseMax()`, `VecMaxPointwiseDivide()`

# External Links
$(_doc_external("Vec/VecPointwiseMaxAbs"))
"""
function VecPointwiseMaxAbs(petsclib::PetscLibType, w::AbstractPetscVec, x::AbstractPetscVec, y::AbstractPetscVec)
    error("VecPointwiseMaxAbs: no generated method for these argument types")
end

@for_petsc function VecPointwiseMaxAbs(petsclib::$UnionPetscLib, w::AbstractPetscVec, x::AbstractPetscVec, y::AbstractPetscVec )

    @chk ccall(
               (:VecPointwiseMaxAbs, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, CVec),
               w, x, y,
              )


	return nothing
end 

"""
	VecPointwiseMin(petsclib::PetscLibType, w::AbstractPetscVec, x::AbstractPetscVec, y::AbstractPetscVec) 
Computes the component-wise minimum `w[i] = min(x[i], y[i])`.

Logically Collective

Input Parameters:
- `x` - the first input vector
- `y` - the second input vector

Output Parameter:
- `w` - the result

Level: advanced

See also: `Vec`, `VecPointwiseDivide()`, `VecPointwiseMult()`, `VecPointwiseMaxAbs()`, `VecMaxPointwiseDivide()`

# External Links
$(_doc_external("Vec/VecPointwiseMin"))
"""
function VecPointwiseMin(petsclib::PetscLibType, w::AbstractPetscVec, x::AbstractPetscVec, y::AbstractPetscVec)
    error("VecPointwiseMin: no generated method for these argument types")
end

@for_petsc function VecPointwiseMin(petsclib::$UnionPetscLib, w::AbstractPetscVec, x::AbstractPetscVec, y::AbstractPetscVec )

    @chk ccall(
               (:VecPointwiseMin, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, CVec),
               w, x, y,
              )


	return nothing
end 

"""
	VecPointwiseMult(petsclib::PetscLibType, w::AbstractPetscVec, x::AbstractPetscVec, y::AbstractPetscVec) 
Computes the component-wise multiplication `w[i] = x[i] * y[i]`.

Logically Collective

Input Parameters:
- `x` - the first vector
- `y` - the second vector

Output Parameter:
- `w` - the result

Level: advanced

See also: `Vec`, `VecPointwiseDivide()`, `VecPointwiseMax()`, `VecPointwiseMin()`, `VecPointwiseMaxAbs()`, `VecMaxPointwiseDivide()`

# External Links
$(_doc_external("Vec/VecPointwiseMult"))
"""
function VecPointwiseMult(petsclib::PetscLibType, w::AbstractPetscVec, x::AbstractPetscVec, y::AbstractPetscVec)
    error("VecPointwiseMult: no generated method for these argument types")
end

@for_petsc function VecPointwiseMult(petsclib::$UnionPetscLib, w::AbstractPetscVec, x::AbstractPetscVec, y::AbstractPetscVec )

    @chk ccall(
               (:VecPointwiseMult, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, CVec),
               w, x, y,
              )


	return nothing
end 

"""
	VecPointwiseSign(petsclib::PetscLibType, y::AbstractPetscVec, x::AbstractPetscVec, sign_type::VecSignMode) 
Computes the component-wise sign `y[i] = sign(x[i])`.

Logically Collective

Input Parameters:
- `x`         - the input vector
- `sign_type` - `VecSignMode` indicating how the function should map zero values.

Output Parameter:
- `y` - the sign vector of `x`

Level: beginner

See also: `Vec`, `VecSignMode`

# External Links
$(_doc_external("Vec/VecPointwiseSign"))
"""
function VecPointwiseSign(petsclib::PetscLibType, y::AbstractPetscVec, x::AbstractPetscVec, sign_type::VecSignMode)
    error("VecPointwiseSign: no generated method for these argument types")
end

@for_petsc function VecPointwiseSign(petsclib::$UnionPetscLib, y::AbstractPetscVec, x::AbstractPetscVec, sign_type::VecSignMode )

    @chk ccall(
               (:VecPointwiseSign, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, VecSignMode),
               y, x, sign_type,
              )


	return nothing
end 

"""
	VecPow(petsclib::PetscLibType, v::AbstractPetscVec, p::PetscScalar) 
Replaces each component of a vector by  x_i^p 

Logically Collective

Input Parameters:
- `v` - the vector
- `p` - the exponent to use on each element

Level: intermediate

See also: `Vec`

# External Links
$(_doc_external("Vec/VecPow"))
"""
function VecPow(petsclib::PetscLibType, v::AbstractPetscVec, p::Number)
    error("VecPow: no generated method for these argument types")
end

@for_petsc function VecPow(petsclib::$UnionPetscLib, v::AbstractPetscVec, p::$PetscScalar )

    @chk ccall(
               (:VecPow, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscScalar),
               v, p,
              )


	return nothing
end 

"""
	VecRealPart(petsclib::PetscLibType, v::AbstractPetscVec) 
Replaces a complex vector with its real part

Collective

Input Parameter:
- `v` - the vector

Level: beginner

See also: `Vec`, `VecNorm()`, `VecImaginaryPart()`

# External Links
$(_doc_external("Vec/VecRealPart"))
"""
function VecRealPart(petsclib::PetscLibType, v::AbstractPetscVec)
    error("VecRealPart: no generated method for these argument types")
end

@for_petsc function VecRealPart(petsclib::$UnionPetscLib, v::AbstractPetscVec )

    @chk ccall(
               (:VecRealPart, $petsc_library),
               PetscErrorCode,
               (CVec,),
               v,
              )


	return nothing
end 

"""
	VecReciprocal(petsclib::PetscLibType, vec::AbstractPetscVec) 
Replaces each component of a vector by its reciprocal.

Logically Collective

Input Parameter:
- `vec` - the vector

Output Parameter:
- `vec` - the vector reciprocal

Level: intermediate

See also: `Vec`, `VecLog()`, `VecExp()`, `VecSqrtAbs()`

# External Links
$(_doc_external("Vec/VecReciprocal"))
"""
function VecReciprocal(petsclib::PetscLibType, vec::AbstractPetscVec)
    error("VecReciprocal: no generated method for these argument types")
end

@for_petsc function VecReciprocal(petsclib::$UnionPetscLib, vec::AbstractPetscVec )

    @chk ccall(
               (:VecReciprocal, $petsc_library),
               PetscErrorCode,
               (CVec,),
               vec,
              )


	return nothing
end 

"""
	VecRegister(petsclib::PetscLibType, sname::String, fnc::external) 
Adds a new vector component implementation

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - The name of a new user-defined creation routine
- `function` - The creation routine

See also: `VecRegisterAll()`, `VecRegisterDestroy()`

# External Links
$(_doc_external("Vec/VecRegister"))
"""
function VecRegister(petsclib::PetscLibType, sname::String, fnc::external)
    error("VecRegister: no generated method for these argument types")
end

@for_petsc function VecRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:VecRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	VecRegisterAll(petsclib::PetscLibType) 
Registers all of the vector types in the `Vec` package.

Not Collective

Level: advanced

See also: `Vec`, `VecType`, `VecRegister()`, `VecRegisterDestroy()`

# External Links
$(_doc_external("Vec/VecRegisterAll"))
"""
function VecRegisterAll(petsclib::PetscLibType)
    error("VecRegisterAll: no generated method for these argument types")
end

@for_petsc function VecRegisterAll(petsclib::$UnionPetscLib)

    @chk ccall(
               (:VecRegisterAll, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	VecReplaceArray(petsclib::PetscLibType, vec::AbstractPetscVec, array::Vector{PetscScalar}) 
Allows one to replace the array in a vector with an
array provided by the user. This is useful to avoid copying an array
into a vector.

Logically Collective; No Fortran Support

Input Parameters:
- `vec`   - the vector
- `array` - the array

Level: developer

See also: `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecPlaceArray()`, `VecResetArray()`

# External Links
$(_doc_external("Vec/VecReplaceArray"))
"""
function VecReplaceArray(petsclib::PetscLibType, vec::AbstractPetscVec, array::AbstractVector{<:Number})
    error("VecReplaceArray: no generated method for these argument types")
end

@for_petsc function VecReplaceArray(petsclib::$UnionPetscLib, vec::AbstractPetscVec, array::Vector{$PetscScalar} )

    @chk ccall(
               (:VecReplaceArray, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{$PetscScalar}),
               vec, array,
              )


	return nothing
end 

"""
	VecResetArray(petsclib::PetscLibType, vec::AbstractPetscVec) 
Resets a vector to use its default memory. Call this
after the use of `VecPlaceArray()`.

Not Collective

Input Parameter:
- `vec` - the vector

Level: developer

See also: `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecReplaceArray()`, `VecPlaceArray()`

# External Links
$(_doc_external("Vec/VecResetArray"))
"""
function VecResetArray(petsclib::PetscLibType, vec::AbstractPetscVec)
    error("VecResetArray: no generated method for these argument types")
end

@for_petsc function VecResetArray(petsclib::$UnionPetscLib, vec::AbstractPetscVec )

    @chk ccall(
               (:VecResetArray, $petsc_library),
               PetscErrorCode,
               (CVec,),
               vec,
              )


	return nothing
end 

"""
	VecRestoreArray(petsclib::PetscLibType, x::AbstractPetscVec, a::Union{Ptr, AbstractArray{PetscScalar}}) 
Restores a vector after `VecGetArray()` has been called and the array is no longer needed

Logically Collective

Input Parameters:
- `x` - the vector
- `a` - location of pointer to array obtained from `VecGetArray()`

Level: beginner

See also: `Vec`, `VecGetArray()`, `VecRestoreArrayRead()`, `VecRestoreArrays()`, `VecPlaceArray()`, `VecRestoreArray2d()`,
`VecGetArrayPair()`, `VecRestoreArrayPair()`

# External Links
$(_doc_external("Vec/VecRestoreArray"))
"""
function VecRestoreArray(petsclib::PetscLibType, x::AbstractPetscVec, a::Union{Ptr, AbstractArray{<:Number}})
    error("VecRestoreArray: no generated method for these argument types")
end

@for_petsc function VecRestoreArray(petsclib::$UnionPetscLib, x::AbstractPetscVec, a::Union{Ptr, AbstractArray{$PetscScalar}} )
	a_ = Ref{Ptr{$PetscScalar}}(a isa Ptr ? a : pointer(a))

    @chk ccall(
               (:VecRestoreArray, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Ptr{$PetscScalar}}),
               x, a_,
              )


	return nothing
end 

# override for VecRestoreArray1d; C signature: VecRestoreArray1d(Vec x, PetscInt m, PetscInt mstart, PetscScalar* a[])
"""
	VecRestoreArray1d(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, mstart::PetscInt, a::PetscArray{PetscScalar, 1}) 
Restores a vector after `VecGetArray1d()` has been called.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of two dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `a`      - location of pointer to array obtained from `VecGetArray1d()`

Level: developer

See also: `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecRestoreArrays()`, `VecPlaceArray()`,
`VecGetArray2d()`, `VecGetArray3d()`, `VecRestoreArray3d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`
`VecGetArray1d()`, `VecRestoreArray2d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecRestoreArray1d"))
"""
function VecRestoreArray1d(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, mstart::PetscInt, a::PetscArray{PetscScalar, 1}) end

@for_petsc function VecRestoreArray1d(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, mstart::$PetscInt, a::PetscArray{$PetscScalar, 1} )
	if a.ptr[]  != C_NULL 

    @chk ccall(
               (:VecRestoreArray1d, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, Ref{Ptr{$PetscScalar}}),
               x, m, mstart, a.ptr,
              )

		
	else
		error("The input array is already restored")
	end


	return nothing
end

# override for VecRestoreArray1dRead; C signature: VecRestoreArray1dRead(Vec x, PetscInt m, PetscInt mstart, PetscScalar* a[])
"""
	VecRestoreArray1dRead(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, mstart::PetscInt, a::PetscArray{PetscScalar, 1}) 
Restores a vector after `VecGetArray1dRead()` has been called.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of two dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `a`      - location of pointer to array obtained from `VecGetArray1dRead()`

Level: developer

See also: `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecRestoreArrays()`, `VecPlaceArray()`,
`VecGetArray2d()`, `VecGetArray3d()`, `VecRestoreArray3d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`
`VecGetArray1d()`, `VecRestoreArray2d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecRestoreArray1dRead"))
"""
function VecRestoreArray1dRead(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, mstart::PetscInt, a::PetscArray{PetscScalar, 1}) end

@for_petsc function VecRestoreArray1dRead(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, mstart::$PetscInt, a::PetscArray{$PetscScalar, 1} )
	if a.ptr[]  != C_NULL 
        @chk ccall(
                (:VecRestoreArray1dRead, $petsc_library),
                PetscErrorCode,
                (CVec, $PetscInt, $PetscInt, Ptr{Ptr{$PetscScalar}}),
                x, m, mstart, a.ptr,
                )
	else
		error("The input array is already restored")
	end


	return nothing
end

# override for VecRestoreArray1dWrite; C signature: VecRestoreArray1dWrite(Vec x, PetscInt m, PetscInt mstart, PetscScalar* a[])
"""
	VecRestoreArray1dWrite(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, mstart::PetscInt, a::PetscArray{PetscScalar, 1}) 
Restores a vector after `VecGetArray1dWrite()` has been called.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of two dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `a`      - location of pointer to array obtained from `VecGetArray1d()`

Level: developer

See also: `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecRestoreArrays()`, `VecPlaceArray()`,
`VecGetArray2d()`, `VecGetArray3d()`, `VecRestoreArray3d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`
`VecGetArray1d()`, `VecRestoreArray2d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecRestoreArray1dWrite"))
"""
function VecRestoreArray1dWrite(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, mstart::PetscInt, a::PetscArray{PetscScalar, 1}) end

@for_petsc function VecRestoreArray1dWrite(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, mstart::$PetscInt, a::PetscArray{$PetscScalar, 1} )
	if a.ptr[]  != C_NULL 

    @chk ccall(
               (:VecRestoreArray1dWrite, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, Ref{Ptr{$PetscScalar}}),
               x, m, mstart, a.ptr,
              )

		
	else
		error("The input array is already restored")
	end


	return nothing
end

# override for VecRestoreArray2d; C signature: VecRestoreArray2d(Vec x, PetscInt m, PetscInt n, PetscInt mstart, PetscInt nstart, PetscScalar** a[])
"""
	VecRestoreArray2d(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, mstart::PetscInt, nstart::PetscInt, a::PetscArray{PetscScalar, 2}) 
Restores a vector after `VecGetArray2d()` has been called.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of two dimensional array
- `n`      - second dimension of the two dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `nstart` - first index in the second coordinate direction (often 0)
- `a`      - location of pointer to array obtained from `VecGetArray2d()`

Level: developer

See also: `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecRestoreArrays()`, `VecPlaceArray()`,
`VecGetArray2d()`, `VecGetArray3d()`, `VecRestoreArray3d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecRestoreArray2d"))
"""
function VecRestoreArray2d(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, mstart::PetscInt, nstart::PetscInt, a::PetscArray{PetscScalar, 2}) end

@for_petsc function VecRestoreArray2d(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, n::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, a::PetscArray{$PetscScalar, 2} )
	if a.ptr[] != C_NULL 

        @chk ccall(
                (:VecRestoreArray2d, $petsc_library),
                PetscErrorCode,
                (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{$PetscScalar}}}),
                x, m, n, mstart, nstart, a.ptr,
                )

	else
		error("The input array is already restored")
	end


	return nothing
end

# override for VecRestoreArray2dRead; C signature: VecRestoreArray2dRead(Vec x, PetscInt m, PetscInt n, PetscInt mstart, PetscInt nstart, PetscScalar** a[])
"""
	VecRestoreArray2dRead(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, mstart::PetscInt, nstart::PetscInt, a::PetscArray{PetscScalar, 2}) 
Restores a vector after `VecGetArray2dRead()` has been called.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of two dimensional array
- `n`      - second dimension of the two dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `nstart` - first index in the second coordinate direction (often 0)
- `a`      - location of pointer to array obtained from VecGetArray2d()

Level: developer

See also: `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecRestoreArrays()`, `VecPlaceArray()`,
`VecGetArray2d()`, `VecGetArray3d()`, `VecRestoreArray3d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecRestoreArray2dRead"))
"""
function VecRestoreArray2dRead(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, mstart::PetscInt, nstart::PetscInt, a::PetscArray{PetscScalar, 2}) end

@for_petsc function VecRestoreArray2dRead(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, n::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, a::PetscArray{$PetscScalar, 2} )
	if a.ptr[]  != C_NULL  

        @chk ccall(
                (:VecRestoreArray2dRead, $petsc_library),
                PetscErrorCode,
                (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{$PetscScalar}}}),
                x, m, n, mstart, nstart, a.ptr,
                )

		
	else
		error("The input array is already restored")
	end


	return nothing
end

# override for VecRestoreArray2dWrite; C signature: VecRestoreArray2dWrite(Vec x, PetscInt m, PetscInt n, PetscInt mstart, PetscInt nstart, PetscScalar** a[])
"""
	VecRestoreArray2dWrite(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, mstart::PetscInt, nstart::PetscInt, a::PetscArray{PetscScalar, 2}) 
Restores a vector after `VecGetArray2dWrite()` has been called.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of two dimensional array
- `n`      - second dimension of the two dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `nstart` - first index in the second coordinate direction (often 0)
- `a`      - location of pointer to array obtained from `VecGetArray2d()`

Level: developer

See also: `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecRestoreArrays()`, `VecPlaceArray()`,
`VecGetArray2d()`, `VecGetArray3d()`, `VecRestoreArray3d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecRestoreArray2dWrite"))
"""
function VecRestoreArray2dWrite(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, mstart::PetscInt, nstart::PetscInt, a::PetscArray{PetscScalar, 2}) end

@for_petsc function VecRestoreArray2dWrite(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, n::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, a::PetscArray{$PetscScalar, 2} )
	if a.ptr[]  != C_NULL  

        @chk ccall(
                (:VecRestoreArray2dWrite, $petsc_library),
                PetscErrorCode,
                (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{$PetscScalar}}}),
                x, m, n, mstart, nstart, a.ptr,
                )

		
	else
		error("The input array is already restored")
	end


	return nothing
end

# override for VecRestoreArray3d; C signature: VecRestoreArray3d(Vec x, PetscInt m, PetscInt n, PetscInt p, PetscInt mstart, PetscInt nstart, PetscInt pstart, PetscScalar** a[])
"""
	VecRestoreArray3d(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, a::PetscArray{PetscScalar, 3}) 
Restores a vector after `VecGetArray3d()` has been called.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of three dimensional array
- `n`      - second dimension of the three dimensional array
- `p`      - third dimension of the three dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `nstart` - first index in the second coordinate direction (often 0)
- `pstart` - first index in the third coordinate direction (often 0)
- `a`      - location of pointer to array obtained from VecGetArray3d()

Level: developer

See also: `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecRestoreArrays()`, `VecPlaceArray()`,
`VecGetArray2d()`, `VecGetArray3d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecRestoreArray3d"))
"""
function VecRestoreArray3d(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, a::PetscArray{PetscScalar, 3}) end

@for_petsc function VecRestoreArray3d(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, n::$PetscInt, p::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, pstart::$PetscInt, a::PetscArray{$PetscScalar, 3} )
	if a.ptr[]  != C_NULL  

    @chk ccall(
               (:VecRestoreArray3d, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{$PetscScalar}}}),
               x, m, n, p, mstart, nstart, pstart, a.ptr,
              )

		
	else
		error("The input array is already restored")
	end


	return nothing
end

# override for VecRestoreArray3dRead; C signature: VecRestoreArray3dRead(Vec x, PetscInt m, PetscInt n, PetscInt p, PetscInt mstart, PetscInt nstart, PetscInt pstart, PetscScalar** a[])
"""
	VecRestoreArray3dRead(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, a::PetscArray{PetscScalar, 3}) 
Restores a vector after `VecGetArray3dRead()` has been called.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of three dimensional array
- `n`      - second dimension of the three dimensional array
- `p`      - third dimension of the three dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `nstart` - first index in the second coordinate direction (often 0)
- `pstart` - first index in the third coordinate direction (often 0)
- `a`      - location of pointer to array obtained from `VecGetArray3dRead()`

Level: developer

See also: `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecRestoreArrays()`, `VecPlaceArray()`,
`VecGetArray2d()`, `VecGetArray3d()`, `VecRestoreArray3d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecRestoreArray3dRead"))
"""
function VecRestoreArray3dRead(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, a::PetscArray{PetscScalar, 3}) end

@for_petsc function VecRestoreArray3dRead(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, n::$PetscInt, p::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, pstart::$PetscInt, a::PetscArray{$PetscScalar, 3} )
	if a.ptr[]  != C_NULL 

    @chk ccall(
               (:VecRestoreArray3dRead, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{$PetscScalar}}}),
               x, m, n, p, mstart, nstart, pstart, a.ptr,
              )

		
	else
		error("The input array is already restored")
	end


	return nothing
end

# override for VecRestoreArray3dWrite; C signature: VecRestoreArray3dWrite(Vec x, PetscInt m, PetscInt n, PetscInt p, PetscInt mstart, PetscInt nstart, PetscInt pstart, PetscScalar** a[])
"""
	VecRestoreArray3dWrite(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, a::PetscArray{PetscScalar, 3}) 
Restores a vector after `VecGetArray3dWrite()` has been called.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of three dimensional array
- `n`      - second dimension of the three dimensional array
- `p`      - third dimension of the three dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `nstart` - first index in the second coordinate direction (often 0)
- `pstart` - first index in the third coordinate direction (often 0)
- `a`      - location of pointer to array obtained from VecGetArray3d()

Level: developer

See also: `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecRestoreArrays()`, `VecPlaceArray()`,
`VecGetArray2d()`, `VecGetArray3d()`, `VecRestoreArray3d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecRestoreArray3dWrite"))
"""
function VecRestoreArray3dWrite(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, a::PetscArray{PetscScalar, 3}) end

@for_petsc function VecRestoreArray3dWrite(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, n::$PetscInt, p::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, pstart::$PetscInt, a::PetscArray{$PetscScalar, 3} )
	if a.ptr[]  != C_NULL  

    @chk ccall(
               (:VecRestoreArray3dWrite, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{$PetscScalar}}}),
               x, m, n, p, mstart, nstart, pstart, a.ptr,
              )

		
	else
		error("The input array is already restored")
	end


	return nothing
end

# override for VecRestoreArray4d; C signature: VecRestoreArray4d(Vec x, PetscInt m, PetscInt n, PetscInt p, PetscInt q, PetscInt mstart, PetscInt nstart, PetscInt pstart, PetscInt qstart, PetscScalar** a[])
"""
	VecRestoreArray4d(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, q::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, qstart::PetscInt, a::PetscArray{PetscScalar, 4}) 
Restores a vector after `VecGetArray4d()` has been called.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of four dimensional array
- `n`      - second dimension of the four dimensional array
- `p`      - third dimension of the four dimensional array
- `q`      - fourth dimension of the four dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `nstart` - first index in the second coordinate direction (often 0)
- `pstart` - first index in the third coordinate direction (often 0)
- `qstart` - first index in the fourth coordinate direction (often 0)
- `a`      - location of pointer to array obtained from VecGetArray4d()

Level: developer

See also: `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecRestoreArrays()`, `VecPlaceArray()`,
`VecGetArray2d()`, `VecGetArray3d()`, `VecRestoreArray3d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`

# External Links
$(_doc_external("Vec/VecRestoreArray4d"))
"""
function VecRestoreArray4d(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, q::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, qstart::PetscInt, a::PetscArray{PetscScalar, 4}) end

@for_petsc function VecRestoreArray4d(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, n::$PetscInt, p::$PetscInt, q::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, pstart::$PetscInt, qstart::$PetscInt, a::PetscArray{$PetscScalar, 4} )
	if a.ptr[]  != C_NULL  

        @chk ccall(
                (:VecRestoreArray4d, $petsc_library),
                PetscErrorCode,
                (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{$PetscScalar}}}),
                x, m, n, p, q, mstart, nstart, pstart, qstart, a.ptr,
                )

		
	else
		error("The input array is already restored")
	end

	return nothing
end

# override for VecRestoreArray4dRead; C signature: VecRestoreArray4dRead(Vec x, PetscInt m, PetscInt n, PetscInt p, PetscInt q, PetscInt mstart, PetscInt nstart, PetscInt pstart, PetscInt qstart, PetscScalar** a[])
"""
	VecRestoreArray4dRead(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, q::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, qstart::PetscInt, a::PetscArray{PetscScalar, 4}) 
Restores a vector after `VecGetArray4d()` has been called.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of four dimensional array
- `n`      - second dimension of the four dimensional array
- `p`      - third dimension of the four dimensional array
- `q`      - fourth dimension of the four dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `nstart` - first index in the second coordinate direction (often 0)
- `pstart` - first index in the third coordinate direction (often 0)
- `qstart` - first index in the fourth coordinate direction (often 0)
- `a`      - location of pointer to array obtained from `VecGetArray4dRead()`

Level: beginner

See also: `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecRestoreArrays()`, `VecPlaceArray()`,
`VecGetArray2d()`, `VecGetArray3d()`, `VecRestoreArray3d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecRestoreArray4dRead"))
"""
function VecRestoreArray4dRead(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, q::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, qstart::PetscInt, a::PetscArray{PetscScalar, 4}) end

@for_petsc function VecRestoreArray4dRead(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, n::$PetscInt, p::$PetscInt, q::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, pstart::$PetscInt, qstart::$PetscInt, a::PetscArray{$PetscScalar, 4} )
	if a.ptr[]  != C_NULL  

    @chk ccall(
               (:VecRestoreArray4dRead, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{$PetscScalar}}}),
               x, m, n, p, q, mstart, nstart, pstart, qstart, a.ptr,
              )

		
	else
		error("The input array is already restored")
	end


	return nothing
end

# override for VecRestoreArray4dWrite; C signature: VecRestoreArray4dWrite(Vec x, PetscInt m, PetscInt n, PetscInt p, PetscInt q, PetscInt mstart, PetscInt nstart, PetscInt pstart, PetscInt qstart, PetscScalar** a[])
"""
	VecRestoreArray4dWrite(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, q::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, qstart::PetscInt, a::PetscArray{PetscScalar, 4}) 
Restores a vector after `VecGetArray4dWrite()` has been called.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of four dimensional array
- `n`      - second dimension of the four dimensional array
- `p`      - third dimension of the four dimensional array
- `q`      - fourth dimension of the four dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `nstart` - first index in the second coordinate direction (often 0)
- `pstart` - first index in the third coordinate direction (often 0)
- `qstart` - first index in the fourth coordinate direction (often 0)
- `a`      - location of pointer to array obtained from `VecGetArray4d()`

Level: developer

See also: `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecRestoreArrays()`, `VecPlaceArray()`,
`VecGetArray2d()`, `VecGetArray3d()`, `VecRestoreArray3d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecRestoreArray4dWrite"))
"""
function VecRestoreArray4dWrite(petsclib::PetscLibType, x::AbstractPetscVec, m::PetscInt, n::PetscInt, p::PetscInt, q::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, qstart::PetscInt, a::PetscArray{PetscScalar, 4}) end

@for_petsc function VecRestoreArray4dWrite(petsclib::$UnionPetscLib, x::AbstractPetscVec, m::$PetscInt, n::$PetscInt, p::$PetscInt, q::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, pstart::$PetscInt, qstart::$PetscInt, a::PetscArray{$PetscScalar, 4} )
	if a.ptr[]  != C_NULL  

    @chk ccall(
               (:VecRestoreArray4dWrite, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{Ptr{Ptr{$PetscScalar}}}}}),
               x, m, n, p, q, mstart, nstart, pstart, qstart, a.ptr,
              )
	else
		error("The input array is already restored")
	end


	return nothing
end

"""
	VecRestoreArrayAndMemType(petsclib::PetscLibType, x::AbstractPetscVec, a::Union{Ptr, AbstractArray{PetscScalar}}) 
Restores a vector after `VecGetArrayAndMemType()` has been called.

Logically Collective; No Fortran Support

Input Parameters:
- `x` - the vector
- `a` - location of pointer to array obtained from `VecGetArrayAndMemType()`

Level: beginner

See also: `Vec`, `VecGetArrayAndMemType()`, `VecGetArray()`, `VecRestoreArrayRead()`, `VecRestoreArrays()`,
`VecPlaceArray()`, `VecRestoreArray2d()`, `VecGetArrayPair()`, `VecRestoreArrayPair()`

# External Links
$(_doc_external("Vec/VecRestoreArrayAndMemType"))
"""
function VecRestoreArrayAndMemType(petsclib::PetscLibType, x::AbstractPetscVec, a::Union{Ptr, AbstractArray{<:Number}})
    error("VecRestoreArrayAndMemType: no generated method for these argument types")
end

@for_petsc function VecRestoreArrayAndMemType(petsclib::$UnionPetscLib, x::AbstractPetscVec, a::Union{Ptr, AbstractArray{$PetscScalar}} )
	a_ = Ref{Ptr{$PetscScalar}}(a isa Ptr ? a : pointer(a))

    @chk ccall(
               (:VecRestoreArrayAndMemType, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Ptr{$PetscScalar}}),
               x, a_,
              )


	return nothing
end 

"""
	VecRestoreArrayPair(petsclib::PetscLibType, x::AbstractPetscVec, y::AbstractPetscVec, xv::Union{Ptr, AbstractArray{PetscScalar}}, yv::Union{Ptr, AbstractArray{PetscScalar}}) 

# External Links
$(_doc_external("Vec/VecRestoreArrayPair"))
"""
function VecRestoreArrayPair(petsclib::PetscLibType, x::AbstractPetscVec, y::AbstractPetscVec, xv::Union{Ptr, AbstractArray{<:Number}}, yv::Union{Ptr, AbstractArray{<:Number}})
    error("VecRestoreArrayPair: no generated method for these argument types")
end

@for_petsc function VecRestoreArrayPair(petsclib::$UnionPetscLib, x::AbstractPetscVec, y::AbstractPetscVec, xv::Union{Ptr, AbstractArray{$PetscScalar}}, yv::Union{Ptr, AbstractArray{$PetscScalar}} )
	xv_ = Ref{Ptr{$PetscScalar}}(xv isa Ptr ? xv : pointer(xv))
	yv_ = Ref{Ptr{$PetscScalar}}(yv isa Ptr ? yv : pointer(yv))

    @chk ccall(
               (:VecRestoreArrayPair, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, Ptr{Ptr{$PetscScalar}}, Ptr{Ptr{$PetscScalar}}),
               x, y, xv_, yv_,
              )


	return nothing
end 

"""
	VecRestoreArrayRead(petsclib::PetscLibType, x::AbstractPetscVec, a::Union{Ptr, AbstractArray{PetscScalar}}) 
Restore array obtained with `VecGetArrayRead()`

Not Collective

Input Parameters:
- `x` - the vector
- `a` - the array

Level: beginner

See also: `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrayPair()`, `VecRestoreArrayPair()`

# External Links
$(_doc_external("Vec/VecRestoreArrayRead"))
"""
function VecRestoreArrayRead(petsclib::PetscLibType, x::AbstractPetscVec, a::Union{Ptr, AbstractArray{<:Number}})
    error("VecRestoreArrayRead: no generated method for these argument types")
end

@for_petsc function VecRestoreArrayRead(petsclib::$UnionPetscLib, x::AbstractPetscVec, a::Union{Ptr, AbstractArray{$PetscScalar}} )
	a_ = Ref{Ptr{$PetscScalar}}(a isa Ptr ? a : pointer(a))

    @chk ccall(
               (:VecRestoreArrayRead, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Ptr{$PetscScalar}}),
               x, a_,
              )


	return nothing
end 

"""
	VecRestoreArrayReadAndMemType(petsclib::PetscLibType, x::AbstractPetscVec, a::Union{Ptr, AbstractArray{PetscScalar}}) 
Restore array obtained with `VecGetArrayReadAndMemType()`

Not Collective; No Fortran Support

Input Parameters:
- `x` - the vector
- `a` - the array

Level: beginner

See also: `Vec`, `VecGetArrayReadAndMemType()`, `VecRestoreArrayAndMemType()`, `VecRestoreArrayWriteAndMemType()`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrayPair()`, `VecRestoreArrayPair()`

# External Links
$(_doc_external("Vec/VecRestoreArrayReadAndMemType"))
"""
function VecRestoreArrayReadAndMemType(petsclib::PetscLibType, x::AbstractPetscVec, a::Union{Ptr, AbstractArray{<:Number}})
    error("VecRestoreArrayReadAndMemType: no generated method for these argument types")
end

@for_petsc function VecRestoreArrayReadAndMemType(petsclib::$UnionPetscLib, x::AbstractPetscVec, a::Union{Ptr, AbstractArray{$PetscScalar}} )
	a_ = Ref{Ptr{$PetscScalar}}(a isa Ptr ? a : pointer(a))

    @chk ccall(
               (:VecRestoreArrayReadAndMemType, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Ptr{$PetscScalar}}),
               x, a_,
              )


	return nothing
end 

"""
	VecRestoreArrayWrite(petsclib::PetscLibType, x::AbstractPetscVec, a::Union{Ptr, AbstractArray{PetscScalar}}) 
Restores a vector after `VecGetArrayWrite()` has been called.

Logically Collective

Input Parameters:
- `x` - the vector
- `a` - location of pointer to array obtained from `VecGetArray()`

Level: beginner

See also: `Vec`, `VecGetArray()`, `VecRestoreArrayRead()`, `VecRestoreArrays()`, `VecPlaceArray()`, `VecRestoreArray2d()`,
`VecGetArrayPair()`, `VecRestoreArrayPair()`, `VecGetArrayWrite()`

# External Links
$(_doc_external("Vec/VecRestoreArrayWrite"))
"""
function VecRestoreArrayWrite(petsclib::PetscLibType, x::AbstractPetscVec, a::Union{Ptr, AbstractArray{<:Number}})
    error("VecRestoreArrayWrite: no generated method for these argument types")
end

@for_petsc function VecRestoreArrayWrite(petsclib::$UnionPetscLib, x::AbstractPetscVec, a::Union{Ptr, AbstractArray{$PetscScalar}} )
	a_ = Ref{Ptr{$PetscScalar}}(a isa Ptr ? a : pointer(a))

    @chk ccall(
               (:VecRestoreArrayWrite, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Ptr{$PetscScalar}}),
               x, a_,
              )


	return nothing
end 

"""
	VecRestoreArrayWriteAndMemType(petsclib::PetscLibType, x::AbstractPetscVec, a::Union{Ptr, AbstractArray{PetscScalar}}) 
Restore array obtained with `VecGetArrayWriteAndMemType()`

Logically Collective; No Fortran Support

Input Parameters:
- `x` - the vector
- `a` - the array

Level: beginner

See also: `Vec`, `VecGetArrayWriteAndMemType()`, `VecRestoreArrayAndMemType()`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrayPair()`, `VecRestoreArrayPair()`

# External Links
$(_doc_external("Vec/VecRestoreArrayWriteAndMemType"))
"""
function VecRestoreArrayWriteAndMemType(petsclib::PetscLibType, x::AbstractPetscVec, a::Union{Ptr, AbstractArray{<:Number}})
    error("VecRestoreArrayWriteAndMemType: no generated method for these argument types")
end

@for_petsc function VecRestoreArrayWriteAndMemType(petsclib::$UnionPetscLib, x::AbstractPetscVec, a::Union{Ptr, AbstractArray{$PetscScalar}} )
	a_ = Ref{Ptr{$PetscScalar}}(a isa Ptr ? a : pointer(a))

    @chk ccall(
               (:VecRestoreArrayWriteAndMemType, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Ptr{$PetscScalar}}),
               x, a_,
              )


	return nothing
end 

"""
	VecRestoreArrays(petsclib::PetscLibType, x::Vector{<:AbstractPetscVec}, n::PetscInt, a::Vector{PetscScalar}) 
Restores a group of vectors after `VecGetArrays()`
has been called.

Logically Collective; No Fortran Support

Input Parameters:
- `x` - the vector
- `n` - the number of vectors
- `a` - location of pointer to arrays obtained from `VecGetArrays()`

See also: `Vec`, `VecGetArrays()`, `VecRestoreArray()`

# External Links
$(_doc_external("Vec/VecRestoreArrays"))
"""
function VecRestoreArrays(petsclib::PetscLibType, x::Vector{<:AbstractPetscVec}, n::Integer, a::AbstractVector{<:Number})
    error("VecRestoreArrays: no generated method for these argument types")
end

@for_petsc function VecRestoreArrays(petsclib::$UnionPetscLib, x::Vector{<:AbstractPetscVec}, n::$PetscInt, a::Vector{$PetscScalar} )

    @chk ccall(
               (:VecRestoreArrays, $petsc_library),
               PetscErrorCode,
               (Ptr{CVec}, $PetscInt, Ptr{Ptr{Ptr{$PetscScalar}}}),
               x, n, a,
              )


	return nothing
end 

"""
	VecRestoreLocalVector(petsclib::PetscLibType, v::AbstractPetscVec, w::AbstractPetscVec) 
Unmaps the local portion of a vector
previously mapped into a vector using `VecGetLocalVector()`.

Logically Collective.

Input Parameters:
- `v` - The local portion of this vector was previously mapped into `w` using `VecGetLocalVector()`.
- `w` - The vector into which the local portion of `v` was mapped.

Level: beginner

See also: `Vec`, `VecCreateLocalVector()`, `VecGetLocalVector()`, `VecGetLocalVectorRead()`, `VecRestoreLocalVectorRead()`, `LocalVectorRead()`, `VecGetArrayRead()`, `VecGetArray()`

# External Links
$(_doc_external("Vec/VecRestoreLocalVector"))
"""
function VecRestoreLocalVector(petsclib::PetscLibType, v::AbstractPetscVec, w::AbstractPetscVec)
    error("VecRestoreLocalVector: no generated method for these argument types")
end

@for_petsc function VecRestoreLocalVector(petsclib::$UnionPetscLib, v::AbstractPetscVec, w::AbstractPetscVec )

    @chk ccall(
               (:VecRestoreLocalVector, $petsc_library),
               PetscErrorCode,
               (CVec, CVec),
               v, w,
              )


	return nothing
end 

"""
	VecRestoreLocalVectorRead(petsclib::PetscLibType, v::AbstractPetscVec, w::AbstractPetscVec) 
Unmaps the local portion of a vector
previously mapped into a vector using `VecGetLocalVectorRead()`.

Not Collective.

Input Parameters:
- `v` - The local portion of this vector was previously mapped into `w` using `VecGetLocalVectorRead()`.
- `w` - The vector into which the local portion of `v` was mapped.

Level: beginner

See also: `Vec`, `VecCreateLocalVector()`, `VecGetLocalVectorRead()`, `VecGetLocalVector()`, `VecGetArrayRead()`, `VecGetArray()`

# External Links
$(_doc_external("Vec/VecRestoreLocalVectorRead"))
"""
function VecRestoreLocalVectorRead(petsclib::PetscLibType, v::AbstractPetscVec, w::AbstractPetscVec)
    error("VecRestoreLocalVectorRead: no generated method for these argument types")
end

@for_petsc function VecRestoreLocalVectorRead(petsclib::$UnionPetscLib, v::AbstractPetscVec, w::AbstractPetscVec )

    @chk ccall(
               (:VecRestoreLocalVectorRead, $petsc_library),
               PetscErrorCode,
               (CVec, CVec),
               v, w,
              )


	return nothing
end 

"""
	VecRestoreSubVector(petsclib::PetscLibType, X::AbstractPetscVec, is::AbstractIS, Y::AbstractPetscVec) 
Restores a subvector extracted using `VecGetSubVector()`

Collective

Input Parameters:
- `X`  - vector from which subvector was obtained
- `is` - index set representing the subset of `X`
- `Y`  - subvector being restored

Level: advanced

See also: `Vec`, `IS`, `VecGetSubVector()`

# External Links
$(_doc_external("Vec/VecRestoreSubVector"))
"""
function VecRestoreSubVector(petsclib::PetscLibType, X::AbstractPetscVec, is::AbstractIS, Y::AbstractPetscVec)
    error("VecRestoreSubVector: no generated method for these argument types")
end

@for_petsc function VecRestoreSubVector(petsclib::$UnionPetscLib, X::AbstractPetscVec, is::AbstractIS, Y::AbstractPetscVec )
	Y_ = Ref(Y.ptr)

    @chk ccall(
               (:VecRestoreSubVector, $petsc_library),
               PetscErrorCode,
               (CVec, CIS, Ptr{CVec}),
               X, is, Y_,
              )

	Y.ptr = Y_[]

	return nothing
end 

"""
	val::PetscScalar = VecTDot(petsclib::PetscLibType, x::AbstractPetscVec, y::AbstractPetscVec) 
Computes an indefinite vector dot product. That is, this
routine does NOT use the complex conjugate.

Collective

Input Parameters:
- `x` - first vector
- `y` - second vector

Output Parameter:
- `val` - the dot product

Level: intermediate

Notes for Users of Complex Numbers:
For complex vectors, `VecTDot()` computes the indefinite form
``
val = (x,y) = y^T x,
``
where y^T denotes the transpose of y.

Use `VecDot()` for the inner product
``
val = (x,y) = y^H x,
``
where y^H denotes the conjugate transpose of y.

See also: `Vec`, `VecDot()`, `VecMTDot()`

# External Links
$(_doc_external("Vec/VecTDot"))
"""
function VecTDot(petsclib::PetscLibType, x::AbstractPetscVec, y::AbstractPetscVec)
    error("VecTDot: no generated method for these argument types")
end

@for_petsc function VecTDot(petsclib::$UnionPetscLib, x::AbstractPetscVec, y::AbstractPetscVec )
	val_ = Ref{$PetscScalar}()

    @chk ccall(
               (:VecTDot, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, Ptr{$PetscScalar}),
               x, y, val_,
              )

	val = val_[]

	return val
end 

"""
	result::PetscScalar = VecTDotBegin(petsclib::PetscLibType, x::AbstractPetscVec, y::AbstractPetscVec) 
Starts a split phase transpose dot product computation.

Input Parameters:
- `x`      - the first vector
- `y`      - the second vector
- `result` - where the result will go (can be `NULL`)

Level: advanced

See also: `VecTDotEnd()`, `VecNormBegin()`, `VecNormEnd()`, `VecNorm()`, `VecDot()`, `VecMDot()`,
`VecDotBegin()`, `VecDotEnd()`, `PetscCommSplitReductionBegin()`

# External Links
$(_doc_external("Vec/VecTDotBegin"))
"""
function VecTDotBegin(petsclib::PetscLibType, x::AbstractPetscVec, y::AbstractPetscVec)
    error("VecTDotBegin: no generated method for these argument types")
end

@for_petsc function VecTDotBegin(petsclib::$UnionPetscLib, x::AbstractPetscVec, y::AbstractPetscVec )
	result_ = Ref{$PetscScalar}()

    @chk ccall(
               (:VecTDotBegin, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, Ptr{$PetscScalar}),
               x, y, result_,
              )

	result = result_[]

	return result
end 

"""
	result::PetscScalar = VecTDotEnd(petsclib::PetscLibType, x::AbstractPetscVec, y::AbstractPetscVec) 
Ends a split phase transpose dot product computation.

Input Parameters:
- `x`      - the first vector (can be `NULL`)
- `y`      - the second vector (can be `NULL`)
- `result` - where the result will go

Level: advanced

See also: `VecTDotBegin()`, `VecNormBegin()`, `VecNormEnd()`, `VecNorm()`, `VecDot()`, `VecMDot()`,
`VecDotBegin()`, `VecDotEnd()`

# External Links
$(_doc_external("Vec/VecTDotEnd"))
"""
function VecTDotEnd(petsclib::PetscLibType, x::AbstractPetscVec, y::AbstractPetscVec)
    error("VecTDotEnd: no generated method for these argument types")
end

@for_petsc function VecTDotEnd(petsclib::$UnionPetscLib, x::AbstractPetscVec, y::AbstractPetscVec )
	result_ = Ref{$PetscScalar}()

    @chk ccall(
               (:VecTDotEnd, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, Ptr{$PetscScalar}),
               x, y, result_,
              )

	result = result_[]

	return result
end 

"""
	n::PetscInt,e::Vector{PetscScalar} = VecUniqueEntries(petsclib::PetscLibType, vec::AbstractPetscVec) 
Compute the number of unique entries, and those entries

Collective

Input Parameter:
- `vec` - the vector

Output Parameters:
- `n` - The number of unique entries
- `e` - The entries, each MPI process receives all the unique entries

Level: intermediate

See also: `Vec`

# External Links
$(_doc_external("Vec/VecUniqueEntries"))
"""
function VecUniqueEntries(petsclib::PetscLibType, vec::AbstractPetscVec)
    error("VecUniqueEntries: no generated method for these argument types")
end

@for_petsc function VecUniqueEntries(petsclib::$UnionPetscLib, vec::AbstractPetscVec )
	n_ = Ref{$PetscInt}()
	e_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:VecUniqueEntries, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{$PetscInt}, Ptr{Ptr{$PetscScalar}}),
               vec, n_, e_,
              )

	n = n_[]
	e = e_[] == C_NULL ? $PetscScalar[] : unsafe_wrap(Array, e_[], n; own = false)

	return n,e
end 

"""
	VecView(petsclib::PetscLibType, vec::AbstractPetscVec, viewer::PetscViewer) 
Views a vector object.

Collective

Input Parameters:
- `vec`    - the vector
- `viewer` - an optional `PetscViewer` visualization context

Level: beginner

See also: `Vec`, `VecViewFromOptions()`, `PetscViewerASCIIOpen()`, `PetscViewerDrawOpen()`, `PetscDrawLGCreate()`,
`PetscViewerSocketOpen()`, `PetscViewerBinaryOpen()`, `VecLoad()`, `PetscViewerCreate()`,
`PetscRealView()`, `PetscScalarView()`, `PetscIntView()`, `PetscViewerHDF5SetTimestep()`

# External Links
$(_doc_external("Vec/VecView"))
"""
function VecView(petsclib::PetscLibType, vec::AbstractPetscVec, viewer::PetscViewer)
    error("VecView: no generated method for these argument types")
end

@for_petsc function VecView(petsclib::$UnionPetscLib, vec::AbstractPetscVec, viewer::PetscViewer )

    @chk ccall(
               (:VecView, $petsc_library),
               PetscErrorCode,
               (CVec, PetscViewer),
               vec, viewer,
              )


	return nothing
end 

"""
	VecViewFromOptions(petsclib::PetscLibType, A::AbstractPetscVec, obj, name::String) 
View a vector based on values in the options database

Collective

Input Parameters:
- `A`    - the vector
- `obj`  - optional object that provides the options prefix for this viewing, use `NULL` to use the prefix of `A`
- `name` - command line option

Options Database Key:
- `-name [viewertype][:...]` - option name and values. See `PetscObjectViewFromOptions()` for the possible arguments

Level: intermediate

See also: `Vec`, `VecView`, `PetscObjectViewFromOptions()`, `VecCreate()`

# External Links
$(_doc_external("Vec/VecViewFromOptions"))
"""
function VecViewFromOptions(petsclib::PetscLibType, A::AbstractPetscVec, obj, name::String)
    error("VecViewFromOptions: no generated method for these argument types")
end

@for_petsc function VecViewFromOptions(petsclib::$UnionPetscLib, A::AbstractPetscVec, obj, name::String )

    @chk ccall(
               (:VecViewFromOptions, $petsc_library),
               PetscErrorCode,
               (CVec, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

"""
	VecViewNative(petsclib::PetscLibType, vec::AbstractPetscVec, viewer::PetscViewer) 
Views a vector object with the original type specific viewer

Collective

Input Parameters:
- `vec`    - the vector
- `viewer` - an optional `PetscViewer` visualization context

Level: developer

See also: `Vec`, `PetscViewerASCIIOpen()`, `PetscViewerDrawOpen()`, `PetscDrawLGCreate()`, `VecView()`,
`PetscViewerSocketOpen()`, `PetscViewerBinaryOpen()`, `VecLoad()`, `PetscViewerCreate()`,
`PetscRealView()`, `PetscScalarView()`, `PetscIntView()`, `PetscViewerHDF5SetTimestep()`

# External Links
$(_doc_external("Vec/VecViewNative"))
"""
function VecViewNative(petsclib::PetscLibType, vec::AbstractPetscVec, viewer::PetscViewer)
    error("VecViewNative: no generated method for these argument types")
end

@for_petsc function VecViewNative(petsclib::$UnionPetscLib, vec::AbstractPetscVec, viewer::PetscViewer )

    @chk ccall(
               (:VecViewNative, $petsc_library),
               PetscErrorCode,
               (CVec, PetscViewer),
               vec, viewer,
              )


	return nothing
end 

"""
	VecWAXPY(petsclib::PetscLibType, w::AbstractPetscVec, alpha::PetscScalar, x::AbstractPetscVec, y::AbstractPetscVec) 
Computes `w = alpha x + y`.

Logically Collective

Input Parameters:
- `alpha` - the scalar
- `x`     - first vector, multiplied by `alpha`
- `y`     - second vector

Output Parameter:
- `w` - the result

Level: intermediate

See also: `Vec`, `VecAXPY()`, `VecAYPX()`, `VecAXPBY()`, `VecMAXPY()`, `VecAXPBYPCZ()`

# External Links
$(_doc_external("Vec/VecWAXPY"))
"""
function VecWAXPY(petsclib::PetscLibType, w::AbstractPetscVec, alpha::Number, x::AbstractPetscVec, y::AbstractPetscVec)
    error("VecWAXPY: no generated method for these argument types")
end

@for_petsc function VecWAXPY(petsclib::$UnionPetscLib, w::AbstractPetscVec, alpha::$PetscScalar, x::AbstractPetscVec, y::AbstractPetscVec )

    @chk ccall(
               (:VecWAXPY, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscScalar, CVec, CVec),
               w, alpha, x, y,
              )


	return nothing
end 

"""
	S::IS = VecWhichBetween(petsclib::PetscLibType, VecLow::AbstractPetscVec, V::AbstractPetscVec, VecHigh::AbstractPetscVec) 
Creates an index set containing the indices
where  `VecLow` < `V` < `VecHigh`

Collective

Input Parameters:
- `VecLow`  - lower bound
- `V`       - Vector to compare
- `VecHigh` - higher bound

Output Parameter:
- `S` - The index set containing the indices i where veclow[i] < v[i] < vechigh[i]

Level: advanced

See also: `Vec`

# External Links
$(_doc_external("Vec/VecWhichBetween"))
"""
function VecWhichBetween(petsclib::PetscLibType, VecLow::AbstractPetscVec, V::AbstractPetscVec, VecHigh::AbstractPetscVec)
    error("VecWhichBetween: no generated method for these argument types")
end

@for_petsc function VecWhichBetween(petsclib::$UnionPetscLib, VecLow::AbstractPetscVec, V::AbstractPetscVec, VecHigh::AbstractPetscVec )
	S_ = Ref{CIS}()

    @chk ccall(
               (:VecWhichBetween, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, CVec, Ptr{CIS}),
               VecLow, V, VecHigh, S_,
              )

	S = IS(S_[], petsclib)

	return S
end 

"""
	S::IS = VecWhichBetweenOrEqual(petsclib::PetscLibType, VecLow::AbstractPetscVec, V::AbstractPetscVec, VecHigh::AbstractPetscVec) 
Creates an index set containing the indices
where  `VecLow` <= `V` <= `VecHigh`

Collective

Input Parameters:
- `VecLow`  - lower bound
- `V`       - Vector to compare
- `VecHigh` - higher bound

Output Parameter:
- `S` - The index set containing the indices i where veclow[i] <= v[i] <= vechigh[i]

Level: advanced

See also: `Vec`

# External Links
$(_doc_external("Vec/VecWhichBetweenOrEqual"))
"""
function VecWhichBetweenOrEqual(petsclib::PetscLibType, VecLow::AbstractPetscVec, V::AbstractPetscVec, VecHigh::AbstractPetscVec)
    error("VecWhichBetweenOrEqual: no generated method for these argument types")
end

@for_petsc function VecWhichBetweenOrEqual(petsclib::$UnionPetscLib, VecLow::AbstractPetscVec, V::AbstractPetscVec, VecHigh::AbstractPetscVec )
	S_ = Ref{CIS}()

    @chk ccall(
               (:VecWhichBetweenOrEqual, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, CVec, Ptr{CIS}),
               VecLow, V, VecHigh, S_,
              )

	S = IS(S_[], petsclib)

	return S
end 

"""
	S::IS = VecWhichEqual(petsclib::PetscLibType, Vec1::AbstractPetscVec, Vec2::AbstractPetscVec) 
Creates an index set containing the indices
where the vectors `Vec1` and `Vec2` have identical elements.

Collective

Input Parameters:
- `Vec1` - the first vector to compare
- `Vec2` - the second two vector to compare

Output Parameter:
- `S` - The index set containing the indices i where vec1[i] == vec2[i]

Level: advanced

See also: `Vec`

# External Links
$(_doc_external("Vec/VecWhichEqual"))
"""
function VecWhichEqual(petsclib::PetscLibType, Vec1::AbstractPetscVec, Vec2::AbstractPetscVec)
    error("VecWhichEqual: no generated method for these argument types")
end

@for_petsc function VecWhichEqual(petsclib::$UnionPetscLib, Vec1::AbstractPetscVec, Vec2::AbstractPetscVec )
	S_ = Ref{CIS}()

    @chk ccall(
               (:VecWhichEqual, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, Ptr{CIS}),
               Vec1, Vec2, S_,
              )

	S = IS(S_[], petsclib)

	return S
end 

"""
	S::IS = VecWhichGreaterThan(petsclib::PetscLibType, Vec1::AbstractPetscVec, Vec2::AbstractPetscVec) 
Creates an index set containing the indices
where the vectors `Vec1` > `Vec2`

Collective

Input Parameters:
- `Vec1` - the first vector to compare
- `Vec2` - the second vector to compare

Output Parameter:
- `S` - The index set containing the indices i where vec1[i] > vec2[i]

Level: advanced

See also: `Vec`

# External Links
$(_doc_external("Vec/VecWhichGreaterThan"))
"""
function VecWhichGreaterThan(petsclib::PetscLibType, Vec1::AbstractPetscVec, Vec2::AbstractPetscVec)
    error("VecWhichGreaterThan: no generated method for these argument types")
end

@for_petsc function VecWhichGreaterThan(petsclib::$UnionPetscLib, Vec1::AbstractPetscVec, Vec2::AbstractPetscVec )
	S_ = Ref{CIS}()

    @chk ccall(
               (:VecWhichGreaterThan, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, Ptr{CIS}),
               Vec1, Vec2, S_,
              )

	S = IS(S_[], petsclib)

	return S
end 

"""
	S::IS = VecWhichInactive(petsclib::PetscLibType, VecLow::AbstractPetscVec, V::AbstractPetscVec, D::AbstractPetscVec, VecHigh::AbstractPetscVec, Strong::PetscBool) 
Creates an `IS` based on a set of vectors

Collective

Input Parameters:
- `VecLow`  - lower bound
- `V`       - Vector to compare
- `D`       - Direction to compare
- `VecHigh` - higher bound
- `Strong`  - indicator for applying strongly inactive test

Output Parameter:
- `S` - The index set containing the indices i where the bound is inactive

Level: advanced

See also: `Vec`

# External Links
$(_doc_external("Vec/VecWhichInactive"))
"""
function VecWhichInactive(petsclib::PetscLibType, VecLow::AbstractPetscVec, V::AbstractPetscVec, D::AbstractPetscVec, VecHigh::AbstractPetscVec, Strong::PetscBool)
    error("VecWhichInactive: no generated method for these argument types")
end

@for_petsc function VecWhichInactive(petsclib::$UnionPetscLib, VecLow::AbstractPetscVec, V::AbstractPetscVec, D::AbstractPetscVec, VecHigh::AbstractPetscVec, Strong::PetscBool )
	S_ = Ref{CIS}()

    @chk ccall(
               (:VecWhichInactive, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, CVec, CVec, PetscBool, Ptr{CIS}),
               VecLow, V, D, VecHigh, Strong, S_,
              )

	S = IS(S_[], petsclib)

	return S
end 

"""
	S::IS = VecWhichLessThan(petsclib::PetscLibType, Vec1::AbstractPetscVec, Vec2::AbstractPetscVec) 
Creates an index set containing the indices
where the vectors `Vec1` < `Vec2`

Collective

Input Parameters:
- `Vec1` - the first vector to compare
- `Vec2` - the second vector to compare

Output Parameter:
- `S` - The index set containing the indices i where vec1[i] < vec2[i]

Level: advanced

See also: `Vec`

# External Links
$(_doc_external("Vec/VecWhichLessThan"))
"""
function VecWhichLessThan(petsclib::PetscLibType, Vec1::AbstractPetscVec, Vec2::AbstractPetscVec)
    error("VecWhichLessThan: no generated method for these argument types")
end

@for_petsc function VecWhichLessThan(petsclib::$UnionPetscLib, Vec1::AbstractPetscVec, Vec2::AbstractPetscVec )
	S_ = Ref{CIS}()

    @chk ccall(
               (:VecWhichLessThan, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, Ptr{CIS}),
               Vec1, Vec2, S_,
              )

	S = IS(S_[], petsclib)

	return S
end 

"""
	VecZeroEntries(petsclib::PetscLibType, vec::AbstractPetscVec) 
puts a `0.0` in each element of a vector

Logically Collective

Input Parameter:
- `vec` - The vector

Level: beginner

See also: `Vec`, `VecCreate()`, `VecSetOptionsPrefix()`, `VecSet()`, `VecSetValues()`

# External Links
$(_doc_external("Vec/VecZeroEntries"))
"""
function VecZeroEntries(petsclib::PetscLibType, vec::AbstractPetscVec)
    error("VecZeroEntries: no generated method for these argument types")
end

@for_petsc function VecZeroEntries(petsclib::$UnionPetscLib, vec::AbstractPetscVec )

    @chk ccall(
               (:VecZeroEntries, $petsc_library),
               PetscErrorCode,
               (CVec,),
               vec,
              )


	return nothing
end 

