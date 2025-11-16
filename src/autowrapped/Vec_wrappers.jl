# autodefined type arguments for class Vec ------
mutable struct _n_ViennaCLVector end
const ViennaCLVector = Ptr{_n_ViennaCLVector}

mutable struct _p_PetscSection end
const PetscSection = Ptr{_p_PetscSection}

mutable struct _p_ISLocalToGlobalMapping end
const ISLocalToGlobalMapping = Ptr{_p_ISLocalToGlobalMapping}

# -------------------------------------------------------
"""
	max::PetscReal = VecMaxPointwiseDivide(petsclib::PetscLibType,x::PetscVec, y::PetscVec) 
Computes the maximum of the componentwise division `max = max_i abs(x[i]/y[i])`.

Logically Collective

Input Parameters:
- `x` - the numerators
- `y` - the denominators

Output Parameter:
- `max` - the result

Level: advanced

-seealso: [](ch_vectors), `Vec`, `VecPointwiseDivide()`, `VecPointwiseMult()`, `VecPointwiseMax()`, `VecPointwiseMin()`, `VecPointwiseMaxAbs()`

# External Links
$(_doc_external("Vec/VecMaxPointwiseDivide"))
"""
function VecMaxPointwiseDivide(petsclib::PetscLibType, x::PetscVec, y::PetscVec) end

@for_petsc function VecMaxPointwiseDivide(petsclib::$UnionPetscLib, x::PetscVec, y::PetscVec )
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
	val::PetscScalar = VecDot(petsclib::PetscLibType,x::PetscVec, y::PetscVec) 
Computes the vector dot product.

Collective

Input Parameters:
- `x` - first vector
- `y` - second vector

Output Parameter:
- `val` - the dot product

Level: intermediate

Notes for Users of Complex Numbers:
For complex vectors, `VecDot()` computes
-seealso: [](ch_vectors), `Vec`, `VecMDot()`, `VecTDot()`, `VecNorm()`, `VecDotBegin()`, `VecDotEnd()`, `VecDotRealPart()`

# External Links
$(_doc_external("Vec/VecDot"))
"""
function VecDot(petsclib::PetscLibType, x::PetscVec, y::PetscVec) end

@for_petsc function VecDot(petsclib::$UnionPetscLib, x::PetscVec, y::PetscVec )
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
	val::PetscReal = VecDotRealPart(petsclib::PetscLibType,x::PetscVec, y::PetscVec) 
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

-seealso: [](ch_vectors), `Vec`, `VecMDot()`, `VecTDot()`, `VecNorm()`, `VecDotBegin()`, `VecDotEnd()`, `VecDot()`, `VecDotNorm2()`

# External Links
$(_doc_external("Vec/VecDotRealPart"))
"""
function VecDotRealPart(petsclib::PetscLibType, x::PetscVec, y::PetscVec) end

@for_petsc function VecDotRealPart(petsclib::$UnionPetscLib, x::PetscVec, y::PetscVec )
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
	val::PetscReal = VecNorm(petsclib::PetscLibType,x::PetscVec, type::NormType) 
Computes the vector norm.

Collective

Input Parameters:
- `x`    - the vector
- `type` - the type of the norm requested

Output Parameter:
- `val` - the norm

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `NormType`, `VecDot()`, `VecTDot()`, `VecDotBegin()`, `VecDotEnd()`, `VecNormAvailable()`,
`VecNormBegin()`, `VecNormEnd()`, `NormType()`

# External Links
$(_doc_external("Vec/VecNorm"))
"""
function VecNorm(petsclib::PetscLibType, x::PetscVec, type::NormType) end

@for_petsc function VecNorm(petsclib::$UnionPetscLib, x::PetscVec, type::NormType )
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
	available::PetscBool,val::PetscReal = VecNormAvailable(petsclib::PetscLibType,x::PetscVec, type::NormType) 
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

-seealso: [](ch_vectors), `Vec`, `VecDot()`, `VecTDot()`, `VecNorm()`, `VecDotBegin()`, `VecDotEnd()`,
`VecNormBegin()`, `VecNormEnd()`

# External Links
$(_doc_external("Vec/VecNormAvailable"))
"""
function VecNormAvailable(petsclib::PetscLibType, x::PetscVec, type::NormType) end

@for_petsc function VecNormAvailable(petsclib::$UnionPetscLib, x::PetscVec, type::NormType )
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
	val::PetscReal = VecNormalize(petsclib::PetscLibType,x::PetscVec) 
Normalizes a vector by its 2

Collective

Input Parameter:
- `x` - the vector

Output Parameter:
- `val` - the vector norm before normalization. May be `NULL` if the value is not needed.

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecNorm()`, `NORM_2`, `NormType`

# External Links
$(_doc_external("Vec/VecNormalize"))
"""
function VecNormalize(petsclib::PetscLibType, x::PetscVec) end

@for_petsc function VecNormalize(petsclib::$UnionPetscLib, x::PetscVec )
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
	p::PetscInt,val::PetscReal = VecMax(petsclib::PetscLibType,x::PetscVec) 
Determines the vector component with maximum real part and its location.

Collective

Input Parameter:
- `x` - the vector

Output Parameters:
- `p`   - the index of `val` (pass `NULL` if you don't want this) in the vector
- `val` - the maximum component

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecNorm()`, `VecMin()`

# External Links
$(_doc_external("Vec/VecMax"))
"""
function VecMax(petsclib::PetscLibType, x::PetscVec) end

@for_petsc function VecMax(petsclib::$UnionPetscLib, x::PetscVec )
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
	p::PetscInt,val::PetscReal = VecMin(petsclib::PetscLibType,x::PetscVec) 
Determines the vector component with minimum real part and its location.

Collective

Input Parameter:
- `x` - the vector

Output Parameters:
- `p`   - the index of `val` (pass `NULL` if you don't want this location) in the vector
- `val` - the minimum component

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecMax()`

# External Links
$(_doc_external("Vec/VecMin"))
"""
function VecMin(petsclib::PetscLibType, x::PetscVec) end

@for_petsc function VecMin(petsclib::$UnionPetscLib, x::PetscVec )
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
	val::PetscScalar = VecTDot(petsclib::PetscLibType,x::PetscVec, y::PetscVec) 
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
-seealso: [](ch_vectors), `Vec`, `VecDot()`, `VecMTDot()`

# External Links
$(_doc_external("Vec/VecTDot"))
"""
function VecTDot(petsclib::PetscLibType, x::PetscVec, y::PetscVec) end

@for_petsc function VecTDot(petsclib::$UnionPetscLib, x::PetscVec, y::PetscVec )
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
	VecAXPY(petsclib::PetscLibType,y::PetscVec, alpha::PetscScalar, x::PetscVec) 
Computes `y = alpha x + y`.

Logically Collective

Input Parameters:
- `alpha` - the scalar
- `x`     - vector scale by `alpha`
- `y`     - vector accumulated into

Output Parameter:
- `y` - output vector

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecAYPX()`, `VecMAXPY()`, `VecWAXPY()`, `VecAXPBYPCZ()`, `VecAXPBY()`

# External Links
$(_doc_external("Vec/VecAXPY"))
"""
function VecAXPY(petsclib::PetscLibType, y::PetscVec, alpha::PetscScalar, x::PetscVec) end

@for_petsc function VecAXPY(petsclib::$UnionPetscLib, y::PetscVec, alpha::$PetscScalar, x::PetscVec )

    @chk ccall(
               (:VecAXPY, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscScalar, CVec),
               y, alpha, x,
              )


	return nothing
end 

"""
	VecAYPX(petsclib::PetscLibType,y::PetscVec, beta::PetscScalar, x::PetscVec) 
Computes `y = x + beta y`.

Logically Collective

Input Parameters:
- `beta` - the scalar
- `x`    - the unscaled vector
- `y`    - the vector to be scaled

Output Parameter:
- `y` - output vector

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecMAXPY()`, `VecWAXPY()`, `VecAXPY()`, `VecAXPBYPCZ()`, `VecAXPBY()`

# External Links
$(_doc_external("Vec/VecAYPX"))
"""
function VecAYPX(petsclib::PetscLibType, y::PetscVec, beta::PetscScalar, x::PetscVec) end

@for_petsc function VecAYPX(petsclib::$UnionPetscLib, y::PetscVec, beta::$PetscScalar, x::PetscVec )

    @chk ccall(
               (:VecAYPX, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscScalar, CVec),
               y, beta, x,
              )


	return nothing
end 

"""
	VecAXPBY(petsclib::PetscLibType,y::PetscVec, alpha::PetscScalar, beta::PetscScalar, x::PetscVec) 
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

-seealso: [](ch_vectors), `Vec`, `VecAYPX()`, `VecMAXPY()`, `VecWAXPY()`, `VecAXPY()`, `VecAXPBYPCZ()`

# External Links
$(_doc_external("Vec/VecAXPBY"))
"""
function VecAXPBY(petsclib::PetscLibType, y::PetscVec, alpha::PetscScalar, beta::PetscScalar, x::PetscVec) end

@for_petsc function VecAXPBY(petsclib::$UnionPetscLib, y::PetscVec, alpha::$PetscScalar, beta::$PetscScalar, x::PetscVec )

    @chk ccall(
               (:VecAXPBY, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscScalar, $PetscScalar, CVec),
               y, alpha, beta, x,
              )


	return nothing
end 

"""
	VecAXPBYPCZ(petsclib::PetscLibType,z::PetscVec, alpha::PetscScalar, beta::PetscScalar, gamma::PetscScalar, x::PetscVec, y::PetscVec) 
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

-seealso: [](ch_vectors), `Vec`, `VecAYPX()`, `VecMAXPY()`, `VecWAXPY()`, `VecAXPY()`, `VecAXPBY()`

# External Links
$(_doc_external("Vec/VecAXPBYPCZ"))
"""
function VecAXPBYPCZ(petsclib::PetscLibType, z::PetscVec, alpha::PetscScalar, beta::PetscScalar, gamma::PetscScalar, x::PetscVec, y::PetscVec) end

@for_petsc function VecAXPBYPCZ(petsclib::$UnionPetscLib, z::PetscVec, alpha::$PetscScalar, beta::$PetscScalar, gamma::$PetscScalar, x::PetscVec, y::PetscVec )

    @chk ccall(
               (:VecAXPBYPCZ, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscScalar, $PetscScalar, $PetscScalar, CVec, CVec),
               z, alpha, beta, gamma, x, y,
              )


	return nothing
end 

"""
	VecWAXPY(petsclib::PetscLibType,w::PetscVec, alpha::PetscScalar, x::PetscVec, y::PetscVec) 
Computes `w = alpha x + y`.

Logically Collective

Input Parameters:
- `alpha` - the scalar
- `x`     - first vector, multiplied by `alpha`
- `y`     - second vector

Output Parameter:
- `w` - the result

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecAXPY()`, `VecAYPX()`, `VecAXPBY()`, `VecMAXPY()`, `VecAXPBYPCZ()`

# External Links
$(_doc_external("Vec/VecWAXPY"))
"""
function VecWAXPY(petsclib::PetscLibType, w::PetscVec, alpha::PetscScalar, x::PetscVec, y::PetscVec) end

@for_petsc function VecWAXPY(petsclib::$UnionPetscLib, w::PetscVec, alpha::$PetscScalar, x::PetscVec, y::PetscVec )

    @chk ccall(
               (:VecWAXPY, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscScalar, CVec, CVec),
               w, alpha, x, y,
              )


	return nothing
end 

"""
	y::Vector{PetscScalar} = VecGetValues(petsclib::PetscLibType,x::PetscVec, ni::PetscInt, ix::Vector{PetscInt}) 
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

-seealso: [](ch_vectors), `Vec`, `VecAssemblyBegin()`, `VecAssemblyEnd()`, `VecSetValues()`

# External Links
$(_doc_external("Vec/VecGetValues"))
"""
function VecGetValues(petsclib::PetscLibType, x::AbstractPetscVec, ni::PetscInt, ix::Vector{PetscInt}) end

@for_petsc function VecGetValues(petsclib::$UnionPetscLib, x::AbstractPetscVec, ni::$PetscInt, ix::Vector{$PetscInt} )
	y = Vector{$PetscScalar}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:VecGetValues, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}),
               x, ni, ix, y,
              )


	return y
end 

"""
	val::Vector{PetscScalar} = VecMTDot(petsclib::PetscLibType,x::PetscVec, nv::PetscInt, y::Vector{PetscVec}) 
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
-seealso: [](ch_vectors), `Vec`, `VecMDot()`, `VecTDot()`

# External Links
$(_doc_external("Vec/VecMTDot"))
"""
function VecMTDot(petsclib::PetscLibType, x::PetscVec, nv::PetscInt, y::Vector{PetscVec}) end

@for_petsc function VecMTDot(petsclib::$UnionPetscLib, x::PetscVec, nv::$PetscInt, y::Vector{PetscVec} )
	val = Vector{$PetscScalar}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:VecMTDot, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{CVec}, Ptr{$PetscScalar}),
               x, nv, y, val,
              )


	return val
end 

"""
	val::Vector{PetscScalar} = VecMDot(petsclib::PetscLibType,x::PetscVec, nv::PetscInt, y::Vector{PetscVec}) 
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
-seealso: [](ch_vectors), `Vec`, `VecMTDot()`, `VecDot()`, `VecDuplicateVecs()`

# External Links
$(_doc_external("Vec/VecMDot"))
"""
function VecMDot(petsclib::PetscLibType, x::PetscVec, nv::PetscInt, y::Vector{PetscVec}) end

@for_petsc function VecMDot(petsclib::$UnionPetscLib, x::PetscVec, nv::$PetscInt, y::Vector{PetscVec} )
	val = Vector{$PetscScalar}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:VecMDot, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{CVec}, Ptr{$PetscScalar}),
               x, nv, y, val,
              )


	return val
end 

"""
	VecMAXPY(petsclib::PetscLibType,y::PetscVec, nv::PetscInt, alpha::Vector{PetscScalar}, x::Vector{PetscVec}) 
Computes `y = y + sum alpha[i] x[i]`

Logically Collective

Input Parameters:
- `nv`    - number of scalars and `x` vectors
- `alpha` - array of scalars
- `y`     - one vector
- `x`     - array of vectors

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecMAXPBY()`,`VecAYPX()`, `VecWAXPY()`, `VecAXPY()`, `VecAXPBYPCZ()`, `VecAXPBY()`, `VecDuplicateVecs()`

# External Links
$(_doc_external("Vec/VecMAXPY"))
"""
function VecMAXPY(petsclib::PetscLibType, y::PetscVec, nv::PetscInt, alpha::Vector{PetscScalar}, x::Vector{PetscVec}) end

@for_petsc function VecMAXPY(petsclib::$UnionPetscLib, y::PetscVec, nv::$PetscInt, alpha::Vector{$PetscScalar}, x::Vector{PetscVec} )

    @chk ccall(
               (:VecMAXPY, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{$PetscScalar}, Ptr{CVec}),
               y, nv, alpha, x,
              )


	return nothing
end 

"""
	VecMAXPBY(petsclib::PetscLibType,y::PetscVec, nv::PetscInt, alpha::Vector{PetscScalar}, beta::PetscScalar, x::Vector{PetscVec}) 
Computes `y = beta y + sum alpha[i] x[i]`

Logically Collective

Input Parameters:
- `nv`    - number of scalars and `x` vectors
- `alpha` - array of scalars
- `beta`  - scalar
- `y`     - one vector
- `x`     - array of vectors

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecMAXPY()`, `VecAYPX()`, `VecWAXPY()`, `VecAXPY()`, `VecAXPBYPCZ()`, `VecAXPBY()`

# External Links
$(_doc_external("Vec/VecMAXPBY"))
"""
function VecMAXPBY(petsclib::PetscLibType, y::PetscVec, nv::PetscInt, alpha::Vector{PetscScalar}, beta::PetscScalar, x::Vector{PetscVec}) end

@for_petsc function VecMAXPBY(petsclib::$UnionPetscLib, y::PetscVec, nv::$PetscInt, alpha::Vector{$PetscScalar}, beta::$PetscScalar, x::Vector{PetscVec} )

    @chk ccall(
               (:VecMAXPBY, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{$PetscScalar}, $PetscScalar, Ptr{CVec}),
               y, nv, alpha, beta, x,
              )


	return nothing
end 

"""
	VecConcatenate(petsclib::PetscLibType,nx::PetscInt, X::Vector{PetscVec}, Y::PetscVec, x_is::Vector{IS}) 
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

-seealso: [](ch_vectors), `Vec`, `VECNEST`, `VECSCATTER`, `VecScatterCreate()`

# External Links
$(_doc_external("Vec/VecConcatenate"))
"""
function VecConcatenate(petsclib::PetscLibType, nx::PetscInt, X::Vector{PetscVec}, Y::PetscVec, x_is::Vector{IS}) end

@for_petsc function VecConcatenate(petsclib::$UnionPetscLib, nx::$PetscInt, X::Vector{PetscVec}, Y::PetscVec, x_is::Vector{IS} )
	Y_ = Ref(Y.ptr)
	x_is_ = Ref(pointer(x_is))

    @chk ccall(
               (:VecConcatenate, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{CVec}, Ptr{CVec}, Ptr{Ptr{IS}}),
               nx, X, Y_, x_is_,
              )

	Y.ptr = C_NULL

	return nothing
end 

"""
	VecGetSubVector(petsclib::PetscLibType,X::PetscVec, is::IS, Y::PetscVec) 
Gets a vector representing part of another vector

Collective

Input Parameters:
- `X`  - vector from which to extract a subvector
- `is` - index set representing portion of `X` to extract

Output Parameter:
- `Y` - subvector corresponding to `is`

Level: advanced

-seealso: [](ch_vectors), `Vec`, `IS`, `VECNEST`, `MatCreateSubMatrix()`

# External Links
$(_doc_external("Vec/VecGetSubVector"))
"""
function VecGetSubVector(petsclib::PetscLibType, X::PetscVec, is::IS, Y::PetscVec) end

@for_petsc function VecGetSubVector(petsclib::$UnionPetscLib, X::PetscVec, is::IS, Y::PetscVec )
	Y_ = Ref(Y.ptr)

    @chk ccall(
               (:VecGetSubVector, $petsc_library),
               PetscErrorCode,
               (CVec, CIS, Ptr{CVec}),
               X, is, Y_,
              )

	Y.ptr = C_NULL

	return nothing
end 

"""
	VecRestoreSubVector(petsclib::PetscLibType,X::PetscVec, is::IS, Y::PetscVec) 
Restores a subvector extracted using `VecGetSubVector()`

Collective

Input Parameters:
- `X`  - vector from which subvector was obtained
- `is` - index set representing the subset of `X`
- `Y`  - subvector being restored

Level: advanced

-seealso: [](ch_vectors), `Vec`, `IS`, `VecGetSubVector()`

# External Links
$(_doc_external("Vec/VecRestoreSubVector"))
"""
function VecRestoreSubVector(petsclib::PetscLibType, X::PetscVec, is::IS, Y::PetscVec) end

@for_petsc function VecRestoreSubVector(petsclib::$UnionPetscLib, X::PetscVec, is::IS, Y::PetscVec )
	Y_ = Ref(Y.ptr)

    @chk ccall(
               (:VecRestoreSubVector, $petsc_library),
               PetscErrorCode,
               (CVec, CIS, Ptr{CVec}),
               X, is, Y_,
              )

	Y.ptr = C_NULL

	return nothing
end 

"""
	w::PetscVec = VecCreateLocalVector(petsclib::PetscLibType,v::PetscVec) 
Creates a vector object suitable for use with `VecGetLocalVector()` and friends. You must call `VecDestroy()` when the
vector is no longer needed.

Not Collective.

Input Parameter:
- `v` - The vector for which the local vector is desired.

Output Parameter:
- `w` - Upon exit this contains the local vector.

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecGetLocalVectorRead()`, `VecRestoreLocalVectorRead()`, `VecGetLocalVector()`, `VecRestoreLocalVector()`

# External Links
$(_doc_external("Vec/VecCreateLocalVector"))
"""
function VecCreateLocalVector(petsclib::PetscLibType, v::PetscVec) end

@for_petsc function VecCreateLocalVector(petsclib::$UnionPetscLib, v::PetscVec )
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
	VecGetLocalVectorRead(petsclib::PetscLibType,v::PetscVec, w::PetscVec) 
Maps the local portion of a vector into a
vector.

Not Collective.

Input Parameter:
- `v` - The vector for which the local vector is desired.

Output Parameter:
- `w` - Upon exit this contains the local vector.

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecCreateLocalVector()`, `VecRestoreLocalVectorRead()`, `VecGetLocalVector()`, `VecGetArrayRead()`, `VecGetArray()`

# External Links
$(_doc_external("Vec/VecGetLocalVectorRead"))
"""
function VecGetLocalVectorRead(petsclib::PetscLibType, v::PetscVec, w::PetscVec) end

@for_petsc function VecGetLocalVectorRead(petsclib::$UnionPetscLib, v::PetscVec, w::PetscVec )

    @chk ccall(
               (:VecGetLocalVectorRead, $petsc_library),
               PetscErrorCode,
               (CVec, CVec),
               v, w,
              )


	return nothing
end 

"""
	VecRestoreLocalVectorRead(petsclib::PetscLibType,v::PetscVec, w::PetscVec) 
Unmaps the local portion of a vector
previously mapped into a vector using `VecGetLocalVectorRead()`.

Not Collective.

Input Parameters:
- `v` - The local portion of this vector was previously mapped into `w` using `VecGetLocalVectorRead()`.
- `w` - The vector into which the local portion of `v` was mapped.

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecCreateLocalVector()`, `VecGetLocalVectorRead()`, `VecGetLocalVector()`, `VecGetArrayRead()`, `VecGetArray()`

# External Links
$(_doc_external("Vec/VecRestoreLocalVectorRead"))
"""
function VecRestoreLocalVectorRead(petsclib::PetscLibType, v::PetscVec, w::PetscVec) end

@for_petsc function VecRestoreLocalVectorRead(petsclib::$UnionPetscLib, v::PetscVec, w::PetscVec )

    @chk ccall(
               (:VecRestoreLocalVectorRead, $petsc_library),
               PetscErrorCode,
               (CVec, CVec),
               v, w,
              )


	return nothing
end 

"""
	VecGetLocalVector(petsclib::PetscLibType,v::PetscVec, w::PetscVec) 
Maps the local portion of a vector into a
vector.

Collective

Input Parameter:
- `v` - The vector for which the local vector is desired.

Output Parameter:
- `w` - Upon exit this contains the local vector.

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecCreateLocalVector()`, `VecRestoreLocalVector()`, `VecGetLocalVectorRead()`, `VecGetArrayRead()`, `VecGetArray()`

# External Links
$(_doc_external("Vec/VecGetLocalVector"))
"""
function VecGetLocalVector(petsclib::PetscLibType, v::PetscVec, w::PetscVec) end

@for_petsc function VecGetLocalVector(petsclib::$UnionPetscLib, v::PetscVec, w::PetscVec )

    @chk ccall(
               (:VecGetLocalVector, $petsc_library),
               PetscErrorCode,
               (CVec, CVec),
               v, w,
              )


	return nothing
end 

"""
	VecRestoreLocalVector(petsclib::PetscLibType,v::PetscVec, w::PetscVec) 
Unmaps the local portion of a vector
previously mapped into a vector using `VecGetLocalVector()`.

Logically Collective.

Input Parameters:
- `v` - The local portion of this vector was previously mapped into `w` using `VecGetLocalVector()`.
- `w` - The vector into which the local portion of `v` was mapped.

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecCreateLocalVector()`, `VecGetLocalVector()`, `VecGetLocalVectorRead()`, `VecRestoreLocalVectorRead()`, `LocalVectorRead()`, `VecGetArrayRead()`, `VecGetArray()`

# External Links
$(_doc_external("Vec/VecRestoreLocalVector"))
"""
function VecRestoreLocalVector(petsclib::PetscLibType, v::PetscVec, w::PetscVec) end

@for_petsc function VecRestoreLocalVector(petsclib::$UnionPetscLib, v::PetscVec, w::PetscVec )

    @chk ccall(
               (:VecRestoreLocalVector, $petsc_library),
               PetscErrorCode,
               (CVec, CVec),
               v, w,
              )


	return nothing
end 

"""
	a::Vector{PetscScalar} = VecGetArray(petsclib::PetscLibType,x::PetscVec) 
Returns a pointer to a contiguous array that contains this
MPI processes's portion of the vector data

Logically Collective

Input Parameter:
- `x` - the vector

Output Parameter:
- `a` - location to put pointer to the array

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecRestoreArray()`, `VecGetArrayRead()`, `VecGetArrays()`, `VecPlaceArray()`, `VecGetArray2d()`,
`VecGetArrayPair()`, `VecRestoreArrayPair()`, `VecGetArrayWrite()`, `VecRestoreArrayWrite()`, `VecGetArrayAndMemType()`

# External Links
$(_doc_external("Vec/VecGetArray"))
"""
function VecGetArray(petsclib::PetscLibType, x::PetscVec) end

@for_petsc function VecGetArray(petsclib::$UnionPetscLib, x::PetscVec )
	a_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:VecGetArray, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Ptr{$PetscScalar}}),
               x, a_,
              )

	a = unsafe_wrap(Array, a_[], VecGetLocalSize(petsclib, x); own = false)

	return a
end 

"""
	VecRestoreArray(petsclib::PetscLibType,x::PetscVec, a::Vector{PetscScalar}) 
Restores a vector after `VecGetArray()` has been called and the array is no longer needed

Logically Collective

Input Parameters:
- `x` - the vector
- `a` - location of pointer to array obtained from `VecGetArray()`

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArrayRead()`, `VecRestoreArrays()`, `VecPlaceArray()`, `VecRestoreArray2d()`,
`VecGetArrayPair()`, `VecRestoreArrayPair()`

# External Links
$(_doc_external("Vec/VecRestoreArray"))
"""
function VecRestoreArray(petsclib::PetscLibType, x::PetscVec, a::Vector{PetscScalar}) end

@for_petsc function VecRestoreArray(petsclib::$UnionPetscLib, x::PetscVec, a::Vector{$PetscScalar} )
	a_ = Ref(pointer(a))

    @chk ccall(
               (:VecRestoreArray, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Ptr{$PetscScalar}}),
               x, a_,
              )


	return nothing
end 

"""
	a::Vector{PetscScalar} = VecGetArrayRead(petsclib::PetscLibType,x::PetscVec) 
Get read

Not Collective

Input Parameter:
- `x` - the vector

Output Parameter:
- `a` - the array

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrayPair()`, `VecRestoreArrayPair()`,
`VecGetArrayAndMemType()`

# External Links
$(_doc_external("Vec/VecGetArrayRead"))
"""
function VecGetArrayRead(petsclib::PetscLibType, x::AbstractPetscVec) end

@for_petsc function VecGetArrayRead(petsclib::$UnionPetscLib, x::AbstractPetscVec )
	a_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:VecGetArrayRead, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Ptr{$PetscScalar}}),
               x, a_,
              )

	a = unsafe_wrap(Array, a_[], VecGetLocalSize(petsclib, x); own = false)

	return a
end 

"""
	VecRestoreArrayRead(petsclib::PetscLibType,x::PetscVec, a::Vector{PetscScalar}) 
Restore array obtained with `VecGetArrayRead()`

Not Collective

Input Parameters:
- `x` - the vector
- `a` - the array

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrayPair()`, `VecRestoreArrayPair()`

# External Links
$(_doc_external("Vec/VecRestoreArrayRead"))
"""
function VecRestoreArrayRead(petsclib::PetscLibType, x::AbstractPetscVec, a::Vector{PetscScalar}) end

@for_petsc function VecRestoreArrayRead(petsclib::$UnionPetscLib, x::AbstractPetscVec, a::Vector{$PetscScalar} )
	a_ = Ref(pointer(a))

    @chk ccall(
               (:VecRestoreArrayRead, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Ptr{$PetscScalar}}),
               x, a_,
              )


	return nothing
end 

"""
	a::Vector{PetscScalar} = VecGetArrayWrite(petsclib::PetscLibType,x::PetscVec) 
Returns a pointer to a contiguous array that WILL contain this
MPI processes's portion of the vector data.

Logically Collective

Input Parameter:
- `x` - the vector

Output Parameter:
- `a` - location to put pointer to the array

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecRestoreArray()`, `VecGetArrayRead()`, `VecGetArrays()`, `VecPlaceArray()`, `VecGetArray2d()`,
`VecGetArrayPair()`, `VecRestoreArrayPair()`, `VecGetArray()`, `VecRestoreArrayWrite()`, `VecGetArrayAndMemType()`

# External Links
$(_doc_external("Vec/VecGetArrayWrite"))
"""
function VecGetArrayWrite(petsclib::PetscLibType, x::AbstractPetscVec) end

@for_petsc function VecGetArrayWrite(petsclib::$UnionPetscLib, x::AbstractPetscVec )
	a_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:VecGetArrayWrite, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Ptr{$PetscScalar}}),
               x, a_,
              )

	a = unsafe_wrap(Array, a_[], VecGetLocalSize(petsclib, x); own = false)

	return a
end 

"""
	VecRestoreArrayWrite(petsclib::PetscLibType,x::PetscVec, a::Vector{PetscScalar}) 
Restores a vector after `VecGetArrayWrite()` has been called.

Logically Collective

Input Parameters:
- `x` - the vector
- `a` - location of pointer to array obtained from `VecGetArray()`

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArrayRead()`, `VecRestoreArrays()`, `VecPlaceArray()`, `VecRestoreArray2d()`,
`VecGetArrayPair()`, `VecRestoreArrayPair()`, `VecGetArrayWrite()`

# External Links
$(_doc_external("Vec/VecRestoreArrayWrite"))
"""
function VecRestoreArrayWrite(petsclib::PetscLibType, x::AbstractPetscVec, a::Vector{PetscScalar}) end

@for_petsc function VecRestoreArrayWrite(petsclib::$UnionPetscLib, x::AbstractPetscVec, a::Vector{$PetscScalar} )
	a_ = Ref(pointer(a))

    @chk ccall(
               (:VecRestoreArrayWrite, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Ptr{$PetscScalar}}),
               x, a_,
              )


	return nothing
end 

"""
	a::Vector{PetscScalar} = VecGetArrays(petsclib::PetscLibType,x::Vector{PetscVec}, n::PetscInt) 
Returns a pointer to the arrays in a set of vectors
that were created by a call to `VecDuplicateVecs()`.

Logically Collective; No Fortran Support

Input Parameters:
- `x` - the vectors
- `n` - the number of vectors

Output Parameter:
- `a` - location to put pointer to the array

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArrays()`

# External Links
$(_doc_external("Vec/VecGetArrays"))
"""
function VecGetArrays(petsclib::PetscLibType, x::Vector{PetscVec}, n::PetscInt) end

@for_petsc function VecGetArrays(petsclib::$UnionPetscLib, x::Vector{PetscVec}, n::$PetscInt )

    @chk ccall(
               (:VecGetArrays, $petsc_library),
               PetscErrorCode,
               (Ptr{CVec}, $PetscInt, Ptr{$PetscScalar}),
               x, n, a,
              )


	return a
end 

"""
	VecRestoreArrays(petsclib::PetscLibType,x::Vector{PetscVec}, n::PetscInt, a::Vector{PetscScalar}) 
Restores a group of vectors after `VecGetArrays()`
has been called.

Logically Collective; No Fortran Support

Input Parameters:
- `x` - the vector
- `n` - the number of vectors
- `a` - location of pointer to arrays obtained from `VecGetArrays()`

-seealso: [](ch_vectors), `Vec`, `VecGetArrays()`, `VecRestoreArray()`

# External Links
$(_doc_external("Vec/VecRestoreArrays"))
"""
function VecRestoreArrays(petsclib::PetscLibType, x::Vector{PetscVec}, n::PetscInt, a::Vector{PetscScalar}) end

@for_petsc function VecRestoreArrays(petsclib::$UnionPetscLib, x::Vector{PetscVec}, n::$PetscInt, a::Vector{$PetscScalar} )

    @chk ccall(
               (:VecRestoreArrays, $petsc_library),
               PetscErrorCode,
               (Ptr{CVec}, $PetscInt, Ptr{$PetscScalar}),
               x, n, a,
              )


	return nothing
end 

"""
	a::Vector{PetscScalar},mtype::PetscMemType = VecGetArrayAndMemType(petsclib::PetscLibType,x::PetscVec) 
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

-seealso: [](ch_vectors), `Vec`, `VecRestoreArrayAndMemType()`, `VecGetArrayReadAndMemType()`, `VecGetArrayWriteAndMemType()`, `VecRestoreArray()`, `VecGetArrayRead()`, `VecGetArrays()`,
`VecPlaceArray()`, `VecGetArray2d()`, `VecGetArrayPair()`, `VecRestoreArrayPair()`, `VecGetArrayWrite()`, `VecRestoreArrayWrite()`

# External Links
$(_doc_external("Vec/VecGetArrayAndMemType"))
"""
function VecGetArrayAndMemType(petsclib::PetscLibType, x::PetscVec) end

@for_petsc function VecGetArrayAndMemType(petsclib::$UnionPetscLib, x::PetscVec )
	a_ = Ref{Ptr{$PetscScalar}}()
	mtype_ = Ref{PetscMemType}()

    @chk ccall(
               (:VecGetArrayAndMemType, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Ptr{$PetscScalar}}, Ptr{PetscMemType}),
               x, a_, mtype_,
              )

	a = unsafe_wrap(Array, a_[], VecGetLocalSize(petsclib, x); own = false)
	mtype = unsafe_string(mtype_[])

	return a,mtype
end 

"""
	VecRestoreArrayAndMemType(petsclib::PetscLibType,x::PetscVec, a::Vector{PetscScalar}) 
Restores a vector after `VecGetArrayAndMemType()` has been called.

Logically Collective; No Fortran Support

Input Parameters:
- `x` - the vector
- `a` - location of pointer to array obtained from `VecGetArrayAndMemType()`

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecGetArrayAndMemType()`, `VecGetArray()`, `VecRestoreArrayRead()`, `VecRestoreArrays()`,
`VecPlaceArray()`, `VecRestoreArray2d()`, `VecGetArrayPair()`, `VecRestoreArrayPair()`

# External Links
$(_doc_external("Vec/VecRestoreArrayAndMemType"))
"""
function VecRestoreArrayAndMemType(petsclib::PetscLibType, x::PetscVec, a::Vector{PetscScalar}) end

@for_petsc function VecRestoreArrayAndMemType(petsclib::$UnionPetscLib, x::PetscVec, a::Vector{$PetscScalar} )
	a_ = Ref(pointer(a))

    @chk ccall(
               (:VecRestoreArrayAndMemType, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Ptr{$PetscScalar}}),
               x, a_,
              )


	return nothing
end 

"""
	a::Vector{PetscScalar},mtype::PetscMemType = VecGetArrayReadAndMemType(petsclib::PetscLibType,x::PetscVec) 
Like `VecGetArrayRead()`, but if the input vector is a device vector, it will return a read
The returned pointer is guaranteed to point to up-to-date data. For host vectors, it functions as `VecGetArrayRead()`.

Not Collective; No Fortran Support

Input Parameter:
- `x` - the vector

Output Parameters:
- `a`     - the array
- `mtype` - memory type of the array

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecRestoreArrayReadAndMemType()`, `VecGetArrayAndMemType()`, `VecGetArrayWriteAndMemType()`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrayPair()`, `VecRestoreArrayPair()`

# External Links
$(_doc_external("Vec/VecGetArrayReadAndMemType"))
"""
function VecGetArrayReadAndMemType(petsclib::PetscLibType, x::PetscVec) end

@for_petsc function VecGetArrayReadAndMemType(petsclib::$UnionPetscLib, x::PetscVec )
	a_ = Ref{Ptr{$PetscScalar}}()
	mtype_ = Ref{PetscMemType}()

    @chk ccall(
               (:VecGetArrayReadAndMemType, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Ptr{$PetscScalar}}, Ptr{PetscMemType}),
               x, a_, mtype_,
              )

	a = unsafe_wrap(Array, a_[], VecGetLocalSize(petsclib, x); own = false)
	mtype = unsafe_string(mtype_[])

	return a,mtype
end 

"""
	VecRestoreArrayReadAndMemType(petsclib::PetscLibType,x::PetscVec, a::Vector{PetscScalar}) 
Restore array obtained with `VecGetArrayReadAndMemType()`

Not Collective; No Fortran Support

Input Parameters:
- `x` - the vector
- `a` - the array

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecGetArrayReadAndMemType()`, `VecRestoreArrayAndMemType()`, `VecRestoreArrayWriteAndMemType()`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrayPair()`, `VecRestoreArrayPair()`

# External Links
$(_doc_external("Vec/VecRestoreArrayReadAndMemType"))
"""
function VecRestoreArrayReadAndMemType(petsclib::PetscLibType, x::PetscVec, a::Vector{PetscScalar}) end

@for_petsc function VecRestoreArrayReadAndMemType(petsclib::$UnionPetscLib, x::PetscVec, a::Vector{$PetscScalar} )
	a_ = Ref(pointer(a))

    @chk ccall(
               (:VecRestoreArrayReadAndMemType, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Ptr{$PetscScalar}}),
               x, a_,
              )


	return nothing
end 

"""
	a::Vector{PetscScalar},mtype::PetscMemType = VecGetArrayWriteAndMemType(petsclib::PetscLibType,x::PetscVec) 
Like `VecGetArrayWrite()`, but if this is a device vector it will always return
a device pointer to the device memory that contains this processor's portion of the vector data.

Logically Collective; No Fortran Support

Input Parameter:
- `x` - the vector

Output Parameters:
- `a`     - the array
- `mtype` - memory type of the array

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecRestoreArrayWriteAndMemType()`, `VecGetArrayReadAndMemType()`, `VecGetArrayAndMemType()`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrayPair()`, `VecRestoreArrayPair()`,

# External Links
$(_doc_external("Vec/VecGetArrayWriteAndMemType"))
"""
function VecGetArrayWriteAndMemType(petsclib::PetscLibType, x::PetscVec) end

@for_petsc function VecGetArrayWriteAndMemType(petsclib::$UnionPetscLib, x::PetscVec )
	a_ = Ref{Ptr{$PetscScalar}}()
	mtype_ = Ref{PetscMemType}()

    @chk ccall(
               (:VecGetArrayWriteAndMemType, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Ptr{$PetscScalar}}, Ptr{PetscMemType}),
               x, a_, mtype_,
              )

	a = unsafe_wrap(Array, a_[], VecGetLocalSize(petsclib, x); own = false)
	mtype = unsafe_string(mtype_[])

	return a,mtype
end 

"""
	VecRestoreArrayWriteAndMemType(petsclib::PetscLibType,x::PetscVec, a::Vector{PetscScalar}) 
Restore array obtained with `VecGetArrayWriteAndMemType()`

Logically Collective; No Fortran Support

Input Parameters:
- `x` - the vector
- `a` - the array

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecGetArrayWriteAndMemType()`, `VecRestoreArrayAndMemType()`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrayPair()`, `VecRestoreArrayPair()`

# External Links
$(_doc_external("Vec/VecRestoreArrayWriteAndMemType"))
"""
function VecRestoreArrayWriteAndMemType(petsclib::PetscLibType, x::PetscVec, a::Vector{PetscScalar}) end

@for_petsc function VecRestoreArrayWriteAndMemType(petsclib::$UnionPetscLib, x::PetscVec, a::Vector{$PetscScalar} )
	a_ = Ref(pointer(a))

    @chk ccall(
               (:VecRestoreArrayWriteAndMemType, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Ptr{$PetscScalar}}),
               x, a_,
              )


	return nothing
end 

"""
	VecPlaceArray(petsclib::PetscLibType,vec::PetscVec, array::Vector{PetscScalar}) 
Allows one to replace the array in a vector with an
array provided by the user. This is useful to avoid copying an array
into a vector.

Logically Collective; No Fortran Support

Input Parameters:
- `vec`   - the vector
- `array` - the array

Level: developer

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecReplaceArray()`, `VecResetArray()`

# External Links
$(_doc_external("Vec/VecPlaceArray"))
"""
function VecPlaceArray(petsclib::PetscLibType, vec::PetscVec, array::Vector{PetscScalar}) end

@for_petsc function VecPlaceArray(petsclib::$UnionPetscLib, vec::PetscVec, array::Vector{$PetscScalar} )

    @chk ccall(
               (:VecPlaceArray, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{$PetscScalar}),
               vec, array,
              )


	return nothing
end 

"""
	VecReplaceArray(petsclib::PetscLibType,vec::PetscVec, array::Vector{PetscScalar}) 
Allows one to replace the array in a vector with an
array provided by the user. This is useful to avoid copying an array
into a vector.

Logically Collective; No Fortran Support

Input Parameters:
- `vec`   - the vector
- `array` - the array

Level: developer

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecPlaceArray()`, `VecResetArray()`

# External Links
$(_doc_external("Vec/VecReplaceArray"))
"""
function VecReplaceArray(petsclib::PetscLibType, vec::PetscVec, array::Vector{PetscScalar}) end

@for_petsc function VecReplaceArray(petsclib::$UnionPetscLib, vec::PetscVec, array::Vector{$PetscScalar} )

    @chk ccall(
               (:VecReplaceArray, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{$PetscScalar}),
               vec, array,
              )


	return nothing
end 

"""
	a::Vector{PetscScalar} = VecGetArray2d(petsclib::PetscLibType,x::PetscVec, m::PetscInt, n::PetscInt, mstart::PetscInt, nstart::PetscInt) 
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

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrays()`, `VecPlaceArray()`,
`VecRestoreArray2d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`, `VecGetArray3d()`, `VecRestoreArray3d()`,
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecGetArray2d"))
"""
function VecGetArray2d(petsclib::PetscLibType, x::PetscVec, m::PetscInt, n::PetscInt, mstart::PetscInt, nstart::PetscInt) end

@for_petsc function VecGetArray2d(petsclib::$UnionPetscLib, x::PetscVec, m::$PetscInt, n::$PetscInt, mstart::$PetscInt, nstart::$PetscInt )

    arr_ptr = Ref{Ptr{Ptr{$PetscScalar}}}()

    @chk ccall(
               (:VecGetArray2d, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{$PetscScalar}}}),
               x, m, n, mstart, nstart, arr_ptr,
              )

    # Assume contiguous storage, use first row pointer
    data_ptr = unsafe_load(arr_ptr[])
    mat = unsafe_wrap(Array, data_ptr, (m, n))

	return PetscArray(mat,[data_ptr])
end 

"""
	a::Vector{PetscScalar} = VecGetArray2dWrite(petsclib::PetscLibType,x::PetscVec, m::PetscInt, n::PetscInt, mstart::PetscInt, nstart::PetscInt) 
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

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrays()`, `VecPlaceArray()`,
`VecRestoreArray2d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`, `VecGetArray3d()`, `VecRestoreArray3d()`,
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecGetArray2dWrite"))
"""
function VecGetArray2dWrite(petsclib::PetscLibType, x::PetscVec, m::PetscInt, n::PetscInt, mstart::PetscInt, nstart::PetscInt) end

@for_petsc function VecGetArray2dWrite(petsclib::$UnionPetscLib, x::PetscVec, m::$PetscInt, n::$PetscInt, mstart::$PetscInt, nstart::$PetscInt )

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

	return PetscArray(mat,[data_ptr])            
end 

"""
	VecRestoreArray2d(petsclib::PetscLibType,x::PetscVec, m::PetscInt, n::PetscInt, mstart::PetscInt, nstart::PetscInt, a::PetscArray{PetscScalar, 2}) 
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

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecRestoreArrays()`, `VecPlaceArray()`,
`VecGetArray2d()`, `VecGetArray3d()`, `VecRestoreArray3d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecRestoreArray2d"))
"""
function VecRestoreArray2d(petsclib::PetscLibType, x::PetscVec, m::PetscInt, n::PetscInt, mstart::PetscInt, nstart::PetscInt, a::PetscArray{PetscScalar, 2}) end

@for_petsc function VecRestoreArray2d(petsclib::$UnionPetscLib, x::PetscVec, m::$PetscInt, n::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, a::PetscArray{$PetscScalar, 2} )
	if a.ptr[]  != C_NULL 

        @chk ccall(
                (:VecRestoreArray2d, $petsc_library),
                PetscErrorCode,
                (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{$PetscScalar}}}),
                x, m, n, mstart, nstart, a,
                )

		a.ptr[] = C_NULL
	else
		error("The input array is already restored")
	end


	return nothing
end 



"""
	VecRestoreArray2dWrite(petsclib::PetscLibType,x::PetscVec, m::PetscInt, n::PetscInt, mstart::PetscInt, nstart::PetscInt, a::PetscArray{PetscScalar, 2}) 
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

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecRestoreArrays()`, `VecPlaceArray()`,
`VecGetArray2d()`, `VecGetArray3d()`, `VecRestoreArray3d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecRestoreArray2dWrite"))
"""
function VecRestoreArray2dWrite(petsclib::PetscLibType, x::PetscVec, m::PetscInt, n::PetscInt, mstart::PetscInt, nstart::PetscInt, a::PetscArray{PetscScalar, 2}) end

@for_petsc function VecRestoreArray2dWrite(petsclib::$UnionPetscLib, x::PetscVec, m::$PetscInt, n::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, a::PetscArray{$PetscScalar, 2} )
	if a.ptr[]  != C_NULL  

        @chk ccall(
                (:VecRestoreArray2dWrite, $petsc_library),
                PetscErrorCode,
                (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{$PetscScalar}}}),
                x, m, n, mstart, nstart, a,
                )

		a.ptr[] = C_NULL
	else
		error("The input array is already restored")
	end


	return nothing
end 



"""
	a::PetscArray = VecGetArray1d(petsclib::PetscLibType,x::PetscVec, m::PetscInt, mstart::PetscInt) 
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

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrays()`, `VecPlaceArray()`,
`VecRestoreArray2d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`, `VecGetArray3d()`, `VecRestoreArray3d()`,
`VecGetArray2d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecGetArray1d"))
"""
function VecGetArray1d(petsclib::PetscLibType, x::PetscVec, m::PetscInt, mstart::PetscInt) end

@for_petsc function VecGetArray1d(petsclib::$UnionPetscLib, x::PetscVec, m::$PetscInt, mstart::$PetscInt )
	a_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:VecGetArray1d, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, Ref{Ptr{$PetscScalar}}),
               x, m, mstart, a_,
              )

	data_ptr = unsafe_load(a_[])
	mat = unsafe_wrap(Array, data_ptr, m) 
	a = PetscArray(mat,[data_ptr]) 

	return a
end 


"""
	a::PetscArray = VecGetArray1dWrite(petsclib::PetscLibType,x::PetscVec, m::PetscInt, mstart::PetscInt) 
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

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrays()`, `VecPlaceArray()`,
`VecRestoreArray2d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`, `VecGetArray3d()`, `VecRestoreArray3d()`,
`VecGetArray2d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecGetArray1dWrite"))
"""
function VecGetArray1dWrite(petsclib::PetscLibType, x::PetscVec, m::PetscInt, mstart::PetscInt) end

@for_petsc function VecGetArray1dWrite(petsclib::$UnionPetscLib, x::PetscVec, m::$PetscInt, mstart::$PetscInt )
	a_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:VecGetArray1dWrite, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, Ref{Ptr{$PetscScalar}}),
               x, m, mstart, a_,
              )

	data_ptr = unsafe_load(a_[])
	mat = unsafe_wrap(Array, data_ptr, m) 
	a = PetscArray(mat,[data_ptr]) 

	return a
end 



"""
	VecRestoreArray1d(petsclib::PetscLibType,x::PetscVec, m::PetscInt, mstart::PetscInt, a::PetscArray{PetscScalar, 1}) 
Restores a vector after `VecGetArray1d()` has been called.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of two dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `a`      - location of pointer to array obtained from `VecGetArray1d()`

Level: developer

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecRestoreArrays()`, `VecPlaceArray()`,
`VecGetArray2d()`, `VecGetArray3d()`, `VecRestoreArray3d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`
`VecGetArray1d()`, `VecRestoreArray2d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecRestoreArray1d"))
"""
function VecRestoreArray1d(petsclib::PetscLibType, x::PetscVec, m::PetscInt, mstart::PetscInt, a::PetscArray{PetscScalar, 1}) end

@for_petsc function VecRestoreArray1d(petsclib::$UnionPetscLib, x::PetscVec, m::$PetscInt, mstart::$PetscInt, a::PetscArray{$PetscScalar, 1} )
	if a.ptr[]  != C_NULL 

    @chk ccall(
               (:VecRestoreArray1d, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, Ref{Ptr{$PetscScalar}}),
               x, m, mstart, a_,
              )

		a.ptr[] = C_NULL
	else
		error("The input array is already restored")
	end


	return nothing
end 


"""
	VecRestoreArray1dWrite(petsclib::PetscLibType,x::PetscVec, m::PetscInt, mstart::PetscInt, a::PetscArray{PetscScalar, 1}) 
Restores a vector after `VecGetArray1dWrite()` has been called.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of two dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `a`      - location of pointer to array obtained from `VecGetArray1d()`

Level: developer

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecRestoreArrays()`, `VecPlaceArray()`,
`VecGetArray2d()`, `VecGetArray3d()`, `VecRestoreArray3d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`
`VecGetArray1d()`, `VecRestoreArray2d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecRestoreArray1dWrite"))
"""
function VecRestoreArray1dWrite(petsclib::PetscLibType, x::PetscVec, m::PetscInt, mstart::PetscInt, a::PetscArray{PetscScalar, 1}) end

@for_petsc function VecRestoreArray1dWrite(petsclib::$UnionPetscLib, x::PetscVec, m::$PetscInt, mstart::$PetscInt, a::PetscArray{$PetscScalar, 1} )
	if a.ptr[]  != C_NULL 

    @chk ccall(
               (:VecRestoreArray1dWrite, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, Ref{Ptr{$PetscScalar}}),
               x, m, mstart, a_,
              )

		a.ptr[] = C_NULL
	else
		error("The input array is already restored")
	end


	return nothing
end 

 

"""
	a::PetscArray = VecGetArray3d(petsclib::PetscLibType,x::PetscVec, m::PetscInt, n::PetscInt, p::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt) 
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

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrays()`, `VecPlaceArray()`,
`VecRestoreArray2d()`, `DMDAVecGetarray()`, `DMDAVecRestoreArray()`, `VecRestoreArray3d()`,
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecGetArray3d"))
"""
function VecGetArray3d(petsclib::PetscLibType, x::PetscVec, m::PetscInt, n::PetscInt, p::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt) end

@for_petsc function VecGetArray3d(petsclib::$UnionPetscLib, x::PetscVec, m::$PetscInt, n::$PetscInt, p::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, pstart::$PetscInt )
	a_ = Ref{Ptr{Ptr{Ptr{$PetscScalar}}}}()

    @chk ccall(
               (:VecGetArray3d, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{Ptr{$PetscScalar}}}}),
               x, m, n, p, mstart, nstart, pstart, a_,
              )

	data_ptr = unsafe_load(a_[])
	mat = unsafe_wrap(Array, data_ptr, (m,n,p)) 
	a = PetscArray(mat,[data_ptr]) 

	return a
end 



"""
	a::PetscArray = VecGetArray3dWrite(petsclib::PetscLibType,x::PetscVec, m::PetscInt, n::PetscInt, p::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt) 
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

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrays()`, `VecPlaceArray()`,
`VecRestoreArray2d()`, `DMDAVecGetarray()`, `DMDAVecRestoreArray()`, `VecGetArray3d()`, `VecRestoreArray3d()`,
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecGetArray3dWrite"))
"""
function VecGetArray3dWrite(petsclib::PetscLibType, x::PetscVec, m::PetscInt, n::PetscInt, p::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt) end

@for_petsc function VecGetArray3dWrite(petsclib::$UnionPetscLib, x::PetscVec, m::$PetscInt, n::$PetscInt, p::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, pstart::$PetscInt )
	a_ = Ref{Ptr{Ptr{Ptr{$PetscScalar}}}}()

    @chk ccall(
               (:VecGetArray3dWrite, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{Ptr{$PetscScalar}}}}),
               x, m, n, p, mstart, nstart, pstart, a_,
              )

	data_ptr = unsafe_load(a_[])
	mat = unsafe_wrap(Array, data_ptr, (m,n,p)) 
	a = PetscArray(mat,[data_ptr]) 

	return a
end 



"""
	VecRestoreArray3d(petsclib::PetscLibType,x::PetscVec, m::PetscInt, n::PetscInt, p::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, a::PetscArray{PetscScalar, 3}) 
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

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecRestoreArrays()`, `VecPlaceArray()`,
`VecGetArray2d()`, `VecGetArray3d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecRestoreArray3d"))
"""
function VecRestoreArray3d(petsclib::PetscLibType, x::PetscVec, m::PetscInt, n::PetscInt, p::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, a::PetscArray{PetscScalar, 3}) end

@for_petsc function VecRestoreArray3d(petsclib::$UnionPetscLib, x::PetscVec, m::$PetscInt, n::$PetscInt, p::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, pstart::$PetscInt, a::PetscArray{$PetscScalar, 3} )
	if a.ptr[]  != C_NULL  

    @chk ccall(
               (:VecRestoreArray3d, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{$PetscScalar}}}),
               x, m, n, p, mstart, nstart, pstart, a,
              )

		a.ptr[] = C_NULL
	else
		error("The input array is already restored")
	end


	return nothing
end 

"""
	VecRestoreArray3dWrite(petsclib::PetscLibType,x::PetscVec, m::PetscInt, n::PetscInt, p::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, a::PetscArray{PetscScalar, 3}) 
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

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecRestoreArrays()`, `VecPlaceArray()`,
`VecGetArray2d()`, `VecGetArray3d()`, `VecRestoreArray3d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecRestoreArray3dWrite"))
"""
function VecRestoreArray3dWrite(petsclib::PetscLibType, x::PetscVec, m::PetscInt, n::PetscInt, p::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, a::PetscArray{PetscScalar, 3}) end

@for_petsc function VecRestoreArray3dWrite(petsclib::$UnionPetscLib, x::PetscVec, m::$PetscInt, n::$PetscInt, p::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, pstart::$PetscInt, a::PetscArray{$PetscScalar, 3} )
	if a.ptr[]  != C_NULL  

    @chk ccall(
               (:VecRestoreArray3dWrite, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{$PetscScalar}}}),
               x, m, n, p, mstart, nstart, pstart, a,
              )

		a.ptr[] = C_NULL
	else
		error("The input array is already restored")
	end


	return nothing
end 



"""
	a::PetscArray = VecGetArray4d(petsclib::PetscLibType,x::PetscVec, m::PetscInt, n::PetscInt, p::PetscInt, q::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, qstart::PetscInt) 
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

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrays()`, `VecPlaceArray()`,
`VecRestoreArray2d()`, `DMDAVecGetarray()`, `DMDAVecRestoreArray()`, `VecGetArray3d()`, `VecRestoreArray3d()`,
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecGetArray4d"))
"""
function VecGetArray4d(petsclib::PetscLibType, x::PetscVec, m::PetscInt, n::PetscInt, p::PetscInt, q::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, qstart::PetscInt) end

@for_petsc function VecGetArray4d(petsclib::$UnionPetscLib, x::PetscVec, m::$PetscInt, n::$PetscInt, p::$PetscInt, q::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, pstart::$PetscInt, qstart::$PetscInt )
	a_ = Ref{Ptr{Ptr{Ptr{Ptr{$PetscScalar}}}}}()

    @chk ccall(
               (:VecGetArray4d, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{Ptr{Ptr{$PetscScalar}}}}}),
               x, m, n, p, q, mstart, nstart, pstart, qstart, a_,
              )

	data_ptr = unsafe_load(a_[])
	mat = unsafe_wrap(Array, data_ptr, (m,n,p,q)) 
	a = PetscArray(mat,[data_ptr]) 

	return a
end 


"""
	a::PetscArray = VecGetArray4dWrite(petsclib::PetscLibType,x::PetscVec, m::PetscInt, n::PetscInt, p::PetscInt, q::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, qstart::PetscInt) 
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

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrays()`, `VecPlaceArray()`,
`VecRestoreArray2d()`, `DMDAVecGetarray()`, `DMDAVecRestoreArray()`, `VecGetArray3d()`, `VecRestoreArray3d()`,
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecGetArray4dWrite"))
"""
function VecGetArray4dWrite(petsclib::PetscLibType, x::PetscVec, m::PetscInt, n::PetscInt, p::PetscInt, q::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, qstart::PetscInt) end

@for_petsc function VecGetArray4dWrite(petsclib::$UnionPetscLib, x::PetscVec, m::$PetscInt, n::$PetscInt, p::$PetscInt, q::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, pstart::$PetscInt, qstart::$PetscInt )
	a_ = Ref{Ptr{Ptr{Ptr{Ptr{$PetscScalar}}}}}()

    @chk ccall(
               (:VecGetArray4dWrite, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{Ptr{Ptr{$PetscScalar}}}}}),
               x, m, n, p, q, mstart, nstart, pstart, qstart, a_,
              )

	data_ptr = unsafe_load(a_[])
	mat = unsafe_wrap(Array, data_ptr, (m,n,p,q)) 
	a = PetscArray(mat,[data_ptr]) 

	return a
end 


"""
	VecRestoreArray4d(petsclib::PetscLibType,x::PetscVec, m::PetscInt, n::PetscInt, p::PetscInt, q::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, qstart::PetscInt, a::PetscArray{PetscScalar, 4}) 
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

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecRestoreArrays()`, `VecPlaceArray()`,
`VecGetArray2d()`, `VecGetArray3d()`, `VecRestoreArray3d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`

# External Links
$(_doc_external("Vec/VecRestoreArray4d"))
"""
function VecRestoreArray4d(petsclib::PetscLibType, x::PetscVec, m::PetscInt, n::PetscInt, p::PetscInt, q::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, qstart::PetscInt, a::PetscArray{PetscScalar, 4}) end

@for_petsc function VecRestoreArray4d(petsclib::$UnionPetscLib, x::PetscVec, m::$PetscInt, n::$PetscInt, p::$PetscInt, q::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, pstart::$PetscInt, qstart::$PetscInt, a::PetscArray{$PetscScalar, 4} )
	if a.ptr[]  != C_NULL  

    @chk ccall(
               (:VecRestoreArray4d, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{$PetscScalar}}}),
               x, m, n, p, q, mstart, nstart, pstart, qstart, a,
              )

		a.ptr[] = C_NULL
	else
		error("The input array is already restored")
	end

	return nothing
end 


"""
	VecRestoreArray4dWrite(petsclib::PetscLibType,x::PetscVec, m::PetscInt, n::PetscInt, p::PetscInt, q::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, qstart::PetscInt, a::PetscArray{PetscScalar, 4}) 
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

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecRestoreArrays()`, `VecPlaceArray()`,
`VecGetArray2d()`, `VecGetArray3d()`, `VecRestoreArray3d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecRestoreArray4dWrite"))
"""
function VecRestoreArray4dWrite(petsclib::PetscLibType, x::PetscVec, m::PetscInt, n::PetscInt, p::PetscInt, q::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, qstart::PetscInt, a::PetscArray{PetscScalar, 4}) end

@for_petsc function VecRestoreArray4dWrite(petsclib::$UnionPetscLib, x::PetscVec, m::$PetscInt, n::$PetscInt, p::$PetscInt, q::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, pstart::$PetscInt, qstart::$PetscInt, a::PetscArray{$PetscScalar, 4} )
	if a.ptr[]  != C_NULL  

    @chk ccall(
               (:VecRestoreArray4dWrite, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{$PetscScalar}}}),
               x, m, n, p, q, mstart, nstart, pstart, qstart, a,
              )

		a.ptr[] = C_NULL
	else
		error("The input array is already restored")
	end


	return nothing
end 

"""
	a::Vector{PetscScalar} = VecGetArray2dRead(petsclib::PetscLibType,x::PetscVec, m::PetscInt, n::PetscInt, mstart::PetscInt, nstart::PetscInt) 
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

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrays()`, `VecPlaceArray()`,
`VecRestoreArray2d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`, `VecGetArray3d()`, `VecRestoreArray3d()`,
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecGetArray2dRead"))
"""
function VecGetArray2dRead(petsclib::PetscLibType, x::PetscVec, m::PetscInt, n::PetscInt, mstart::PetscInt, nstart::PetscInt) end

@for_petsc function VecGetArray2dRead(petsclib::$UnionPetscLib, x::PetscVec, m::$PetscInt, n::$PetscInt, mstart::$PetscInt, nstart::$PetscInt )

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

	return PetscArray(mat,[data_ptr])
end 

"""
	VecRestoreArray2dRead(petsclib::PetscLibType,x::PetscVec, m::PetscInt, n::PetscInt, mstart::PetscInt, nstart::PetscInt, a::PetscArray{PetscScalar, 2}) 
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

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecRestoreArrays()`, `VecPlaceArray()`,
`VecGetArray2d()`, `VecGetArray3d()`, `VecRestoreArray3d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecRestoreArray2dRead"))
"""
function VecRestoreArray2dRead(petsclib::PetscLibType, x::PetscVec, m::PetscInt, n::PetscInt, mstart::PetscInt, nstart::PetscInt, a::PetscArray{PetscScalar, 2}) end

@for_petsc function VecRestoreArray2dRead(petsclib::$UnionPetscLib, x::PetscVec, m::$PetscInt, n::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, a::PetscArray{$PetscScalar, 2} )
	if a.ptr[]  != C_NULL  

    @chk ccall(
               (:VecRestoreArray2dRead, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{$PetscScalar}}}),
               x, m, n, mstart, nstart, a,
              )

		a.ptr[] = C_NULL
	else
		error("The input array is already restored")
	end


	return nothing
end 


"""
	VecRestoreArray1dRead(petsclib::PetscLibType,x::PetscVec, m::PetscInt, mstart::PetscInt, a::PetscArray{PetscScalar, 1}) 
Restores a vector after `VecGetArray1dRead()` has been called.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of two dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `a`      - location of pointer to array obtained from `VecGetArray1dRead()`

Level: developer

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecRestoreArrays()`, `VecPlaceArray()`,
`VecGetArray2d()`, `VecGetArray3d()`, `VecRestoreArray3d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`
`VecGetArray1d()`, `VecRestoreArray2d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecRestoreArray1dRead"))
"""
function VecRestoreArray1dRead(petsclib::PetscLibType, x::PetscVec, m::PetscInt, mstart::PetscInt, a::PetscArray{PetscScalar, 1}) end

@for_petsc function VecRestoreArray1dRead(petsclib::$UnionPetscLib, x::PetscVec, m::$PetscInt, mstart::$PetscInt, a::PetscArray{$PetscScalar, 1} )
	if a.ptr[]  != C_NULL 
    @chk ccall(
               (:VecRestoreArray1dRead, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, Ref{Ptr{$PetscScalar}}),
               x, m, mstart, a_,
              )

		a.ptr[] = C_NULL
	else
		error("The input array is already restored")
	end


	return nothing
end 

 



"""
	VecRestoreArray1dRead(petsclib::PetscLibType,x::PetscVec, m::PetscInt, mstart::PetscInt, a::Vector{PetscScalar}) 
Restores a vector after `VecGetArray1dRead()` has been called.

Logically Collective

Input Parameters:
- `x`      - the vector
- `m`      - first dimension of two dimensional array
- `mstart` - first index you will use in first coordinate direction (often 0)
- `a`      - location of pointer to array obtained from `VecGetArray1dRead()`

Level: developer

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecRestoreArrays()`, `VecPlaceArray()`,
`VecGetArray2d()`, `VecGetArray3d()`, `VecRestoreArray3d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`
`VecGetArray1d()`, `VecRestoreArray2d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecRestoreArray1dRead"))
"""
function VecRestoreArray1dRead(petsclib::PetscLibType, x::PetscVec, m::PetscInt, mstart::PetscInt, a::Vector{PetscScalar}) end

@for_petsc function VecRestoreArray1dRead(petsclib::$UnionPetscLib, x::PetscVec, m::$PetscInt, mstart::$PetscInt, a::Vector{$PetscScalar} )
	a_ = Ref(pointer(a))

    @chk ccall(
               (:VecRestoreArray1dRead, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, Ptr{Ptr{$PetscScalar}}),
               x, m, mstart, a_,
              )


	return nothing
end 

"""
	a::PetscArray = VecGetArray3dRead(petsclib::PetscLibType,x::PetscVec, m::PetscInt, n::PetscInt, p::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt) 
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

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrays()`, `VecPlaceArray()`,
`VecRestoreArray2d()`, `DMDAVecGetarray()`, `DMDAVecRestoreArray()`, `VecGetArray3d()`, `VecRestoreArray3d()`,
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecGetArray3dRead"))
"""
function VecGetArray3dRead(petsclib::PetscLibType, x::PetscVec, m::PetscInt, n::PetscInt, p::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt) end

@for_petsc function VecGetArray3dRead(petsclib::$UnionPetscLib, x::PetscVec, m::$PetscInt, n::$PetscInt, p::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, pstart::$PetscInt )
	a_ = Ref{Ptr{Ptr{Ptr{$PetscScalar}}}}()

    @chk ccall(
               (:VecGetArray3dRead, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{Ptr{$PetscScalar}}}}),
               x, m, n, p, mstart, nstart, pstart, a_,
              )

	data_ptr = unsafe_load(a_[])
	mat = unsafe_wrap(Array, data_ptr, (m,n,p)) 
	a = PetscArray(mat,[data_ptr]) 

	return a
end 

 

"""
	VecRestoreArray3dRead(petsclib::PetscLibType,x::PetscVec, m::PetscInt, n::PetscInt, p::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, a::PetscArray{PetscScalar, 3}) 
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

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecRestoreArrays()`, `VecPlaceArray()`,
`VecGetArray2d()`, `VecGetArray3d()`, `VecRestoreArray3d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecRestoreArray3dRead"))
"""
function VecRestoreArray3dRead(petsclib::PetscLibType, x::PetscVec, m::PetscInt, n::PetscInt, p::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, a::PetscArray{PetscScalar, 3}) end

@for_petsc function VecRestoreArray3dRead(petsclib::$UnionPetscLib, x::PetscVec, m::$PetscInt, n::$PetscInt, p::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, pstart::$PetscInt, a::PetscArray{$PetscScalar, 3} )
	if a.ptr[]  != C_NULL 

    @chk ccall(
               (:VecRestoreArray3dRead, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{$PetscScalar}}}),
               x, m, n, p, mstart, nstart, pstart, a,
              )

		a.ptr[] = C_NULL
	else
		error("The input array is already restored")
	end


	return nothing
end 


"""
	a::PetscArray = VecGetArray4dRead(petsclib::PetscLibType,x::PetscVec, m::PetscInt, n::PetscInt, p::PetscInt, q::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, qstart::PetscInt) 
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

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecGetArrays()`, `VecPlaceArray()`,
`VecRestoreArray2d()`, `DMDAVecGetarray()`, `DMDAVecRestoreArray()`, `VecGetArray3d()`, `VecRestoreArray3d()`,
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecGetArray4dRead"))
"""
function VecGetArray4dRead(petsclib::PetscLibType, x::PetscVec, m::PetscInt, n::PetscInt, p::PetscInt, q::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, qstart::PetscInt) end

@for_petsc function VecGetArray4dRead(petsclib::$UnionPetscLib, x::PetscVec, m::$PetscInt, n::$PetscInt, p::$PetscInt, q::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, pstart::$PetscInt, qstart::$PetscInt )
	a_ = Ref{Ptr{Ptr{Ptr{Ptr{$PetscScalar}}}}}()

    @chk ccall(
               (:VecGetArray4dRead, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{Ptr{Ptr{$PetscScalar}}}}}),
               x, m, n, p, q, mstart, nstart, pstart, qstart, a_,
              )

	data_ptr = unsafe_load(a_[])
	mat = unsafe_wrap(Array, data_ptr, (m,n,p,q)) 
	a = PetscArray(mat,[data_ptr]) 

	return a
end 


"""
	VecRestoreArray4dRead(petsclib::PetscLibType,x::PetscVec, m::PetscInt, n::PetscInt, p::PetscInt, q::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, qstart::PetscInt, a::PetscArray{PetscScalar, 4}) 
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

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecRestoreArrays()`, `VecPlaceArray()`,
`VecGetArray2d()`, `VecGetArray3d()`, `VecRestoreArray3d()`, `DMDAVecGetArray()`, `DMDAVecRestoreArray()`
`VecGetArray1d()`, `VecRestoreArray1d()`, `VecGetArray4d()`, `VecRestoreArray4d()`

# External Links
$(_doc_external("Vec/VecRestoreArray4dRead"))
"""
function VecRestoreArray4dRead(petsclib::PetscLibType, x::PetscVec, m::PetscInt, n::PetscInt, p::PetscInt, q::PetscInt, mstart::PetscInt, nstart::PetscInt, pstart::PetscInt, qstart::PetscInt, a::PetscArray{PetscScalar, 4}) end

@for_petsc function VecRestoreArray4dRead(petsclib::$UnionPetscLib, x::PetscVec, m::$PetscInt, n::$PetscInt, p::$PetscInt, q::$PetscInt, mstart::$PetscInt, nstart::$PetscInt, pstart::$PetscInt, qstart::$PetscInt, a::PetscArray{$PetscScalar, 4} )
	if a.ptr[]  != C_NULL  

    @chk ccall(
               (:VecRestoreArray4dRead, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ref{Ptr{Ptr{$PetscScalar}}}),
               x, m, n, p, q, mstart, nstart, pstart, qstart, a,
              )

		a.ptr[] = C_NULL
	else
		error("The input array is already restored")
	end


	return nothing
end 



"""
	state::PetscInt = VecLockGet(petsclib::PetscLibType,x::PetscVec) 
Get the current lock status of a vector

Logically Collective

Input Parameter:
- `x` - the vector

Output Parameter:
- `state` - greater than zero indicates the vector is locked for read; less than zero indicates the vector is
locked for write; equal to zero means the vector is unlocked, that is, it is free to read or write.

Level: advanced

-seealso: [](ch_vectors), `Vec`, `VecRestoreArray()`, `VecGetArrayRead()`, `VecLockReadPush()`, `VecLockReadPop()`

# External Links
$(_doc_external("Vec/VecLockGet"))
"""
function VecLockGet(petsclib::PetscLibType, x::PetscVec) end

@for_petsc function VecLockGet(petsclib::$UnionPetscLib, x::PetscVec )
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
	VecLockGetLocation(petsclib::PetscLibType,x::PetscVec, file::Vector{Cchar}, func::Vector{Cchar}, line::Cint) 

# External Links
$(_doc_external("Vec/VecLockGetLocation"))
"""
function VecLockGetLocation(petsclib::PetscLibType, x::PetscVec, file::Vector{Cchar}, func::Vector{Cchar}, line::Cint) end

@for_petsc function VecLockGetLocation(petsclib::$UnionPetscLib, x::PetscVec, file::Vector{Cchar}, func::Vector{Cchar}, line::Cint )
	file_ = Ref(pointer(file))
	func_ = Ref(pointer(func))

    @chk ccall(
               (:VecLockGetLocation, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Ptr{Cchar}}, Ptr{Ptr{Cchar}}, Ptr{Cint}),
               x, file_, func_, line,
              )


	return nothing
end 

"""
	VecLockReadPush(petsclib::PetscLibType,x::PetscVec) 
Push a read

Logically Collective

Input Parameter:
- `x` - the vector

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecRestoreArray()`, `VecGetArrayRead()`, `VecLockReadPop()`, `VecLockGet()`

# External Links
$(_doc_external("Vec/VecLockReadPush"))
"""
function VecLockReadPush(petsclib::PetscLibType, x::PetscVec) end

@for_petsc function VecLockReadPush(petsclib::$UnionPetscLib, x::PetscVec )

    @chk ccall(
               (:VecLockReadPush, $petsc_library),
               PetscErrorCode,
               (CVec,),
               x,
              )


	return nothing
end 

"""
	VecLockReadPop(petsclib::PetscLibType,x::PetscVec) 
Pop a read

Logically Collective

Input Parameter:
- `x` - the vector

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecRestoreArray()`, `VecGetArrayRead()`, `VecLockReadPush()`, `VecLockGet()`

# External Links
$(_doc_external("Vec/VecLockReadPop"))
"""
function VecLockReadPop(petsclib::PetscLibType, x::PetscVec) end

@for_petsc function VecLockReadPop(petsclib::$UnionPetscLib, x::PetscVec )

    @chk ccall(
               (:VecLockReadPop, $petsc_library),
               PetscErrorCode,
               (CVec,),
               x,
              )


	return nothing
end 

"""
	VecLockWriteSet(petsclib::PetscLibType,x::PetscVec, flg::PetscBool) 
Lock or unlock a vector for exclusive read/write access

Logically Collective

Input Parameters:
- `x`   - the vector
- `flg` - `PETSC_TRUE` to lock the vector for exclusive read/write access; `PETSC_FALSE` to unlock it.

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecRestoreArray()`, `VecGetArrayRead()`, `VecLockReadPush()`, `VecLockReadPop()`, `VecLockGet()`

# External Links
$(_doc_external("Vec/VecLockWriteSet"))
"""
function VecLockWriteSet(petsclib::PetscLibType, x::PetscVec, flg::PetscBool) end

@for_petsc function VecLockWriteSet(petsclib::$UnionPetscLib, x::PetscVec, flg::PetscBool )

    @chk ccall(
               (:VecLockWriteSet, $petsc_library),
               PetscErrorCode,
               (CVec, PetscBool),
               x, flg,
              )


	return nothing
end 

"""
	VecGetLocalToGlobalMapping(petsclib::PetscLibType,X::PetscVec, mapping::ISLocalToGlobalMapping) 
Gets the local

Not Collective

Input Parameter:
- `X` - the vector

Output Parameter:
- `mapping` - the mapping

Level: advanced

-seealso: [](ch_vectors), `Vec`, `VecSetValuesLocal()`, `VecSetLocalToGlobalMapping()`

# External Links
$(_doc_external("Vec/VecGetLocalToGlobalMapping"))
"""
function VecGetLocalToGlobalMapping(petsclib::PetscLibType, X::PetscVec, mapping::ISLocalToGlobalMapping) end

@for_petsc function VecGetLocalToGlobalMapping(petsclib::$UnionPetscLib, X::PetscVec, mapping::ISLocalToGlobalMapping )

    @chk ccall(
               (:VecGetLocalToGlobalMapping, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{ISLocalToGlobalMapping}),
               X, mapping,
              )


	return nothing
end 

"""
	VecAssemblyBegin(petsclib::PetscLibType,vec::PetscVec) 
Begins assembling the vector; that is ensuring all the vector's entries are stored on the correct MPI process. This routine should
be called after completing all calls to `VecSetValues()`.

Collective

Input Parameter:
- `vec` - the vector

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecAssemblyEnd()`, `VecSetValues()`

# External Links
$(_doc_external("Vec/VecAssemblyBegin"))
"""
function VecAssemblyBegin(petsclib::PetscLibType, vec::PetscVec) end

@for_petsc function VecAssemblyBegin(petsclib::$UnionPetscLib, vec::PetscVec )

    @chk ccall(
               (:VecAssemblyBegin, $petsc_library),
               PetscErrorCode,
               (CVec,),
               vec,
              )


	return nothing
end 

"""
	VecAssemblyEnd(petsclib::PetscLibType,vec::PetscVec) 
Completes assembling the vector.  This routine should be called after `VecAssemblyBegin()`.

Collective

Input Parameter:
- `vec` - the vector

Options Database Keys:
- `-vec_view`                 - Prints vector in `PETSC_VIEWER_DEFAULT` format
- `-vec_view ::ascii_matlab`  - Prints vector in `PETSC_VIEWER_ASCII_MATLAB` format to stdout
- `-vec_view matlab:filename` - Prints vector in MATLAB .mat file to filename (requires PETSc configured with --with-matlab)
- `-vec_view draw`            - Activates vector viewing using drawing tools
- `-display <name>`           - Sets display name (default is host)
- `-draw_pause <sec>`         - Sets number of seconds to pause after display
- `-vec_view socket`          - Activates vector viewing using a socket

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecAssemblyBegin()`, `VecSetValues()`

# External Links
$(_doc_external("Vec/VecAssemblyEnd"))
"""
function VecAssemblyEnd(petsclib::PetscLibType, vec::PetscVec) end

@for_petsc function VecAssemblyEnd(petsclib::$UnionPetscLib, vec::PetscVec )

    @chk ccall(
               (:VecAssemblyEnd, $petsc_library),
               PetscErrorCode,
               (CVec,),
               vec,
              )


	return nothing
end 

"""
	VecPointwiseMax(petsclib::PetscLibType,w::PetscVec, x::PetscVec, y::PetscVec) 
Computes the component

Logically Collective

Input Parameters:
- `x` - the first input vector
- `y` - the second input vector

Output Parameter:
- `w` - the result

Level: advanced

-seealso: [](ch_vectors), `Vec`, `VecPointwiseDivide()`, `VecPointwiseMult()`, `VecPointwiseMin()`, `VecPointwiseMaxAbs()`, `VecMaxPointwiseDivide()`

# External Links
$(_doc_external("Vec/VecPointwiseMax"))
"""
function VecPointwiseMax(petsclib::PetscLibType, w::PetscVec, x::PetscVec, y::PetscVec) end

@for_petsc function VecPointwiseMax(petsclib::$UnionPetscLib, w::PetscVec, x::PetscVec, y::PetscVec )

    @chk ccall(
               (:VecPointwiseMax, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, CVec),
               w, x, y,
              )


	return nothing
end 

"""
	VecPointwiseMin(petsclib::PetscLibType,w::PetscVec, x::PetscVec, y::PetscVec) 
Computes the component

Logically Collective

Input Parameters:
- `x` - the first input vector
- `y` - the second input vector

Output Parameter:
- `w` - the result

Level: advanced

-seealso: [](ch_vectors), `Vec`, `VecPointwiseDivide()`, `VecPointwiseMult()`, `VecPointwiseMaxAbs()`, `VecMaxPointwiseDivide()`

# External Links
$(_doc_external("Vec/VecPointwiseMin"))
"""
function VecPointwiseMin(petsclib::PetscLibType, w::PetscVec, x::PetscVec, y::PetscVec) end

@for_petsc function VecPointwiseMin(petsclib::$UnionPetscLib, w::PetscVec, x::PetscVec, y::PetscVec )

    @chk ccall(
               (:VecPointwiseMin, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, CVec),
               w, x, y,
              )


	return nothing
end 

"""
	VecPointwiseMaxAbs(petsclib::PetscLibType,w::PetscVec, x::PetscVec, y::PetscVec) 
Computes the component

Logically Collective

Input Parameters:
- `x` - the first input vector
- `y` - the second input vector

Output Parameter:
- `w` - the result

Level: advanced

-seealso: [](ch_vectors), `Vec`, `VecPointwiseDivide()`, `VecPointwiseMult()`, `VecPointwiseMin()`, `VecPointwiseMax()`, `VecMaxPointwiseDivide()`

# External Links
$(_doc_external("Vec/VecPointwiseMaxAbs"))
"""
function VecPointwiseMaxAbs(petsclib::PetscLibType, w::PetscVec, x::PetscVec, y::PetscVec) end

@for_petsc function VecPointwiseMaxAbs(petsclib::$UnionPetscLib, w::PetscVec, x::PetscVec, y::PetscVec )

    @chk ccall(
               (:VecPointwiseMaxAbs, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, CVec),
               w, x, y,
              )


	return nothing
end 

"""
	VecPointwiseDivide(petsclib::PetscLibType,w::PetscVec, x::PetscVec, y::PetscVec) 
Computes the component

Logically Collective

Input Parameters:
- `x` - the numerator vector
- `y` - the denominator vector

Output Parameter:
- `w` - the result

Level: advanced

-seealso: [](ch_vectors), `Vec`, `VecPointwiseMult()`, `VecPointwiseMax()`, `VecPointwiseMin()`, `VecPointwiseMaxAbs()`, `VecMaxPointwiseDivide()`

# External Links
$(_doc_external("Vec/VecPointwiseDivide"))
"""
function VecPointwiseDivide(petsclib::PetscLibType, w::PetscVec, x::PetscVec, y::PetscVec) end

@for_petsc function VecPointwiseDivide(petsclib::$UnionPetscLib, w::PetscVec, x::PetscVec, y::PetscVec )

    @chk ccall(
               (:VecPointwiseDivide, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, CVec),
               w, x, y,
              )


	return nothing
end 

"""
	VecPointwiseMult(petsclib::PetscLibType,w::PetscVec, x::PetscVec, y::PetscVec) 
Computes the component

Logically Collective

Input Parameters:
- `x` - the first vector
- `y` - the second vector

Output Parameter:
- `w` - the result

Level: advanced

-seealso: [](ch_vectors), `Vec`, `VecPointwiseDivide()`, `VecPointwiseMax()`, `VecPointwiseMin()`, `VecPointwiseMaxAbs()`, `VecMaxPointwiseDivide()`

# External Links
$(_doc_external("Vec/VecPointwiseMult"))
"""
function VecPointwiseMult(petsclib::PetscLibType, w::PetscVec, x::PetscVec, y::PetscVec) end

@for_petsc function VecPointwiseMult(petsclib::$UnionPetscLib, w::PetscVec, x::PetscVec, y::PetscVec )

    @chk ccall(
               (:VecPointwiseMult, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, CVec),
               w, x, y,
              )


	return nothing
end 

"""
	newv::PetscVec = VecDuplicate(petsclib::PetscLibType,v::PetscVec) 
Creates a new vector of the same type as an existing vector.

Collective

Input Parameter:
- `v` - a vector to mimic

Output Parameter:
- `newv` - location to put new vector

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecDestroy()`, `VecDuplicateVecs()`, `VecCreate()`, `VecCopy()`

# External Links
$(_doc_external("Vec/VecDuplicate"))
"""
function VecDuplicate(petsclib::PetscLibType, v::PetscVec) end

@for_petsc function VecDuplicate(petsclib::$UnionPetscLib, v::PetscVec )
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
	VecDestroy(petsclib::PetscLibType,v::PetscVec) 
Destroys a vector.

Collective

Input Parameter:
- `v` - the vector

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecCreate()`, `VecDuplicate()`, `VecDestroyVecs()`

# External Links
$(_doc_external("Vec/VecDestroy"))
"""
function VecDestroy(petsclib::PetscLibType, v::PetscVec) end

@for_petsc function VecDestroy(petsclib::$UnionPetscLib, v::PetscVec )
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
	V::Vector{PetscVec} = VecDuplicateVecs(petsclib::PetscLibType,v::PetscVec, m::PetscInt) 
Creates several vectors of the same type as an existing vector.

Collective

Input Parameters:
- `m` - the number of vectors to obtain
- `v` - a vector to mimic

Output Parameter:
- `V` - location to put pointer to array of vectors

Level: intermediate

-seealso: [](ch_vectors), `Vec`, [](ch_fortran), `VecDestroyVecs()`, `VecDuplicate()`, `VecCreate()`, `VecMDot()`, `VecMAXPY()`, `KSPGMRES`,
`KSPGMRESSetPreAllocateVectors()`

# External Links
$(_doc_external("Vec/VecDuplicateVecs"))
"""
function VecDuplicateVecs(petsclib::PetscLibType, v::PetscVec, m::PetscInt) end

@for_petsc function VecDuplicateVecs(petsclib::$UnionPetscLib, v::PetscVec, m::$PetscInt )
	V_ = Ref{Ptr{PetscVec}}()

    @chk ccall(
               (:VecDuplicateVecs, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{Ptr{CVec}}),
               v, m, V_,
              )

	V = unsafe_wrap(Array, V_[], VecGetLocalSize(petsclib, x); own = false)

	return V
end 

"""
	VecDestroyVecs(petsclib::PetscLibType,m::PetscInt, vv::Vector{PetscVec}) 
Frees a block of vectors obtained with `VecDuplicateVecs()`.

Collective

Input Parameters:
- `m`  - the number of vectors previously obtained, if zero no vectors are destroyed
- `vv` - pointer to pointer to array of vector pointers, if `NULL` no vectors are destroyed

Level: intermediate

-seealso: [](ch_vectors), `Vec`, [](ch_fortran), `VecDuplicateVecs()`, `VecDestroyVecsf90()`

# External Links
$(_doc_external("Vec/VecDestroyVecs"))
"""
function VecDestroyVecs(petsclib::PetscLibType, m::PetscInt, vv::Vector{PetscVec}) end

@for_petsc function VecDestroyVecs(petsclib::$UnionPetscLib, m::$PetscInt, vv::Vector{PetscVec} )
	vv_ = Ref(pointer(vv))

    @chk ccall(
               (:VecDestroyVecs, $petsc_library),
               PetscErrorCode,
               ($PetscInt, Ptr{Ptr{CVec}}),
               m, vv_,
              )


	return nothing
end 

"""
	VecViewFromOptions(petsclib::PetscLibType,A::PetscVec, obj::PetscObject, name::Vector{Cchar}) 
View a vector based on values in the options database

Collective

Input Parameters:
- `A`    - the vector
- `obj`  - optional object that provides the options prefix for this viewing, use 'NULL' to use the prefix of `A`
- `name` - command line option

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecView`, `PetscObjectViewFromOptions()`, `VecCreate()`

# External Links
$(_doc_external("Vec/VecViewFromOptions"))
"""
function VecViewFromOptions(petsclib::PetscLibType, A::PetscVec, obj::PetscObject, name::Vector{Cchar}) end

@for_petsc function VecViewFromOptions(petsclib::$UnionPetscLib, A::PetscVec, obj::PetscObject, name::Vector{Cchar} )

    @chk ccall(
               (:VecViewFromOptions, $petsc_library),
               PetscErrorCode,
               (CVec, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

"""
	VecView(petsclib::PetscLibType,vec::PetscVec, viewer::PetscViewer) 
Views a vector object.

Collective

Input Parameters:
- `vec`    - the vector
- `viewer` - an optional `PetscViewer` visualization context

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecViewFromOptions()`, `PetscViewerASCIIOpen()`, `PetscViewerDrawOpen()`, `PetscDrawLGCreate()`,
`PetscViewerSocketOpen()`, `PetscViewerBinaryOpen()`, `VecLoad()`, `PetscViewerCreate()`,
`PetscRealView()`, `PetscScalarView()`, `PetscIntView()`, `PetscViewerHDF5SetTimestep()`

# External Links
$(_doc_external("Vec/VecView"))
"""
function VecView(petsclib::PetscLibType, vec::PetscVec, viewer::PetscViewer) end

@for_petsc function VecView(petsclib::$UnionPetscLib, vec::PetscVec, viewer::PetscViewer )

    @chk ccall(
               (:VecView, $petsc_library),
               PetscErrorCode,
               (CVec, PetscViewer),
               vec, viewer,
              )


	return nothing
end 

"""
	VecViewNative(petsclib::PetscLibType,vec::PetscVec, viewer::PetscViewer) 
Views a vector object with the original type specific viewer

Collective

Input Parameters:
- `vec`    - the vector
- `viewer` - an optional `PetscViewer` visualization context

Level: developer

-seealso: [](ch_vectors), `Vec`, `PetscViewerASCIIOpen()`, `PetscViewerDrawOpen()`, `PetscDrawLGCreate()`, `VecView()`
`PetscViewerSocketOpen()`, `PetscViewerBinaryOpen()`, `VecLoad()`, `PetscViewerCreate()`,
`PetscRealView()`, `PetscScalarView()`, `PetscIntView()`, `PetscViewerHDF5SetTimestep()`

# External Links
$(_doc_external("Vec/VecViewNative"))
"""
function VecViewNative(petsclib::PetscLibType, vec::PetscVec, viewer::PetscViewer) end

@for_petsc function VecViewNative(petsclib::$UnionPetscLib, vec::PetscVec, viewer::PetscViewer )

    @chk ccall(
               (:VecViewNative, $petsc_library),
               PetscErrorCode,
               (CVec, PetscViewer),
               vec, viewer,
              )


	return nothing
end 

"""
	size::PetscInt = VecGetSize(petsclib::PetscLibType,x::PetscVec) 
Returns the global number of elements of the vector.

Not Collective

Input Parameter:
- `x` - the vector

Output Parameter:
- `size` - the global length of the vector

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecGetLocalSize()`

# External Links
$(_doc_external("Vec/VecGetSize"))
"""
function VecGetSize(petsclib::PetscLibType, x::AbstractPetscVec) end

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
	size::PetscInt = VecGetLocalSize(petsclib::PetscLibType,x::PetscVec) 
Returns the number of elements of the vector stored
in local memory (that is on this MPI process)

Not Collective

Input Parameter:
- `x` - the vector

Output Parameter:
- `size` - the length of the local piece of the vector

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecGetSize()`

# External Links
$(_doc_external("Vec/VecGetLocalSize"))
"""
function VecGetLocalSize(petsclib::PetscLibType, x::AbstractPetscVec) end

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
	low::PetscInt,high::PetscInt = VecGetOwnershipRange(petsclib::PetscLibType,x::PetscVec) 
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

-seealso: [](ch_vectors), `Vec`, `MatGetOwnershipRange()`, `MatGetOwnershipRanges()`, `VecGetOwnershipRanges()`, `PetscSplitOwnership()`,
`VecSetSizes()`, `VecCreateMPI()`, `PetscLayout`, `DMDAGetGhostCorners()`, `DM`

# External Links
$(_doc_external("Vec/VecGetOwnershipRange"))
"""
function VecGetOwnershipRange(petsclib::PetscLibType, x::PetscVec) end

@for_petsc function VecGetOwnershipRange(petsclib::$UnionPetscLib, x::PetscVec )
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
	ranges::Vector{PetscInt} = VecGetOwnershipRanges(petsclib::PetscLibType,x::PetscVec) 
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

-seealso: [](ch_vectors), `Vec`, `MatGetOwnershipRange()`, `MatGetOwnershipRanges()`, `VecGetOwnershipRange()`, `PetscSplitOwnership()`,
`VecSetSizes()`, `VecCreateMPI()`, `PetscLayout`, `DMDAGetGhostCorners()`, `DM`

# External Links
$(_doc_external("Vec/VecGetOwnershipRanges"))
"""
function VecGetOwnershipRanges(petsclib::PetscLibType, x::PetscVec) end

@for_petsc function VecGetOwnershipRanges(petsclib::$UnionPetscLib, x::PetscVec )
	ranges_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:VecGetOwnershipRanges, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Ptr{$PetscInt}}),
               x, ranges_,
              )

	ranges = unsafe_wrap(Array, ranges_[], VecGetLocalSize(petsclib, x); own = false)

	return ranges
end 

"""
	VecResetArray(petsclib::PetscLibType,vec::PetscVec) 
Resets a vector to use its default memory. Call this
after the use of `VecPlaceArray()`.

Not Collective

Input Parameter:
- `vec` - the vector

Level: developer

-seealso: [](ch_vectors), `Vec`, `VecGetArray()`, `VecRestoreArray()`, `VecReplaceArray()`, `VecPlaceArray()`

# External Links
$(_doc_external("Vec/VecResetArray"))
"""
function VecResetArray(petsclib::PetscLibType, vec::PetscVec) end

@for_petsc function VecResetArray(petsclib::$UnionPetscLib, vec::PetscVec )

    @chk ccall(
               (:VecResetArray, $petsc_library),
               PetscErrorCode,
               (CVec,),
               vec,
              )


	return nothing
end 

"""
	VecLoad(petsclib::PetscLibType,vec::PetscVec, viewer::PetscViewer) 
Loads a vector that has been stored in binary or HDF5 format
with `VecView()`.

Collective

Input Parameters:
- `vec`    - the newly loaded vector, this needs to have been created with `VecCreate()` or
some related function before the call to `VecLoad()`.
- `viewer` - binary file viewer, obtained from `PetscViewerBinaryOpen()` or
HDF5 file viewer, obtained from `PetscViewerHDF5Open()`

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `PetscViewerBinaryOpen()`, `VecView()`, `MatLoad()`

# External Links
$(_doc_external("Vec/VecLoad"))
"""
function VecLoad(petsclib::PetscLibType, vec::PetscVec, viewer::PetscViewer) end

@for_petsc function VecLoad(petsclib::$UnionPetscLib, vec::PetscVec, viewer::PetscViewer )

    @chk ccall(
               (:VecLoad, $petsc_library),
               PetscErrorCode,
               (CVec, PetscViewer),
               vec, viewer,
              )


	return nothing
end 

"""
	VecReciprocal(petsclib::PetscLibType,vec::PetscVec) 
Replaces each component of a vector by its reciprocal.

Logically Collective

Input Parameter:
- `vec` - the vector

Output Parameter:
- `vec` - the vector reciprocal

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecLog()`, `VecExp()`, `VecSqrtAbs()`

# External Links
$(_doc_external("Vec/VecReciprocal"))
"""
function VecReciprocal(petsclib::PetscLibType, vec::PetscVec) end

@for_petsc function VecReciprocal(petsclib::$UnionPetscLib, vec::PetscVec )

    @chk ccall(
               (:VecReciprocal, $petsc_library),
               PetscErrorCode,
               (CVec,),
               vec,
              )


	return nothing
end 

"""
	VecZeroEntries(petsclib::PetscLibType,vec::PetscVec) 
puts a `0.0` in each element of a vector

Logically Collective

Input Parameter:
- `vec` - The vector

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecCreate()`, `VecSetOptionsPrefix()`, `VecSet()`, `VecSetValues()`

# External Links
$(_doc_external("Vec/VecZeroEntries"))
"""
function VecZeroEntries(petsclib::PetscLibType, vec::PetscVec) end

@for_petsc function VecZeroEntries(petsclib::$UnionPetscLib, vec::PetscVec )

    @chk ccall(
               (:VecZeroEntries, $petsc_library),
               PetscErrorCode,
               (CVec,),
               vec,
              )


	return nothing
end 

"""
	bs::PetscInt = VecGetBlockSize(petsclib::PetscLibType,v::PetscVec) 
Gets the blocksize for the vector, i.e. what is used for `VecSetValuesBlocked()`
and `VecSetValuesBlockedLocal()`.

Not Collective

Input Parameter:
- `v` - the vector

Output Parameter:
- `bs` - the blocksize

Level: advanced

-seealso: [](ch_vectors), `Vec`, `VecSetValuesBlocked()`, `VecSetLocalToGlobalMapping()`, `VecSetBlockSize()`

# External Links
$(_doc_external("Vec/VecGetBlockSize"))
"""
function VecGetBlockSize(petsclib::PetscLibType, v::PetscVec) end

@for_petsc function VecGetBlockSize(petsclib::$UnionPetscLib, v::PetscVec )
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
	VecAppendOptionsPrefix(petsclib::PetscLibType,v::PetscVec, prefix::Vector{Cchar}) 
Appends to the prefix used for searching for all
`Vec` options in the database.

Logically Collective

Input Parameters:
- `v`      - the `Vec` context
- `prefix` - the prefix to prepend to all option names

Level: advanced

-seealso: [](ch_vectors), `Vec`, `VecGetOptionsPrefix()`

# External Links
$(_doc_external("Vec/VecAppendOptionsPrefix"))
"""
function VecAppendOptionsPrefix(petsclib::PetscLibType, v::PetscVec, prefix::Vector{Cchar}) end

@for_petsc function VecAppendOptionsPrefix(petsclib::$UnionPetscLib, v::PetscVec, prefix::Vector{Cchar} )

    @chk ccall(
               (:VecAppendOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Cchar}),
               v, prefix,
              )


	return nothing
end 

"""
	VecGetOptionsPrefix(petsclib::PetscLibType,v::PetscVec, prefix::Vector{Cchar}) 
Sets the prefix used for searching for all
Vec options in the database.

Not Collective

Input Parameter:
- `v` - the `Vec` context

Output Parameter:
- `prefix` - pointer to the prefix string used

Level: advanced

-seealso: [](ch_vectors), `Vec`, `VecAppendOptionsPrefix()`

# External Links
$(_doc_external("Vec/VecGetOptionsPrefix"))
"""
function VecGetOptionsPrefix(petsclib::PetscLibType, v::PetscVec, prefix::Vector{Cchar}) end

@for_petsc function VecGetOptionsPrefix(petsclib::$UnionPetscLib, v::PetscVec, prefix::Vector{Cchar} )
	prefix_ = Ref(pointer(prefix))

    @chk ccall(
               (:VecGetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Ptr{Cchar}}),
               v, prefix_,
              )


	return nothing
end 

"""
	VecGetState(petsclib::PetscLibType,v::PetscVec, state::PetscObjectState) 
Gets the state of a `Vec`.

Not Collective

Input Parameter:
- `v` - the `Vec` context

Output Parameter:
- `state` - the object state

Level: advanced

-seealso: [](ch_vectors), `Vec`, `VecCreate()`, `PetscObjectStateGet()`

# External Links
$(_doc_external("Vec/VecGetState"))
"""
function VecGetState(petsclib::PetscLibType, v::PetscVec, state::PetscObjectState) end

@for_petsc function VecGetState(petsclib::$UnionPetscLib, v::PetscVec, state::PetscObjectState )

    @chk ccall(
               (:VecGetState, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{PetscObjectState}),
               v, state,
              )


	return nothing
end 

"""
	VecCopy(petsclib::PetscLibType,x::PetscVec, y::PetscVec) 
Copies a vector `y = x`

Logically Collective

Input Parameter:
- `x` - the vector

Output Parameter:
- `y` - the copy

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecDuplicate()`

# External Links
$(_doc_external("Vec/VecCopy"))
"""
function VecCopy(petsclib::PetscLibType, x::PetscVec, y::PetscVec) end

@for_petsc function VecCopy(petsclib::$UnionPetscLib, x::PetscVec, y::PetscVec )

    @chk ccall(
               (:VecCopy, $petsc_library),
               PetscErrorCode,
               (CVec, CVec),
               x, y,
              )


	return nothing
end 

"""
	VecGetLayout(petsclib::PetscLibType,x::PetscVec, map::PetscLayout) 
get `PetscLayout` describing a vector layout

Not Collective

Input Parameter:
- `x` - the vector

Output Parameter:
- `map` - the layout

Level: developer

-seealso: [](ch_vectors), `PetscLayout`, `Vec`, `VecGetSize()`, `VecGetOwnershipRange()`, `VecGetOwnershipRanges()`

# External Links
$(_doc_external("Vec/VecGetLayout"))
"""
function VecGetLayout(petsclib::PetscLibType, x::PetscVec, map::PetscLayout) end

@for_petsc function VecGetLayout(petsclib::$UnionPetscLib, x::PetscVec, map::PetscLayout )

    @chk ccall(
               (:VecGetLayout, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{PetscLayout}),
               x, map,
              )


	return nothing
end 

"""
	VecFlag(petsclib::PetscLibType,xin::PetscVec, flg::PetscInt) 
set infinity into the local part of the vector on any subset of MPI processes

Logically Collective

Input Parameters:
- `xin` - the vector, can be `NULL` but only if on all processes
- `flg` - indicates if this processes portion of the vector should be set to infinity

Level: developer

-seealso: [](ch_vectors), `Vec`, `PetscLayout`, `VecGetLayout()`, `VecGetSize()`, `VecGetOwnershipRange()`, `VecGetOwnershipRanges()`

# External Links
$(_doc_external("Vec/VecFlag"))
"""
function VecFlag(petsclib::PetscLibType, xin::PetscVec, flg::PetscInt) end

@for_petsc function VecFlag(petsclib::$UnionPetscLib, xin::PetscVec, flg::$PetscInt )

    @chk ccall(
               (:VecFlag, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt),
               xin, flg,
              )


	return nothing
end 

"""
	VecBindToCPU(petsclib::PetscLibType,v::PetscVec, flg::PetscBool) 
marks a vector to temporarily stay on the CPU and perform computations on the CPU

Logically collective

Input Parameters:
- `v`   - the vector
- `flg` - bind to the CPU if value of `PETSC_TRUE`

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecBoundToCPU()`

# External Links
$(_doc_external("Vec/VecBindToCPU"))
"""
function VecBindToCPU(petsclib::PetscLibType, v::PetscVec, flg::PetscBool) end

@for_petsc function VecBindToCPU(petsclib::$UnionPetscLib, v::PetscVec, flg::PetscBool )

    @chk ccall(
               (:VecBindToCPU, $petsc_library),
               PetscErrorCode,
               (CVec, PetscBool),
               v, flg,
              )


	return nothing
end 

"""
	flg::PetscBool = VecBoundToCPU(petsclib::PetscLibType,v::PetscVec) 
query if a vector is bound to the CPU

Not collective

Input Parameter:
- `v` - the vector

Output Parameter:
- `flg` - the logical flag

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecBindToCPU()`

# External Links
$(_doc_external("Vec/VecBoundToCPU"))
"""
function VecBoundToCPU(petsclib::PetscLibType, v::PetscVec) end

@for_petsc function VecBoundToCPU(petsclib::$UnionPetscLib, v::PetscVec )
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
	flg::PetscBool = VecGetBindingPropagates(petsclib::PetscLibType,v::PetscVec) 
Gets whether the state of being bound to the CPU for a GPU vector type propagates to child and some other associated objects

Input Parameter:
- `v` - the vector

Output Parameter:
- `flg` - flag indicating whether the boundtocpu flag will be propagated

Level: developer

-seealso: [](ch_vectors), `Vec`, `VecSetBindingPropagates()`

# External Links
$(_doc_external("Vec/VecGetBindingPropagates"))
"""
function VecGetBindingPropagates(petsclib::PetscLibType, v::PetscVec) end

@for_petsc function VecGetBindingPropagates(petsclib::$UnionPetscLib, v::PetscVec )
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
	VecGetPinnedMemoryMin(petsclib::PetscLibType,v::PetscVec, mbytes::Csize_t) 
Get the minimum data size for which pinned memory will be used for host (CPU) allocations.

Logically Collective

Input Parameter:
- `v` - the vector

Output Parameter:
- `mbytes` - minimum data size in bytes

Level: developer

-seealso: [](ch_vectors), `Vec`, `VecSetPinnedMemoryMin()`

# External Links
$(_doc_external("Vec/VecGetPinnedMemoryMin"))
"""
function VecGetPinnedMemoryMin(petsclib::PetscLibType, v::PetscVec, mbytes::Csize_t) end

@for_petsc function VecGetPinnedMemoryMin(petsclib::$UnionPetscLib, v::PetscVec, mbytes::Csize_t )

    @chk ccall(
               (:VecGetPinnedMemoryMin, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{Csize_t}),
               v, mbytes,
              )


	return nothing
end 

"""
	VecGetOffloadMask(petsclib::PetscLibType,v::PetscVec, mask::PetscOffloadMask) 
Get the offload mask of a `Vec`

Not Collective

Input Parameter:
- `v` - the vector

Output Parameter:
- `mask` - corresponding `PetscOffloadMask` enum value.

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecCreateSeqCUDA()`, `VecCreateSeqViennaCL()`, `VecGetArray()`, `VecGetType()`

# External Links
$(_doc_external("Vec/VecGetOffloadMask"))
"""
function VecGetOffloadMask(petsclib::PetscLibType, v::PetscVec, mask::PetscOffloadMask) end

@for_petsc function VecGetOffloadMask(petsclib::$UnionPetscLib, v::PetscVec, mask::PetscOffloadMask )

    @chk ccall(
               (:VecGetOffloadMask, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{PetscOffloadMask}),
               v, mask,
              )


	return nothing
end 

"""
	norm::PetscReal,norm_loc::PetscInt,norma::PetscReal,norma_loc::PetscInt,normr::PetscReal,normr_loc::PetscInt = VecErrorWeightedNorms(petsclib::PetscLibType,U::PetscVec, Y::PetscVec, E::PetscVec, wnormtype::NormType, atol::PetscReal, vatol::PetscVec, rtol::PetscReal, vrtol::PetscVec, ignore_max::PetscReal) 
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

-seealso: [](ch_vectors), `Vec`, `NormType`, `TSErrorWeightedNorm()`, `TSErrorWeightedENorm()`

# External Links
$(_doc_external("Vec/VecErrorWeightedNorms"))
"""
function VecErrorWeightedNorms(petsclib::PetscLibType, U::PetscVec, Y::PetscVec, E::PetscVec, wnormtype::NormType, atol::PetscReal, vatol::PetscVec, rtol::PetscReal, vrtol::PetscVec, ignore_max::PetscReal) end

@for_petsc function VecErrorWeightedNorms(petsclib::$UnionPetscLib, U::PetscVec, Y::PetscVec, E::PetscVec, wnormtype::NormType, atol::$PetscReal, vatol::PetscVec, rtol::$PetscReal, vrtol::PetscVec, ignore_max::$PetscReal )
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
	VecRegisterAll(petsclib::PetscLibType) 
Registers all of the vector types in the `Vec` package.

Not Collective

Level: advanced

-seealso: [](ch_vectors), `Vec`, `VecType`, `VecRegister()`, `VecRegisterDestroy()`

# External Links
$(_doc_external("Vec/VecRegisterAll"))
"""
function VecRegisterAll(petsclib::PetscLibType) end

@for_petsc function VecRegisterAll(petsclib::$UnionPetscLib)

    @chk ccall(
               (:VecRegisterAll, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	type::VecType = VecGetType(petsclib::PetscLibType,vec::PetscVec) 
Gets the vector type name (as a string) from a `Vec`.

Not Collective

Input Parameter:
- `vec` - The vector

Output Parameter:
- `type` - The `VecType` of the vector

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecType`, `VecCreate()`, `VecDuplicate()`, `VecDuplicateVecs()`

# External Links
$(_doc_external("Vec/VecGetType"))
"""
function VecGetType(petsclib::PetscLibType, vec::AbstractPetscVec) end

@for_petsc function VecGetType(petsclib::$UnionPetscLib, vec::AbstractPetscVec )
	type_ = Ref{VecType}()

    @chk ccall(
               (:VecGetType, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{VecType}),
               vec, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	VecRegister(petsclib::PetscLibType,sname::Vector{Cchar}, fnc::external) 
Adds a new vector component implementation

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - The name of a new user-defined creation routine
- `function` - The creation routine

-seealso: `VecRegisterAll()`, `VecRegisterDestroy()`

# External Links
$(_doc_external("Vec/VecRegister"))
"""
function VecRegister(petsclib::PetscLibType, sname::Vector{Cchar}, fnc::external) end

@for_petsc function VecRegister(petsclib::$UnionPetscLib, sname::Vector{Cchar}, fnc::external )

    @chk ccall(
               (:VecRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	VecInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `Vec` package. It is called
from PetscDLLibraryRegister_petscvec() when using dynamic libraries, and on the first call to `VecCreate()`
when using shared or static libraries.

Level: developer

-seealso: `PetscInitialize()`

# External Links
$(_doc_external("Vec/VecInitializePackage"))
"""
function VecInitializePackage(petsclib::PetscLibType) end

@for_petsc function VecInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:VecInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	VecFinalizePackage(petsclib::PetscLibType) 
This function finalizes everything in the Vec package. It is called
from PetscFinalize().

Level: developer

-seealso: `PetscInitialize()`

# External Links
$(_doc_external("Vec/VecFinalizePackage"))
"""
function VecFinalizePackage(petsclib::PetscLibType) end

@for_petsc function VecFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:VecFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	vec::PetscVec = VecCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Creates an empty vector object. The type can then be set with `VecSetType()`,
or `VecSetFromOptions().`

Collective

Input Parameter:
- `comm` - The communicator for the vector object

Output Parameter:
- `vec` - The vector object

Level: beginner

-seealso: [](ch_vectors), `Vec`, `VecSetType()`, `VecSetSizes()`, `VecCreateMPIWithArray()`, `VecCreateMPI()`, `VecDuplicate()`,
`VecDuplicateVecs()`, `VecCreateGhost()`, `VecCreateSeq()`, `VecPlaceArray()`

# External Links
$(_doc_external("Vec/VecCreate"))
"""
function VecCreate(petsclib::PetscLibType, comm::MPI_Comm) end

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
	vec::PetscVec = VecCreateFromOptions(petsclib::PetscLibType,comm::MPI_Comm, prefix::Vector{Cchar}, bs::PetscInt, m::PetscInt, n::PetscInt) 
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

-seealso: [](ch_vectors), `Vec`, `VecSetType()`, `VecSetSizes()`, `VecCreateMPIWithArray()`, `VecCreateMPI()`, `VecDuplicate()`,
`VecDuplicateVecs()`, `VecCreateGhost()`, `VecCreateSeq()`, `VecPlaceArray()`, `VecCreate()`, `VecType`

# External Links
$(_doc_external("Vec/VecCreateFromOptions"))
"""
function VecCreateFromOptions(petsclib::PetscLibType, comm::MPI_Comm, prefix::Vector{Cchar}, bs::PetscInt, m::PetscInt, n::PetscInt) end

@for_petsc function VecCreateFromOptions(petsclib::$UnionPetscLib, comm::MPI_Comm, prefix::Vector{Cchar}, bs::$PetscInt, m::$PetscInt, n::$PetscInt )
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
	v::PetscVec = VecCreateMPI(petsclib::PetscLibType,comm::MPI_Comm, n::PetscInt, N::PetscInt) 
Creates a parallel vector.

Collective

Input Parameters:
- `comm` - the MPI communicator to use
- `n`    - local vector length (or `PETSC_DECIDE` to have calculated if `N` is given)
- `N`    - global vector length (or `PETSC_DETERMINE` to have calculated if `n` is given)

Output Parameter:
- `v` - the vector

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecType`, `VecCreateSeq()`, `VecCreate()`, `VecDuplicate()`, `VecDuplicateVecs()`, `VecCreateGhost()`,
`VecCreateMPIWithArray()`, `VecCreateGhostWithArray()`, `VecMPISetGhost()`, `PetscLayout`,
`VecGetOwnershipRange()`, `VecGetOwnershipRanges()`

# External Links
$(_doc_external("Vec/VecCreateMPI"))
"""
function VecCreateMPI(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt, N::PetscInt) end

@for_petsc function VecCreateMPI(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, N::$PetscInt )
	v_ = Ref{CVec}()

    @chk ccall(
               (:VecCreateMPI, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, Ptr{CVec}),
               comm, n, N, v_,
              )

	v = PetscVec(v_[], petsclib)

	return v
end 

"""
	VecGhostGetLocalForm(petsclib::PetscLibType,g::PetscVec, l::PetscVec) 
Obtains the local ghosted representation of
a parallel vector (obtained with `VecCreateGhost()`, `VecCreateGhostWithArray()` or `VecCreateSeq()`).

Logically Collective

Input Parameter:
- `g` - the global vector

Output Parameter:
- `l` - the local (ghosted) representation,`NULL` if `g` is not ghosted

Level: advanced

-seealso: [](ch_vectors), `VecGhostUpdateBegin()`, `VecGhostUpdateEnd()`, `Vec`, `VecType`, `VecCreateGhost()`, `VecGhostRestoreLocalForm()`, `VecCreateGhostWithArray()`

# External Links
$(_doc_external("Vec/VecGhostGetLocalForm"))
"""
function VecGhostGetLocalForm(petsclib::PetscLibType, g::PetscVec, l::PetscVec) end

@for_petsc function VecGhostGetLocalForm(petsclib::$UnionPetscLib, g::PetscVec, l::PetscVec )
	l_ = Ref(l.ptr)

    @chk ccall(
               (:VecGhostGetLocalForm, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{CVec}),
               g, l_,
              )

	l.ptr = C_NULL

	return nothing
end 

"""
	flg::PetscBool = VecGhostIsLocalForm(petsclib::PetscLibType,g::PetscVec, l::PetscVec) 
Checks if a given vector is the local form of a global vector

Not Collective

Input Parameters:
- `g` - the global vector
- `l` - the local vector

Output Parameter:
- `flg` - `PETSC_TRUE` if `l` is the local form

Level: advanced

-seealso: [](ch_vectors), `Vec`, `VecType`, `VecCreateGhost()`, `VecGhostRestoreLocalForm()`, `VecCreateGhostWithArray()`, `VecGhostGetLocalForm()`

# External Links
$(_doc_external("Vec/VecGhostIsLocalForm"))
"""
function VecGhostIsLocalForm(petsclib::PetscLibType, g::PetscVec, l::PetscVec) end

@for_petsc function VecGhostIsLocalForm(petsclib::$UnionPetscLib, g::PetscVec, l::PetscVec )
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
	VecGhostRestoreLocalForm(petsclib::PetscLibType,g::PetscVec, l::PetscVec) 
Restores the local ghosted representation of
a parallel vector obtained with `VecGhostGetLocalForm()`.

Logically Collective

Input Parameters:
- `g` - the global vector
- `l` - the local (ghosted) representation

Level: advanced

-seealso: [](ch_vectors), `VecGhostUpdateBegin()`, `VecGhostUpdateEnd()`, `Vec`, `VecType`, `VecCreateGhost()`, `VecGhostGetLocalForm()`, `VecCreateGhostWithArray()`

# External Links
$(_doc_external("Vec/VecGhostRestoreLocalForm"))
"""
function VecGhostRestoreLocalForm(petsclib::PetscLibType, g::PetscVec, l::PetscVec) end

@for_petsc function VecGhostRestoreLocalForm(petsclib::$UnionPetscLib, g::PetscVec, l::PetscVec )
	l_ = Ref(l.ptr)

    @chk ccall(
               (:VecGhostRestoreLocalForm, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{CVec}),
               g, l_,
              )

	l.ptr = C_NULL

	return nothing
end 

"""
	VecGhostUpdateBegin(petsclib::PetscLibType,g::PetscVec, insertmode::InsertMode, scattermode::ScatterMode) 
Begins the vector scatter to update the vector from
local representation to global or global representation to local.

Neighbor-wise Collective

Input Parameters:
- `g`           - the vector (obtained with `VecCreateGhost()` or `VecDuplicate()`)
- `insertmode`  - one of `ADD_VALUES`, `MAX_VALUES`, `MIN_VALUES` or `INSERT_VALUES`
- `scattermode` - one of `SCATTER_FORWARD` (update ghosts) or `SCATTER_REVERSE` (update local values from ghosts)

Level: advanced

-seealso: [](ch_vectors), `Vec`, `VecType`, `VecCreateGhost()`, `VecGhostUpdateEnd()`, `VecGhostGetLocalForm()`,
`VecGhostRestoreLocalForm()`, `VecCreateGhostWithArray()`

# External Links
$(_doc_external("Vec/VecGhostUpdateBegin"))
"""
function VecGhostUpdateBegin(petsclib::PetscLibType, g::PetscVec, insertmode::InsertMode, scattermode::ScatterMode) end

@for_petsc function VecGhostUpdateBegin(petsclib::$UnionPetscLib, g::PetscVec, insertmode::InsertMode, scattermode::ScatterMode )

    @chk ccall(
               (:VecGhostUpdateBegin, $petsc_library),
               PetscErrorCode,
               (CVec, InsertMode, ScatterMode),
               g, insertmode, scattermode,
              )


	return nothing
end 

"""
	VecGhostUpdateEnd(petsclib::PetscLibType,g::PetscVec, insertmode::InsertMode, scattermode::ScatterMode) 
End the vector scatter to update the vector from
local representation to global or global representation to local.

Neighbor-wise Collective

Input Parameters:
- `g`           - the vector (obtained with `VecCreateGhost()` or `VecDuplicate()`)
- `insertmode`  - one of `ADD_VALUES`, `MAX_VALUES`, `MIN_VALUES` or `INSERT_VALUES`
- `scattermode` - one of `SCATTER_FORWARD` (update ghosts) or `SCATTER_REVERSE` (update local values from ghosts)

Level: advanced

-seealso: [](ch_vectors), `Vec`, `VecType`, `VecCreateGhost()`, `VecGhostUpdateBegin()`, `VecGhostGetLocalForm()`,
`VecGhostRestoreLocalForm()`, `VecCreateGhostWithArray()`

# External Links
$(_doc_external("Vec/VecGhostUpdateEnd"))
"""
function VecGhostUpdateEnd(petsclib::PetscLibType, g::PetscVec, insertmode::InsertMode, scattermode::ScatterMode) end

@for_petsc function VecGhostUpdateEnd(petsclib::$UnionPetscLib, g::PetscVec, insertmode::InsertMode, scattermode::ScatterMode )

    @chk ccall(
               (:VecGhostUpdateEnd, $petsc_library),
               PetscErrorCode,
               (CVec, InsertMode, ScatterMode),
               g, insertmode, scattermode,
              )


	return nothing
end 

"""
	vv::PetscVec = VecCreateMPIWithArray(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, n::PetscInt, N::PetscInt, array::Vector{PetscScalar}) 
Creates a parallel, array
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

-seealso: [](ch_vectors), `Vec`, `VecType`, `VecCreateSeqWithArray()`, `VecCreate()`, `VecDuplicate()`, `VecDuplicateVecs()`, `VecCreateGhost()`,
`VecCreateMPI()`, `VecCreateGhostWithArray()`, `VecPlaceArray()`

# External Links
$(_doc_external("Vec/VecCreateMPIWithArray"))
"""
function VecCreateMPIWithArray(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, n::PetscInt, N::PetscInt, array::Vector{PetscScalar}) end

@for_petsc function VecCreateMPIWithArray(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, n::$PetscInt, N::$PetscInt, array::Vector{$PetscScalar} )
	vv_ = Ref{CVec}()

    @chk ccall(
               (:VecCreateMPIWithArray, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscScalar}, Ptr{CVec}),
               comm, bs, n, N, array, vv_,
              )

	vv = PetscVec(vv_[], petsclib)

	return vv
end 

"""
	vv::PetscVec = VecCreateGhostWithArray(petsclib::PetscLibType,comm::MPI_Comm, n::PetscInt, N::PetscInt, nghost::PetscInt, ghosts::Vector{PetscInt}, array::Vector{PetscScalar}) 
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

-seealso: [](ch_vectors), `Vec`, `VecType`, `VecCreate()`, `VecGhostGetLocalForm()`, `VecGhostRestoreLocalForm()`,
`VecCreateGhost()`, `VecCreateSeqWithArray()`, `VecCreateMPIWithArray()`,
`VecCreateGhostBlock()`, `VecCreateGhostBlockWithArray()`, `VecMPISetGhost()`, `VecGhostUpdateBegin()`, `VecGhostUpdateEnd()`

# External Links
$(_doc_external("Vec/VecCreateGhostWithArray"))
"""
function VecCreateGhostWithArray(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt, N::PetscInt, nghost::PetscInt, ghosts::Vector{PetscInt}, array::Vector{PetscScalar}) end

@for_petsc function VecCreateGhostWithArray(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, N::$PetscInt, nghost::$PetscInt, ghosts::Vector{$PetscInt}, array::Vector{$PetscScalar} )
	vv_ = Ref{CVec}()

    @chk ccall(
               (:VecCreateGhostWithArray, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}, Ptr{CVec}),
               comm, n, N, nghost, ghosts, array, vv_,
              )

	vv = PetscVec(vv_[], petsclib)

	return vv
end 

"""
	VecGhostGetGhostIS(petsclib::PetscLibType,X::PetscVec, ghost::IS) 
Return ghosting indices of a ghost vector

Input Parameters:
- `X` - ghost vector

Output Parameter:
- `ghost` - ghosting indices

Level: beginner

-seealso: `VecCreateGhostWithArray()`, `VecCreateMPIWithArray()`

# External Links
$(_doc_external("Vec/VecGhostGetGhostIS"))
"""
function VecGhostGetGhostIS(petsclib::PetscLibType, X::PetscVec, ghost::IS) end

@for_petsc function VecGhostGetGhostIS(petsclib::$UnionPetscLib, X::PetscVec, ghost::IS )

    @chk ccall(
               (:VecGhostGetGhostIS, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{IS}),
               X, ghost,
              )


	return nothing
end 

"""
	vv::PetscVec = VecCreateGhost(petsclib::PetscLibType,comm::MPI_Comm, n::PetscInt, N::PetscInt, nghost::PetscInt, ghosts::Vector{PetscInt}) 
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

-seealso: [](ch_vectors), `Vec`, `VecType`, `VecCreateSeq()`, `VecCreate()`, `VecDuplicate()`, `VecDuplicateVecs()`, `VecCreateMPI()`,
`VecGhostGetLocalForm()`, `VecGhostRestoreLocalForm()`, `VecGhostUpdateBegin()`,
`VecCreateGhostWithArray()`, `VecCreateMPIWithArray()`, `VecGhostUpdateEnd()`,
`VecCreateGhostBlock()`, `VecCreateGhostBlockWithArray()`, `VecMPISetGhost()`


# External Links
$(_doc_external("Vec/VecCreateGhost"))
"""
function VecCreateGhost(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt, N::PetscInt, nghost::PetscInt, ghosts::Vector{PetscInt}) end

@for_petsc function VecCreateGhost(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, N::$PetscInt, nghost::$PetscInt, ghosts::Vector{$PetscInt} )
	vv_ = Ref{CVec}()

    @chk ccall(
               (:VecCreateGhost, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{CVec}),
               comm, n, N, nghost, ghosts, vv_,
              )

	vv = PetscVec(vv_[], petsclib)

	return vv
end 

"""
	VecMPISetGhost(petsclib::PetscLibType,vv::PetscVec, nghost::PetscInt, ghosts::Vector{PetscInt}) 
Sets the ghost points for an MPI ghost vector

Collective

Input Parameters:
- `vv`     - the MPI vector
- `nghost` - number of local ghost points
- `ghosts` - global indices of ghost points, these do not need to be in increasing order (sorted)

Level: advanced

-seealso: [](ch_vectors), `Vec`, `VecType`, `VecCreateSeq()`, `VecCreate()`, `VecDuplicate()`, `VecDuplicateVecs()`, `VecCreateMPI()`,
`VecGhostGetLocalForm()`, `VecGhostRestoreLocalForm()`, `VecGhostUpdateBegin()`,
`VecCreateGhostWithArray()`, `VecCreateMPIWithArray()`, `VecGhostUpdateEnd()`,
`VecCreateGhostBlock()`, `VecCreateGhostBlockWithArray()`

# External Links
$(_doc_external("Vec/VecMPISetGhost"))
"""
function VecMPISetGhost(petsclib::PetscLibType, vv::PetscVec, nghost::PetscInt, ghosts::Vector{PetscInt}) end

@for_petsc function VecMPISetGhost(petsclib::$UnionPetscLib, vv::PetscVec, nghost::$PetscInt, ghosts::Vector{$PetscInt} )

    @chk ccall(
               (:VecMPISetGhost, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{$PetscInt}),
               vv, nghost, ghosts,
              )


	return nothing
end 

"""
	vv::PetscVec = VecCreateGhostBlockWithArray(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, n::PetscInt, N::PetscInt, nghost::PetscInt, ghosts::Vector{PetscInt}, array::Vector{PetscScalar}) 
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

-seealso: [](ch_vectors), `Vec`, `VecType`, `VecCreate()`, `VecGhostGetLocalForm()`, `VecGhostRestoreLocalForm()`,
`VecCreateGhost()`, `VecCreateSeqWithArray()`, `VecCreateMPIWithArray()`,
`VecCreateGhostWithArray()`, `VecCreateGhostBlock()`, `VecGhostUpdateBegin()`, `VecGhostUpdateEnd()`

# External Links
$(_doc_external("Vec/VecCreateGhostBlockWithArray"))
"""
function VecCreateGhostBlockWithArray(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, n::PetscInt, N::PetscInt, nghost::PetscInt, ghosts::Vector{PetscInt}, array::Vector{PetscScalar}) end

@for_petsc function VecCreateGhostBlockWithArray(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, n::$PetscInt, N::$PetscInt, nghost::$PetscInt, ghosts::Vector{$PetscInt}, array::Vector{$PetscScalar} )
	vv_ = Ref{CVec}()

    @chk ccall(
               (:VecCreateGhostBlockWithArray, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscScalar}, Ptr{CVec}),
               comm, bs, n, N, nghost, ghosts, array, vv_,
              )

	vv = PetscVec(vv_[], petsclib)

	return vv
end 

"""
	vv::PetscVec = VecCreateGhostBlock(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, n::PetscInt, N::PetscInt, nghost::PetscInt, ghosts::Vector{PetscInt}) 
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

-seealso: [](ch_vectors), `Vec`, `VecType`, `VecCreateSeq()`, `VecCreate()`, `VecDuplicate()`, `VecDuplicateVecs()`, `VecCreateMPI()`,
`VecGhostGetLocalForm()`, `VecGhostRestoreLocalForm()`, `VecGhostUpdateBegin()`, `VecGhostUpdateEnd()`
`VecCreateGhostWithArray()`, `VecCreateMPIWithArray()`, `VecCreateGhostBlockWithArray()`

# External Links
$(_doc_external("Vec/VecCreateGhostBlock"))
"""
function VecCreateGhostBlock(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, n::PetscInt, N::PetscInt, nghost::PetscInt, ghosts::Vector{PetscInt}) end

@for_petsc function VecCreateGhostBlock(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, n::$PetscInt, N::$PetscInt, nghost::$PetscInt, ghosts::Vector{$PetscInt} )
	vv_ = Ref{CVec}()

    @chk ccall(
               (:VecCreateGhostBlock, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscInt}, Ptr{CVec}),
               comm, bs, n, N, nghost, ghosts, vv_,
              )

	vv = PetscVec(vv_[], petsclib)

	return vv
end 

"""
	v::PetscVec = VecCreateMPIKokkosWithArray(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, n::PetscInt, N::PetscInt, darray::Vector{PetscScalar}) 

# External Links
$(_doc_external("Vec/VecCreateMPIKokkosWithArray"))
"""
function VecCreateMPIKokkosWithArray(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, n::PetscInt, N::PetscInt, darray::Vector{PetscScalar}) end

@for_petsc function VecCreateMPIKokkosWithArray(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, n::$PetscInt, N::$PetscInt, darray::Vector{$PetscScalar} )
	v_ = Ref{CVec}()

    @chk ccall(
               (:VecCreateMPIKokkosWithArray, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscScalar}, Ptr{CVec}),
               comm, bs, n, N, darray, v_,
              )

	v = PetscVec(v_[], petsclib)

	return v
end 

"""
	array::ViennaCLVector,vv::PetscVec = VecCreateMPIViennaCLWithArray(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, n::PetscInt, N::PetscInt) 

# External Links
$(_doc_external("Vec/VecCreateMPIViennaCLWithArray"))
"""
function VecCreateMPIViennaCLWithArray(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, n::PetscInt, N::PetscInt) end

@for_petsc function VecCreateMPIViennaCLWithArray(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, n::$PetscInt, N::$PetscInt )
	array_ = Ref{ViennaCLVector}()
	vv_ = Ref{CVec}()

    @chk ccall(
               (:VecCreateMPIViennaCLWithArray, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{ViennaCLVector}, Ptr{CVec}),
               comm, bs, n, N, array_, vv_,
              )

	array = array_[]
	vv = PetscVec(vv_[], petsclib)

	return array,vv
end 

"""
	viennaclvec::ViennaCLVector,vv::PetscVec = VecCreateMPIViennaCLWithArrays(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, n::PetscInt, N::PetscInt, cpuarray::Vector{PetscScalar}) 

# External Links
$(_doc_external("Vec/VecCreateMPIViennaCLWithArrays"))
"""
function VecCreateMPIViennaCLWithArrays(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, n::PetscInt, N::PetscInt, cpuarray::Vector{PetscScalar}) end

@for_petsc function VecCreateMPIViennaCLWithArrays(petsclib::$UnionPetscLib, comm::MPI_Comm, bs::$PetscInt, n::$PetscInt, N::$PetscInt, cpuarray::Vector{$PetscScalar} )
	viennaclvec_ = Ref{ViennaCLVector}()
	vv_ = Ref{CVec}()

    @chk ccall(
               (:VecCreateMPIViennaCLWithArrays, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscScalar}, Ptr{ViennaCLVector}, Ptr{CVec}),
               comm, bs, n, N, cpuarray, viennaclvec_, vv_,
              )

	viennaclvec = viennaclvec_[]
	vv = PetscVec(vv_[], petsclib)

	return viennaclvec,vv
end 

"""
	v::PetscVec = VecCreateSeq(petsclib::PetscLibType,comm::MPI_Comm, n::PetscInt) 
Creates a standard, sequential array

Collective

Input Parameters:
- `comm` - the communicator, should be `PETSC_COMM_SELF`
- `n`    - the vector length

Output Parameter:
- `v` - the vector

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecType`, `VecCreateMPI()`, `VecCreate()`, `VecDuplicate()`, `VecDuplicateVecs()`, `VecCreateGhost()`

# External Links
$(_doc_external("Vec/VecCreateSeq"))
"""
function VecCreateSeq(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt) end

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
	V::PetscVec = VecCreateSeqWithArray(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, n::PetscInt, array::Vector{PetscScalar}) 
Creates a standard,sequential array
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

-seealso: `VecCreateMPIWithArray()`, `VecCreate()`, `VecDuplicate()`, `VecDuplicateVecs()`,
`VecCreateGhost()`, `VecCreateSeq()`, `VecPlaceArray()`

# External Links
$(_doc_external("Vec/VecCreateSeqWithArray"))
"""
function VecCreateSeqWithArray(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, n::PetscInt, array::Vector{PetscScalar}) end

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
	a::PetscScalar = VecKokkosPlaceArray(petsclib::PetscLibType,v::PetscVec) 

# External Links
$(_doc_external("Vec/VecKokkosPlaceArray"))
"""
function VecKokkosPlaceArray(petsclib::PetscLibType, v::PetscVec) end

@for_petsc function VecKokkosPlaceArray(petsclib::$UnionPetscLib, v::PetscVec )
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
	VecKokkosResetArray(petsclib::PetscLibType,v::PetscVec) 

# External Links
$(_doc_external("Vec/VecKokkosResetArray"))
"""
function VecKokkosResetArray(petsclib::PetscLibType, v::PetscVec) end

@for_petsc function VecKokkosResetArray(petsclib::$UnionPetscLib, v::PetscVec )

    @chk ccall(
               (:VecKokkosResetArray, $petsc_library),
               PetscErrorCode,
               (CVec,),
               v,
              )


	return nothing
end 

"""
	v::PetscVec = VecCreateSeqKokkosWithArray(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, n::PetscInt, darray::Vector{PetscScalar}) 

# External Links
$(_doc_external("Vec/VecCreateSeqKokkosWithArray"))
"""
function VecCreateSeqKokkosWithArray(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, n::PetscInt, darray::Vector{PetscScalar}) end

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
	v::PetscVec = VecCreateSeqKokkos(petsclib::PetscLibType,comm::MPI_Comm, n::PetscInt) 

# External Links
$(_doc_external("Vec/VecCreateSeqKokkos"))
"""
function VecCreateSeqKokkos(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt) end

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
	v::PetscVec = VecCreateSeqViennaCL(petsclib::PetscLibType,comm::MPI_Comm, n::PetscInt) 

# External Links
$(_doc_external("Vec/VecCreateSeqViennaCL"))
"""
function VecCreateSeqViennaCL(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt) end

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
	viennaclvec::ViennaCLVector,V::PetscVec = VecCreateSeqViennaCLWithArrays(petsclib::PetscLibType,comm::MPI_Comm, bs::PetscInt, n::PetscInt, cpuarray::Vector{PetscScalar}) 

# External Links
$(_doc_external("Vec/VecCreateSeqViennaCLWithArrays"))
"""
function VecCreateSeqViennaCLWithArrays(petsclib::PetscLibType, comm::MPI_Comm, bs::PetscInt, n::PetscInt, cpuarray::Vector{PetscScalar}) end

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
	v::PetscVec = VecCreateShared(petsclib::PetscLibType,comm::MPI_Comm, n::PetscInt, N::PetscInt) 
Creates a parallel vector that uses shared memory.

Collective

Input Parameters:
- `comm` - the MPI communicator to use
- `n`    - local vector length (or `PETSC_DECIDE` to have calculated if `N` is given)
- `N`    - global vector length (or `PETSC_DECIDE` to have calculated if `n` is given)

Output Parameter:
- `v` - the vector

Level: advanced

-seealso: [](ch_vectors), `Vec`, `VecType`, `VecCreateSeq()`, `VecCreate()`, `VecCreateMPI()`, `VecDuplicate()`, `VecDuplicateVecs()`,
`VecCreateGhost()`, `VecCreateMPIWithArray()`, `VecCreateGhostWithArray()`

# External Links
$(_doc_external("Vec/VecCreateShared"))
"""
function VecCreateShared(petsclib::PetscLibType, comm::MPI_Comm, n::PetscInt, N::PetscInt) end

@for_petsc function VecCreateShared(petsclib::$UnionPetscLib, comm::MPI_Comm, n::$PetscInt, N::$PetscInt )
	v_ = Ref{CVec}()

    @chk ccall(
               (:VecCreateShared, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, $PetscInt, Ptr{CVec}),
               comm, n, N, v_,
              )

	v = PetscVec(v_[], petsclib)

	return v
end 

"""
	VecNestGetSubVec(petsclib::PetscLibType,X::PetscVec, idxm::PetscInt, sx::PetscVec) 
Returns a single, sub

Not Collective

Input Parameters:
- `X`    - nest vector
- `idxm` - index of the vector within the nest

Output Parameter:
- `sx` - vector at index `idxm` within the nest

Level: developer

-seealso: `VECNEST`,  [](ch_vectors), `Vec`, `VecType`, `VecNestGetSize()`, `VecNestGetSubVecs()`

# External Links
$(_doc_external("Vec/VecNestGetSubVec"))
"""
function VecNestGetSubVec(petsclib::PetscLibType, X::PetscVec, idxm::PetscInt, sx::PetscVec) end

@for_petsc function VecNestGetSubVec(petsclib::$UnionPetscLib, X::PetscVec, idxm::$PetscInt, sx::PetscVec )
	sx_ = Ref(sx.ptr)

    @chk ccall(
               (:VecNestGetSubVec, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{CVec}),
               X, idxm, sx_,
              )

	sx.ptr = C_NULL

	return nothing
end 

"""
	N::PetscInt = VecNestGetSubVecs(petsclib::PetscLibType,X::PetscVec, sx::Vector{PetscVec}) 
Returns the entire array of vectors defining a nest vector.

Not Collective

Input Parameter:
- `X` - nest vector

Output Parameters:
- `N`  - number of nested vecs
- `sx` - array of vectors, can pass in `NULL`

Level: developer

-seealso: `VECNEST`,  [](ch_vectors), `Vec`, `VecType`, `VecNestGetSize()`, `VecNestGetSubVec()`

# External Links
$(_doc_external("Vec/VecNestGetSubVecs"))
"""
function VecNestGetSubVecs(petsclib::PetscLibType, X::PetscVec, sx::Vector{PetscVec}) end

@for_petsc function VecNestGetSubVecs(petsclib::$UnionPetscLib, X::PetscVec, sx::Vector{PetscVec} )
	N_ = Ref{$PetscInt}()
	sx_ = Ref(pointer(sx))

    @chk ccall(
               (:VecNestGetSubVecs, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{$PetscInt}, Ptr{Ptr{CVec}}),
               X, N_, sx_,
              )

	N = N_[]

	return N
end 

"""
	VecNestSetSubVec(petsclib::PetscLibType,X::PetscVec, idxm::PetscInt, sx::PetscVec) 
Set a single component vector in a nest vector at specified index.

Not Collective

Input Parameters:
- `X`    - nest vector
- `idxm` - index of the vector within the nest vector
- `sx`   - vector at index `idxm` within the nest vector

Level: developer

-seealso: `VECNEST`,  [](ch_vectors), `Vec`, `VecType`, `VecNestSetSubVecs()`, `VecNestGetSubVec()`

# External Links
$(_doc_external("Vec/VecNestSetSubVec"))
"""
function VecNestSetSubVec(petsclib::PetscLibType, X::PetscVec, idxm::PetscInt, sx::PetscVec) end

@for_petsc function VecNestSetSubVec(petsclib::$UnionPetscLib, X::PetscVec, idxm::$PetscInt, sx::PetscVec )

    @chk ccall(
               (:VecNestSetSubVec, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, CVec),
               X, idxm, sx,
              )


	return nothing
end 

"""
	VecNestSetSubVecs(petsclib::PetscLibType,X::PetscVec, N::PetscInt, idxm::Vector{PetscInt}, sx::Vector{PetscVec}) 
Sets the component vectors at the specified indices in a nest vector.

Not Collective

Input Parameters:
- `X`    - nest vector
- `N`    - number of component vecs in `sx`
- `idxm` - indices of component vectors that are to be replaced
- `sx`   - array of vectors

Level: developer

-seealso: `VECNEST`,  [](ch_vectors), `Vec`, `VecType`, `VecNestGetSize()`, `VecNestGetSubVec()`

# External Links
$(_doc_external("Vec/VecNestSetSubVecs"))
"""
function VecNestSetSubVecs(petsclib::PetscLibType, X::PetscVec, N::PetscInt, idxm::Vector{PetscInt}, sx::Vector{PetscVec}) end

@for_petsc function VecNestSetSubVecs(petsclib::$UnionPetscLib, X::PetscVec, N::$PetscInt, idxm::Vector{$PetscInt}, sx::Vector{PetscVec} )

    @chk ccall(
               (:VecNestSetSubVecs, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{$PetscInt}, Ptr{CVec}),
               X, N, idxm, sx,
              )


	return nothing
end 

"""
	N::PetscInt = VecNestGetSize(petsclib::PetscLibType,X::PetscVec) 
Returns the size of the nest vector.

Not Collective

Input Parameter:
- `X` - nest vector

Output Parameter:
- `N` - number of nested vecs

Level: developer

-seealso: `VECNEST`,  [](ch_vectors), `Vec`, `VecType`, `VecNestGetSubVec()`, `VecNestGetSubVecs()`

# External Links
$(_doc_external("Vec/VecNestGetSize"))
"""
function VecNestGetSize(petsclib::PetscLibType, X::PetscVec) end

@for_petsc function VecNestGetSize(petsclib::$UnionPetscLib, X::PetscVec )
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
	Y::PetscVec = VecCreateNest(petsclib::PetscLibType,comm::MPI_Comm, nb::PetscInt, is::Vector{IS}, x::Vector{PetscVec}) 
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

-seealso: `VECNEST`,  [](ch_vectors), `Vec`, `VecType`, `VecCreate()`, `MatCreateNest()`, `DMSetVecType()`

# External Links
$(_doc_external("Vec/VecCreateNest"))
"""
function VecCreateNest(petsclib::PetscLibType, comm::MPI_Comm, nb::PetscInt, is::Vector{IS}, x::Vector{PetscVec}) end

@for_petsc function VecCreateNest(petsclib::$UnionPetscLib, comm::MPI_Comm, nb::$PetscInt, is::Vector{IS}, x::Vector{PetscVec} )
	Y_ = Ref{CVec}()

    @chk ccall(
               (:VecCreateNest, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, Ptr{IS}, Ptr{CVec}, Ptr{CVec}),
               comm, nb, is, x, Y_,
              )

	Y = PetscVec(Y_[], petsclib)

	return Y
end 

"""
	VecDotBegin(petsclib::PetscLibType,x::PetscVec, y::PetscVec, result::PetscScalar) 
Starts a split phase dot product computation.

Input Parameters:
- `x`      - the first vector
- `y`      - the second vector
- `result` - where the result will go (can be `NULL`)

Level: advanced

-seealso: `VecDotEnd()`, `VecNormBegin()`, `VecNormEnd()`, `VecNorm()`, `VecDot()`, `VecMDot()`,
`VecTDotBegin()`, `VecTDotEnd()`, `PetscCommSplitReductionBegin()`

# External Links
$(_doc_external("Vec/VecDotBegin"))
"""
function VecDotBegin(petsclib::PetscLibType, x::PetscVec, y::PetscVec, result::PetscScalar) end

@for_petsc function VecDotBegin(petsclib::$UnionPetscLib, x::PetscVec, y::PetscVec, result::$PetscScalar )

    @chk ccall(
               (:VecDotBegin, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, Ptr{$PetscScalar}),
               x, y, result,
              )


	return nothing
end 

"""
	VecDotEnd(petsclib::PetscLibType,x::PetscVec, y::PetscVec, result::PetscScalar) 
Ends a split phase dot product computation.

Input Parameters:
- `x`      - the first vector (can be `NULL`)
- `y`      - the second vector (can be `NULL`)
- `result` - where the result will go

Level: advanced

-seealso: `VecDotBegin()`, `VecNormBegin()`, `VecNormEnd()`, `VecNorm()`, `VecDot()`, `VecMDot()`,
`VecTDotBegin()`, `VecTDotEnd()`, `PetscCommSplitReductionBegin()`

# External Links
$(_doc_external("Vec/VecDotEnd"))
"""
function VecDotEnd(petsclib::PetscLibType, x::PetscVec, y::PetscVec, result::PetscScalar) end

@for_petsc function VecDotEnd(petsclib::$UnionPetscLib, x::PetscVec, y::PetscVec, result::$PetscScalar )

    @chk ccall(
               (:VecDotEnd, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, Ptr{$PetscScalar}),
               x, y, result,
              )


	return nothing
end 

"""
	VecTDotBegin(petsclib::PetscLibType,x::PetscVec, y::PetscVec, result::PetscScalar) 
Starts a split phase transpose dot product computation.

Input Parameters:
- `x`      - the first vector
- `y`      - the second vector
- `result` - where the result will go (can be `NULL`)

Level: advanced

-seealso: `VecTDotEnd()`, `VecNormBegin()`, `VecNormEnd()`, `VecNorm()`, `VecDot()`, `VecMDot()`,
`VecDotBegin()`, `VecDotEnd()`, `PetscCommSplitReductionBegin()`

# External Links
$(_doc_external("Vec/VecTDotBegin"))
"""
function VecTDotBegin(petsclib::PetscLibType, x::PetscVec, y::PetscVec, result::PetscScalar) end

@for_petsc function VecTDotBegin(petsclib::$UnionPetscLib, x::PetscVec, y::PetscVec, result::$PetscScalar )

    @chk ccall(
               (:VecTDotBegin, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, Ptr{$PetscScalar}),
               x, y, result,
              )


	return nothing
end 

"""
	VecTDotEnd(petsclib::PetscLibType,x::PetscVec, y::PetscVec, result::PetscScalar) 
Ends a split phase transpose dot product computation.

Input Parameters:
- `x`      - the first vector (can be `NULL`)
- `y`      - the second vector (can be `NULL`)
- `result` - where the result will go

Level: advanced

-seealso: `VecTDotBegin()`, `VecNormBegin()`, `VecNormEnd()`, `VecNorm()`, `VecDot()`, `VecMDot()`,
`VecDotBegin()`, `VecDotEnd()`

# External Links
$(_doc_external("Vec/VecTDotEnd"))
"""
function VecTDotEnd(petsclib::PetscLibType, x::PetscVec, y::PetscVec, result::PetscScalar) end

@for_petsc function VecTDotEnd(petsclib::$UnionPetscLib, x::PetscVec, y::PetscVec, result::$PetscScalar )

    @chk ccall(
               (:VecTDotEnd, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, Ptr{$PetscScalar}),
               x, y, result,
              )


	return nothing
end 

"""
	VecNormBegin(petsclib::PetscLibType,x::PetscVec, ntype::NormType, result::PetscReal) 
Starts a split phase norm computation.

Input Parameters:
- `x`      - the first vector
- `ntype`  - norm type, one of `NORM_1`, `NORM_2`, `NORM_MAX`, `NORM_1_AND_2`
- `result` - where the result will go (can be `NULL`)

Level: advanced

-seealso: `VecNormEnd()`, `VecNorm()`, `VecDot()`, `VecMDot()`, `VecDotBegin()`, `VecDotEnd()`, `PetscCommSplitReductionBegin()`

# External Links
$(_doc_external("Vec/VecNormBegin"))
"""
function VecNormBegin(petsclib::PetscLibType, x::PetscVec, ntype::NormType, result::PetscReal) end

@for_petsc function VecNormBegin(petsclib::$UnionPetscLib, x::PetscVec, ntype::NormType, result::$PetscReal )

    @chk ccall(
               (:VecNormBegin, $petsc_library),
               PetscErrorCode,
               (CVec, NormType, Ptr{$PetscReal}),
               x, ntype, result,
              )


	return nothing
end 

"""
	VecNormEnd(petsclib::PetscLibType,x::PetscVec, ntype::NormType, result::PetscReal) 
Ends a split phase norm computation.

Input Parameters:
- `x`      - the first vector
- `ntype`  - norm type, one of `NORM_1`, `NORM_2`, `NORM_MAX`, `NORM_1_AND_2`
- `result` - where the result will go

Level: advanced

-seealso: `VecNormBegin()`, `VecNorm()`, `VecDot()`, `VecMDot()`, `VecDotBegin()`, `VecDotEnd()`, `PetscCommSplitReductionBegin()`

# External Links
$(_doc_external("Vec/VecNormEnd"))
"""
function VecNormEnd(petsclib::PetscLibType, x::PetscVec, ntype::NormType, result::PetscReal) end

@for_petsc function VecNormEnd(petsclib::$UnionPetscLib, x::PetscVec, ntype::NormType, result::$PetscReal )

    @chk ccall(
               (:VecNormEnd, $petsc_library),
               PetscErrorCode,
               (CVec, NormType, Ptr{$PetscReal}),
               x, ntype, result,
              )


	return nothing
end 

"""
	VecMDotBegin(petsclib::PetscLibType,x::PetscVec, nv::PetscInt, y::Vector{PetscVec}, result::Vector{PetscScalar}) 
Starts a split phase multiple dot product computation.

Input Parameters:
- `x`      - the first vector
- `nv`     - number of vectors
- `y`      - array of vectors
- `result` - where the result will go (can be `NULL`)

Level: advanced

-seealso: `VecMDotEnd()`, `VecNormBegin()`, `VecNormEnd()`, `VecNorm()`, `VecDot()`, `VecMDot()`,
`VecTDotBegin()`, `VecTDotEnd()`, `VecMTDotBegin()`, `VecMTDotEnd()`, `PetscCommSplitReductionBegin()`

# External Links
$(_doc_external("Vec/VecMDotBegin"))
"""
function VecMDotBegin(petsclib::PetscLibType, x::PetscVec, nv::PetscInt, y::Vector{PetscVec}, result::Vector{PetscScalar}) end

@for_petsc function VecMDotBegin(petsclib::$UnionPetscLib, x::PetscVec, nv::$PetscInt, y::Vector{PetscVec}, result::Vector{$PetscScalar} )

    @chk ccall(
               (:VecMDotBegin, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{CVec}, Ptr{$PetscScalar}),
               x, nv, y, result,
              )


	return nothing
end 

"""
	result::Vector{PetscScalar} = VecMDotEnd(petsclib::PetscLibType,x::PetscVec, nv::PetscInt, y::Vector{PetscVec}) 
Ends a split phase multiple dot product computation.

Input Parameters:
- `x`  - the first vector (can be `NULL`)
- `nv` - number of vectors
- `y`  - array of vectors (can be `NULL`)

Output Parameter:
- `result` - where the result will go

Level: advanced

-seealso: `VecMDotBegin()`, `VecNormBegin()`, `VecNormEnd()`, `VecNorm()`, `VecDot()`, `VecMDot()`,
`VecTDotBegin()`, `VecTDotEnd()`, `VecMTDotBegin()`, `VecMTDotEnd()`, `PetscCommSplitReductionBegin()`

# External Links
$(_doc_external("Vec/VecMDotEnd"))
"""
function VecMDotEnd(petsclib::PetscLibType, x::PetscVec, nv::PetscInt, y::Vector{PetscVec}) end

@for_petsc function VecMDotEnd(petsclib::$UnionPetscLib, x::PetscVec, nv::$PetscInt, y::Vector{PetscVec} )
	result = Vector{$PetscScalar}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:VecMDotEnd, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{CVec}, Ptr{$PetscScalar}),
               x, nv, y, result,
              )


	return result
end 

"""
	VecMTDotBegin(petsclib::PetscLibType,x::PetscVec, nv::PetscInt, y::Vector{PetscVec}, result::Vector{PetscScalar}) 
Starts a split phase transpose multiple dot product computation.

Input Parameters:
- `x`      - the first vector
- `nv`     - number of vectors
- `y`      - array of  vectors
- `result` - where the result will go (can be `NULL`)

Level: advanced

-seealso: `VecMTDotEnd()`, `VecNormBegin()`, `VecNormEnd()`, `VecNorm()`, `VecDot()`, `VecMDot()`,
`VecDotBegin()`, `VecDotEnd()`, `VecMDotBegin()`, `VecMDotEnd()`, `PetscCommSplitReductionBegin()`

# External Links
$(_doc_external("Vec/VecMTDotBegin"))
"""
function VecMTDotBegin(petsclib::PetscLibType, x::PetscVec, nv::PetscInt, y::Vector{PetscVec}, result::Vector{PetscScalar}) end

@for_petsc function VecMTDotBegin(petsclib::$UnionPetscLib, x::PetscVec, nv::$PetscInt, y::Vector{PetscVec}, result::Vector{$PetscScalar} )

    @chk ccall(
               (:VecMTDotBegin, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{CVec}, Ptr{$PetscScalar}),
               x, nv, y, result,
              )


	return nothing
end 

"""
	result::Vector{PetscScalar} = VecMTDotEnd(petsclib::PetscLibType,x::PetscVec, nv::PetscInt, y::Vector{PetscVec}) 
Ends a split phase transpose multiple dot product computation.

Input Parameters:
- `x`  - the first vector (can be `NULL`)
- `nv` - number of vectors
- `y`  - array of  vectors (can be `NULL`)

Output Parameter:
- `result` - where the result will go

Level: advanced

-seealso: `VecMTDotBegin()`, `VecNormBegin()`, `VecNormEnd()`, `VecNorm()`, `VecDot()`, `VecMDot()`,
`VecDotBegin()`, `VecDotEnd()`, `VecMDotBegin()`, `VecMDotEnd()`, `PetscCommSplitReductionBegin()`

# External Links
$(_doc_external("Vec/VecMTDotEnd"))
"""
function VecMTDotEnd(petsclib::PetscLibType, x::PetscVec, nv::PetscInt, y::Vector{PetscVec}) end

@for_petsc function VecMTDotEnd(petsclib::$UnionPetscLib, x::PetscVec, nv::$PetscInt, y::Vector{PetscVec} )
	result = Vector{$PetscScalar}(undef, nv);  

    @chk ccall(
               (:VecMTDotEnd, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscInt, Ptr{CVec}, Ptr{$PetscScalar}),
               x, nv, y, result,
              )


	return result
end 

"""
	VecExp(petsclib::PetscLibType,v::PetscVec) 
Replaces each component of a vector by e^x_i

Not Collective

Input Parameter:
- `v` - The vector

Output Parameter:
- `v` - The vector of exponents

Level: beginner

-seealso: `Vec`, `VecLog()`, `VecAbs()`, `VecSqrtAbs()`, `VecReciprocal()`


# External Links
$(_doc_external("Vec/VecExp"))
"""
function VecExp(petsclib::PetscLibType, v::PetscVec) end

@for_petsc function VecExp(petsclib::$UnionPetscLib, v::PetscVec )

    @chk ccall(
               (:VecExp, $petsc_library),
               PetscErrorCode,
               (CVec,),
               v,
              )


	return nothing
end 

"""
	VecLog(petsclib::PetscLibType,v::PetscVec) 
Replaces each component of a vector by log(x_i), the natural logarithm

Not Collective

Input Parameter:
- `v` - The vector

Output Parameter:
- `v` - The vector of logs

Level: beginner

-seealso: `Vec`, `VecExp()`, `VecAbs()`, `VecSqrtAbs()`, `VecReciprocal()`


# External Links
$(_doc_external("Vec/VecLog"))
"""
function VecLog(petsclib::PetscLibType, v::PetscVec) end

@for_petsc function VecLog(petsclib::$UnionPetscLib, v::PetscVec )

    @chk ccall(
               (:VecLog, $petsc_library),
               PetscErrorCode,
               (CVec,),
               v,
              )


	return nothing
end 

"""
	VecAbs(petsclib::PetscLibType,v::PetscVec) 
Replaces every element in a vector with its absolute value.

Logically Collective

Input Parameter:
- `v` - the vector

Level: intermediate

-seealso: `Vec`, `VecExp()`, `VecSqrtAbs()`, `VecReciprocal()`, `VecLog()`

# External Links
$(_doc_external("Vec/VecAbs"))
"""
function VecAbs(petsclib::PetscLibType, v::PetscVec) end

@for_petsc function VecAbs(petsclib::$UnionPetscLib, v::PetscVec )

    @chk ccall(
               (:VecAbs, $petsc_library),
               PetscErrorCode,
               (CVec,),
               v,
              )


	return nothing
end 

"""
	VecConjugate(petsclib::PetscLibType,x::PetscVec) 
Conjugates a vector. That is, replace every entry in a vector with its complex conjugate

Logically Collective

Input Parameter:
- `x` - the vector

Level: intermediate

-seealso: [](ch_vectors), `Vec`, `VecSet()`

# External Links
$(_doc_external("Vec/VecConjugate"))
"""
function VecConjugate(petsclib::PetscLibType, x::PetscVec) end

@for_petsc function VecConjugate(petsclib::$UnionPetscLib, x::PetscVec )

    @chk ccall(
               (:VecConjugate, $petsc_library),
               PetscErrorCode,
               (CVec,),
               x,
              )


	return nothing
end 

"""
	VecImaginaryPart(petsclib::PetscLibType,v::PetscVec) 
Replaces a complex vector with its imginary part

Collective

Input Parameter:
- `v` - the vector

Level: beginner

-seealso: `Vec`, `VecNorm()`, `VecRealPart()`

# External Links
$(_doc_external("Vec/VecImaginaryPart"))
"""
function VecImaginaryPart(petsclib::PetscLibType, v::PetscVec) end

@for_petsc function VecImaginaryPart(petsclib::$UnionPetscLib, v::PetscVec )

    @chk ccall(
               (:VecImaginaryPart, $petsc_library),
               PetscErrorCode,
               (CVec,),
               v,
              )


	return nothing
end 

"""
	VecRealPart(petsclib::PetscLibType,v::PetscVec) 
Replaces a complex vector with its real part

Collective

Input Parameter:
- `v` - the vector

Level: beginner

-seealso: `Vec`, `VecNorm()`, `VecImaginaryPart()`

# External Links
$(_doc_external("Vec/VecRealPart"))
"""
function VecRealPart(petsclib::PetscLibType, v::PetscVec) end

@for_petsc function VecRealPart(petsclib::$UnionPetscLib, v::PetscVec )

    @chk ccall(
               (:VecRealPart, $petsc_library),
               PetscErrorCode,
               (CVec,),
               v,
              )


	return nothing
end 

"""
	dp::PetscScalar,nm::PetscReal = VecDotNorm2(petsclib::PetscLibType,s::PetscVec, t::PetscVec) 
computes the inner product of two vectors and the 2

Collective

Input Parameters:
- `s` - first vector
- `t` - second vector

Output Parameters:
- `dp` - s'conj(t)
- `nm` - t'conj(t)

Level: advanced

-seealso: `Vec`, `VecDot()`, `VecNorm()`, `VecDotBegin()`, `VecNormBegin()`, `VecDotEnd()`, `VecNormEnd()`


# External Links
$(_doc_external("Vec/VecDotNorm2"))
"""
function VecDotNorm2(petsclib::PetscLibType, s::PetscVec, t::PetscVec) end

@for_petsc function VecDotNorm2(petsclib::$UnionPetscLib, s::PetscVec, t::PetscVec )
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
	mean::PetscScalar = VecMean(petsclib::PetscLibType,v::PetscVec) 
Computes the arithmetic mean of all the components of a vector.

Collective

Input Parameter:
- `v` - the vector

Output Parameter:
- `mean` - the result

Level: beginner

-seealso: `Vec`, `VecSum()`, `VecNorm()`

# External Links
$(_doc_external("Vec/VecMean"))
"""
function VecMean(petsclib::PetscLibType, v::PetscVec) end

@for_petsc function VecMean(petsclib::$UnionPetscLib, v::PetscVec )
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
	VecPermute(petsclib::PetscLibType,x::PetscVec, row::IS, inv::PetscBool) 
Permutes a vector in place using the given ordering.

Input Parameters:
- `x`   - The vector
- `row` - The ordering
- `inv` - The flag for inverting the permutation

Level: beginner

-seealso: `Vec`, `MatPermute()`

# External Links
$(_doc_external("Vec/VecPermute"))
"""
function VecPermute(petsclib::PetscLibType, x::PetscVec, row::IS, inv::PetscBool) end

@for_petsc function VecPermute(petsclib::$UnionPetscLib, x::PetscVec, row::IS, inv::PetscBool )

    @chk ccall(
               (:VecPermute, $petsc_library),
               PetscErrorCode,
               (CVec, CIS, PetscBool),
               x, row, inv,
              )


	return nothing
end 

"""
	flg::PetscBool = VecEqual(petsclib::PetscLibType,vec1::PetscVec, vec2::PetscVec) 
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

-seealso: `Vec`

# External Links
$(_doc_external("Vec/VecEqual"))
"""
function VecEqual(petsclib::PetscLibType, vec1::PetscVec, vec2::PetscVec) end

@for_petsc function VecEqual(petsclib::$UnionPetscLib, vec1::PetscVec, vec2::PetscVec )
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
	n::PetscInt,e::Vector{PetscScalar} = VecUniqueEntries(petsclib::PetscLibType,vec::PetscVec) 
Compute the number of unique entries, and those entries

Collective

Input Parameter:
- `vec` - the vector

Output Parameters:
- `n` - The number of unique entries
- `e` - The entries, each MPI process receives all the unique entries

Level: intermediate

-seealso: `Vec`

# External Links
$(_doc_external("Vec/VecUniqueEntries"))
"""
function VecUniqueEntries(petsclib::PetscLibType, vec::PetscVec) end

@for_petsc function VecUniqueEntries(petsclib::$UnionPetscLib, vec::PetscVec )
	n_ = Ref{$PetscInt}()
	e_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:VecUniqueEntries, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{$PetscInt}, Ptr{Ptr{$PetscScalar}}),
               vec, n_, e_,
              )

	n = n_[]
	e = unsafe_wrap(Array, e_[], VecGetLocalSize(petsclib, x); own = false)

	return n,e
end 

"""
	VecFilter(petsclib::PetscLibType,v::PetscVec, tol::PetscReal) 
Set all values in the vector with an absolute value less than or equal to the tolerance to zero

Input Parameters:
- `v`   - The vector
- `tol` - The zero tolerance

Output Parameter:
- `v` - The filtered vector

Level: intermediate

-seealso: `VecCreate()`, `VecSet()`, `MatFilter()`

# External Links
$(_doc_external("Vec/VecFilter"))
"""
function VecFilter(petsclib::PetscLibType, v::PetscVec, tol::PetscReal) end

@for_petsc function VecFilter(petsclib::$UnionPetscLib, v::PetscVec, tol::$PetscReal )

    @chk ccall(
               (:VecFilter, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscReal),
               v, tol,
              )


	return nothing
end 

"""
	VecWhichEqual(petsclib::PetscLibType,Vec1::PetscVec, Vec2::PetscVec, S::IS) 
Creates an index set containing the indices
where the vectors `Vec1` and `Vec2` have identical elements.

Collective

Input Parameters:
- `Vec1` - the first vector to compare
- `Vec2` - the second two vector to compare

Output Parameter:
- `S` - The index set containing the indices i where vec1[i] == vec2[i]

Level: advanced

-seealso: `Vec`

# External Links
$(_doc_external("Vec/VecWhichEqual"))
"""
function VecWhichEqual(petsclib::PetscLibType, Vec1::PetscVec, Vec2::PetscVec, S::IS) end

@for_petsc function VecWhichEqual(petsclib::$UnionPetscLib, Vec1::PetscVec, Vec2::PetscVec, S::IS )

    @chk ccall(
               (:VecWhichEqual, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, Ptr{IS}),
               Vec1, Vec2, S,
              )


	return nothing
end 

"""
	VecWhichLessThan(petsclib::PetscLibType,Vec1::PetscVec, Vec2::PetscVec, S::IS) 
Creates an index set containing the indices
where the vectors `Vec1` < `Vec2`

Collective

Input Parameters:
- `Vec1` - the first vector to compare
- `Vec2` - the second vector to compare

Output Parameter:
- `S` - The index set containing the indices i where vec1[i] < vec2[i]

Level: advanced

-seealso: `Vec`

# External Links
$(_doc_external("Vec/VecWhichLessThan"))
"""
function VecWhichLessThan(petsclib::PetscLibType, Vec1::PetscVec, Vec2::PetscVec, S::IS) end

@for_petsc function VecWhichLessThan(petsclib::$UnionPetscLib, Vec1::PetscVec, Vec2::PetscVec, S::IS )

    @chk ccall(
               (:VecWhichLessThan, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, Ptr{IS}),
               Vec1, Vec2, S,
              )


	return nothing
end 

"""
	VecWhichGreaterThan(petsclib::PetscLibType,Vec1::PetscVec, Vec2::PetscVec, S::IS) 
Creates an index set containing the indices
where the vectors `Vec1` > `Vec2`

Collective

Input Parameters:
- `Vec1` - the first vector to compare
- `Vec2` - the second vector to compare

Output Parameter:
- `S` - The index set containing the indices i where vec1[i] > vec2[i]

Level: advanced

-seealso: `Vec`

# External Links
$(_doc_external("Vec/VecWhichGreaterThan"))
"""
function VecWhichGreaterThan(petsclib::PetscLibType, Vec1::PetscVec, Vec2::PetscVec, S::IS) end

@for_petsc function VecWhichGreaterThan(petsclib::$UnionPetscLib, Vec1::PetscVec, Vec2::PetscVec, S::IS )

    @chk ccall(
               (:VecWhichGreaterThan, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, Ptr{IS}),
               Vec1, Vec2, S,
              )


	return nothing
end 

"""
	VecWhichBetween(petsclib::PetscLibType,VecLow::PetscVec, V::PetscVec, VecHigh::PetscVec, S::IS) 
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

-seealso: `Vec`

# External Links
$(_doc_external("Vec/VecWhichBetween"))
"""
function VecWhichBetween(petsclib::PetscLibType, VecLow::PetscVec, V::PetscVec, VecHigh::PetscVec, S::IS) end

@for_petsc function VecWhichBetween(petsclib::$UnionPetscLib, VecLow::PetscVec, V::PetscVec, VecHigh::PetscVec, S::IS )

    @chk ccall(
               (:VecWhichBetween, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, CVec, Ptr{IS}),
               VecLow, V, VecHigh, S,
              )


	return nothing
end 

"""
	VecWhichBetweenOrEqual(petsclib::PetscLibType,VecLow::PetscVec, V::PetscVec, VecHigh::PetscVec, S::IS) 
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

-seealso: `Vec`

# External Links
$(_doc_external("Vec/VecWhichBetweenOrEqual"))
"""
function VecWhichBetweenOrEqual(petsclib::PetscLibType, VecLow::PetscVec, V::PetscVec, VecHigh::PetscVec, S::IS) end

@for_petsc function VecWhichBetweenOrEqual(petsclib::$UnionPetscLib, VecLow::PetscVec, V::PetscVec, VecHigh::PetscVec, S::IS )

    @chk ccall(
               (:VecWhichBetweenOrEqual, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, CVec, Ptr{IS}),
               VecLow, V, VecHigh, S,
              )


	return nothing
end 

"""
	VecWhichInactive(petsclib::PetscLibType,VecLow::PetscVec, V::PetscVec, D::PetscVec, VecHigh::PetscVec, Strong::PetscBool, S::IS) 
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

-seealso: `Vec`

# External Links
$(_doc_external("Vec/VecWhichInactive"))
"""
function VecWhichInactive(petsclib::PetscLibType, VecLow::PetscVec, V::PetscVec, D::PetscVec, VecHigh::PetscVec, Strong::PetscBool, S::IS) end

@for_petsc function VecWhichInactive(petsclib::$UnionPetscLib, VecLow::PetscVec, V::PetscVec, D::PetscVec, VecHigh::PetscVec, Strong::PetscBool, S::IS )

    @chk ccall(
               (:VecWhichInactive, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, CVec, CVec, PetscBool, Ptr{IS}),
               VecLow, V, D, VecHigh, Strong, S,
              )


	return nothing
end 

"""
	VecISAXPY(petsclib::PetscLibType,vfull::PetscVec, is::IS, alpha::PetscScalar, vreduced::PetscVec) 
Adds a reduced vector to the appropriate elements of a full
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

-seealso: `VecISCopy()`, `VecISSet()`, `VecAXPY()`

# External Links
$(_doc_external("Vec/VecISAXPY"))
"""
function VecISAXPY(petsclib::PetscLibType, vfull::PetscVec, is::IS, alpha::PetscScalar, vreduced::PetscVec) end

@for_petsc function VecISAXPY(petsclib::$UnionPetscLib, vfull::PetscVec, is::IS, alpha::$PetscScalar, vreduced::PetscVec )

    @chk ccall(
               (:VecISAXPY, $petsc_library),
               PetscErrorCode,
               (CVec, CIS, $PetscScalar, CVec),
               vfull, is, alpha, vreduced,
              )


	return nothing
end 

"""
	VecISCopy(petsclib::PetscLibType,vfull::PetscVec, is::IS, mode::ScatterMode, vreduced::PetscVec) 
Copies between a reduced vector and the appropriate elements of a full

Logically Collective

Input Parameters:
- `vfull`    - the full-space vector
- `is`       - the index set for the reduced space
- `mode`     - the direction of copying, `SCATTER_FORWARD` or `SCATTER_REVERSE`
- `vreduced` - the reduced-space vector

Output Parameter:
- `vfull` - the sum of the full-space vector and reduced-space vector

Level: advanced

-seealso: `VecISSet()`, `VecISAXPY()`, `VecCopy()`

# External Links
$(_doc_external("Vec/VecISCopy"))
"""
function VecISCopy(petsclib::PetscLibType, vfull::PetscVec, is::IS, mode::ScatterMode, vreduced::PetscVec) end

@for_petsc function VecISCopy(petsclib::$UnionPetscLib, vfull::PetscVec, is::IS, mode::ScatterMode, vreduced::PetscVec )

    @chk ccall(
               (:VecISCopy, $petsc_library),
               PetscErrorCode,
               (CVec, CIS, ScatterMode, CVec),
               vfull, is, mode, vreduced,
              )


	return nothing
end 

"""
	VecISSet(petsclib::PetscLibType,V::PetscVec, S::IS, c::PetscScalar) 
Sets the elements of a vector, specified by an index set, to a constant

Logically Collective

Input Parameters:
- `V` - the vector
- `S` - index set for the locations in the vector
- `c` - the constant

Level: advanced

-seealso: `VecISCopy()`, `VecISAXPY()`, `VecISShift()`, `VecSet()`

# External Links
$(_doc_external("Vec/VecISSet"))
"""
function VecISSet(petsclib::PetscLibType, V::PetscVec, S::IS, c::PetscScalar) end

@for_petsc function VecISSet(petsclib::$UnionPetscLib, V::PetscVec, S::IS, c::$PetscScalar )

    @chk ccall(
               (:VecISSet, $petsc_library),
               PetscErrorCode,
               (CVec, CIS, $PetscScalar),
               V, S, c,
              )


	return nothing
end 

"""
	VecISShift(petsclib::PetscLibType,V::PetscVec, S::IS, c::PetscScalar) 
Shifts the elements of a vector, specified by an index set, by a constant

Logically Collective

Input Parameters:
- `V` - the vector
- `S` - index set for the locations in the vector
- `c` - the constant

Level: advanced

-seealso: `VecISCopy()`, `VecISAXPY()`, `VecISSet()`, `VecShift()`

# External Links
$(_doc_external("Vec/VecISShift"))
"""
function VecISShift(petsclib::PetscLibType, V::PetscVec, S::IS, c::PetscScalar) end

@for_petsc function VecISShift(petsclib::$UnionPetscLib, V::PetscVec, S::IS, c::$PetscScalar )

    @chk ccall(
               (:VecISShift, $petsc_library),
               PetscErrorCode,
               (CVec, CIS, $PetscScalar),
               V, S, c,
              )


	return nothing
end 

"""
	VecBoundGradientProjection(petsclib::PetscLibType,G::PetscVec, X::PetscVec, XL::PetscVec, XU::PetscVec, GP::PetscVec) 
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

-seealso: `Vec`

# External Links
$(_doc_external("Vec/VecBoundGradientProjection"))
"""
function VecBoundGradientProjection(petsclib::PetscLibType, G::PetscVec, X::PetscVec, XL::PetscVec, XU::PetscVec, GP::PetscVec) end

@for_petsc function VecBoundGradientProjection(petsclib::$UnionPetscLib, G::PetscVec, X::PetscVec, XL::PetscVec, XU::PetscVec, GP::PetscVec )

    @chk ccall(
               (:VecBoundGradientProjection, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, CVec, CVec, CVec),
               G, X, XL, XU, GP,
              )


	return nothing
end 

"""
	VecPow(petsclib::PetscLibType,v::PetscVec, p::PetscScalar) 
Replaces each component of a vector by  x_i^p 

Logically Collective

Input Parameters:
- `v` - the vector
- `p` - the exponent to use on each element

Level: intermediate

-seealso: `Vec`

# External Links
$(_doc_external("Vec/VecPow"))
"""
function VecPow(petsclib::PetscLibType, v::PetscVec, p::PetscScalar) end

@for_petsc function VecPow(petsclib::$UnionPetscLib, v::PetscVec, p::$PetscScalar )

    @chk ccall(
               (:VecPow, $petsc_library),
               PetscErrorCode,
               (CVec, $PetscScalar),
               v, p,
              )


	return nothing
end 

"""
	VecMedian(petsclib::PetscLibType,Vec1::PetscVec, Vec2::PetscVec, Vec3::PetscVec, VMedian::PetscVec) 
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

-seealso: `Vec`

# External Links
$(_doc_external("Vec/VecMedian"))
"""
function VecMedian(petsclib::PetscLibType, Vec1::PetscVec, Vec2::PetscVec, Vec3::PetscVec, VMedian::PetscVec) end

@for_petsc function VecMedian(petsclib::$UnionPetscLib, Vec1::PetscVec, Vec2::PetscVec, Vec3::PetscVec, VMedian::PetscVec )

    @chk ccall(
               (:VecMedian, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, CVec, CVec),
               Vec1, Vec2, Vec3, VMedian,
              )


	return nothing
end 

"""
	values::Vector{PetscScalar} = VecGetValuesSection(petsclib::PetscLibType,v::PetscVec, s::PetscSection, point::PetscInt) 
Gets all the values associated with a given point, according to the section, in the given `Vec`

Not Collective

Input Parameters:
- `v`     - the `Vec`
- `s`     - the organizing `PetscSection`
- `point` - the point

Output Parameter:
- `values` - the array of output values

Level: developer

-seealso: `PetscSection`, `PetscSectionCreate()`, `VecSetValuesSection()`

# External Links
$(_doc_external("Vec/VecGetValuesSection"))
"""
function VecGetValuesSection(petsclib::PetscLibType, v::PetscVec, s::PetscSection, point::PetscInt) end

@for_petsc function VecGetValuesSection(petsclib::$UnionPetscLib, v::PetscVec, s::PetscSection, point::$PetscInt )
	values_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:VecGetValuesSection, $petsc_library),
               PetscErrorCode,
               (CVec, PetscSection, $PetscInt, Ptr{Ptr{$PetscScalar}}),
               v, s, point, values_,
              )

	values = unsafe_wrap(Array, values_[], VecGetLocalSize(petsclib, x); own = false)

	return values
end 

"""
	VecGetDM(petsclib::PetscLibType,v::PetscVec, dm::PetscDM) 
Gets the `DM` defining the data layout of the vector

Not Collective

Input Parameter:
- `v` - The `Vec`

Output Parameter:
- `dm` - The `DM`

Level: intermediate

Note:
A `Vec` may not have a `DM` associated with it.

See also: 
=== 
`DM`, `VecSetDM()`, `DMGetLocalVector()`, `DMGetGlobalVector()`, `DMSetVecType()`

# External Links
$(_doc_external("DM/VecGetDM"))
"""
function VecGetDM(petsclib::PetscLibType, v::PetscVec, dm::PetscDM) end

@for_petsc function VecGetDM(petsclib::$UnionPetscLib, v::PetscVec, dm::PetscDM )
	dm_ = Ref(dm.ptr)

    @chk ccall(
               (:VecGetDM, $petsc_library),
               PetscErrorCode,
               (CVec, Ptr{CDM}),
               v, dm_,
              )

	dm.ptr = C_NULL

	return nothing
end 

"""
	VecFischer(petsclib::PetscLibType,X::PetscVec, F::PetscVec, L::PetscVec, U::PetscVec, FB::PetscVec) 
Evaluates the Fischer
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

-seealso: `Vec`, `VecSFischer()`, `MatDFischer()`, `MatDSFischer()`

# External Links
$(_doc_external("Tao/VecFischer"))
"""
function VecFischer(petsclib::PetscLibType, X::PetscVec, F::PetscVec, L::PetscVec, U::PetscVec, FB::PetscVec) end

@for_petsc function VecFischer(petsclib::$UnionPetscLib, X::PetscVec, F::PetscVec, L::PetscVec, U::PetscVec, FB::PetscVec )

    @chk ccall(
               (:VecFischer, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, CVec, CVec, CVec),
               X, F, L, U, FB,
              )


	return nothing
end 

"""
	xv::Vector{PetscScalar},yv::Vector{PetscScalar} = VecGetArrayPair(petsclib::PetscLibType,x::PetscVec, y::PetscVec) 

# External Links
$(_doc_external("Vec/VecGetArrayPair"))
"""
function VecGetArrayPair(petsclib::PetscLibType, x::PetscVec, y::PetscVec) end

@for_petsc function VecGetArrayPair(petsclib::$UnionPetscLib, x::PetscVec, y::PetscVec )
	xv_ = Ref{Ptr{$PetscScalar}}()
	yv_ = Ref{Ptr{$PetscScalar}}()

    @chk ccall(
               (:VecGetArrayPair, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, Ptr{Ptr{$PetscScalar}}, Ptr{Ptr{$PetscScalar}}),
               x, y, xv_, yv_,
              )

	xv = unsafe_wrap(Array, xv_[], VecGetLocalSize(petsclib, x); own = false)
	yv = unsafe_wrap(Array, yv_[], VecGetLocalSize(petsclib, x); own = false)

	return xv,yv
end 

"""
	VecRestoreArrayPair(petsclib::PetscLibType,x::PetscVec, y::PetscVec, xv::Vector{PetscScalar}, yv::Vector{PetscScalar}) 

# External Links
$(_doc_external("Vec/VecRestoreArrayPair"))
"""
function VecRestoreArrayPair(petsclib::PetscLibType, x::PetscVec, y::PetscVec, xv::Vector{PetscScalar}, yv::Vector{PetscScalar}) end

@for_petsc function VecRestoreArrayPair(petsclib::$UnionPetscLib, x::PetscVec, y::PetscVec, xv::Vector{$PetscScalar}, yv::Vector{$PetscScalar} )
	xv_ = Ref(pointer(xv))
	yv_ = Ref(pointer(yv))

    @chk ccall(
               (:VecRestoreArrayPair, $petsc_library),
               PetscErrorCode,
               (CVec, CVec, Ptr{Ptr{$PetscScalar}}, Ptr{Ptr{$PetscScalar}}),
               x, y, xv_, yv_,
              )


	return nothing
end 

