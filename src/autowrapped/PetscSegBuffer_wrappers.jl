# autodefined type arguments for class ------
mutable struct _n_PetscSegBuffer end
const PetscSegBuffer = Ptr{_n_PetscSegBuffer}
# -------------------------------------------------------
"""
	seg::PetscSegBuffer = PetscSegBufferCreate(petsclib::PetscLibType,unitbytes::Csize_t, expected::PetscCount) 
create a segmented buffer

Not Collective, No Fortran Support

Input Parameters:
- `unitbytes` - number of bytes that each entry will contain
- `expected`  - expected/typical number of entries

Output Parameter:
- `seg` - `PetscSegBuffer` object

Level: developer

-seealso: `PetscSegBufferGet()`, `PetscSegBufferExtractAlloc()`, `PetscSegBufferExtractTo()`, `PetscSegBufferExtractInPlace()`, `PetscSegBufferDestroy()`,
`PetscSegBuffer`

# External Links
$(_doc_external("Sys/PetscSegBufferCreate"))
"""
function PetscSegBufferCreate(petsclib::PetscLibType, unitbytes::Csize_t, expected::PetscCount) end

@for_petsc function PetscSegBufferCreate(petsclib::$UnionPetscLib, unitbytes::Csize_t, expected::PetscCount )
	seg_ = Ref{PetscSegBuffer}()

    @chk ccall(
               (:PetscSegBufferCreate, $petsc_library),
               PetscErrorCode,
               (Csize_t, PetscCount, Ptr{PetscSegBuffer}),
               unitbytes, expected, seg_,
              )

	seg = seg_[]

	return seg
end 

"""
	PetscSegBufferGet(petsclib::PetscLibType,seg::PetscSegBuffer, count::PetscCount, buf::Cvoid) 
get new buffer space from a segmented buffer

Not Collective, No Fortran Support

Input Parameters:
- `seg`   - `PetscSegBuffer` buffer
- `count` - number of entries needed

Output Parameter:
- `buf` - address of new buffer for contiguous data

Level: developer

-seealso: `PetscSegBufferCreate()`, `PetscSegBufferExtractAlloc()`, `PetscSegBufferExtractTo()`, `PetscSegBufferExtractInPlace()`, `PetscSegBufferDestroy()`,
`PetscSegBuffer`, `PetscSegBufferGetInts()`

# External Links
$(_doc_external("Sys/PetscSegBufferGet"))
"""
function PetscSegBufferGet(petsclib::PetscLibType, seg::PetscSegBuffer, count::PetscCount, buf::Cvoid) end

@for_petsc function PetscSegBufferGet(petsclib::$UnionPetscLib, seg::PetscSegBuffer, count::PetscCount, buf::Cvoid )

    @chk ccall(
               (:PetscSegBufferGet, $petsc_library),
               PetscErrorCode,
               (PetscSegBuffer, PetscCount, Ptr{Cvoid}),
               seg, count, buf,
              )


	return nothing
end 

"""
	PetscSegBufferDestroy(petsclib::PetscLibType,seg::PetscSegBuffer) 
destroy segmented buffer

Not Collective, No Fortran Support

Input Parameter:
- `seg` - address of segmented buffer object

Level: developer

-seealso: `PetscSegBuffer`, `PetscSegBufferCreate()`

# External Links
$(_doc_external("Sys/PetscSegBufferDestroy"))
"""
function PetscSegBufferDestroy(petsclib::PetscLibType, seg::PetscSegBuffer) end

@for_petsc function PetscSegBufferDestroy(petsclib::$UnionPetscLib, seg::PetscSegBuffer )

    @chk ccall(
               (:PetscSegBufferDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscSegBuffer},),
               seg,
              )


	return nothing
end 

"""
	PetscSegBufferExtractTo(petsclib::PetscLibType,seg::PetscSegBuffer, contig::Cvoid) 
extract contiguous data to provided buffer and reset segmented buffer

Not Collective, No Fortran Support

Input Parameters:
- `seg`    - segmented buffer
- `contig` - allocated buffer to hold contiguous data

Level: developer

-seealso: `PetscSegBufferCreate()`, `PetscSegBufferGet()`, `PetscSegBufferDestroy()`, `PetscSegBufferExtractAlloc()`, `PetscSegBufferExtractInPlace()`,
`PetscSegBuffer`

# External Links
$(_doc_external("Sys/PetscSegBufferExtractTo"))
"""
function PetscSegBufferExtractTo(petsclib::PetscLibType, seg::PetscSegBuffer, contig::Cvoid) end

@for_petsc function PetscSegBufferExtractTo(petsclib::$UnionPetscLib, seg::PetscSegBuffer, contig::Cvoid )

    @chk ccall(
               (:PetscSegBufferExtractTo, $petsc_library),
               PetscErrorCode,
               (PetscSegBuffer, Ptr{Cvoid}),
               seg, contig,
              )


	return nothing
end 

"""
	PetscSegBufferExtractAlloc(petsclib::PetscLibType,seg::PetscSegBuffer, contiguous::Cvoid) 
extract contiguous data to new allocation and reset segmented buffer

Not Collective, No Fortran Support

Input Parameter:
- `seg` - `PetscSegBuffer` buffer

Output Parameter:
- `contiguous` - address of new array containing contiguous data, caller frees with `PetscFree()`

Level: developer

-seealso: `PetscSegBufferCreate()`, `PetscSegBufferGet()`, `PetscSegBufferDestroy()`, `PetscSegBufferExtractTo()`, `PetscSegBufferExtractInPlace()`,
`PetscSegBuffer`

# External Links
$(_doc_external("Sys/PetscSegBufferExtractAlloc"))
"""
function PetscSegBufferExtractAlloc(petsclib::PetscLibType, seg::PetscSegBuffer, contiguous::Cvoid) end

@for_petsc function PetscSegBufferExtractAlloc(petsclib::$UnionPetscLib, seg::PetscSegBuffer, contiguous::Cvoid )

    @chk ccall(
               (:PetscSegBufferExtractAlloc, $petsc_library),
               PetscErrorCode,
               (PetscSegBuffer, Ptr{Cvoid}),
               seg, contiguous,
              )


	return nothing
end 

"""
	PetscSegBufferExtractInPlace(petsclib::PetscLibType,seg::PetscSegBuffer, contig::Cvoid) 
extract in

Not Collective, No Fortran Support

Input Parameter:
- `seg` - `PetscSegBuffer` object

Output Parameter:
- `contig` - address of pointer to contiguous memory, may be `NULL`

Level: developer

-seealso: `PetscSegBuffer`, `PetscSegBufferExtractAlloc()`, `PetscSegBufferExtractTo()`

# External Links
$(_doc_external("Sys/PetscSegBufferExtractInPlace"))
"""
function PetscSegBufferExtractInPlace(petsclib::PetscLibType, seg::PetscSegBuffer, contig::Cvoid) end

@for_petsc function PetscSegBufferExtractInPlace(petsclib::$UnionPetscLib, seg::PetscSegBuffer, contig::Cvoid )

    @chk ccall(
               (:PetscSegBufferExtractInPlace, $petsc_library),
               PetscErrorCode,
               (PetscSegBuffer, Ptr{Cvoid}),
               seg, contig,
              )


	return nothing
end 

"""
	PetscSegBufferGetSize(petsclib::PetscLibType,seg::PetscSegBuffer, usedsize::PetscCount) 
get currently used number of entries of a `PetscSegBuffer`

Not Collective, No Fortran Support

Input Parameter:
- `seg` - `PetscSegBuffer` object

Output Parameter:
- `usedsize` - number of used units

Level: developer

-seealso: `PetscSegBuffer`, `PetscSegBufferExtractAlloc()`, `PetscSegBufferExtractTo()`, `PetscSegBufferCreate()`, `PetscSegBufferGet()`

# External Links
$(_doc_external("Sys/PetscSegBufferGetSize"))
"""
function PetscSegBufferGetSize(petsclib::PetscLibType, seg::PetscSegBuffer, usedsize::PetscCount) end

@for_petsc function PetscSegBufferGetSize(petsclib::$UnionPetscLib, seg::PetscSegBuffer, usedsize::PetscCount )

    @chk ccall(
               (:PetscSegBufferGetSize, $petsc_library),
               PetscErrorCode,
               (PetscSegBuffer, Ptr{PetscCount}),
               seg, usedsize,
              )


	return nothing
end 

"""
	PetscSegBufferUnuse(petsclib::PetscLibType,seg::PetscSegBuffer, unused::PetscCount) 
return some unused entries obtained with an overzealous `PetscSegBufferGet()`

Not Collective, No Fortran Support

Input Parameters:
- `seg`    - `PetscSegBuffer` object
- `unused` - number of unused units to return

Level: developer

-seealso: `PetscSegBuffer`, `PetscSegBufferCreate()`, `PetscSegBufferGet()`

# External Links
$(_doc_external("Sys/PetscSegBufferUnuse"))
"""
function PetscSegBufferUnuse(petsclib::PetscLibType, seg::PetscSegBuffer, unused::PetscCount) end

@for_petsc function PetscSegBufferUnuse(petsclib::$UnionPetscLib, seg::PetscSegBuffer, unused::PetscCount )

    @chk ccall(
               (:PetscSegBufferUnuse, $petsc_library),
               PetscErrorCode,
               (PetscSegBuffer, PetscCount),
               seg, unused,
              )


	return nothing
end 

"""
	PetscSegBufferGetInts(petsclib::PetscLibType,seg::PetscSegBuffer, count::PetscCount, slot::PetscInt) 

# External Links
$(_doc_external("Sys/PetscSegBufferGetInts"))
"""
function PetscSegBufferGetInts(petsclib::PetscLibType, seg::PetscSegBuffer, count::PetscCount, slot::PetscInt) end

@for_petsc function PetscSegBufferGetInts(petsclib::$UnionPetscLib, seg::PetscSegBuffer, count::PetscCount, slot::$PetscInt )

    @chk ccall(
               (:PetscSegBufferGetInts, $petsc_library),
               PetscErrorCode,
               (PetscSegBuffer, PetscCount, $PetscInt),
               seg, count, slot,
              )


	return nothing
end 

