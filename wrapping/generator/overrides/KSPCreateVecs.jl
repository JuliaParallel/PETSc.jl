# override for KSPCreateVecs; C signature: KSPCreateVecs(KSP ksp, PetscInt rightn, Vec* right[], PetscInt leftn, Vec* left[])
"""
	right::Vector{PetscVec},left::Vector{PetscVec} = KSPCreateVecs(petsclib::PetscLibType,ksp::AbstractPetscKSP, rightn::PetscInt, leftn::PetscInt) 
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

-seealso: [](ch_ksp), `MatCreateVecs()`, `VecDestroyVecs()`, `KSPSetWorkVecs()`

# External Links
$(_doc_external("KSP/KSPCreateVecs"))
"""
function KSPCreateVecs(petsclib::PetscLibType, ksp::AbstractPetscKSP, rightn::PetscInt, leftn::PetscInt) end

@for_petsc function KSPCreateVecs(petsclib::$UnionPetscLib, ksp::AbstractPetscKSP, rightn::$PetscInt, leftn::$PetscInt )
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

