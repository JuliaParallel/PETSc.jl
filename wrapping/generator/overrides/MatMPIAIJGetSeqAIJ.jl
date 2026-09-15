# override for MatMPIAIJGetSeqAIJ; C signature: MatMPIAIJGetSeqAIJ(Mat A, Mat* Ad, Mat* Ao, PetscInt* colmap[])
"""
	Ad::PetscMat,Ao::PetscMat,colmap::Vector{PetscInt} = MatMPIAIJGetSeqAIJ(petsclib::PetscLibType,A::AbstractPetscMat) 
Returns the local pieces of this distributed matrix

Not Collective

Input Parameter:
- `A` - The `MATMPIAIJ` matrix

Output Parameters:
- `Ad`     - The local diagonal block as a `MATSEQAIJ` matrix
- `Ao`     - The local off-diagonal block as a `MATSEQAIJ` matrix
- `colmap` - An array mapping local column numbers of `Ao` to global column numbers of the parallel matrix

Level: intermediate

-seealso: [](ch_matrices), `Mat`, `MATMPIAIJ`, `MatMPIAIJGetLocalMat()`, `MatMPIAIJGetLocalMatCondensed()`, `MatCreateAIJ()`, `MATSEQAIJ`

# External Links
$(_doc_external("Mat/MatMPIAIJGetSeqAIJ"))
"""
function MatMPIAIJGetSeqAIJ(petsclib::PetscLibType, A::AbstractPetscMat) end

# `Ad`, `Ao` and `colmap` are owned by `A` (do not destroy them); `colmap` has one entry per
# column of `Ao`
@for_petsc function MatMPIAIJGetSeqAIJ(petsclib::$UnionPetscLib, A::AbstractPetscMat )
	Ad_ = Ref{CMat}()
	Ao_ = Ref{CMat}()
	colmap_ = Ref{Ptr{$PetscInt}}(C_NULL)

    @chk ccall(
               (:MatMPIAIJGetSeqAIJ, $petsc_library),
               PetscErrorCode,
               (CMat, Ptr{CMat}, Ptr{CMat}, Ptr{Ptr{$PetscInt}}),
               A, Ad_, Ao_, colmap_,
              )

	Ad = PetscMat(Ad_[], petsclib)
	Ao = PetscMat(Ao_[], petsclib)
	_, ncols_o = MatGetLocalSize(petsclib, Ao)
	colmap = colmap_[] == C_NULL ? $PetscInt[] : unsafe_wrap(Array, colmap_[], Int(ncols_o); own = false)

	return Ad,Ao,colmap
end

