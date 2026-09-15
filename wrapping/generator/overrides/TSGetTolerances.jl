# override for TSGetTolerances; C signature: TSGetTolerances(TS ts, PetscReal* atol, Vec* vatol, PetscReal* rtol, Vec* vrtol)
"""
	atol::PetscReal,vatol::PetscVec,rtol::PetscReal,vrtol::PetscVec = TSGetTolerances(petsclib::PetscLibType,ts::AbstractTS) 
Get tolerances for local truncation error when using adaptive controller

Logically Collective

Input Parameter:
- `ts` - time integration context

Output Parameters:
- `atol`  - scalar absolute tolerances, `NULL` to ignore
- `vatol` - vector of absolute tolerances, `NULL` to ignore
- `rtol`  - scalar relative tolerances, `NULL` to ignore
- `vrtol` - vector of relative tolerances, `NULL` to ignore

Level: beginner

-seealso: [](ch_ts), `TS`, `TSAdapt`, `TSErrorWeightedNorm()`, `TSSetTolerances()`

# External Links
$(_doc_external("Ts/TSGetTolerances"))
"""
function TSGetTolerances(petsclib::PetscLibType, ts::AbstractTS) end

@for_petsc function TSGetTolerances(petsclib::$UnionPetscLib, ts::AbstractTS)
	atol_ = Ref{$PetscReal}()
	vatol_ = Ref{CVec}(C_NULL)
	rtol_ = Ref{$PetscReal}()
	vrtol_ = Ref{CVec}(C_NULL)

    @chk ccall(
               (:TSGetTolerances, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscReal}, Ptr{CVec}, Ptr{$PetscReal}, Ptr{CVec}),
               ts, atol_, vatol_, rtol_, vrtol_,
              )

	# The per-component vectors belong to the TS, so they get no finalizer.
	# They come back NULL when only scalar tolerances are set.
	vatol = PetscVec(vatol_[], petsclib)
	vrtol = PetscVec(vrtol_[], petsclib)

	return atol_[],vatol,rtol_[],vrtol
end

