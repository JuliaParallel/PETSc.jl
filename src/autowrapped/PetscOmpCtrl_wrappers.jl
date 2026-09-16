"""
	PetscOmpCtrlBarrier(petsclib::PetscLibType, ctrl::PetscOmpCtrl) 
Do barrier on MPI ranks in omp_comm contained by the PETSc OMP controller (to let slave ranks free their CPU)

Input Parameter:
- `ctrl` - a PETSc OMP controller

See also: `PetscOmpCtrlOmpRegionOnMasterBegin()`, `PetscOmpCtrlOmpRegionOnMasterEnd()`, `PetscOmpCtrlCreate()`, `PetscOmpCtrlDestroy()`

# External Links
$(_doc_external("Sys/PetscOmpCtrlBarrier"))
"""
function PetscOmpCtrlBarrier(petsclib::PetscLibType, ctrl::PetscOmpCtrl)
    error("PetscOmpCtrlBarrier: no generated method for these argument types")
end

@for_petsc function PetscOmpCtrlBarrier(petsclib::$UnionPetscLib, ctrl::PetscOmpCtrl )

    @chk ccall(
               (:PetscOmpCtrlBarrier, $petsc_library),
               PetscErrorCode,
               (PetscOmpCtrl,),
               ctrl,
              )


	return nothing
end 

"""
	pctrl::PetscOmpCtrl = PetscOmpCtrlCreate(petsclib::PetscLibType, petsc_comm::MPI_Comm, nthreads::PetscInt) 
create a PETSc OpenMP controller, which manages PETSc's interaction with third party libraries that use OpenMP

Input Parameters:
- `petsc_comm` - a communicator some PETSc object (for example, a matrix) lives in
- `nthreads`   - number of threads per MPI rank to spawn in a library using OpenMP. If nthreads = -1, let PETSc decide a suitable value

Output Parameter:
- `pctrl` - a PETSc OpenMP controller

Level: developer

See also: `PetscOmpCtrlDestroy()`, `PetscOmpCtrlGetOmpComms()`, `PetscOmpCtrlBarrier()`, `PetscOmpCtrlOmpRegionOnMasterBegin()`, `PetscOmpCtrlOmpRegionOnMasterEnd()`

# External Links
$(_doc_external("Sys/PetscOmpCtrlCreate"))
"""
function PetscOmpCtrlCreate(petsclib::PetscLibType, petsc_comm::MPI_Comm, nthreads::Integer)
    error("PetscOmpCtrlCreate: no generated method for these argument types")
end

@for_petsc function PetscOmpCtrlCreate(petsclib::$UnionPetscLib, petsc_comm::MPI_Comm, nthreads::$PetscInt )
	pctrl_ = Ref{PetscOmpCtrl}()

    @chk ccall(
               (:PetscOmpCtrlCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, $PetscInt, Ptr{PetscOmpCtrl}),
               petsc_comm, nthreads, pctrl_,
              )

	pctrl = pctrl_[]

	return pctrl
end 

"""
	PetscOmpCtrlDestroy(petsclib::PetscLibType, pctrl::Union{PetscOmpCtrl, Ref{PetscOmpCtrl}}) 
destroy the PETSc OpenMP controller

Input Parameter:
- `pctrl` - a PETSc OpenMP controller

Level: developer

See also: `PetscOmpCtrlCreate()`, `PetscOmpCtrlGetOmpComms()`, `PetscOmpCtrlBarrier()`, `PetscOmpCtrlOmpRegionOnMasterBegin()`, `PetscOmpCtrlOmpRegionOnMasterEnd()`

# External Links
$(_doc_external("Sys/PetscOmpCtrlDestroy"))
"""
function PetscOmpCtrlDestroy(petsclib::PetscLibType, pctrl::Union{PetscOmpCtrl, Ref{PetscOmpCtrl}})
    error("PetscOmpCtrlDestroy: no generated method for these argument types")
end

@for_petsc function PetscOmpCtrlDestroy(petsclib::$UnionPetscLib, pctrl::Union{PetscOmpCtrl, Ref{PetscOmpCtrl}} )
	pctrl_ = pctrl isa Base.RefValue ? pctrl : Ref{PetscOmpCtrl}(pctrl)

    @chk ccall(
               (:PetscOmpCtrlDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscOmpCtrl},),
               pctrl_,
              )


	return nothing
end 

"""
	omp_comm::MPI_Comm,omp_master_comm::MPI_Comm,is_omp_master::PetscBool = PetscOmpCtrlGetOmpComms(petsclib::PetscLibType, ctrl::PetscOmpCtrl) 
Get MPI communicators from a PETSc OMP controller

Input Parameter:
- `ctrl` - a PETSc OMP controller

Output Parameters:
- `omp_comm`        - a communicator that includes a master rank and slave ranks where master spawns threads
- `omp_master_comm` - on master ranks, return a communicator that include master ranks of each omp_comm;
on slave ranks, `MPI_COMM_NULL` will be return in reality.
- `is_omp_master`   - true if the calling process is an OMP master rank.

See also: `PetscOmpCtrlCreate()`, `PetscOmpCtrlDestroy()`, `PetscOmpCtrlBarrier()`, `PetscOmpCtrlOmpRegionOnMasterBegin()`, `PetscOmpCtrlOmpRegionOnMasterEnd()`

# External Links
$(_doc_external("Sys/PetscOmpCtrlGetOmpComms"))
"""
function PetscOmpCtrlGetOmpComms(petsclib::PetscLibType, ctrl::PetscOmpCtrl)
    error("PetscOmpCtrlGetOmpComms: no generated method for these argument types")
end

@for_petsc function PetscOmpCtrlGetOmpComms(petsclib::$UnionPetscLib, ctrl::PetscOmpCtrl )
	omp_comm_ = Ref{MPI.MPI_Comm}()
	omp_master_comm_ = Ref{MPI.MPI_Comm}()
	is_omp_master_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscOmpCtrlGetOmpComms, $petsc_library),
               PetscErrorCode,
               (PetscOmpCtrl, Ptr{MPI.MPI_Comm}, Ptr{MPI.MPI_Comm}, Ptr{PetscBool}),
               ctrl, omp_comm_, omp_master_comm_, is_omp_master_,
              )

	omp_comm = MPI.Comm(omp_comm_[])
	omp_master_comm = MPI.Comm(omp_master_comm_[])
	is_omp_master = is_omp_master_[]

	return omp_comm,omp_master_comm,is_omp_master
end 

"""
	PetscOmpCtrlOmpRegionOnMasterBegin(petsclib::PetscLibType, ctrl::PetscOmpCtrl) 
Mark the beginning of an OpenMP library call on master ranks

Input Parameter:
- `ctrl` - a PETSc OMP controller

See also: `PetscOmpCtrlOmpRegionOnMasterEnd()`, `PetscOmpCtrlCreate()`, `PetscOmpCtrlDestroy()`, `PetscOmpCtrlBarrier()`

# External Links
$(_doc_external("Sys/PetscOmpCtrlOmpRegionOnMasterBegin"))
"""
function PetscOmpCtrlOmpRegionOnMasterBegin(petsclib::PetscLibType, ctrl::PetscOmpCtrl)
    error("PetscOmpCtrlOmpRegionOnMasterBegin: no generated method for these argument types")
end

@for_petsc function PetscOmpCtrlOmpRegionOnMasterBegin(petsclib::$UnionPetscLib, ctrl::PetscOmpCtrl )

    @chk ccall(
               (:PetscOmpCtrlOmpRegionOnMasterBegin, $petsc_library),
               PetscErrorCode,
               (PetscOmpCtrl,),
               ctrl,
              )


	return nothing
end 

"""
	PetscOmpCtrlOmpRegionOnMasterEnd(petsclib::PetscLibType, ctrl::PetscOmpCtrl) 
Mark the end of an OpenMP library call on master ranks

Input Parameter:
- `ctrl` - a PETSc OMP controller

See also: `PetscOmpCtrlOmpRegionOnMasterBegin()`, `PetscOmpCtrlCreate()`, `PetscOmpCtrlDestroy()`, `PetscOmpCtrlBarrier()`

# External Links
$(_doc_external("Sys/PetscOmpCtrlOmpRegionOnMasterEnd"))
"""
function PetscOmpCtrlOmpRegionOnMasterEnd(petsclib::PetscLibType, ctrl::PetscOmpCtrl)
    error("PetscOmpCtrlOmpRegionOnMasterEnd: no generated method for these argument types")
end

@for_petsc function PetscOmpCtrlOmpRegionOnMasterEnd(petsclib::$UnionPetscLib, ctrl::PetscOmpCtrl )

    @chk ccall(
               (:PetscOmpCtrlOmpRegionOnMasterEnd, $petsc_library),
               PetscErrorCode,
               (PetscOmpCtrl,),
               ctrl,
              )


	return nothing
end 

