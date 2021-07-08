const CDM = Ptr{Cvoid}

abstract type AbstractDM{PetscLib} end
mutable struct DM{PetscLib} <: AbstractDM{PetscLib}
    ptr::CDM
    opts::Options{PetscLib}
    DM{PetscLib}(ptr, opts = Options(PetscLib)) where {PetscLib} =
        new{PetscLib}(ptr, opts)
end

"""
    DMSetUp!(da::DM)

see [PETSc manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DM/DMSetUp.html)
"""
function DMSetUp! end

@for_petsc function DMSetUp!(da::DM{$PetscLib})
    with(da.opts) do
        @chk ccall(
            (:DMSetFromOptions, $petsc_library),
            PetscErrorCode,
            (CDM,),
            da,
        )

        @chk ccall((:DMSetUp, $petsc_library), PetscErrorCode, (CDM,), da)
    end
end

@for_petsc begin
    function destroy(da::DM{$PetscLib})
        finalized($PetscScalar) || @chk ccall(
            (:DMDestroy, $petsc_library),
            PetscErrorCode,
            (Ptr{CDM},),
            da,
        )
        da.ptr = C_NULL
        return nothing
    end
end
