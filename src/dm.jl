const CDM = Ptr{Cvoid}

abstract type AbstractDM{T, PetscLib} end
mutable struct DM{T, PetscLib} <: AbstractDM{T, PetscLib}
    ptr::CDM
    opts::Options{PetscLib}
    DM{T, PetscLib}(
        ptr,
        opts = Options(PetscLib),
    ) where {T, PetscLib} = new{T, PetscLib}(ptr, opts)
end

"""
    DMSetUp!(da::DM)

see [PETSc manual](https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DM/DMSetUp.html)
"""
function DMSetUp! end

@for_libpetsc begin
    function destroy(da::DM{$PetscScalar})
        finalized($PetscScalar) ||
            @chk ccall((:DMDestroy, $libpetsc), PetscErrorCode, (Ptr{CDM},), da)
        da.ptr = C_NULL
        return nothing
    end

    function DMSetUp!(da::DM{$PetscScalar})
        with(da.opts) do
            @chk ccall(
                (:DMSetFromOptions, $libpetsc),
                PetscErrorCode,
                (CDM,),
                da,
            )

            @chk ccall((:DMSetUp, $libpetsc), PetscErrorCode, (CDM,), da)
        end
    end
end
