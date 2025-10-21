# define common PETSc constants
# this excludes configurable constants (e.g. PetscScalar) which are set in lib.jl

const PetscErrorCode = Cint

struct PetscError <: Exception
    code::PetscErrorCode
end

macro chk(expr)
    :((errcode = $(esc(expr))) == 0 || throw(PetscError(errcode)))
end