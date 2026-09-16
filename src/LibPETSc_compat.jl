# ============================================================================
#   Version-dependent bindings
# ============================================================================
#
# src/autowrapped is generated from a single PETSc release, but the 0.4.x line
# binds both 3.22.x and 3.25.x. The few calls whose C signature differs between
# them ask the library which release it is instead of trusting the generated
# signature. Everything else is ABI-identical across the two.

"""
    petsc_version(petsclib) -> VersionNumber

The release of the PETSc library behind `petsclib`, as the library reports it.
"""
function petsc_version(petsclib)
    major, minor, subminor, _release = PetscGetVersionNumber(petsclib)
    return VersionNumber(Int(major), Int(minor), Int(subminor))
end

# PETSc 3.25 inserted a KSPDMActive mask before the flag of KSPSetDMActive.
# KSP_DMACTIVE_ALL covers operator, right-hand side and initial guess, 
# which is what the older two-argument call meant.
const KSP_DMACTIVE_ALL_MASK = Cint(1 + 2 + 4)

# The first release whose signatures differ from the generated ones.
const PETSC_SIGNATURE_BREAK = v"3.25"
