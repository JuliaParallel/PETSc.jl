# define common PETSc constants
# this excludes configurable constants (e.g. PetscScalar) which are set in lib.jl

const PetscErrorCode = Cint

"""
    PetscError(code)

A PETSc call returned the nonzero error code `code`. `func` is the C function
the failing call went to (`Symbol("")` when not known), and `message` is PETSc's
text for the error: its generic description followed by the specific message
set where the error started, deepest in the call chain.

# External Links
$(_doc_external("Sys/PetscErrorMessage"))
"""
struct PetscError <: Exception
    code::PetscErrorCode
    func::Symbol
    message::String
end
PetscError(code) = PetscError(code, Symbol(""), "")

function Base.showerror(io::IO, err::PetscError)
    print(io, "PetscError(", err.code, ")")
    isempty(String(err.func)) || print(io, " from ", err.func)
    isempty(err.message)      || print(io, ": ", err.message)
    return nothing
end

# PETSc's text for error `code`, read right after the failing call: PETSc keeps
# the specific message only until the next error. `library` is what the failing
# `ccall` named. Building an exception must not throw, so a library that can't
# be reached gives an empty message.
@noinline function petsc_error(code, func::Symbol, library)
    message = try
        petsc_error_message(code, library)
    catch
        ""
    end
    return PetscError(code, func, message)
end

function petsc_error_message(code, library::AbstractString)
    # The library is already loaded, so this only looks up its handle
    handle = Libdl.dlopen(library)
    fptr   = Libdl.dlsym(handle, :PetscErrorMessage)
    text, specific = Ref{Ptr{Cchar}}(C_NULL), Ref{Ptr{Cchar}}(C_NULL)
    ccall(fptr, PetscErrorCode, (PetscErrorCode, Ptr{Ptr{Cchar}}, Ptr{Ptr{Cchar}}),
          code, text, specific) == 0 || return ""
    parts = String[unsafe_string(p[]) for p in (text, specific) if p[] != C_NULL]
    return join(filter(!isempty, strip.(parts)), ". ")
end
petsc_error_message(_, _) = ""

# The function name and library of `ccall((:Name, library), ...)`, or nothing
# when the call goes through a function pointer
function ccall_target(expr)
    Meta.isexpr(expr, :call)    && expr.args[1] === :ccall  || return nothing
    target = expr.args[2]
    Meta.isexpr(target, :tuple) && length(target.args) == 2 || return nothing
    name = target.args[1]
    name isa QuoteNode          && name.value isa Symbol    || return nothing
    return name.value, target.args[2]
end

# Throw a `PetscError` when the PETSc call `expr` returns a nonzero code
macro chk(expr)
    target = ccall_target(expr)
    if target === nothing
        return :((errcode = $(esc(expr))) == 0 || throw(PetscError(errcode)))
    end
    name, library = target
    return :((errcode = $(esc(expr))) == 0 ||
        throw(petsc_error(errcode, $(QuoteNode(name)), $(esc(library)))))
end
