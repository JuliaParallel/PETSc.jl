
# we need to parametrise by type so we know which library to call
abstract type Viewer{T} end

const CPetscViewer         = Ptr{Cvoid}

Base.cconvert(::Type{CPetscViewer}, obj::Viewer) = obj.ptr
Base.unsafe_convert(::Type{Ptr{CPetscViewer}}, obj::Viewer) = 
    convert(Ptr{CPetscViewer}, pointer_from_objref(obj))


mutable struct ViewerStdout{T} <: Viewer{T}
    ptr::CPetscViewer
    comm::MPI.Comm
end

@for_libpetsc begin
    function ViewerStdout{$PetscScalar}(comm::MPI.Comm=MPI.COMM_SELF)
        ptr = ccall((:PETSC_VIEWER_STDOUT_, $libpetsc), CPetscViewer, (MPI.MPI_Comm,), comm)
        return ViewerStdout{$PetscScalar}(ptr, comm)
    end
    function Base.push!(viewer::Viewer{$PetscScalar}, format::PetscViewerFormat)
        @chk ccall((:PetscViewerPushFormat, $libpetsc), PetscErrorCode, (CPetscViewer,PetscViewerFormat), viewer, format)
        return nothing
    end
    function Base.pop!(viewer::Viewer{$PetscScalar})
        @chk ccall((:PetscViewerPopFormat, $libpetsc), PetscErrorCode, (CPetscViewer,), viewer)
        return nothing
    end

end

function with(f, viewer::Viewer, format::PetscViewerFormat)
    push!(viewer, format)
    try
        f()
    finally
        pop!(viewer)
    end
end

# ideally we would capture the output directly, but this looks difficult
# easiest option is to redirect stdout
# based on suggestion from https://github.com/JuliaLang/julia/issues/32567
function _show(io::IO, obj)
    old_stdout = stdout
    try
        rd, = redirect_stdout()
        view(obj)
        Libc.flush_cstdio()
        flush(stdout)
        write(io, readavailable(rd))
    finally
        redirect_stdout(old_stdout)
    end
    return nothing
end



#=
# PETSc_jll isn't built with X support
mutable struct ViewerDraw <: Viewer
    ptr::CPetscViewer
    comm::MPI.Comm
end
function ViewerDraw(comm::MPI.Comm)
    ptr = ccall((:PETSC_VIEWER_DRAW_, libpetsc), CPetscViewer, (MPI.MPI_Comm,), comm)
    return ViewerDraw(ptr, comm)
end
=#