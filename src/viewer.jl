
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