#
# A few utility functions to help wrapping PETSc


export PETSc_unsafe_wrap, PETSc_RefPtr

MPI_Comm = MPI.Comm


"""
    PETSc_unsafe_wrap(r_array, dims::NTuple{N,Int}; own=false)

`unsafe_wrap` applied to a PETSc array which can change in size from 1D to 4D
"""
function PETSc_unsafe_wrap(r_array, dims::NTuple{N,IT}; own=false) where {N,IT}

    if N==1
        array = unsafe_wrap(Array, r_array[], dims; own = own)
    elseif N==2
        array = unsafe_wrap(Array, unsafe_load(r_array[]), dims; own = own)
    elseif N==3
        array = unsafe_wrap(Array, unsafe_load(unsafe_load(r_array[])), dims; own = own)
    elseif N==4
        array = unsafe_wrap(Array, unsafe_load(unsafe_load(unsafe_load(r_array[]))), dims; own = own)
    else
        error("not implemented for N=$N")
    end

    return array
end

"""
    PETSc_RefPtr(N::Int, PetscType::DataType=Float64)

Creates a reference to a pointer for different dimensions `N` 
"""
function PETSc_RefPtr(N::IT, PetscType::DataType=Float64) where {IT}

    if N==1
        r_array = Ref{Ptr{PetscType}}()
    elseif N==2
        r_array = Ref{Ptr{Ptr{PetscType}}}()
    elseif N==3
        r_array = Ref{Ptr{Ptr{Ptr{PetscType}}}}()
    elseif N==4
        r_array = Ref{Ptr{Ptr{Ptr{Ptr{PetscType}}}}}()
    else
        error("not implemented for N=$N")
    end

    return r_array
end


"""
    PETSc_RefPtr(dims::NTuple{N,Int}, PetscType::DataType=Float64) 

Creates a reference to a pointer of a PETSc array of dimensions `dims`
"""
PETSc_RefPtr(dims::NTuple{N,IT}, PetscType::DataType=Float64) where {N,IT} = PETSc_RefPtr(N, PetscType)