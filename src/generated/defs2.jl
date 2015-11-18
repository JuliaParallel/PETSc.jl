# hacky definitions that Clang did not pick up
typealias PetscCopyMode Int32
const PETSC_COPY_VALUES = Int32(0)

function ISCreateGeneral(arg1::MPI_Comm,arg2::Integer,arg3::Union{Ptr{Int64},StridedArray{Int64},Ptr{Int64},Ref{Int64}}, arg3a::Integer,arg4::Union{Ptr{IS{Float64}},StridedArray{IS{Float64}},Ptr{IS{Float64}},Ref{IS{Float64}}})
  err = ccall((:ISCreateGeneral,petscRealDouble),PetscErrorCode,(comm_type,Int64,Ptr{Int64},Cint,Ptr{IS{Float64}}),arg1,arg2,arg3,arg3a,arg4)
  return err
end
