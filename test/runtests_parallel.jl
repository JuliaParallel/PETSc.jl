# run tests in parallel
include("runtests_setup.jl")

comm = MPI.COMM_WORLD
comm_size = MPI.Comm_size(MPI.COMM_WORLD)
comm_rank = MPI.Comm_rank(MPI.COMM_WORLD)



# size of the system owned by this process (ie. local sys size)
sys_size = PETSc.C.PetscInt(3)

# create a 3x3 block for each process
# create these with smallest precision, so they can be promoted
tmp3 = convert(Array{Complex64, 1}, [1.0 + 0im; 2.0 + 1.0im; 3.0 + 2.0im])
tmp4 = convert( Array{Complex64, 2}, [1.0 + 1im   2 + 2im  3 + 3im; 4 + 4im  5 + 5im 7 + 7im; 7 + 7im 8 + 8im 9 + 9im])




for ST in PETSc.C.petsc_type
  # @testset "Scalar type $ST" begin # uncomment when nested test results can be printed

  # Have both local and global versions of system for testing
  global A_julia = zeros(ST, sys_size, sys_size)
  global rhs = zeros(ST, sys_size)
  global A_julia_global = zeros(ST, comm_size*sys_size, comm_size*sys_size)
  global rhs_global = zeros(ST, comm_size*sys_size)

  # create right hand side
  for i=1:sys_size
    rhs[i] = convert(ST, RC(tmp3[i]))

    for j=0:(comm_size - 1)  # global 
      idx = i + j*sys_size
      rhs_global[idx] = convert(ST, RC(tmp3[i]))
    end    


  end

  # create matrix A
  for i=1:sys_size
    for j=1:sys_size
      A_julia[i,j] = convert(ST, RC(tmp4[i,j]))

      for k=0:(comm_size - 1)
        idxm = i + k*sys_size
        idxn = j + k*sys_size
        A_julia_global[idxm, idxn] = convert(ST, RC(tmp4[i,j]) )
      end
 
    end
  end

  if comm_rank == 0
    println("rhs_global = ", rhs_global)
    println("A_julia_global = ", A_julia_global)
  end

  #println("rhs = ", rhs)
  #println("A_julia = ", A_julia)
  x_julia = A_julia\rhs
  #println("x_julia = ", x_julia)



  include("vecp.jl")
  include("c_funcs.jl")
  # end
end


