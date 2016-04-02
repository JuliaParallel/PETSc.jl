
@testset "Parallel Vec{$ST}" begin

  comm_rank = MPI.Comm_rank(MPI.COMM_WORLD)
  comm_size = MPI.Comm_size(MPI.COMM_WORLD)
  r = comm_rank:(comm_rank+2)
  r2 = (comm_rank+1):(comm_rank+3)
  vec1j = ST[r;]
  vec1j2 = ST[r2;]
  vec1 = PETSc.Vec(vec1j, comm=MPI.COMM_WORLD)

  @test length(vec1) == 3*comm_size
  @test vec1 == vec1j
  @test vec1 != vec1j2

  # ghost vectors
  # figure out global indices of ghost points
  nlocal = 4  # number of local values
  start_idx = 1 + comm_rank*nlocal
  end_idx = (comm_rank + 1)*nlocal
  local_range = start_idx:end_idx
  ghost_idx = [end_idx+1, end_idx+2]
  if comm_rank == (comm_size-1)  # last process
    ghost_idx = [1, 2]
  end

  v = VecGhost(Float64, nlocal, ghost_idx)
  # populate local values
  vlocal = VecLocal(v)
  for i=1:nlocal
    vlocal[i] = local_range[i]
  end

  restore(vlocal)
  ghost_update!(v)

  # check the ghost values
  vlocal = VecLocal(v)
  if comm_rank == (comm_size - 1)
    @test vlocal[nlocal + 1] == 1
    @test vlocal[nlocal + 2] == 2
  else
    @test vlocal[nlocal + 1] == end_idx + 1
    @test vlocal[nlocal + 2] == end_idx + 2
  end

  scatter!(v)
  if comm_rank == (comm_size - 1)
    @test vlocal[nlocal + 1] == 1
    @test vlocal[nlocal + 2] == 2
  else
    @test vlocal[nlocal + 1] == end_idx + 1
    @test vlocal[nlocal + 2] == end_idx + 2
  end


  vec = Vec([1.0, 2.0, 3.0])
  @test length(vec) == comm_size*3
  @test lengthlocal(vec) == 3

  vec_arr = LocalArray(vec)
  @test length(vec_arr) == 3
  LocalArrayRestore(vec_arr)

  vec_arr = LocalArrayRead(vec)
  @test length(vec_arr) == 3
  LocalArrayRestore(vec_arr)

  @testset "Application Ordering{$ST}" begin

    # create an application ordering that reverses each segment a
    # vector
    vec1 = Vec(ST, vtype=PETSc.C.VECMPI, mlocal=3, comm=MPI.COMM_WORLD)
    local_range = localpart(vec1)

    start_idx = local_range[1]
    end_idx = local_range[end]
    app_idx = collect(PetscInt, end_idx:-1:start_idx)
    petsc_idx = collect(PetscInt, local_range)

    ao = AO(ST, app_idx, petsc_idx)
    
    # transform back
    idx_arr = copy(petsc_idx)
    map_petsc_to_app!(ao, idx_arr)
    @test idx_arr == app_idx

    # transform back
    map_app_to_petsc!(ao, idx_arr)
    @test idx_arr == petsc_idx

    # do the same for index sets
    is_petsc = IS(ST, petsc_idx, comm=MPI.COMM_WORLD)
    is_app = IS(ST, app_idx, comm=MPI.COMM_WORLD)
    ao2 = AO(is_app, is_petsc)

    working_is = copy(is_petsc)
    map_petsc_to_app!(ao2, working_is)
    @test working_is == is_app

    map_app_to_petsc!(ao2, working_is)
    @test working_is == is_petsc

  end
end 
    
