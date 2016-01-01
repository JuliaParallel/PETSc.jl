# used to modify function signatures
type_dict = Dict{Any, Any}(
  :PetscScalar => :Float32,
  :PetscReal => :Float32,
  :PetscInt => :Int64,
  :PetscLogDouble => :Float64,
  :_PetscCDIntNd => :Void,
  :_PetscCDArrNd => :Void,
)
