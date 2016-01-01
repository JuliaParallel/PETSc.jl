# used to modify function signatures
type_dict = Dict{Any, Any}(
  :PetscScalar => :Complex128,
  :PetscReal => :Float64,
  :PetscInt => :Int64,
  :PetscLogDouble => :Float64,
  :_PetscCDIntNd => :Void,
  :_PetscCDArrNd => :Void,
)
