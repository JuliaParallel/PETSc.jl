# used to modify function signatures
type_dict = Dict{Any, Any}(
  :PetscScalar => :Float64,
  :PetscReal => :Float64,
  :PetscInt => :Int64,
  :PetscLogDouble => :Float64,
  :_PetscCDIntNd => :Void,
  :_PetscCDArrNd => :Void,
  :NLF_DAAD => :Void,
)
