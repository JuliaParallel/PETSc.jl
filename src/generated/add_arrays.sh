typenames=(Int32, Int64, Bool, Float32, Float64, PetscInt, PetscScalar, PetscReal, PetscBool)

for val in "${typenames[@]}";
do
  echo $val
  sed -i "/^function/s/Ptr{$val}/Union(Ptr{$val}, AbstractArray{$val}, Ptr{Void})/g" ./PETSc.jl
done
