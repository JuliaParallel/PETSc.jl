
# create Vec
facts("\n --- Testing Vector Function ---") do
vtype = PETSc.C.VECMPI
vec = PETSc.Vec(Float64, vtype)
#PETSc.settype!(vec, vtype)
PETSc.setsizes!(vec, 4)
len_ret = length(vec)

@fact len_ret => 4
size_ret = size(vec)
@fact size_ret => (4,) 
len_local = lengthlocal(vec)
@fact len_local => 4
size_local = sizelocal(vec)
@fact size_local => (4,)

vectype_ret = PETSc.gettype(vec)

@fact vectype_ret => PETSc.C.VECMPI



vec[1] = 2
val_ret = vec[1]
@fact val_ret => 2

vec2 = similar(vec)
PETSc.AssemblyBegin(vec2)
PETSc.AssemblyEnd(vec2)
val2_ret = vec2[1]

@fact val2_ret => not(val_ret)

#vec3 = similar(vec, Float64, 5)
#len3_ret = length(vec3)
#@fact len3_ret => 5


end
