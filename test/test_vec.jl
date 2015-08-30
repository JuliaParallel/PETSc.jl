
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

vec3 = similar(vec, Float64, 5)
len3_ret = length(vec3)
@fact len3_ret => 5

vec4 = copy(vec)

for i=1:length(vec)
  @fact vec4[i] => roughly(vec[i])
end

idx = [1,3, 4]
vec4[idx] = 2
vals_ret = vec4[idx]
for i=1:length(idx)
  @fact vals_ret[i] => 2
end

fill!(vec4, 3)

for i=1:length(vec4)
  @fact vec4[i] => roughly(3)
end

vec4[1:2] = 4

@fact vec4[1:2] => [4.0, 4.0]


vals = [1, 3., 4]
vec4[idx] = vals

for i=1:length(idx)
  @fact vec4[idx[i]] => vals[i]
end

logicals = Array(Bool, length(vec4))
for i=1:length(logicals)
  logicals[i] = false
end
logicals[2] = true

vec4[logicals] = 5

@fact vec4[2] => roughly(5)
@fact vec4[1] => not(5)

vals = rand(1)
vec4[logicals] = vals
println("vals = ", vals)
println("logicals = ", logicals)
println("vec4 = ", vec4)
@fact vec4[2] => roughly(vals[1])
@fact vec4[1] => not(vals[1])



end
