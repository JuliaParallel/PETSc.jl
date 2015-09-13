
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


println("testing logical indexing")
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

# reset vec4
 vec4_j = zeros(length(vec4))
for i=1:length(vec4)
  vec4[i] = -i
  vec4_j[i] = -i
end

println("testing math functions")

vec4_j = abs(vec4_j)
abs!(vec4)

for i=1:length(vec4)
  @fact vec4[i] => vec4_j[i]
end

vec4_j = exp(vec4_j)
exp!(vec4)

for i=1:length(vec4)
  @fact vec4[i] => vec4_j[i]
end


vec4_j = log(vec4_j)
log!(vec4)

for i=1:length(vec4)
  @fact vec4[i] => vec4_j[i]
end

#=
(max_val, max_index) = findmax(vec4)
max_val_j, max_index_j = findmax(vec4_j)

@fact max_val => roughly(maxval_j)
@fact max_index => max_index_j

(min_val, min_index) = findmin(vec4)
min_val_j, min_index_j = findmin(vec4_j)

@fact min_val => roughly(min_val_j)
@fact min_index => min_index_j
=#

#=
val = norm(vec4, 1)
val_j = norm(vec4_j, 1)

@fact val => val_j

val = norm(vec4, 2)
val_j = norm(vec4_j, 2)

@fact val => val_j

val = norm(vec4, Inf)
val_j = norm(vec4_j, Inf)

@fact val => val_j
=#
#=
normalize!(vec4)
vec4_j = vec4_j/norm(vec4_j, 2)

for i=1:length(vec4)
  @fact vec4[i] => vec4_j[i]
end
=#

println("testing dot product")

val = dot(vec4, vec)
#val_j = vec4.'*vec
val_j = dot(vec4, vec)
println("val = ", val)
println("val_j = ", val_j)

@fact val => val_j

# make copies of vecs 1 2 4

println("testing level 1 Blas")

vecj = zeros(length(vec))
vec2j = zeros(length(vec))
vec4j = zeros(length(vec))

for i=1:length(vec)
  vecj[i] = vec[i]
  vec2j[i] = vec2[i]
  vec4j[i] = vec4[i]
end

axpy!(2, vec, vec2)
vec2j = 2*vecj + vec2j

for i=1:length(vec)
  @fact vec2j[i] => vec2[i]
end

axpy!(2, vec, vec2, vec4)
vec4j = 2*vecj + vec2j 

for i=1:length(vec)
  @fact vec2j[i] => vec2[i]
end


aypx!(vec, 2, vec2)
vec2j = 2*vec2j + vec

for i=1:length(vec)
  @fact vec2j[i] => vec2[i]
end


axpby!(2, vec, 3, vec2)
vec2j = 2*vecj + 3*vec2j

for i=1:length(vec)
  @fact vec2j[i] => vec2[i]
end


axpbypcz!(2, vec, 3, vec2, 4, vec4)
vec4j = 2*vecj + 3*vec2j + 4*vec4j

for i=1:length(vec)
  @fact vec4j[i] => vec4[i]
end

vecs = Array(typeof(vec), 2)
vecs[1] = vec
vecs[2] = vec2
#vecs = [vec; vec2]
alphas = [2.0, 3.0]
println("vecs = ", vecs)
println("typeof(vecs) = ", typeof(vecs))

PETSc.maxpy!(vec4, alphas, vecs)
vec4j = vec4j + 2.0*vecj + 3.0*vec2j
println("vec4 = ", vec4)
println("vec4j = ", vec4j)
for i=1:length(vec)
  @fact vec4j[i] => vec4[i]
end


vec5 = Vec(Float64, 3, PETSc.C.VECMPI)
vec6 = similar(vec5)
vec5j = zeros(3)
vec6j = zeros(3)

for i=1:3
  vec5[i] = i
  vec6[i] = i+3
  vec5j[i] = i
  vec6j[i] = i+3
end

println("vec5 = ", vec5)
println("vec6 = ", vec6)

vec7 = vec5.*vec6
vec7j = vec5j.*vec6j
println("vec7j = ", vec7j)
println("vec7 = ", vec7)
for i=1:3
  @fact vec7[i] => roughly(vec7j[i])
end

vec8 = vec5./vec6
vec8j = vec5j./vec6j

for i=1:3
  @fact vec8[i] => roughly(vec8j[i])
end

vec9 = vec5.^3
vec9j = vec5j.^3

for i=1:3
  @fact vec9[i] => roughly(vec9j[i])
end

vec10 = vec5 + vec6
vec10j = vec5j + vec6j
for i=1:3
  @fact vec10[i] => roughly(vec10j[i])
end

vec11 = vec5 - vec6
vec11j = vec5j - vec6j

for i=1:3
  @fact vec11[i] => roughly(vec11j[i])
end








end
