
# create Vec
facts("\n --- Testing Vector Function ---") do
vtype = PETSc.C.VECMPI
vec = PETSc.Vec(ST, vtype)
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

println("point 1")

vt = complex(2.,2)  # use vt to hold temporary values
vec[1] = RC(vt)
val_ret = vec[1]
@fact val_ret => RC(vt)

vec2 = similar(vec)
PETSc.AssemblyBegin(vec2)
PETSc.AssemblyEnd(vec2)
val2_ret = vec2[1]

@fact val2_ret => not(val_ret)

vec3 = similar(vec, ST, 5)
len3_ret = length(vec3)
@fact len3_ret => 5

println("\n\n\n Copying vec to vec4")
vec4 = copy(vec)
println("vec = ", vec)
println("vec4 = ", vec4)

println("\n\n\n")
for i=1:length(vec)
  @fact vec4[i] => roughly(vec[i])
end

println("point 2")

idx = [1,3, 4]
vt = RC(complex(2.,2))
println("idx = ",idx)
println("typeof(idx) = ", typeof(idx))
println("size(vec4) = ", size(vec4))
vec4[idx] = vt
println("set vec4 values")
println("idx = ", idx)
println("size(vec4) = ", size(vec4))
vals_ret = vec4[idx]
println("retrieved vec4 values")
println("vals_ret = ", vals_ret)
println("point 2.4")
println("idx = ", idx)
println("typeof(idx) = ", typeof(idx))
println("length(idx) = ", length(idx))
for i=1:length(idx)
  println("i = ", i)
#  println("vals_ret[i] = ", vals_ret[i])
  @fact vals_ret[i] => vt
end

println("point 2.5")
vt = RC(complex(3.,3))
fill!(vec4, vt)

for i=1:length(vec4)
  @fact vec4[i] => roughly(vt)
end

vt = RC(complex( 4.,4))
vec4[1:2] = vt


println("point 3")
@fact vec4[1:2] => [vt, vt]


vals = [RC(complex(1,1.)), RC(complex(3.,3)), RC(complex(4., 3))]
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

vt = RC(complex(5,5.))
vec4[logicals] = vt

@fact vec4[2] => roughly(vt)
@fact vec4[1] => not(vt)

vt = RC(complex(rand(), rand()))
vals = [vt]
vec4[logicals] = vals
println("vals = ", vals)
println("logicals = ", logicals)
println("vec4 = ", vec4)
@fact vec4[2] => roughly(vals[1])
@fact vec4[1] => not(vals[1])

# reset vec4
vec4_j = zeros(ST, length(vec4))
for i=1:length(vec4)
  vec4[i] = RC(complex(Float64(-i), Float64(-i))) 
  vec4_j[i] = RC(complex(Float64(-i), Float64(-i))) 

end

println("testing math functions")

println("testing abs")
vec4_j = abs(vec4_j)
abs!(vec4)

for i=1:length(vec4)
  @fact vec4[i] => vec4_j[i]
end

println("testing exp")

println("before, vec and vec4 = ")
println("vec4 = ", vec4)
println("vec4_j = ", vec4_j)

vec4_j = exp(vec4_j)
exp!(vec4)

println("vec4 = ", vec4)
println("vec4_j = ", vec4_j)


for i=1:length(vec4)
  @fact vec4[i] => roughly(vec4_j[i], atol=1e-4)
end

println("testing log")
vec4_j = log(vec4_j)
log!(vec4)

for i=1:length(vec4)
  @fact vec4[i] => roughly(vec4_j[i], atol=1e-4)
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

vecj = zeros(ST, length(vec))
vec2j = zeros(ST, length(vec))
vec4j = zeros(ST, length(vec))

for i=1:length(vec)
  vecj[i] = vec[i]
  vec2j[i] = vec2[i]
  vec4j[i] = vec4[i]
end

println("testing axpy")
vt = RC(complex(2.,2))
axpy!(vt, vec, vec2)
vec2j = vt*vecj + vec2j
println("vec2j = ", vec2j)
println("vec2 = ", vec2)
for i=1:length(vec)
  println("vec2j[i] = ", vec2j[i], ", vec2[i] = ", vec2[i])
  @fact vec2j[i] => vec2[i]
end

println("testing 4 argument axpy")
axpy!(vt, vec, vec2, vec4)
vec4j = vt*vecj + vec2j 

for i=1:length(vec)
  @fact vec2j[i] => vec2[i]
end


println("testing aypx")
aypx!(vec, vt, vec2)
vec2j = vt*vec2j + vec

for i=1:length(vec)
  @fact vec2j[i] => vec2[i]
end

println("testing axpby")
println("before operation:")
println("vec = ", vec)
println("vecj = ", vecj)
println("vec2 = ", vec2)
println("vec2j = ", vec2j)

vt2 = RC(complex(3.,3))
vt3 = RC(complex(4.,4))
axpby!(vt, vec, vt2, vec2)
vec2j = vt*vecj + vt2*vec2j
println("after operation:")
println("vec = ", vec)
println("vecj = ", vecj)
println("vec2 = ", vec2)
println("vec2j = ", vec2j)
for i=1:length(vec)
  @fact vec2j[i] => vec2[i]
end


axpbypcz!(vt, vec, vt2, vec2, vt3, vec4)
vec4j = vt*vecj + vt2*vec2j + vt3*vec4j

for i=1:length(vec)
  @fact vec4j[i] => vec4[i]
end

vecs = Array(typeof(vec), 2)
vecs[1] = vec
vecs[2] = vec2
#vecs = [vec; vec2]
alphas = [vt2, vt3]
println("vecs = ", vecs)
println("typeof(vecs) = ", typeof(vecs))

PETSc.maxpy!(vec4, alphas, vecs)
vec4j = vec4j + vt2*vecj + vt3*vec2j
println("vec4 = ", vec4)
println("vec4j = ", vec4j)
for i=1:length(vec)
  @fact vec4j[i] => vec4[i]
end


vec5 = Vec(ST, 3, PETSc.C.VECMPI)
vec6 = similar(vec5)
vec5j = zeros(ST, 3)
vec6j = zeros(ST, 3)

for i=1:3
  i_float = Float64(i)
  
  vec5[i] = RC(complex(i_float, i_float))
  vec6[i] = RC(complex(i_float+3, i_float+3))
  vec5j[i] = RC(complex(i_float, i_float))
  vec6j[i] = RC(complex(i_float +3, i_float+3))
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
