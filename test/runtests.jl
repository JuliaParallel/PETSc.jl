using PETSc
using FactCheck

# determine scalar type of current run
global ST = Float64  # scalar type

function RC(x::Number)
# used to test real, complex
  if ST == Float64
    return Float64(real(x))
  elseif ST == Float32
    return Float32(real(x))
  else  # scalar_type == 3
    return complex(x)
  end
end

function RC(x::AbstractArray)
# used to test real, complex
  if ST == Float64
    tmp = similar(x, ST)
    for i=1:length(x)
      tmp[i] = Float64(real(x[i]))
    end
    return tmp
  elseif ST == Float32
    tmp = similar(x, ST)
    for i=1:length(x)
      tmp[i] = ST(real(x[i]))
    end
    return tmp
  else  # scalar_type == 3
    return x
  end
end

for ST in PETSc.C.petsc_type
  println("\n\nTesting ", ST)
  include("error.jl")
  include("vec.jl")
  include("mat.jl")
  include("ksp.jl")
  include("is.jl")
end

facts("\ntesting options") do
    OPTIONS["foo"]=true
    for ST in PETSc.C.petsc_type
        @fact haskey(OPTIONS[ST], :foo) --> true
        @fact OPTIONS[ST][:foo] --> "true"
    end
    OPTIONS["foo"]=nothing
    for ST in PETSc.C.petsc_type
        @fact haskey(OPTIONS[ST], :foo) --> false
        OPTIONS[ST]["bar"] = 17
        @fact pop!(OPTIONS[ST], "bar") --> "17"
        @fact pop!(OPTIONS[ST], "bar", 23) --> 23
        @fact haskey(OPTIONS[ST], :bar) --> false
        withoptions(ST, "baz"=>"aloha") do
            @fact OPTIONS[ST][:baz] --> "aloha"
        end
        @fact haskey(OPTIONS[ST], :baz) --> false
    end
end
