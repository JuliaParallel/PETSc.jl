using PETSc
if VERSION >= v"0.5.0-dev+7720"
  using Base.Test
else
  using BaseTestNext
end

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
  # @testset "Scalar type $ST" begin # uncomment when nested test results can be printed
    include("error.jl")
    include("ksp.jl")
    include("vec.jl")
    include("is.jl")
    include("mat.jl")
  # end
end

@test PETSc.petsc_sizeof(PETSc.C.PETSC_BOOL) == 4

@testset "Options" begin
    OPTIONS["foo"]=true
    for ST in PETSc.C.petsc_type
        @test haskey(OPTIONS[ST], :foo)
        @test OPTIONS[ST][:foo] == "true"
    end
    OPTIONS["foo"]=nothing
    for ST in PETSc.C.petsc_type
        @test !haskey(OPTIONS[ST], :foo)
        OPTIONS[ST]["bar"] = 17
        @test pop!(OPTIONS[ST], "bar") == "17"
        @test pop!(OPTIONS[ST], "bar", 23) == 23
        @test !haskey(OPTIONS[ST], :bar)
        withoptions(ST, "baz"=>"aloha") do
            @test OPTIONS[ST][:baz] == "aloha"
        end
        @test !haskey(OPTIONS[ST], :baz)
        OPTIONS[ST]["whee"] = 12
        delete!(OPTIONS[ST], "whee", 13)
        @test !haskey(OPTIONS[ST], :whee)
        @test isempty(similar(OPTIONS[ST]))
    end
end
