include("runtests_setup.jl")

for ST in PETSc.C.petsc_type
  # @testset "Scalar type $ST" begin # uncomment when nested test results can be printed
#    include("error.jl")
#    include("ksp.jl")
#    include("vec.jl")
#    include("is.jl")
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
