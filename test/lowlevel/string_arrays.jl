using Test
using Libdl
using PETSc, MPI
using PETSc: LibPETSc

MPI.Initialized() || MPI.Init()

# Wrappers whose C argument is a `char ***` return a Vector{String}.

@testset "string-array out-parameters" begin
    for petsclib in PETSc.petsclibs[1:min(2, length(PETSc.petsclibs))]
        PETSc.initialize(petsclib)

        @testset "PetscStrToArray" begin
            @test LibPETSc.PetscStrToArray(petsclib, "a,bc,,def", Cchar(',')) == ["a", "bc", "def"]   # PETSc skips empty pieces
            @test LibPETSc.PetscStrToArray(petsclib, "one", Cchar(' ')) == ["one"]
            # repeated calls free what they allocate
            for _ in 1:100
                LibPETSc.PetscStrToArray(petsclib, "x y z", Cchar(' '))
            end
        end

        @testset "PetscGetArgs and PetscGetArguments" begin
            argc, args = LibPETSc.PetscGetArgs(petsclib)
            @test argc == length(args)
            @test all(a -> a isa String, args)
            rest = LibPETSc.PetscGetArguments(petsclib)
            @test rest isa Vector{String}
            @test rest == args[2:end]
        end

        @testset "PetscFunctionListGet" begin
            fl = Ref{LibPETSc.PetscFunctionList}(C_NULL)
            noop = @cfunction(() -> nothing, Cvoid, ())   # a NULL function is not added
            for name in ("beta", "alpha", "gamma")
                # PetscFunctionListAdd is a C macro over this function
                add = Libdl.dlsym(Libdl.dlopen(petsclib.petsc_library), :PetscFunctionListAdd_Private)
                err = ccall(add, Cint, (Ptr{LibPETSc.PetscFunctionList}, Cstring, Ptr{Cvoid}), fl, name, noop)
                @test err == 0
            end
            @test sort(LibPETSc.PetscFunctionListGet(petsclib, fl[])) == ["alpha", "beta", "gamma"]
            @test LibPETSc.PetscFunctionListGet(petsclib, fl[]) isa Vector{String}
            LibPETSc.PetscFunctionListDestroy(petsclib, fl)
        end
    end
end
