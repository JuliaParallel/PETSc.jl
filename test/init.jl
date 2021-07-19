using Test
using PETSc

@testset "init" begin
    for petsclib in PETSc.petsclibs
        # The first time through finalize should be false since we have never
        # initialized petsclib yet...
        initial_finalized_value = false

        # since we haven't called anything these should be false!
        @test !(PETSc.initialized(petsclib))
        @test PETSc.finalized(petsclib) == initial_finalized_value

        # The second time through  time through finalize should be true since
        # we have initialized petsclib yet...
        initial_finalized_value = true

        # initialize PETSc
        PETSc.initialize(petsclib)

        # Check values again
        @test PETSc.initialized(petsclib)
        @test !(PETSc.finalized(petsclib))

        PETSc.finalize(petsclib)

        # Check values again
        @test !(PETSc.initialized(petsclib))
        @test PETSc.finalized(petsclib)
    end
end

