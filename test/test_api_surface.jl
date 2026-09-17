# The API surface is fully accounted for by the register (naming.md §1.1).
#
# This is `scripts/api_surface.jl --check`, run from the test suite so that a new
# high-level binding cannot land without an entry in `scripts/renames.jl`. The
# fix for a failure is to add the name to RENAMES, INTERNAL or UNCHANGED_PUBLIC
# there, never to loosen the check.

using Test
using PETSc

include(joinpath(@__DIR__, "..", "scripts", "api_surface.jl"))

@testset "api_surface --check" begin
    absent = APISurface.unregistered()
    if !isempty(absent)
        @info "bindings in no set of scripts/renames.jl" absent
    end
    @test isempty(absent)
    @test APISurface.check(devnull)
    # The surface is non-empty, so that a filter bug cannot make the check pass
    # by finding nothing at all.
    @test length(APISurface.surface()) > 100
end
