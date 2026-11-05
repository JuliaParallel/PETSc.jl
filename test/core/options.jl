using Test
using PETSc
using SparseArrays: sparse

isdefined(Main, :PETScTestUtils) || include(joinpath(@__DIR__, "..", "testutils.jl"))
using .PETScTestUtils: without_petsc_traceback

@testset "options" begin
    kw_opts = (
        ksp_monitor = nothing,
        ksp_view = true,
        false_opt = false,
        da_grid_x = 100,
        da_grid_y = 100,
        pc_type = "mg",
        pc_mg_levels = 1,
        mg_levels_0_pc_type = "ilu",
    )
    for petsclib in PETSc.petsclibs
        #petsclib = PETSc.petsclibs[1]
        PETSc.initialize(petsclib)

        opts = PETSc.PetscOptions(petsclib; kw_opts...)
        
        # Check that all the keys got added
        for (key, val) in pairs(kw_opts)
            key = string(key)
            if val === false
                @test_throws KeyError opts[key]
            else
                val = isnothing(val) || val === true ? "" : string(val)
                @test val == opts[key]
            end
        end

        # Check that we throw when a bad key is asked for
        @test_throws KeyError opts["bad key"]

        # try to insert some new keys
        opts["new_opt"] = 1
        opts["nothing_opt"] = nothing
        @test "" == opts["nothing_opt"]
        @test "1" == opts["new_opt"]
        
        # Check that viewer is working
        _stdout = stdout
        (rd, wr) = redirect_stdout()
        @show opts
        #@test readline(rd) == "PETSc Options database:"
        @test readline(rd) == "#PETSc Option Table entries:"
        @test readline(rd) == "-da_grid_x 100 # (source: code)"
        @test readline(rd) == "-da_grid_y 100 # (source: code)"
        @test readline(rd) == "-ksp_monitor # (source: code)"
        @test readline(rd) == "-ksp_view # (source: code)"
        @test readline(rd) == "-mg_levels_0_pc_type ilu # (source: code)"
        @test readline(rd) == "-new_opt 1 # (source: code)"
        @test readline(rd) == "-nothing_opt # (source: code)"
        @test readline(rd) == "-pc_mg_levels 1 # (source: code)"
        @test readline(rd) == "-pc_type mg # (source: code)"
        @test readline(rd) == "#End of PETSc Option Table entries"

        redirect_stdout(_stdout)
        

        PETSc.finalize(petsclib)
    end
end

@testset "setting an option to false clears it" begin
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
    opts = PETSc.PetscOptions(petsclib)
    has(key) = PETSc.LibPETSc.PetscOptionsHasName(petsclib, opts, "", key)

    opts[:x] = true
    @test has("-x")
    opts[:x] = false
    @test !has("-x")
    opts[:pc_type] = "lu"
    opts[:pc_type] = false
    @test !has("-pc_type")
    # a key never set is a no-op
    opts[:never_set] = false
    @test !has("-never_set")

    # a viewer option PETSc can't parse as a value (`-ksp_monitor false` is an error)
    PetscScalar = petsclib.PetscScalar
    S = sparse(PetscScalar[2 0; 0 4])
    ksp = PETSc.KSP(petsclib, PETSc.MPI.COMM_SELF, S; ksp_monitor = false, pc_type = "none")
    b = PETSc.PetscVec(petsclib, PetscScalar[2, 4])
    x = PETSc.solve!(similar(b), ksp, b)
    @test x[:] ≈ [1, 1]
    PETSc.destroy!.((x, b, ksp, opts))
    PETSc.finalize(petsclib)
end

@testset "solver options layer under the options database" begin
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
    PetscLib = typeof(petsclib)
    LP = PETSc.LibPETSc
    T = petsclib.PetscScalar
    comm = PETSc.MPI.COMM_SELF
    global_opts = LP.PetscOptions{PetscLib}(C_NULL, petsclib.age; own = false)
    has(key) = LP.PetscOptionsHasName(petsclib, global_opts, "", key)
    S = sparse(T[2 0; 0 4])
    b = PETSc.PetscVec(petsclib, T[2, 4])
    x = similar(b)
    rtol(ksp) = PETSc.tolerances(ksp).rtol

    LP.PetscOptionsSetValue(petsclib, global_opts, "-ksp_rtol", "1e-3")
    try
        # the database wins on a shared key; the other constructor keys apply
        ksp = PETSc.KSP(petsclib, comm, S; ksp_rtol = 1e-6, pc_type = "none")
        PETSc.solve!(x, ksp, b)
        @test rtol(ksp) ≈ 1e-3
        @test PETSc.type_name(PETSc.pc(ksp)) == :none
        @test !has("-pc_type")
        @test LP.PetscOptionsGetString(petsclib, global_opts, "", "-ksp_rtol") == "1e-3"

        # a solver without options reads the database at its first solve
        plain = PETSc.KSP(petsclib, comm, S)
        PETSc.solve!(x, plain, b)
        @test rtol(plain) ≈ 1e-3
        PETSc.destroy!.((ksp, plain))
    finally
        LP.PetscOptionsClearValue(petsclib, global_opts, "-ksp_rtol")
    end

    # a nested solver's key beats the one an enclosing block set, which comes back
    outer = PETSc.KSP(petsclib, comm, S; ksp_rtol = 1e-2)
    inner = PETSc.KSP(petsclib, comm, S; ksp_rtol = 1e-8)
    PETSc.with_options(outer) do
        PETSc.solve!(x, inner, b)
        @test LP.PetscOptionsGetString(petsclib, global_opts, "", "-ksp_rtol") == "0.01"
    end
    @test rtol(inner) ≈ 1e-8
    @test !has("-ksp_rtol")

    # a database pushed through LibPETSc inside a block: its keys are the user's
    for (user_kwargs, expected) in (((ksp_rtol = 1e-4,), 1e-4), ((;), 1e-8))
        userdb = PETSc.PetscOptions(petsclib; user_kwargs...)
        user_has(key) = LP.PetscOptionsHasName(petsclib, userdb, "", key)
        fresh = PETSc.KSP(petsclib, comm, S; ksp_rtol = 1e-8)
        PETSc.with_options(outer) do
            LP.PetscOptionsPush(petsclib, userdb)
            try
                PETSc.solve!(x, fresh, b)
            finally
                LP.PetscOptionsPop(petsclib)
            end
            @test LP.PetscOptionsGetString(petsclib, global_opts, "", "-ksp_rtol") == "0.01"
        end
        @test rtol(fresh) ≈ expected
        if isempty(user_kwargs)
            @test !user_has("-ksp_rtol")
        else
            @test userdb["ksp_rtol"] == "0.0001"
        end
        @test !has("-ksp_rtol")
        PETSc.destroy!.((fresh, userdb))
    end

    # a throw inside the block leaves no key behind
    bad = PETSc.KSP(petsclib, comm, S; ksp_type = "no_such_type")
    @test_throws LP.PetscError without_petsc_traceback(() -> PETSc.solve!(x, bad, b), petsclib)
    @test !has("-ksp_type")
    PETSc.destroy!.((outer, inner, bad, x, b))
    PETSc.finalize(petsclib)
end

@testset "DM options layer under the options database" begin
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
    PetscLib = typeof(petsclib)
    LP = PETSc.LibPETSc
    comm = PETSc.MPI.COMM_SELF
    nb = PETSc.DM_BOUNDARY_NONE
    global_opts = LP.PetscOptions{PetscLib}(C_NULL, petsclib.age; own = false)
    npoints(dm) = (v = PETSc.global_vec(dm); n = length(v); PETSc.destroy!(v); n)
    da(; kw...) = PETSc.DMDA(petsclib, comm, (nb,), (10,), 1; stencil_width = 1, kw...)

    # a constructor key absent from the database applies, and is gone afterwards
    dm = da(; da_refine = 2)
    @test npoints(dm) == 37
    @test !LP.PetscOptionsHasName(petsclib, global_opts, "", "-da_refine")
    PETSc.destroy!(dm)

    # the database wins on a shared key, and reaches a DM built without options
    LP.PetscOptionsSetValue(petsclib, global_opts, "-da_refine", "1")
    try
        dm, plain = da(; da_refine = 2), da()
        @test npoints(dm) == 19
        @test npoints(plain) == 19
        PETSc.destroy!.((dm, plain))
    finally
        LP.PetscOptionsClearValue(petsclib, global_opts, "-da_refine")
    end
    PETSc.finalize(petsclib)
end

@testset "parse_options" begin
    @test begin

        julia = joinpath(Sys.BINDIR, Base.julia_exename())
        # Force packages to be installed 
        run(`$(julia) --startup-file=no --project -e "using Pkg; Pkg.instantiate()" `)

        run(`$(julia) --startup-file=no --project -e "using PETSc
                                 using Test
                                 opts = PETSc.parse_options(ARGS)
                                 @test length(opts) == 4
                                 @test opts.ksp_monitor === nothing
                                 @test opts.ksp_view === nothing
                                 @test opts.pc_type === \"mg\"
                                 @test opts.da_grid_x === \"100\"
                                 @test_throws ErrorException opts.bad_opt
                                 "
                                 --
                                 -ksp_monitor
                                 -da_grid_x 100
                                 -ksp_view
                                 -pc_type mg`)
        true
    end
end

@testset "parse_option - options" begin
    opt = (tup = (1, 2, 3), string_tup = "1,2,3", string_int = "4", int = 4)
    @test PETSc.parse_option(opt, :bad_key, (1, 1, 1)) === (1, 1, 1)
    @test PETSc.parse_option(opt, :string_tup, (1, 1, 1)) === opt.tup
    @test PETSc.parse_option(opt, :string_tup, "a") === opt.string_tup
    @test PETSc.parse_option(opt, :bad_key, "a") === "a"
    @test PETSc.parse_option(opt, :int, Float64(7)) === Float64(opt.int)
    @test PETSc.parse_option(opt, :bad_key, Float64(7)) === Float64(7)
    @test PETSc.parse_option(opt, :tup, Float64.((1, 1, 1))) === Float64.(opt.tup)
    @test PETSc.parse_option(opt, :bad_key, Float64.((1, 1, 1))) ===
          Float64.((1, 1, 1))
    @test PETSc.parse_option(opt, :string_tup, Float64.((1, 1, 1))) ===
          Float64.(opt.tup)
    @test PETSc.parse_option(opt, :tup, (1, 1, 1)) === opt.tup
    @test PETSc.parse_option(opt, :bad_key, (1, 1, 1)) === (1, 1, 1)
    @test PETSc.parse_option(opt, :string_tup, (1, 1, 1)) === opt.tup
    @test PETSc.parse_option(opt, :string_tup, "a") === opt.string_tup
    @test PETSc.parse_option(opt, :bad_key, "a") === "a"
    @test PETSc.parse_option(opt, :int, 7) === opt.int
    @test PETSc.parse_option(opt, :bad_key, 7) === 7
    @test PETSc.parse_option(opt, :int, Float64(7)) === Float64(opt.int)
    @test PETSc.parse_option(opt, :bad_key, Float64(7)) === Float64(7)
    @test PETSc.parse_option(opt, :tup, Float64.((1, 1, 1))) === Float64.(opt.tup)
    @test PETSc.parse_option(opt, :bad_key, Float64.((1, 1, 1))) ===
          Float64.((1, 1, 1))
    @test PETSc.parse_option(opt, :string_tup, Float64.((1, 1, 1))) ===
          Float64.(opt.tup)
end

# a malformed command line is user input, so `ArgumentError`.
@testset "parse_options rejects malformed arguments" begin
    @test_throws ArgumentError PETSc.parse_options(["ksp_monitor"])
    @test_throws ArgumentError PETSc.parse_options(["-a=b=c"])
    @test PETSc.parse_options(["-pc_type", "mg"]) === (pc_type = "mg",)
end

# A value that starts with '-' is a negative number, as PETSc reads it
@testset "parse_options keeps negative and empty values" begin
    o = PETSc.parse_options(["-mat_shift", "-0.5", "-ksp_monitor", "-x0", "-.5"])
    @test Dict(pairs(o)) == Dict(:mat_shift => "-0.5", :ksp_monitor => nothing, :x0 => "-.5")
    @test PETSc.parse_options(["-prefix", ""]) === (prefix = "",)
end

@testset "a DM constructor that throws leaves the options stack as it was" begin
    petsclib = PETSc.petsclibs[1]
    PETSc.initialize(petsclib)
    comm = PETSc.MPI.COMM_SELF
    PetscLib = typeof(petsclib)
    nb = PETSc.DM_BOUNDARY_NONE
    # With nothing pushed, PETSc reads the global database
    global_opts = PETSc.LibPETSc.PetscOptions{PetscLib}(C_NULL, petsclib.age; own = false)
    visible(key) = PETSc.LibPETSc.PetscOptionsHasName(petsclib, global_opts, "", key)

    base = PETSc.DMStag(petsclib, comm, (nb,), (10,), (1, 1); stencil_width = 1)
    # Each call fails inside PETSc, on an option value that isn't an integer
    cases = [
        "-da_refine" => () -> PETSc.DMDA(petsclib, comm, (nb,), (10,), 1; stencil_width = 1, da_refine = "abc"),
        "-stag_grid_x" => () -> PETSc.DMStag(petsclib, comm, (nb,), (10,), (1, 1); stencil_width = 1, stag_grid_x = "abc"),
        "-stag_grid_x" => () -> PETSc.DMStag(base, (1, 1); stag_grid_x = "abc"),
        "-dm_refine" => () -> PETSc.DMPlex(petsclib, comm; dm_refine = "abc"),
        "-dm_refine" => () -> PETSc.DMPlex(petsclib, comm, 2, false, [2, 2]; dm_refine = "abc"),
        "-petscspace_degree" => () -> PETSc.fe_create_default(petsclib, comm, 2, 1, false; degree = -3),
    ]
    for (key, construct) in cases
        @test_throws PETSc.LibPETSc.PetscError without_petsc_traceback(construct, petsclib)
        @test isempty(PETSc.options_stack(PetscLib))
        @test !visible(key)
    end
    PETSc.destroy!(base)
    PETSc.finalize(petsclib)
end

@testset "with_options has its docstring" begin
    # the method for `opts` has a docstring of its own, next to the one for `obj`
    docs = Base.Docs.meta(PETSc)[Base.Docs.Binding(PETSc, :with_options)].docs
    @test any(sig -> occursin("AbstractPetscOptions", string(sig)), keys(docs))
end
