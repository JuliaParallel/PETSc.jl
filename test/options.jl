using Test
using PETSc

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

        opts = PETSc.Options(petsclib; kw_opts...)
        
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

@testset "typedget - options" begin
    opt = (tup = (1, 2, 3), string_tup = "1,2,3", string_int = "4", int = 4)
    @test PETSc.typedget(opt, :bad_key, (1, 1, 1)) === (1, 1, 1)
    @test PETSc.typedget(opt, :string_tup, (1, 1, 1)) === opt.tup
    @test PETSc.typedget(opt, :string_tup, "a") === opt.string_tup
    @test PETSc.typedget(opt, :bad_key, "a") === "a"
    @test PETSc.typedget(opt, :int, Float64(7)) === Float64(opt.int)
    @test PETSc.typedget(opt, :bad_key, Float64(7)) === Float64(7)
    @test PETSc.typedget(opt, :tup, Float64.((1, 1, 1))) === Float64.(opt.tup)
    @test PETSc.typedget(opt, :bad_key, Float64.((1, 1, 1))) ===
          Float64.((1, 1, 1))
    @test PETSc.typedget(opt, :string_tup, Float64.((1, 1, 1))) ===
          Float64.(opt.tup)
    @test PETSc.typedget(opt, :tup, (1, 1, 1)) === opt.tup
    @test PETSc.typedget(opt, :bad_key, (1, 1, 1)) === (1, 1, 1)
    @test PETSc.typedget(opt, :string_tup, (1, 1, 1)) === opt.tup
    @test PETSc.typedget(opt, :string_tup, "a") === opt.string_tup
    @test PETSc.typedget(opt, :bad_key, "a") === "a"
    @test PETSc.typedget(opt, :int, 7) === opt.int
    @test PETSc.typedget(opt, :bad_key, 7) === 7
    @test PETSc.typedget(opt, :int, Float64(7)) === Float64(opt.int)
    @test PETSc.typedget(opt, :bad_key, Float64(7)) === Float64(7)
    @test PETSc.typedget(opt, :tup, Float64.((1, 1, 1))) === Float64.(opt.tup)
    @test PETSc.typedget(opt, :bad_key, Float64.((1, 1, 1))) ===
          Float64.((1, 1, 1))
    @test PETSc.typedget(opt, :string_tup, Float64.((1, 1, 1))) ===
          Float64.(opt.tup)
end