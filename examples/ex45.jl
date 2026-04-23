# INCLUDE IN MPI TEST
# Translation of PETSc's ex45.c to PETSc.jl
using PETSc, MPI, Printf

# boundary data and forcing  
exact(x, y, z) = sin(2π * x) * sin(2π * y) * sin(2π * z)
# The discrete operator is -∇², so forcing = -∇²(exact)
# For exact(x,y,z) = sin(2πx)sin(2πy)sin(2πz), we have ∇²u = -3(2π)²sin(2πx)sin(2πy)sin(2πz)
# So forcing = -∇²u = 12π²sin(2πx)sin(2πy)sin(2πz)
forcing(x, y, z) = 12 * π^2 * sin(2π * x) * sin(2π * y) * sin(2π * z)

# Get the PETSc lib with our chosen `PetscScalar` type
petsclib = PETSc.getlib(; PetscScalar = Float64)

# If we want to use a custom PETSc library, we can set it with `PETSc.set_library!` before calling `PETSc.getlib`. For example:
#=
using PETSc
PETSc.set_library!(
    "/path/to/custom/libpetsc.so";
    PetscScalar = Float64,
    PetscInt    = Int64,
)
# Restart Julia — PETSc_jll is not loaded and your library is used automatically.
=#
# This creates a LocalPreferences.toml file. 
#petsclib = PETSc.getlib(; PetscScalar = Float64, PetscInt = Int64)


# Initialize PETSc
PETSc.initialize(petsclib, log_view=true)

function solve_ex45(N=7; da_grid_x=7, da_grid_y=7, da_grid_z=7, kwargs...)
    comm = MPI.COMM_WORLD

    cli_opts = PETSc.parse_options(ARGS)
    cli = Dict{Symbol,Any}()
    for (k, v) in pairs(cli_opts)
        cli[k] = v
    end
    for (k, v) in pairs(kwargs)
        cli[k] = v
    end

    function _as_int(x, default)
        if x === nothing
            return default
        elseif x isa Integer
            return Int(x)
        else
            return parse(Int, string(x))
        end
    end

    N   = _as_int(get(cli, :N,   N),   N)
    dim = _as_int(get(cli, :dim, 3),   3)
    Nx  = _as_int(get(cli, :Nx,  N),   N)
    Ny  = _as_int(get(cli, :Ny,  N),   N)
    Nz  = _as_int(get(cli, :Nz,  N),   N)

    opts = (; [(k => cli[k]) for k in keys(cli)]...)

    boundary    = Tuple(fill(PETSc.DM_BOUNDARY_NONE, dim))
    stencil     = PETSc.DMDA_STENCIL_STAR
    global_size = dim == 3 ? (Nx, Ny, Nz) : (Nx, Ny, 1)
    da = PETSc.DMDA(petsclib, comm, boundary, global_size, 1, 1, stencil; opts...)

    final_grid = PETSc.getinfo(da).global_size
    if MPI.Comm_rank(comm) == 0
        @printf("Solving on %d × %d × %d grid\n", final_grid[1], final_grid[2], final_grid[3])
    end

    ksp = PETSc.KSP(da; opts...)

    PETSc.setcomputeoperators!(ksp) do J, jac, ksp
        dm      = PETSc.getDM(ksp)
        corners = PETSc.getcorners(dm)
        info    = PETSc.getinfo(dm)
        N       = info.global_size

        # Grid spacings (uniform)
        Hx, Hy, Hz = 1.0 ./ (N .- 1)

        # Scale entire system by Hx*Hy*Hz (= h³ for uniform grid) so that
        # matrix entries are O(h) ~ O(1/N) rather than O(1/h²) ~ O(N²).
        # This keeps MG smoothers well-conditioned regardless of resolution.
        #
        # Scaled operator: (Hx*Hy*Hz) * (-∇²_h)
        # For uniform h: scaled diagonal  = 6h³/h² = 6h
        #                scaled off-diag  = -h³/h² = -h
        # All entries O(h) → matrix is well-scaled for MG.
        scale = Hx * Hy * Hz

        # Scaled coefficients
        cx = scale / (Hx * Hx)   # = Hy*Hz/Hx
        cy = scale / (Hy * Hy)   # = Hx*Hz/Hy
        cz = scale / (Hz * Hz)   # = Hx*Hy/Hz
        diag3d = 2.0 * (cx + cy + cz)
        diag2d = 2.0 * (cx + cy)

        for k in corners.lower[3]:corners.upper[3]
            for j in corners.lower[2]:corners.upper[2]
                for i in corners.lower[1]:corners.upper[1]
                    idx = CartesianIndex(i, j, k)
                    is_boundary = (i == 1 || j == 1 || k == 1 ||
                                   i == N[1] || j == N[2] || k == N[3])
                    if dim == 3
                        if is_boundary
                            jac[idx, idx] = 1.0
                        else
                            jac[idx, idx + CartesianIndex(-1,  0,  0)] = -cx
                            jac[idx, idx + CartesianIndex( 1,  0,  0)] = -cx
                            jac[idx, idx + CartesianIndex( 0, -1,  0)] = -cy
                            jac[idx, idx + CartesianIndex( 0,  1,  0)] = -cy
                            jac[idx, idx + CartesianIndex( 0,  0, -1)] = -cz
                            jac[idx, idx + CartesianIndex( 0,  0,  1)] = -cz
                            jac[idx, idx]                               =  diag3d
                        end
                    elseif dim == 2
                        if is_boundary
                            jac[idx, idx] = 1.0
                        else
                            jac[idx, idx + CartesianIndex(-1, 0, 0)] = -cx
                            jac[idx, idx + CartesianIndex( 1, 0, 0)] = -cx
                            jac[idx, idx + CartesianIndex( 0,-1, 0)] = -cy
                            jac[idx, idx + CartesianIndex( 0, 1, 0)] = -cy
                            jac[idx, idx]                             =  diag2d
                        end
                    end
                end
            end
        end
        PETSc.assemble!(jac)
        return 0
    end

    PETSc.setcomputerhs!(ksp) do b_vec, ksp
        dm      = PETSc.getDM(ksp)
        corners = PETSc.getcorners(dm)
        info    = PETSc.getinfo(dm)
        N       = info.global_size

        Hx, Hy, Hz = 1.0 ./ (N .- 1)
        scale = Hx * Hy * Hz   # same scaling as matrix

        g_x = range(0, length = N[1], stop = 1)
        g_y = range(0, length = N[2], stop = 1)
        g_z = range(0, length = N[3], stop = 1)

        l_x = g_x[(corners.lower[1]):(corners.upper[1])]
        l_y = g_y[(corners.lower[2]):(corners.upper[2])]
        l_z = g_z[(corners.lower[3]):(corners.upper[3])]

        PETSc.withlocalarray!(b_vec; read=false) do b
            sz = corners.size

            if dim == 3
                b3 = reshape(b, sz...)
                for kk = 1:sz[3], jj = 1:sz[2], ii = 1:sz[1]
                    x, y, z = l_x[ii], l_y[jj], l_z[kk]
                    gi = corners.lower[1] + (ii - 1)
                    gj = corners.lower[2] + (jj - 1)
                    gk = corners.lower[3] + (kk - 1)
                    if gi == 1 || gj == 1 || gk == 1 ||
                       gi == N[1] || gj == N[2] || gk == N[3]
                        # Dirichlet BC: identity row in matrix → RHS = exact value
                        b3[ii, jj, kk] = exact(x, y, z)
                    else
                        # Scale forcing by same factor as matrix
                        b3[ii, jj, kk] = forcing(x, y, z) * scale
                    end
                end
            elseif dim == 2
                b2 = reshape(b, sz...)
                for jj = 1:sz[2], ii = 1:sz[1]
                    x, y = l_x[ii], l_y[jj]
                    gi = corners.lower[1] + (ii - 1)
                    gj = corners.lower[2] + (jj - 1)
                    if gi == 1 || gj == 1 || gi == N[1] || gj == N[2]
                        b2[ii, jj, 1] = exact(x, y, 0.0)
                    else
                        b2[ii, jj, 1] = forcing(x, y, 0.0) * scale
                    end
                end
            end
        end
        PETSc.assemble!(b_vec)
        return 0
    end

    solve_time = @elapsed PETSc.solve!(ksp)
    niter      = LibPETSc.KSPGetIterationNumber(petsclib, ksp)
    final_grid = PETSc.getinfo(da).global_size

    x    = LibPETSc.KSPGetSolution(petsclib, ksp)
    norm = LibPETSc.KSPGetResidualNorm(petsclib, ksp)

    # Compute L2 and max error against analytic solution
    info  = PETSc.getinfo(da)
    Nglob = info.global_size
    Hx, Hy, Hz = 1.0 ./ (Nglob .- 1)
    vol = Hx * Hy * Hz

    local_err2 = 0.0
    local_max  = 0.0
    PETSc.withlocalarray!(x; read=true) do xu
        corners = PETSc.getcorners(da)
        dims    = corners.size
        uas     = reshape(xu, dims...)
        for kk = 1:dims[3], jj = 1:dims[2], ii = 1:dims[1]
            gi = corners.lower[1] + (ii - 1)
            gj = corners.lower[2] + (jj - 1)
            gk = corners.lower[3] + (kk - 1)
            xc = (gi - 1) * Hx
            yc = (gj - 1) * Hy
            zc = (gk - 1) * Hz
            u_num = uas[ii, jj, kk]
            u_ex  = exact(xc, yc, zc)
            err   = u_num - u_ex
            local_err2 += err^2 * vol
            local_max   = max(local_max, abs(err))
        end
    end

    global_err2 = MPI.Allreduce(local_err2, MPI.SUM, comm)
    global_max  = MPI.Allreduce(local_max,  MPI.MAX, comm)
    L2 = sqrt(global_err2)

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        @printf("L2 error: %.6e\n", L2)
        @printf("Max error: %.6e\n", global_max)
        @printf("Residual norm %g\n", norm)
        @printf("Solve time: %.6f seconds\n", solve_time)
    end

    PETSc.destroy(ksp)
    return (; norm, final_grid, niter, solve_time, L2, global_max)
end

if !isinteractive() && abspath(PROGRAM_FILE) == @__FILE__
    res = solve_ex45()
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        @printf("Final residual norm %g\n", res.norm)
        @printf("Fine grid resolution: %d × %d × %d\n",
                res.final_grid[1], res.final_grid[2], res.final_grid[3])
        @printf("Outer KSP iterations: %d\n", res.niter)
        @printf("Solve time: %.6f seconds\n", res.solve_time)
        @printf("L2 error: %.6e\n", res.L2)
        @printf("Max error: %.6e\n", res.global_max)
    end
    PETSc.finalize(petsclib)
end