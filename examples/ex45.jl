# INCLUDE IN MPI TEST
# Translation of PETSc's ex45.c to PETSc.jl
using PETSc, MPI, Printf

# boundary data and forcing  
exact(x, y, z) = sin(2π * x) * sin(2π * y) * sin(2π * z)
# The discrete operator is -∇², so forcing = -∇²(exact)
# For exact(x,y,z) = sin(2πx)sin(2πy)sin(2πz), we have ∇²u = -3(2π)²sin(2πx)sin(2πy)sin(2πz)
# So forcing = -∇²u = 12π²sin(2πx)sin(2πy)sin(2πz)
forcing(x, y, z) = 12 * π^2 * sin(2π * x) * sin(2π * y)  * sin(2π * z)

# Get the PETSc lib with our chosen `PetscScalar` type
petsclib = PETSc.getlib(; PetscScalar = Float64)

# Initialize PETSc
PETSc.initialize(petsclib)

function solve_ex45(N=7; da_grid_x=7, da_grid_y=7, da_grid_z=7, kwargs...)
    comm = MPI.COMM_WORLD

    # Parse command-line options and merge with keyword arguments passed
    # programmatically. Keyword args take precedence over CLI.
    cli_opts = PETSc.parse_options(ARGS)           # NamedTuple
    # Build a flexible Dict{Symbol,Any} to allow merging CLI and kwargs
    cli = Dict{Symbol,Any}()
    for (k, v) in pairs(cli_opts)
        cli[k] = v
    end
    for (k, v) in pairs(kwargs)
        cli[k] = v
    end

    # Default RHS: continuous analytic forcing (scaled by cell volume)
    # We keep CLI/kwargs simple — any specific forcing-mode option was
    # removed: continuous mode is now the only behavior.

    # Determine N (keyword/arg default overridden by CLI/kwargs)
    function _as_int(x, default)
        if x === nothing
            return default
        elseif x isa Integer
            return Int(x)
        else
            return parse(Int, string(x))
        end
    end

    N = _as_int(get(cli, :N, N), N)
    dim = _as_int(get(cli, :dim, 3), 3)

    # Re-pack merged options into a NamedTuple for keyword splatting
    opts = (; [(k => cli[k]) for k in keys(cli)]...)

    # Create a 3D DMDA (defaults 7x7x7 like the C example)
    boundary = Tuple(fill(PETSc.DM_BOUNDARY_NONE,dim))
    stencil = PETSc.DMDA_STENCIL_STAR
    global_size = Tuple(fill(N,dim))
    # Pass parsed opts into DMDA so other DM options are applied
    da = PETSc.DMDA(petsclib, comm, boundary, global_size, 1, 1, stencil; opts...)

    # Print fine grid resolution at start (only on rank 0)
    final_grid = PETSc.getinfo(da).global_size
    if MPI.Comm_rank(comm) == 0
        @printf("Solving on %d × %d × %d grid\n", final_grid[1], final_grid[2], final_grid[3])
    end

    # Create KSP tied to DM and pass parsed PETSc options so
    # command-line flags (e.g. -ksp_monitor) are applied like in `ex50.jl`
    ksp = PETSc.KSP(da; opts...)

    # Set compute operators (matrix assembly) and RHS
    PETSc.setcomputeoperators!(ksp) do J, jac, ksp
        dm = PETSc.getDMDA(ksp)
        corners = PETSc.getcorners(dm)
        info = PETSc.getinfo(dm)
        N = info.global_size

        Hx,Hy,Hz = 1.0 ./ (N .- 1)
        HxHydHz  = Hx * Hy / Hz
        HxHzdHy  = Hx * Hz / Hy
        HyHzdHx  = Hy * Hz / Hx

        HxdHy = Hx / Hy
        HydHx = Hy / Hx
        for k in corners.lower[3]:corners.upper[3]
            for j in corners.lower[2]:corners.upper[2]
                for i in corners.lower[1]:corners.upper[1]
                    idx = CartesianIndex(i, j, k)
                    # boundary check (global indices start at 1 here)
                    is_boundary = (i == 1 || j == 1 || k == 1 || i == N[1] || j == N[2] || k == N[3])
                    if dim==3
                        if is_boundary
                            # Dirichlet BC: set row to identity
                            jac[idx, idx] = 1.0
                        else
                            jac[idx, idx + CartesianIndex(0,-1,0)] = -HxHzdHy
                            jac[idx, idx + CartesianIndex(-1,0,0)] = -HyHzdHx
                            jac[idx, idx] = 2.0 * (HxHydHz + HxHzdHy + HyHzdHx)
                            jac[idx, idx + CartesianIndex(1,0,0)] = -HyHzdHx
                            jac[idx, idx + CartesianIndex(0,1,0)] = -HxHzdHy
                            jac[idx, idx + CartesianIndex(0,0,-1)] = -HxHydHz
                            jac[idx, idx + CartesianIndex(0,0,1)] = -HxHydHz
                        end
                    elseif dim==2
                        if is_boundary
                            # Dirichlet BC: set row to identity
                            jac[idx, idx] = 1.0
                        else
                            # Interior point - full 5-point stencil
                            jac[idx, idx + CartesianIndex(0, -1, 0)] = -HxdHy
                            jac[idx, idx + CartesianIndex(-1, 0, 0)] = -HydHx
                            jac[idx, idx] = 2.0 * (HxdHy + HydHx)
                            jac[idx, idx + CartesianIndex(1, 0, 0)] = -HydHx
                            jac[idx, idx + CartesianIndex(0, 1, 0)] = -HxdHy
                        end
                    end

                end
            end
        end
        PETSc.assemble!(jac)
        return 0
    end

    PETSc.setcomputerhs!(ksp) do b_vec, ksp
        dm = PETSc.getDMDA(ksp)
        corners = PETSc.getcorners(dm)
        info = PETSc.getinfo(dm)
        N = info.global_size

        Hx,Hy,Hz = 1.0 ./ (N .- 1)
        HxHydHz  = Hx * Hy / Hz
        HxHzdHy  = Hx * Hz / Hy
        HyHzdHx  = Hy * Hz / Hx

        g_x = range(0, length = N[1], stop = 1)
        g_y = range(0, length = N[2], stop = 1)
        g_z = range(0, length = N[3], stop = 1)
    
        l_x = g_x[(corners.lower[1]):(corners.upper[1])]
        l_y = g_y[(corners.lower[2]):(corners.upper[2])]
        l_z = g_z[(corners.lower[3]):(corners.upper[3])]
    
        PETSc.withlocalarray!(b_vec; read=false) do b
            # reshape into local block (xm, ym, zm)
            N = corners.size;
            #nx = Int(corners.size[1]); ny = Int(corners.size[2]); nz = Int(corners.size[3])

            if dim==3
                #is_boundary = (x == 1 || y == 1 || x == global_size[1] || y == global_size[2])
                b3 = reshape(b, N...)
                for kk = 1:N[3], jj = 1:N[2], ii = 1:N[1]
                    # global indices
                    x,y,z = l_x[ii], l_y[jj], l_z[kk]

                    gi = corners.lower[1] + (ii - 1)
                    gj = corners.lower[2] + (jj - 1)
                    gk = corners.lower[3] + (kk - 1)
                    if gi == 1 || gj == 1 || gk == 1 || gi == N[1] || gj == N[2] || gk == N[3]
                        # Dirichlet BC: RHS is analytic boundary value
                        b3[ii, jj, kk] = exact(x, y, z)
                    else
                        # Continuous forcing (default): analytic -∇² evaluated
                        # at the point and scaled by the cell volume so the
                        # RHS uses the same quadrature scaling as the matrix.
                        vol_cell = Hx * Hy * Hz
                        b3[ii, jj, kk] = forcing(x, y, z) * vol_cell
                    end
                end
            else

            end

        end
        PETSc.assemble!(b_vec)
        return 0
    end

    # Time the solve so callers can inspect runtime
    solve_time = @elapsed PETSc.solve!(ksp)

    # Number of outer KSP iterations
    niter = LibPETSc.KSPGetIterationNumber(petsclib, ksp)

    # Final (fine) grid resolution
    final_grid = PETSc.getinfo(da).global_size

    # Get the residual norm
    x = LibPETSc.KSPGetSolution(petsclib, ksp)
    norm = LibPETSc.KSPGetResidualNorm(petsclib, ksp)

    # Compute L2 and max error against analytic solution `exact(x,y,z)`
    info = PETSc.getinfo(da)
    Nglob = info.global_size
    Hx,Hy,Hz = 1.0 ./ (Nglob .- 1)
    vol = Hx * Hy * Hz

    local_err2 = 0.0
    local_max = 0.0
    PETSc.withlocalarray!(x; read=true) do xu
        corners = PETSc.getcorners(da)
        dims = corners.size
        uas = reshape(xu, dims...)
        for kk = 1:dims[3], jj = 1:dims[2], ii = 1:dims[1]
            gi = corners.lower[1] + (ii - 1)
            gj = corners.lower[2] + (jj - 1)
            gk = corners.lower[3] + (kk - 1)
            xc = (gi - 1) * Hx
            yc = (gj - 1) * Hy
            zc = (gk - 1) * Hz
            u_num = uas[ii, jj, kk]
            u_ex = exact(xc, yc, zc)
            err = u_num - u_ex
            local_err2 += err^2 * vol
            local_max = max(local_max, abs(err))
        end
    end

    global_err2 = MPI.Allreduce(local_err2, MPI.SUM, comm)
    global_max = MPI.Allreduce(local_max, MPI.MAX, comm)
    L2 = sqrt(global_err2)

    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        @printf("L2 error: %.6e\n", L2)
        @printf("Max error: %.6e\n", global_max)
        @printf("Residual norm %g\n", norm)
    end

    # Clean up
    PETSc.destroy(ksp)
    return (;norm,final_grid,niter,solve_time,L2,global_max)
end

if !isinteractive() && abspath(PROGRAM_FILE) == @__FILE__
    res = solve_ex45()
    # res is a NamedTuple: (norm, final_grid, niter, solve_time, L2, max)
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        @printf("Final residual norm %g\n", res.norm)
        @printf("Fine grid resolution: %d × %d × %d\n", res.final_grid[1], res.final_grid[2], res.final_grid[3])
        @printf("Outer KSP iterations: %d\n", res.niter)
        @printf("Solve time: %.6f seconds\n", res.solve_time)
        @printf("L2 error: %.6e\n", res.L2)
        @printf("Max error: %.6e\n", res.global_max)
    end
    PETSc.finalize(petsclib)
end
