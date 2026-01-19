# Getting started


- [Getting started](#getting-started)
    - [1a. Installation using pre-built libraries](#1a-installation-using-pre-built-libraries)
    - [1b. Installation using pre-built libraries](#1b-installation-using-pre-built-libraries)
    - [2. Solving a linear system of equations](#2-solving-a-linear-system-of-equations)
    - [3. Nonlinear example](#3-nonlinear-example)
    - [4. Next steps](#4-next-steps)

### 1a. Installation using pre-built libraries 
The easiest way to install the package is: 
```julia
julia> ]
(@v1.12) pkg> add PETSc
```
which will install a pre-built PETSc library (`PETSc_jll`) as well as `MPI.jl` on your system. This will work both in serial and in parallel on your machine.

!!! warning "Windows Users"
    The prebuild binaries currently do not work on Windows as we had to build `PETSc_jll` without MPI due to compatibility issues with `MicrosoftMPI_jll`.
    
    **Windows users are therefore advised to install the [Windows Subsystem for Linux](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux) (WSL) and run PETSc.jl from within WSL.** This will provide full functionality with both serial and parallel (MPI) support.

### 1b. Installation using pre-built libraries 
On many high-performance clusters, you will have to use the provided `MPI` installation for that cluster and the default download above will not be sufficient. Alternatively, you may be interested in a PETSc installation that comes with additional external packages. Ensure that this PETSc installation is compiled as a dynamic (and not a static) library, after which you need to specify the correct library with:

```julia
# Create custom library instance
petsclib = set_petsclib("/path/to/custom/libpetsc.so"; 
                             PetscScalar=Float64, PetscInt=Int64)
# Use it like any precompiled library
PETSc.initialize(petsclib, log_view=true)
# ... your code ...
PETSc.finalize(petsclib)
```


### 2. Solving a linear system of equations

Lets consider the following elliptic equation:
```math
\begin{aligned}
{\partial^2 T \over \partial x^2}  = 0, T(0) = 1, T(1) = 11
\end{aligned}
```

```julia
julia> using PETSc
julia> petsclib = PETSc.petsclibs[1]
julia> PETSc.initialize(petsclib, log_view=false)
```
Note that if you initialize the PETSc library with the option `log_view=true` you will get a detailed overview of your simulation once you close the library with `PETSc.finalize(petsclib)` (note that it's set to `false` by default). 
Next, lets define the number of gridpoints and the spacing between them with:
```julia
julia> n   =  11
julia> Δx  =  1. / (n - 1)
```

Let's first define the matrix with coefficients:
```julia
julia> nnz =  ones(Int64,n); nnz[2:n-1] .= 3;
julia> A   =  PETSc.MatSeqAIJ(petsclib,n,n,nnz);
julia> for i=2:n-1
            A[i,i-1] =  1/Δx^2
            A[i,i  ] = -2/Δx^2
            A[i,i+1] =  1/Δx^2
       end;
julia> A[1,1]=1; A[n,n]=1;  # boundary conditions (Dirichlet)
julia> PETSc.assemble!(A)   # Assemble the matrix
```
This creates a sequential matrix (on one processor):
```julia
julia> A
PETSc seqaij Mat of size (11, 11)
```
Now, lets define the right-hand-size vector `rhs` as a julia vector:
```julia
julia> rhs = zeros(n); rhs[1]=1; rhs[11]=11;
```

Next, we define the linear solver for the matrix `A`, which is done by setting a `KSP` solver: 
```julia
julia> ksp = PETSc.KSP(A; ksp_rtol=1e-8, pc_type="jacobi", ksp_monitor=true)
```
Note that you can specify all PETSc command-line options as keywords here.

Solving the system of equations is simple:
```julia
julia> sol = ksp\rhs
  0 KSP Residual norm 1.104536101719e+01 
  1 KSP Residual norm 4.939635614091e+00 
  2 KSP Residual norm 2.410295378065e+00 
  3 KSP Residual norm 1.462993806273e+00 
  4 KSP Residual norm 1.004123728835e+00 
  5 KSP Residual norm 7.700861485629e-01 
  6 KSP Residual norm 6.165623662013e-01 
  7 KSP Residual norm 4.972507567923e-01 
  8 KSP Residual norm 4.074986825669e-01 
  9 KSP Residual norm 3.398492183940e-01 
 10 KSP Residual norm 3.283015493450e-15 
11-element Vector{Float64}:
  1.0
  2.000000000000001
  3.0000000000000013
  4.000000000000001
  5.000000000000002
  6.0
  7.0000000000000036
  8.000000000000002
  9.000000000000004
 10.000000000000002
 11.0
```
This solves the system of equations with the default (iterative) solver of PETSc. You can find out what is being used by adding the `ksp_view=true` option to the `ksp` object:

```julia
julia> ksp = PETSc.KSP(A; ksp_rtol=1e-8, pc_type="jacobi", ksp_monitor=true, ksp_view=true)
```
which will tell you:
```julia
KSP Object: 1 MPI process
  type: gmres
    restart=30, using Classical (unmodified) Gram-Schmidt Orthogonalization with no iterative refinement
    happy breakdown tolerance 1e-30
  maximum iterations=10000, initial guess is zero
  tolerances: relative=1e-08, absolute=1e-50, divergence=10000.
  left preconditioning
  using PRECONDITIONED norm type for convergence test
PC Object: 1 MPI process
  type: jacobi
    type DIAGONAL
  linear system matrix = precond matrix:
  Mat Object: 1 MPI process
    type: seqaij
    rows=11, cols=11
    total: nonzeros=29, allocated nonzeros=29
    total number of mallocs used during MatSetValues calls=0
      not using I-node routines
```      
So we are using a `gmres` solver with a `jacobi` preconditioner.
The power of PETSc is that you can change the solver on the fly.
We can solve the same system with a direct solver, using:
```julia
julia> ksp = PETSc.KSP(A; ksp_rtol=1e-8, pc_type="lu", ksp_monitor=true, ksp_view=true, ksp_type="preonly");
julia> sol = ksp\rhs;
  0 KSP Residual norm 1.104536101719e+01
  1 KSP Residual norm 2.846442092393e-13
KSP Object: 1 MPI process
  type: preonly
  maximum iterations=10000, initial guess is zero
  tolerances: relative=1e-05, absolute=1e-50, divergence=10000.
  left preconditioning
  using NONE norm type for convergence test
PC Object: 1 MPI process
  type: lu
    out-of-place factorization
    tolerance for zero pivot 2.22045e-14
    matrix ordering: nd
    factor fill ratio given 5., needed 1.31034
      Factored matrix follows:
        Mat Object: 1 MPI process
          type: seqaij
          rows=11, cols=11
          package used to perform factorization: petsc
          total: nonzeros=38, allocated nonzeros=38
            not using I-node routines
  linear system matrix = precond matrix:
  Mat Object: 1 MPI process
    type: seqaij
    rows=11, cols=11
    total: nonzeros=29, allocated nonzeros=29
    total number of mallocs used during MatSetValues calls=0
      not using I-node routines
```
which converges, as expected, in 1 iteration.

And since we are using julia, plotting the solution can be done with
```julia
julia> using Plots
julia> plot(0:Δx:1,sol, ylabel="solution",xlabel="x")
```

![linear_solution](../assets/img/linear_KSP_solution.png)


Note that in general, one knows more about the equations than just the matrix entries. In this case, for example, we use a finite difference discretization to solve the equation. We could also have used a finite element code to do the same. 
This is important information that helps PETSc to distribute the problem over several processors, or to setup multigrid preconditioners etc.
PETSc uses the `DM` infrastructure for such cases. `DMDA` is for (collocated) finite differences, `DMStag` for staggered finite differences, `DMPlex` for finite element/finite volume discretisations and `DMForest` for adaptive mesh refinement problems.
If possible, use this infrastructure as it simplifies your life. Have a look at the [examples](https://github.com/JuliaParallel/PETSc.jl/tree/main/examples) or [tests](https://github.com/JuliaParallel/PETSc.jl/tree/main/test) of `PETSc.jl`.

### 3. Nonlinear example
Let's solve the coupled system of nonlinear equations: 
```math
\begin{aligned}
x^2 + x y  &= 3 \\
x y + y^2  &= 6
\end{aligned}
```
for ``x`` and ``y``, which can be written in terms of a residual vector ``r``:
```math
r = \binom{ x^2 + x y  - 3} {x y + y^2  - 6}
```

We start by initializing PETSc:
```julia
julia> using PETSc, MPI
julia> petsclib = PETSc.petsclibs[1]
julia> PETSc.initialize(petsclib, log_view=false)
julia> PetscScalar = petsclib.PetscScalar
julia> PetscInt = petsclib.PetscInt
```
In order to solve this, we need to provide a residual function that computes ``r``:
```julia
julia> function Residual!(rx, snes, x)
            rx[1] = x[1]^2 + x[1] * x[2] - 3
            rx[2] = x[1] * x[2] + x[2]^2 - 6
            return PetscInt(0)  # petsc success
        end
```

In addition, we need to provide the Jacobian:
```math
J = 
\begin{pmatrix}
\frac{\partial r_1}{ \partial x} & \frac{\partial r_1}{ \partial y}  \\
\frac{\partial r_2}{ \partial x} & \frac{\partial r_2}{ \partial y}  \\
\end{pmatrix}
= 
\begin{pmatrix}
2x + y & x  \\
y & x + 2y  \\
\end{pmatrix}
```
In Julia, this is:
```julia
julia> function updateJ!(J, snes, x)
            J[1, 1] = 2x[1] + x[2]
            J[1, 2] = x[1]
            J[2, 1] = x[2]
            J[2, 2] = x[1] + 2x[2]

            PETSc.assemble!(J)
            return PetscInt(0)
        end
```
In order to solve this using the PETSc nonlinear equation solvers, you first define the `SNES` solver together with the jacobian and residual functions as 
```julia
julia> snes = PETSc.SNES(petsclib,MPI.COMM_SELF; ksp_rtol=1e-4, pc_type="none")
julia> r = PETSc.VecSeq(petsclib, zeros(PetscScalar, 2))
julia> PETSc.setfunction!(snes, Residual!, r)
julia> J = zeros(2,2)
julia> PJ = PETSc.MatSeqDense(petsclib,J)
julia> PETSc.setjacobian!(updateJ!, snes, PJ)
```

You can solve this as:
```julia
julia> x = PETSc.VecSeq(petsclib, [2.0, 3.0])
julia> b = PETSc.VecSeq(petsclib, [0.0, 0.0])
julia> PETSc.solve!(x, snes, b)
julia> x[:]
2-element Vector{Float64}:
 1.000000003259629
 1.999999998137355
```
which indeed recovers the analytical solution ``(x=1, y=2)``.

If we are done, finalize it with:
```julia
julia> PETSc.finalize(petsclib)
```


### 4. Next steps 
Now that you have the basics, ypu can start playing with some more complete examples.
Here some suggestions:
1. [laplacian.jl](https://github.com/JuliaParallel/PETSc.jl/blob/main/examples/laplacian.jl) - simple (non-parallel) laplacian example using Julia sparse matrixes
2. [dmda_laplacian.jl](https://github.com/JuliaParallel/PETSc.jl/blob/main/examples/dmda_laplacian.jl) - 2D laplacian example with Dirichlet boundary conditions using the `DMDA` framework. Examples are given on how to run it with various (multigrid) solvers.
3. [ex50.jl](https://github.com/JuliaParallel/PETSc.jl/blob/main/examples/ex50.jl) - 2D parallel `DMDA` laplacian example with Neumann boundary conditions, which requires the nullspace to be removed.  
4. [SNES_ex2.jl](https://github.com/JuliaParallel/PETSc.jl/blob/main/examples/SNES_ex2.jl) - 1D laplacian with nonlinear terms where the hand-derived jacobian is hardcoded
5. [SNES_ex2b.jl](https://github.com/JuliaParallel/PETSc.jl/blob/main/examples/SNES_ex2b.jl) - as `SNES_ex2.jl` but using automatic differentiation to derive the jacobian.
6. [Liouville_Bratu_Gelfand.jl](https://github.com/JuliaParallel/PETSc.jl/blob/main/examples/Liouville\_Bratu\_Gelfand.jl) - 1D/2D/3D poisson equation with nonlinear terms. Shows how to combine `DMDA` with `SNES` solvers and solve them in parallel, if you have a jacobian.
7. [porosity_waves.jl](https://github.com/JuliaParallel/PETSc.jl/blob/main/examples/porosity_waves.jl) nD MPI-parallel example of 2 coupled nonlinear PDE's using the `DMDA` framework along with automatic differentiation to derive the jacobian.

Working through those examples should give you a fair idea of how to use PETSc. 

If you have further questions, please have a look at the [test](https://github.com/JuliaParallel/PETSc.jl/tree/main/test) directory; this is run automatically every time we make a new commit, and we do our best to keep it working.

Furthermore, the left menu gives additional instructions on how to use the low-level functions.