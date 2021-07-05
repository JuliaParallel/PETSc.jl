# Getting started


### 1a. Installation using pre-build libraries 
The easiest way to install the package is: 
```julia
julia> ]
(@v1.6) pkg> add https://github.com/JuliaParallel/PETSc.jl
```
which will install a pre-build PETSc library (`PETSc_jll`) as well as `MPI.jl` on your system. This will work both in serial and in parallel on your machine.

### 1b. Installation using pre-build libraries 
On many high-performance clusters, you will have to use the provided `MPI` installation for that cluster and the default download above will not be sufficient. Alternatively, you may be interested in a PETSc installation that comes with additional external packages. Ensure that this PETSc installation is compiled as a dynamic (and not a static) library, after which you need to set the environmental variable `JULIA_PETSC_LIBRARY` to link to your PETSc installation: 

```
$export JULIA_PETSC_LIBRARY = /path/to/your/petsc/installation:
```

Now rebuild the package:
```julia
julia> ]
pkg> build PETSc
```

### 2. Linear example
To be added 

### 3. Nonlinear example
Let's solve the coupled system of nonlinear equations: 
```math
\begin{aligned}
x^2 + x y  &= 3 \\
x y + y^2  &= 6
\end{aligned}
```
for ``x`` and ``y``, which can be written in terms of a residual vector ``f``:
```math
f = \binom{ x^2 + x y  - 3} {x y + y^2  - 6}
```

In order to solve this, we need to provide a function that computes the residual vector `f`:
```julia
julia> using PETSc, MPI
julia> function F!(fx, x)
         fx[1] = x[1]^2 + x[1]*x[2] - 3
         fx[2] = x[1]*x[2] + x[2]^2 - 6
       end
```

In addition, we need to provide the Jacobian:
```math
J = 
\begin{pmatrix}
\frac{\partial f_1}{ \partial x} & \frac{\partial f_1}{ \partial y}  \\
\frac{\partial f_2}{ \partial x} & \frac{\partial f_2}{ \partial y}  \\
\end{pmatrix}
= 
\begin{pmatrix}
2x + y & x  \\
y & x + 2y  \\
\end{pmatrix}
```
As julia function, this is coded as:
```julia
julia> J = zeros(2,2)
julia> PJ = PETSc.MatSeqDense(J)
julia> function updateJ!(x, args...)
            J[1,1] = 2x[1] + x[2]
            J[1,2] = x[1]
            J[2,1] = x[2]
            J[2,2] = x[1] + 2x[2]
        end
```
In order to solve this using the PETSc nonlinear equation solvers, you first define the `SNES` solver together with the jacobian and residual functions as 
```julia
julia> S = PETSc.SNES{Float64}(MPI.COMM_SELF; ksp_rtol=1e-4, pc_type="none")
julia> PETSc.setfunction!(S, F!, PETSc.VecSeq(zeros(2)))
julia> PETSc.setjacobian!(S, updateJ!, PJ, PJ)
```

You can solve this as:
```julia
julia> PETSc.solve!([2.0,3.0], S)
2-element Vector{Float64}:
 1.000000003259629
 1.999999998137355
```
which indeed recovers the analytical solution ``(x=1, y=2)``.