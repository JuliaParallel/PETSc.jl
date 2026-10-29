# 1D staggered FD example
#
# Solves  u' = p,  p' = c  on [0, 1] with u(0) = a and u(1) = b, with u on the
# vertices and p on the elements of a DMStag, and compares with the analytical
# solution
#
#   p(x) = b - a - c/2 + c x,   u(x) = a + (b - a - c/2) x + c x²/2.
#
# The discretization is exact for this quadratic, so the error is roundoff.
#
# Run it on one rank. The matrix has a zero diagonal on the vertex rows, so the
# default per-rank ILU fails with a zero pivot on more ranks; the assembly itself
# is written with global element numbers and works on any number.
using PETSc, MPI, LinearAlgebra

MPI.Initialized() || MPI.Init()

petsclib = PETSc.petsclibs[1]
PETSc.initialize(petsclib)
PetscScalar = PETSc.scalartype(petsclib)

comm = MPI.COMM_WORLD

nx   = 10
x0   = 0.0
xend = 1.0

# One DOF per vertex (u) and one per element (p)
bnd = (PETSc.DM_BOUNDARY_NONE,)
dm  = PETSc.DMStag(petsclib, comm, bnd, (nx,), (1, 1);
                   stencil_width = 1, stencil_type = :box)
PETSc.set_uniform_coordinates!(dm, (x0,), (xend,))

a = 1.0; b = 2.0; c = 1.0
if bnd[1] == PETSc.DM_BOUNDARY_PERIODIC
    b = a
    c = 0.0
end

# In 1D the vertices are the faces of the elements
vertex       = PETSc.face_location(dm, 1)           # the left vertex of an element
right_vertex = PETSc.LibPETSc.DMSTAG_RIGHT
element      = PETSc.element_location(dm)

corners = PETSc.corners(dm)
owned   = corners.lower[1]:corners.upper[1]
# The last rank also owns the vertex at the right end, which belongs to a
# partial element past the last one
owned_vertices = corners.lower[1]:(corners.upper[1] + corners.nextra[1])

# ── Reference solution ───────────────────────────────────────────────────────
xa       = PETSc.global_vec(dm)
xa_local = PETSc.local_vec(dm)
PETSc.with_product_coordinates(dm) do X
    PETSc.with_field_views!(dm, xa_local; fields = (vertex => 0, element => 0), read = false) do (u, p)
        for i in owned_vertices
            u[i] = a + (b - a - c / 2) * X[i, 1] + c / 2 * X[i, 1]^2
        end
        for i in owned
            p[i] = b - a - c / 2 + c * X[i, 2]
        end
    end
end
PETSc.local_to_global!(xa, dm, xa_local)

# Element width
h = PetscScalar((xend - x0) / nx)

# ── Right-hand side ──────────────────────────────────────────────────────────
rhs = PETSc.global_vec(dm)
for i in owned
    left  = PETSc.stencil(dm, vertex, (i,))
    elem  = PETSc.stencil(dm, element, (i,))
    # u(0) = a at the first vertex, p' = c at the others; p - u' = 0 in the element
    PETSc.set_values!(rhs, dm, [left, elem], PetscScalar[i == 1 ? a : c, 0])
end
# u(1) = b at the right end, which the element nx owns
nx in owned && PETSc.set_values!(rhs, dm, [PETSc.stencil(dm, right_vertex, (nx,))], PetscScalar[b])
PETSc.assemble!(rhs)

# ── Matrix ───────────────────────────────────────────────────────────────────
A = PETSc.PetscMat(dm)
for i in owned
    left  = PETSc.stencil(dm, vertex, (i,))
    elem  = PETSc.stencil(dm, element, (i,))
    right = PETSc.stencil(dm, right_vertex, (i,))

    if i == 1
        # u(0) = a
        PETSc.set_values!(A, dm, [left], [left], PetscScalar[1])
    else
        # (p_i - p_{i-1}) / h = c at the left vertex of element i
        prev = PETSc.stencil(dm, element, (i - 1,))
        PETSc.set_values!(A, dm, [left], [elem, prev, left], PetscScalar[1 / h, -1 / h, 0])
    end
    if i == nx
        # u(1) = b
        PETSc.set_values!(A, dm, [right], [right], PetscScalar[1])
    end

    # p_i - (u_right - u_left) / h = 0 in the element
    PETSc.set_values!(A, dm, [elem], [right, left, elem], PetscScalar[-1 / h, 1 / h, 1])
end
PETSc.assemble!(A)

# ── Solve and compare ────────────────────────────────────────────────────────
x   = PETSc.global_vec(dm)
ksp = PETSc.KSP(A)
PETSc.solve!(x, ksp, rhs)

err         = xa - x
xa_norm     = norm(xa)
error_norm  = norm(err)
errRel      = error_norm / xa_norm

print("Error (abs): ", error_norm, "\nError (rel): ", errRel, "\n")

foreach(PETSc.destroy!, (err, ksp, A, rhs, x, xa_local, xa, dm))
PETSc.finalize(petsclib)

nothing
