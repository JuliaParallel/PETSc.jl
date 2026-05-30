import Gmsh: gmsh
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("test")
rect = gmsh.model.occ.addRectangle(-1.0, -1.0, 0.0, 2.0, 2.0)
gmsh.model.occ.synchronize()
bnd = gmsh.model.getBoundary([(2, rect)], false, false, false)
outer_tags = [abs(t) for (d, t) in bnd]
gmsh.model.addPhysicalGroup(1, outer_tags, 1, "marker")
gmsh.model.addPhysicalGroup(2, [rect], 1, "matrix")
gmsh.model.mesh.setSize(gmsh.model.getEntities(0), 0.5)
gmsh.model.mesh.generate(2)
gmsh.write("/tmp/test_labels.msh")
gmsh.finalize()
println("Mesh written")

using MPI, PETSc
petsclib = PETSc.getlib(PetscScalar = Float64)
MPI.Initialized() || MPI.Init()
PETSc.initialize(petsclib)

# Load mesh directly
dm = LibPETSc.DMPlexCreateFromFile(petsclib, MPI.COMM_WORLD, "/tmp/test_labels.msh", "test", LibPETSc.PetscBool(true))

# Check which labels exist
for lname in ("marker", "Face Sets", "Cell Sets", "boundary", "matrix")
    lbl = PETSc.getlabel(dm, lname)
    exists = convert(Ptr{Cvoid}, lbl) != C_NULL
    println("Label '$lname': ", exists ? "EXISTS" : "not found")
end

println("Done")

