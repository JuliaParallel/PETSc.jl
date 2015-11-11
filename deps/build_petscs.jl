# build all three version of Petsc

petsc_name = "petsc-3.6.0"
fmt = ".tar.gz"

# download Petsc if needed
file_name = string(petsc_name, fmt)

if !isfile(file_name)
  println("DOWNLOADING $file_name")
  download("http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/$file_name", file_name)
end

const verbose=false # change to true for verbose compilation output

# the PETSc configure script does not support Python 3, so we need to
# find Python 2 if we can:
pyconfigvar(python::AbstractString, var::AbstractString) = chomp(readall(`$python -c "import distutils.sysconfig; print(distutils.sysconfig.get_config_var('$var'))"`))
pythonvers(python) = VersionNumber(pyconfigvar(python, "VERSION"))
pythonok(python) = try
  return pythonvers(python) < v"3"
catch
  return false
end
const PYTHON = get(ENV, "PYTHON", "python")
const python = pythonok(PYTHON) ? PYTHON : pythonok("python") ? "python" :
  pythonok("python2") ? "python2" : pythonok("python2.7") ? "python2.7" :
  error("could not find Python version < 3 for PETSc configure script")

arches = Dict() # dict of dirname => PETSC_ARCH

for scalar_type in ("real", "complex"), precision in ("double", "single")
  # extract the tarfile into a new directory:
  name_i = ucfirst(scalar_type) * ucfirst(precision)
  name_i == "ComplexSingle" && continue # PETSc doesn't support this (weird!)
  println("****************************************************************")
  if isdir(name_i)
    verbose && println("deleting existing installation $name_i")
    rm("$name_i", recursive=true)
  end
  mkdir(name_i)
  println("EXTRACTING PETSc to $name_i...")
  run(`tar xz -C $name_i -f $file_name`)

  cd(joinpath(name_i, petsc_name)) do
    withenv("PETSC_DIR"=>nothing, "PETSC_ARCH"=>nothing) do
      println("CONFIGURING PETSc $name_i...")
      PETSC_ARCH=""
      for line in eachline(`$python configure --with-64-bit-indices=true --with-scalar-type=$scalar_type --with-precision=$precision`)
        verbose && print("    $line")
        if ismatch(r"^\s*PETSC_ARCH:", line)
          PETSC_ARCH = chomp(split(line)[2])
        end
      end
      isempty(PETSC_ARCH) && error("couldn't determine PETSC_ARCH")
      arches[name_i] = PETSC_ARCH
      println("BUILDING PETSc $name_i...")
      for line in eachline(`make MAKE_NP=$CPU_CORES PETSC_DIR=$(pwd()) PETSC_ARCH=$PETSC_ARCH`)
        verbose && print("    $line")
      end
      println("TESTING PETSc $name_i...")
      for line in eachline(`make MAKE_NP=$CPU_CORES PETSC_DIR=$(pwd()) PETSC_ARCH=$PETSC_ARCH test`)
        verbose && print("    $line")
      end
    end
  end
end

# create deps.jl file with library locations
open("deps.jl", "w") do f
  for name in ("RealDouble", "RealSingle", "ComplexDouble")
    path = abspath(name, petsc_name, arches[name], "lib", "libpetsc")
    println(f, "const petsc$name = \"", escape_string(path), "\"")
  end
  println(f, "const PetscInt = Int64")
end
