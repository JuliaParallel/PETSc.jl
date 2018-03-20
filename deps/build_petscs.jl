# build all three version of Petsc

petsc_name = "petsc-3.6.0"
fmt = ".tar.gz"

# download Petsc if needed
file_name = string(petsc_name, fmt)

const buildopts_fname = "buildopts.conf"

# load build options from last build
const fpath = joinpath(dirname(@__FILE__), buildopts_fname)
if isfile(fpath)
  println("USING build options in file ", fpath)
  envvars = readdlm(fpath)
  for entry = 1:size(envvars, 1)
    key = envvars[entry, 1]
    val =  envvars[entry, 2]
    if haskey(ENV, key)
      println("WARNING: build option $key cannot be specified as both environmental variable and buildopts.conf option, ignoring value from file...")
    else
      ENV[key] = val
    end
  end
end


# figure out what to build
const build_names = ["RealDouble", "RealSingle", "ComplexDouble"]
build_control = Dict() # dirname => whether or not to build it
arches = Dict() # dict of dirname => (PETSC_DIR, PETSC_ARCH)
build_any = false  # whether or not to build anything
have_petsc = Array{Bool}(3)  # whether or not each version of Petsc is usable
for (i, name) in enumerate(build_names)
  if haskey(ENV, string("JULIA_PETSC_", name, "_DIR")) || haskey(ENV, string("JULIA_PETSC_", name, "_ARCH"))
    if !(haskey(ENV, string("JULIA_PETSC_", name, "_DIR")) && haskey(ENV, string("JULIA_PETSC_", name, "_DIR")))
      error("Must have either both or neither DIR and ARCH for JULIA_PETSC_$name")
    else
      build_control[name] = false
      have_petsc[i] = true
      arches[name] = (ENV[string("JULIA_PETSC_",name, "_DIR")], ENV[string("JULIA_PETSC_", name, "_ARCH")])
    end
  elseif haskey(ENV, string("JULIA_PETSC_", name, "_NOBUILD"))
    build_control[name] = false
    have_petsc[i] = false
    arches[name] = ("NODIR", "42")
  else
    build_control[name] = true
    have_petsc[i] = true
    build_any = true
  end
end

# write options to file used for rebuild
suffixes_name = ["_DIR", "_ARCH", "_NOBUILD"]
suffixes_noname = ["_OPT", "_FLAGS"]
envvars = Dict{Any, Any}()
for (i, name) in enumerate(build_names)
  for suf in suffixes_name
    key_i = string("JULIA_PETSC_", name, suf)
    if haskey(ENV, key_i)
      envvars[key_i] = ENV[key_i]
    end
  end

  for suf in suffixes_noname
    key_i = string("JULIA_PETSC_", suf)
    if haskey(ENV, key_i)
      envvars[key_i] = ENV[key_i]
    end
  end
end

f = open(fpath, "w")
for (key, val) in envvars
  println(f, key, " ", val)
end
close(f)

# proceed with the build

if !isfile(file_name) && build_any
  println("DOWNLOADING $file_name")
  download("http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/$file_name", file_name)
end

const verbose=false # change to true for verbose compilation output

# flags the user cannot specify
forbidden_flags = ["--with-64-bit-indices", "with-scalar-type", "--with-precision"]
# optimize or debug build
if haskey(ENV, "JULIA_PETSC_OPT")
  println("Performing optimized build")
  opt_flags = "--with-debugging=0 COPTFLAGS='-O3 -march=native -mtune=native' -CXXOPTFLAGS='-O3 -march=native -mtune=native' FOPTFLAGS='-O3 -march=native -mtune=native'"
elseif haskey(ENV, "JULIA_PETSC_FLAGS")
  println("Performing user specified build")
  opt_flags = ENV["JULIA_PETSC_FLAGS"]
  for i in forbidden_flags
    if contains(opt_flags, i)
      throw(ArgumentError("JULIA_PETSC_FLAGS cannot contain $i"))
    end
  end
else
  println("Performing standard (debug) build")
  opt_flags = "--with-debugging=1"
end

# the PETSc configure script does not support Python 3, so we need to
# find Python 2 if we can:
if build_any
  pyconfigvar(python::AbstractString, var::AbstractString) = chomp(readstring(`$python -c "import distutils.sysconfig; print(distutils.sysconfig.get_config_var('$var'))"`))
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
end


for scalar_type in ("real", "complex"), precision in ("double", "single")
  # extract the tarfile into a new directory:
  name_i = ucfirst(scalar_type) * ucfirst(precision)
  name_i == "ComplexSingle" && continue # PETSc doesn't support this (weird!)
  !build_control[name_i] && continue # skip if the user said not to build this one
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
      for line in eachline(`$python configure --with-fc=0 --with-cxx=0 --with-64-bit-indices=true --with-scalar-type=$scalar_type --with-precision=$precision $opt_flags`)
        verbose && print("    $line")
        if ismatch(r"^\s*PETSC_ARCH:", line)
          PETSC_ARCH = chomp(split(line)[2])
        end
      end
      isempty(PETSC_ARCH) && error("couldn't determine PETSC_ARCH")
#      println("PETSC_ARCH = ", PETSC_ARCH)
      PETSC_DIR = pwd()  # because we cd into this directory
#      println("PETSC_DIR = ", PETSC_DIR)
      arches[name_i] = (PETSC_DIR, PETSC_ARCH)
      println("BUILDING PETSc $name_i...")
      for line in eachline(`make MAKE_NP=$(Sys.CPU_CORES) PETSC_DIR=$(pwd()) PETSC_ARCH=$PETSC_ARCH`)
        verbose && print("    $line")
      end
      println("TESTING PETSc $name_i...")
      for line in eachline(`make MAKE_NP=$(Sys.CPU_CORES) PETSC_DIR=$(pwd()) PETSC_ARCH=$PETSC_ARCH test`)
        verbose && print("    $line")
      end
    end
  end
end

# create deps.jl file with library locations

open("deps.jl", "w") do f
  for (i, name) in enumerate(build_names)
    if haskey(arches, name)
      PETSC_DIR, PETSC_ARCH = arches[name]
      path = abspath(PETSC_DIR, PETSC_ARCH, "lib", "libpetsc")
  #    path = abspath(name, petsc_name, arches[name], "lib", "libpetsc")
      println(f, "const petsc$name = \"", escape_string(path), "\"")
    end
  end
  println(f, "const PetscInt = Int64")
  println(f, "const have_petsc = [", have_petsc[1], " ", have_petsc[2], " ", have_petsc[3], "]")
end
