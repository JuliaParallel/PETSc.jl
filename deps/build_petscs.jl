# build all three version of Petsc

petsc_name = "petsc-3.6.0"
fmt = ".tar.gz"

# download Petsc if needed
file_name = string(petsc_name, fmt)

if !isfile(file_name)
  println("Downloading $file_name")
  download("http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/$file_name")
end

start_dir = pwd()
if false

# delete existing installation
name_i = "RealDouble"
if isdir(name_i)
  println("deleting existing installation $name_i")
  rm("$name_i", recursive=true)
end
# create new one
mkdir(name_i)
cd(name_i)
run(`tar xfvz ../$file_name`)
cd(petsc_name)
run(`../../install_petsc.sh --with-scalar-type=real --with-precision=double`)

cd(start_dir)

# delete existing installation
name_i = "RealSingle"
if isdir(name_i)
  println("deleting existing installation $name_i")
  rm("$name_i", recursive=true)
end
# create new one
mkdir(name_i)
cd(name_i)
run(`tar xfvz ../$file_name`)
cd(petsc_name)
run(`../../install_petsc.sh --with-scalar-type=real --with-precision=single`)


cd(start_dir)

# delete existing installation
name_i = "ComplexDouble"
if isdir(name_i)
  println("deleting existing installation $name_i")
  rm("$name_i", recursive=true)
end
# create new one
mkdir(name_i)
cd(name_i)
run(`tar xfvz ../$file_name`)
cd(petsc_name)
run(`../../install_petsc.sh --with-scalar-type=complex --with-precision=double`)


cd(start_dir)

end
# create library locations

  # get PETSC_ARCH
  f = open("./RealDouble/petsc_evars")
  evars = readlines(f)
  close(f)
  petsc_arch = evars[2][1:(end-1)]
  petscRealDouble = joinpath(pwd(), "RealDouble/$petsc_name/$petsc_arch/lib/libpetsc")
  petscRealSingle = joinpath(pwd(), "RealSingle/$petsc_name/$petsc_arch/lib/libpetsc")
  petscComplexDouble = joinpath(pwd(), "ComplexDouble/$petsc_name/$petsc_arch/lib/libpetsc")




  f = open("lib_locations.jl", "w")
  write(f, "global const petscRealDouble = \"$petscRealDouble\"\n")
  write(f, "global const petscRealSingle = \"$petscRealSingle\"\n")
  write(f, "global const petscComplexDouble = \"$petscComplexDouble\"\n")
  close(f)


  if isfile("../src/generated/lib_locations.jl")
    rm("../src/generated/lib_locations.jl")
  end

  cp("lib_locations.jl", "../src/generated/lib_locations.jl")
