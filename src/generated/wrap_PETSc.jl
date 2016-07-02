# Script to wrap PETSc header files with Julia code
#=
using Clang_orig.cindex
using Clang_orig.wrap_c
using Compat

import Clang_orig.wrap_c.repr_jl
=#

using Clang.cindex
using Clang.wrap_c
using Compat

import Clang.wrap_c.repr_jl

include("rewriter.jl")
PETSC_INCLUDE = "../../deps/RealDouble/petsc-3.6.0/include"
ARCH_INCLUDE = "../../deps/RealDouble/petsc-3.6.0/arch-linux2-c-debug/include"  # TEMPORARY
petsc_header = [joinpath(PETSC_INCLUDE, "petsc.h")]
h1 = "/usr/lib/gcc/x86_64-linux-gnu/4.8/include"  # get some C datatype definitions like size_t
h2 = joinpath(PETSC_INCLUDE, "petscsys.h")
h3 = joinpath(ARCH_INCLUDE, "petscconf.h")
# Set up include paths
clang_includes = ASCIIString[]
push!(clang_includes, PETSC_INCLUDE)
push!(clang_includes, ARCH_INCLUDE)
push!(clang_includes, h1)
#push!(clang_includes, h2)
#push!(clang_includes, h3)


println("clang_includes = ", clang_includes)

# Clang arguments
clang_extraargs = [ "-std=c99", "-D", "__STDC_LIMIT_MACROS", "-D", "__STDC_CONSTANT_MACROS", "-D", "PETSC_USE_REAL_SINGLE" ]

# Callback to test if a header should actually be wrapped (for exclusion)
function wrap_header(top_hdr::ASCIIString, cursor_header::ASCIIString)
  return startswith(dirname(cursor_header), PETSC_INCLUDE)
end

lib_file(hdr::ASCIIString) = "petsc"
output_file(hdr::ASCIIString) = "PETSc.jl"

function wrap_cursor(name::ASCIIString, cursor)
  println("\nwrapping name = ", name)
  exc = true
#= debug: wrap everything
  if length(name) > 0
    if name[1] == '_'
      exc = false
    end
  end
=#
  return exc
end

const wc = wrap_c.init(;
                        headers = petsc_header,
                        output_file = "libPETSc_h.jl",
                        common_file = "libPETSc_common.jl",
                        clang_includes      = clang_includes,
                        clang_args          = clang_extraargs,
                        header_wrapped      = wrap_header,
                        header_library      = lib_file,
                        header_outputfile   = output_file,
                        cursor_wrapped      = wrap_cursor,
                        clang_diagnostics = true,
                        rewriter = petsc_rewriter)

run(wc)
