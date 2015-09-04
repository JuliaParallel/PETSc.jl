#
# Example script to wrap PETSc
#

using Clang.cindex
using Clang.wrap_c
using Compat


import Clang.wrap_c.repr_jl

include("rewriter.jl")
#=
function repr_jl(ptr::cindex.Pointer)
  ptee = pointee_type(ptr)
  ex1 = Expr(:curly, :Ptr, repr_jl(ptee))
  ex2 = Expr(:curly, :AbstractArray, repr_jl(ptee))
  Expr(:curly, Union, ex1, ex2)
end
=#
PETSC_INCLUDE = "/home/jared/build/petsc-3.6.0/include"
#PETSC_INCLUDE = "/home/jared/.julia/v0.4/PETSc/deps/petsc-3.6.0/arch-linux2-c-opt/include"
MPI_INCLUDE = "/usr/include/mpich"

petsc_header = [joinpath(PETSC_INCLUDE, "petscksp.h")]
#petsc_header = [joinpath(PETSC_INCLUDE, "petsc.h")]
h1 = "/usr/include"
h2 = joinpath(PETSC_INCLUDE, "petsc/mpiuni")
h3 = joinpath(PETSC_INCLUDE, "petsc/private")
h4 = joinpath(PETSC_INCLUDE, "petsc/private/kernels")

# Set up include paths
clang_includes = ASCIIString[]
push!(clang_includes, PETSC_INCLUDE)
#push!(clang_includes, MPI_INCLUDE)
#push!(clang_includes, h1)
#push!(clang_includes, h2)
#push!(clang_includes, h3)
#push!(clang_includes, h4)


println("clang_includes = ", clang_includes)
# Clang arguments
#clang_extraargs = ["-v"]
 clang_extraargs = [ "-std=c99", "-D", "__STDC_LIMIT_MACROS", "-D", "__STDC_CONSTANT_MACROS"]

# Callback to test if a header should actually be wrapped (for exclusion)
function wrap_header(top_hdr::ASCIIString, cursor_header::ASCIIString)
    return startswith(dirname(cursor_header), PETSC_INCLUDE)
end

lib_file(hdr::ASCIIString) = "petsc"
output_file(hdr::ASCIIString) = "PETSc.jl"

function wrap_cursor(name::ASCIIString, cursor)
    println("name = ", name)
    exc = true
    if length(name) > 0
      if name[1] == '_'
        exc = false
      end
#    exc = false
#    exc |= contains(name, "MPI")
    end
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
                        clang_diagnostics = true
                        rewriter = petsc_rewriter)

run(wc)
