# Script to wrap PETSc header files with Julia code

using Clang.cindex
using Clang.wrap_c
using Compat

import Clang.wrap_c.repr_jl

for (ts,arg) in ((:"RealSingle",:"PETSC_USE_REAL_SINGLE"),(:"RealDouble",:"PETSC_USE_REAL_DOUBLE"),(:"ComplexDouble",:"PETSC_USE_REAL_DOUBLE"))
    @eval begin
        str = string($ts)
        petsc_libname = :petscComplexDouble
        if $ts == "RealSingle"
          petsc_libname = :petscRealSingle
        elseif $ts == "RealDouble"
          petsc_libname = :petscRealDouble
        end
        println(petsc_libname)
        include("$str.jl")
        include("rewriter.jl")
        PETSC_INCLUDE = "../../deps/$str/petsc-3.6.0/include"
        petsc_header = [joinpath(PETSC_INCLUDE, "petsc.h")]
        h1 = "/usr/lib/gcc/x86_64-linux-gnu/4.8/include"  # get some C datatype definitions like size_t
        h2 = joinpath(PETSC_INCLUDE, "petscsys.h")
        h3 = "/home/kshyatt/mpi/openmpi/include"
        # Set up include paths
        clang_includes = ASCIIString[]
        push!(clang_includes, h2)
        push!(clang_includes, PETSC_INCLUDE)
        push!(clang_includes, h1)
        push!(clang_includes, h3)


        println("clang_includes = ", clang_includes)

        # Clang arguments
        clang_extraargs = [ "-std=c99", "-D", "__STDC_LIMIT_MACROS", "-D", "__STDC_CONSTANT_MACROS", "-D", $arg ]

        # Callback to test if a header should actually be wrapped (for exclusion)
        function wrap_header(top_hdr::ASCIIString, cursor_header::ASCIIString)
          return startswith(dirname(cursor_header), PETSC_INCLUDE)
        end

        lib_file(hdr::ASCIIString) = "petsc"
        output_file(hdr::ASCIIString) = "PETSc$str.jl"

        function wrap_cursor(name::ASCIIString, cursor)
          #println("name = ", name)
          exc = true
          if length(name) > 0
            if name[1] == '_'
              exc = false
            end
          end
          return exc
        end

        const wc = wrap_c.init(;
                    headers = petsc_header,
                    output_file = "libPETSc_h.jl",
                    common_file = "libPETSc_common$str.jl",
                    clang_includes      = clang_includes,
                    clang_args          = clang_extraargs,
                    header_wrapped      = wrap_header,
                    header_library      = lib_file,
                    header_outputfile   = output_file,
                    cursor_wrapped      = wrap_cursor,
                    clang_diagnostics = true,
                    rewriter = petsc_rewriter)

        run(wc)
    end
end
