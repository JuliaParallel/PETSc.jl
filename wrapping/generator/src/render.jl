# Rendering of one wrapper function (docstring, stub, @for_petsc method)

function doc_header(r::Rules, args::Vector{FArg}, fname::String)
    ins, outs, outs_doc = String[], String[], String[]
    for a in args
        if a.output
            push!(outs, a.name)
            push!(outs_doc, "$(a.name)::$(a.typename)")
        else
            push!(ins, "$(a.name)::$(abstract_arg_type(r, a.typename))")
        end
    end
    str_in = join(ins, ", ")
    str_out = join(outs, ",")
    call = isempty(ins) ? "$fname(petsclib::PetscLibType)" : "$fname(petsclib::PetscLibType,$str_in)"
    hdr = isempty(outs) ? call : "$(join(outs_doc, ",")) = $call"
    return hdr, str_in, str_out, length(ins), length(outs)
end

function ccall_header(args::Vector{FArg})
    types = join((a.ccall_str for a in args), ", ")
    names = join((a.name_ccall for a in args), ", ")
    tuple = length(args) == 1 ? "($types,)" : "($types)"
    return tuple, names
end

_titlecase(s) = isempty(s) ? s : uppercase(s[1:1]) * lowercase(s[2:end])

"""
    render_function(io, r, fn, args, doc_lines)

Writes the docstring, the untyped stub and the `@for_petsc` method for `fn`.
"""
function render_function(io::IO, r::Rules, fn::Fn, args::Vector{FArg}, doc_lines, mansec::String)
    name = fn.name
    hdr, str_in, str_out, num_in, num_out = doc_header(r, args, name)
    ctypes, cnames = ccall_header(args)
    println(io, "\"\"\"")
    println(io, "\t$hdr ")
    if doc_lines !== nothing
        for c in doc_lines
            println(io, replace(c, "\\" => "\\\\"))   # backslashes would be escape sequences in the docstring
        end
    end
    println(io, "")
    println(io, "# External Links")
    println(io, "\$(_doc_external(\"$mansec/$name\"))")
    println(io, "\"\"\"")
    if num_in > 0
        println(io, "function $name(petsclib::PetscLibType, $str_in) end")
    else
        println(io, "function $name(petsclib::PetscLibType) end")
    end
    println(io, "")
    if num_in > 0
        println(io, "@for_petsc function $name(petsclib::\$UnionPetscLib, $(dispatch(r, str_in)) )")
    else
        println(io, "@for_petsc function $name(petsclib::\$UnionPetscLib)")
    end
    for a in args
        isempty(a.init) || println(io, "\t$(dispatch(r, a.init))")
    end
    println(io, "")
    println(io, "    @chk ccall(")
    println(io, "               (:$name, \$petsc_library),")
    println(io, "               PetscErrorCode,")
    println(io, "               $(dispatch(r, ctypes)),")
    isempty(args) || println(io, "               $cnames,")
    println(io, "              )")
    println(io, "")
    for a in args
        isempty(a.extract) || println(io, "\t$(a.extract)")
    end
    if num_out > 0
        println(io, "\n\treturn $str_out")
    else
        println(io, "\n\treturn nothing")
    end
    println(io, "end \n")
    return nothing
end
