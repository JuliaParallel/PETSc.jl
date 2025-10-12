# This generates Julia bindings for the PETSc library
#
# It uses the getAPI.py python code
#
# It replaces the Clang.jl infrastructure and re-uses the python binding tools.
# Should be run from petsc/config/utils/; please specify the petsc directory @ the beginning

using PythonCall
import Base: contains, String

String(x::Py) =  pyconvert(String, x)
contains(a::Py, b::String) = occursin(b, String(a))
isopaque(a::Py) = pyconvert(Bool,a.opaque)

# directory where PETSc is installed on your system:
petsc_dir = "/Users/kausb/Downloads/petsc"
start_dir = pwd()

if  !@isdefined classes
    # for some reason, if we run this twice it returns empty classes

    # Set the Python path to include the directory of getAPI.py
    cd(petsc_dir*"/config/utils")
    PythonCall.pyimport("sys").path.append(".")

    # Import the getAPI module
    getAPI = pyimport("getAPI")

    # Run the getAPI python function
    curdir = pwd()
    cd("../../")
    classes, enums, senums, typedefs, structs, funcs, files, mansecs, submansecs = getAPI.getAPI()
    cd(start_dir)
end


"""
Retrieves a record from a python struct and reshapes it in a julia way
"""
function struct_record(b::Py)
    str  = replace(String(b.type), "," => " ")  # deal with multiple args
    if contains(str,"[")
        # We need to make an NTuple out of the parameter
        le = split(replace(str,"]"=>""),"[")[2] # length of parameter

        str  = split(replace(str,"]"=>""),"[")[1]
        spl  = strip.(split(str))
        spl[1] = "NTuple{$le, $(spl[1])}"
    elseif contains(str,"*")
        # We need to make a Tuple out of the parameter
        str  = replace(replace(str,")"=>"") ,"("=>"")

        spl  = strip.(split(str))
        spl[1] = "Tuple{$(join(spl[2:end],", "))}"
    else


        spl  = strip.(split(str))
    end

    type, args = spl[1],  spl[2:end]

    return args.*"::$type"      # add type to each argument
end

# contains a function argument
struct f_args
    name::String
    typename::String
    array::Bool
    cons::Bool
    optional::Bool
    stars::Int64
    stringlen::Bool
end

# returns 
func_args(a::Py) = f_args(String(a.name), String(a.typename), Bool(a.array), Bool(a.const),Bool(a.optional), pyconvert(Int64, a.stars), Bool(a.stringlen))

"""
    func_arg_string(a::f_args)
returns a string representation of the function argument, to be used in the julia type signature
"""
func_arg_string(a::f_args) = "$(a.name)::$(a.typename)"
func_arg_string(a::Py) = func_arg_string(func_args(a))

function ccall_arg_string(a::f_args) 
    type_str = a.typename
    if a.stars==1
        type_str = "Ptr{$(type_str)}"
    end
    name_str = a.name
    type_str = replace_function_string(type_str)
    name_str = replace_function_string(name_str)


    return type_str, name_str
end

function replace_function_string(type_str::String)
    type_str = replace(type_str, "void"=>"Cvoid",
                                 "char"=>"Cchar",
                                 "function"=>"fnc")
    return type_str
end

ccall_arg_string(a::Py) = ccall_arg_string(func_args(a))

# returns a string with all variables of the function, to be used to create a julia header
function julia_function_header(a::Py)
    str = ""
    for (i,arg) in enumerate(a.arguments)
        if i>1
            str *= ", "
        end
        str *= func_arg_string(arg)
    end
    return str
end

# returns 2 strings needed for the ccall routine
function julia_ccall_header(a::Py)
    n_arg = length(a.arguments)
    type_str, name_str = "(",""
    for (i,arg) in enumerate(a.arguments)
        if i>1
            type_str *= ", "
            name_str *= ", "
        end

        type, name = ccall_arg_string(arg)

        type_str *= type
        name_str *= name    
        if n_arg==1
            type_str *= ","
        end

    end
    type_str *= ")"
    return type_str, name_str, n_arg
end



"""
    write_enum(enum_val::Py, io = stdout)
Writes an enum to either file or REPL
"""
function write_enum(enum_val::Py, io = stdout)

    println(io,"@enum $(enum_val.name) begin")
    for (i,name) in enumerate(enum_val.values)
        if contains(String(name),"DEPRECATED")
            break
        end
        if  contains(name,"=")
            println(io,"    $(name)")
        else
            println(io,"    $(name) = $(i-1)")
        end
    end
    println(io,"end \n")
    return nothing
end

"""
Writes a function to the specified IO stream
"""
function write_funcs(funcs_val::Py, io = stdout)
    if contains(funcs_val.name,"_") || isopaque(funcs_val)
        return nothing
    end

    # print function
    name          = String(funcs_val.name)
    julia_fct_str = julia_function_header(funcs_val)
    julia_fct_str = replace_function_string(julia_fct_str)
    ccall_type, ccall_name, n_arg = julia_ccall_header(funcs_val)
    
    println(io,"\"\"\"");
    println(io,"\t$(funcs_val.name)($julia_fct_str) ");
    println(io,"\"\"\"");
    println(io,"function $(funcs_val.name)($julia_fct_str) end");
    println(io,"");
    println(io,"@for_petsc function $(funcs_val.name)(::\$UnionPetscLib, $julia_fct_str)");
    println(io,"    @chk ccall(");
    println(io,"               (:$(name), \$petsc_library),");
    println(io,"               PetscErrorCode,");
    println(io,"               $ccall_type,");
    if n_arg>0
    println(io,"               $ccall_name,");
    end
    println(io,"              )");
    println(io,"end \n");

    return nothing
end

function write_struct(struct_val::Py, io = stdout)

    println(io,"mutable struct $(struct_val.name)")

    for (i,s) in enumerate(struct_val.records)
        str = extract_struct_entry(s)
        if !isnothing(str)   
            println(io,"    $str")
        else
            # stop when we find the firtst #ifdef
            break
        end
    end
    println(io,"    $(struct_val.name)() = new()")
    println(io,"end \n")
    return nothing
end

function extract_struct_entry(s::Py)
    str = String(s.type)
    str = replace(str,"const "=>"")
    type = split(str)[1]
    name = split(str)[2:end]
    if contains(str,"[")
        name = split(replace(str,"]"=>""),"[")
        le   = name[2] # length of parameter
        name = name[1] 
        type = split(name)[1]
        name = split(name)[2:end]
        type = "NTuple{$le, $type}"
    end
    if length(name)==1
        name = name[1]
    end

    if contains(type,"#if") || contains(type,"#else") || contains(type,"#endif")
        return
    end
    if any(contains.(name,"*"))
        type = "Ptr{$type}"
        name = replace.(name,"*"=>"","("=>"",")"=>"")
    end
    type = replace_types(type)
    name = replace(name,"function"=>"_function")

    if any(contains.(name,"["))
        name = split(replace.(name,"]"=>""),"[")
        le   = name[2] # length of parameter
        name = name[1] 
        type = "NTuple{$le, $type}"
    end
    str_out = ""
    if isa(name, Vector)
        name = replace.(name,","=>"")
        for na in name
            str_out *= "$na::$type \n    "
        end
    else
        str_out *= "$name::$type"
    end
    return str_out
end

replace_types(type::AbstractString) = replace(type, "size_t"=>"Csize_t","ptrdiff_t"=>"Cptrdiff_t",
                                            "short"=>"Cshort","int32_t"=>"Int32",
                                            "float"=>"Cfloat",
                                            "char"=>"Cchar", "void"=>"Cvoid","int"=>"Cint",
                                            "double"=>"Cdouble");


function write_typedefs(typedef_val::Py, io = stdout)
    type = String(typedef_val.value)
    type = replace_types(type)

    println(io,"const $(typedef_val.name) = $type")
    return nothing
end


function write_keys_to_file(filename::String, start_dir::String, object::Py, fn::Function)
    open(joinpath(start_dir, filename), "w") do file
        # Call the write_enum function and pass the file as the io argument
        for val in object.keys()
            fn(object[val], file)
        end
    end
end

function write_skeys_to_file(filename::String, start_dir::String, object::Py)
    open(joinpath(start_dir, filename), "w") do io
        # Call the write_enum function and pass the file as the io argument
        println(io, "# not quite sure yet how to deal with this")
        for val in object.keys()
            println(io, "$(String(object[val].name))=Ptr{Cchar}")
        end
    end
end

function write_functions_to_file(filename::String, start_dir::String, classes::Py, function_name::String)
    open(joinpath(start_dir, filename), "w") do file
        # Call the write_enum function and pass the file as the io argument
        for f in classes[function_name].functions
            name = String(f)
            @info name
            write_funcs(classes[function_name].functions[name], file)
            #write_funcs(classes[function_name].functions[String(f)])
            
        end
    end
end

function write_structs_to_file(filename::String, start_dir::String, structs::Py)
    open(joinpath(start_dir, filename), "w") do file
        # Call the write_enum function and pass the file as the io argument
        for struct_val in structs
            write_struct(structs[String(struct_val)], file)
        end
    end
end

function write_typedefs_to_file(filename::String, start_dir::String, typedefs::Py)
    open(joinpath(start_dir, filename), "w") do file
        # Call the write_enum function and pass the file as the io argument
        for typedef_val in typedefs
            write_typedefs(typedefs[String(typedef_val)], file)
        end
    end 
end



#write_keys_to_file("enums_wrappers.jl",  start_dir,  enums, write_enum)  # Write enums to file
#write_skeys_to_file("senums_wrappers.jl",start_dir, senums)              # Write string enums to file
#write_structs_to_file("struct_wrappers.jl", start_dir, structs)          # Write all structs to file
#write_typedefs_to_file("typedefs_wrappers.jl", start_dir, typedefs)      # Write all typedefs to file

# Write KSP functions to file (this should be expanded to all other classes)
write_functions_to_file("KSP_wrappers.jl",start_dir, classes, "KSP")     
    

# open question:
#   why are structs such as _p_PetscSF not included in the python structs above, even though 
#   they are define in petscsftypes.h?  


# if all is well, we should be able to say:
#