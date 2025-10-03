# This generates Julia bindings for the PETSc library
#
# It uses the getAPI.py python code
#
# It replaces the Clang.jl infrastructure and re-uses the python binding tools.
# Should be run from petsc/config/utils/

using PythonCall
import Base: contains, String

String(x::Py) =  pyconvert(String, x)
contains(a::Py, b::String) = occursin(b, String(a))
isopaque(a::Py) = pyconvert(Bool,a.opaque)

# directory where PETSc is installed on your system:
petsc_dir = "/Users/kausb/Downloads/petsc"
start_dir = pwd()

# Set the Python path to include the directory of getAPI.py
cd(petsc_dir*"/config/utils")
PythonCall.pyimport("sys").path.append(".")

# Import the getAPI module
getAPI = pyimport("getAPI")

# Run the getAPI python function
curdir = pwd()
cd("../../")
classes, enums, senums, typedefs, structs, funcs, files, mansecs, submansecs = getAPI.getAPI()
cd(curdir)

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
    return type_str, name_str
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
    type_str, name_str = "(",""
    for (i,arg) in enumerate(a.arguments)
        if i>1
            type_str *= ", "
            name_str *= ", "
        end
        type, name = ccall_arg_string(arg)

        type_str *= type
        name_str *= name    
    end
    type_str *= ")"
    return type_str, name_str
end



"""
    write_enum(enum_val::Py, io = stdout)
Writes an enum to either file or REPL
"""
function write_enum(enum_val::Py, io = stdout)

    println(io,"@enum $(enum_val.name) begin")
    for (i,name) in enumerate(enum_val.values)
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
    ccall_type, ccall_name = julia_ccall_header(funcs_val)
    
    println(io,"@for_petsc $(funcs_val.name)(::\$UnionPetscLib, $julia_fct_str)");
    println(io,"    @chk ccall(");
    println(io,"              (:$(name), \$petsc_library),");
    println(io,"              PetscErrorCode,");
    println(io,"              $ccall_type,");
    println(io,"              $ccall_name,");
    println(io,"    )");
    println(io,"end \n");

    return nothing
end

function write_struct(struct_val::Py, io = stdout)

    println(io,"mutable struct $(struct_val.name)")

    for (i,name) in enumerate(struct_val.records)
        if  contains(name,"=")
            println(io,"    $(name)")
        else
            println(io,"    $(name) = $(i-1)")
        end
    end
    println(io,"end \n")
    return nothing
end


function write_keys_to_file(filename::String, start_dir::String, object::Py, fn::Function)
    open(start_dir*filename, "w") do file
        # Call the write_enum function and pass the file as the io argument
        for val in object.keys()
            fn(object[val], file)
        end
    end
end

function write_skeys_to_file(filename::String, start_dir::String, object::Py)
    open(start_dir*filename, "w") do io
        # Call the write_enum function and pass the file as the io argument
        println(io, "# not quite sure yet how to deal with this")
        for val in object.keys()
            println(io, "$(String(object[val].name))=Ptr{Cchar}")
        end
    end
end





#write_keys_to_file("enums_wrappers.jl",  start_dir,  enums, write_enum)  # Write enums to file
#write_skeys_to_file("senums_wrappers.jl",start_dir, senums)              # Write string enums to file



funcs_val = classes["KSP"].functions["KSPBuildSolution"]
