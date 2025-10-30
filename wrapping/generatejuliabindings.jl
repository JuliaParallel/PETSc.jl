# This generates Julia bindings for the PETSc library
#
# It uses the getAPI.py python code.
# We also read the documentation of each of the functions.
#
# It replaces the Clang.jl infrastructure and re-uses the python binding tools.
# Should be run from petsc/config/utils/; please specify the petsc directory @ the beginning

using PythonCall
import Base: contains, String

include("find_doc_strings.jl")

String(x::Py) =  pyconvert(String, x)
contains(a::Py, b::String) = occursin(b, String(a))
isopaque(a::Py) = pyconvert(Bool,a.opaque)

# directory where PETSc is installed on your system:
petsc_dir = "/Users/kausb/Downloads/petsc"
start_dir = pwd()

const CustomTypes = ["PetscVec","PetscMat","PetscDM", "PetscKSP", "PetscSNES","PetscOptions"]

# Remove entry from vector of strings
remove_entry(x::Vector{String}, entry::String) = filter!(y -> y != entry, x)

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

"""
    mutable struct f_args
Contains PETSc function arguments, which are used 
"""
mutable struct f_args
    name::String            # variable name
    name_ccall::String      # name in ccall array (sometimes different from output name)
    typename::String        # type
    array::Bool             # is it an array?
    _const::Bool            # is it a const?
    optional::Bool          # optional?
    stars::Int64            # how many stars?
    isfunction::Bool        # is this a function?
    stringlen::Bool
    output::Bool             # is this an output argument?
    init_arg::String
    extract_arg::String
    ccall_str :: String
end

replace_types(type::AbstractString) = replace(type, 
                                            "size_t"=>"Csize_t",
                                            "ptrdiff_t"=>"Cptrdiff_t",
                                            "short"=>"Cshort",
                                            "int32_t"=>"Int32",
                                            "float"=>"Cfloat",
                                            "bool"=>"Bool",
                                            "char"=>"Cchar", 
                                            "void"=>"Cvoid",
                                            "double"=>"Cdouble",
                                            r"int"=>"Cint",
                                            "FILE"=>"Libc.FILE",

                                            # handle custom-defined types:
                                            r"\bVec\b"=>"PetscVec",
                                            r"\bMat\b"=>"PetscMat",
                                            r"\bDM\b"=>"PetscDM",
                                            r"\bKSP\b"=>"PetscKSP",
                                            r"\bSNES\b"=>"PetscSNES",
                                            );


replace_names(type_str::String)     =  replace(type_str, 
                                            "function"=>"fnc",
                                            "local"=>"loc");

# Since we use the @for_petsc macro to generate multiple dispatch versions of our functions,
# we have to add dollar signs to these types, before writing them     
replace_dispatch_types(type::AbstractString) = replace(type, 
                                            "PetscScalar"=>"\$PetscScalar",
                                            "PetscReal"=>"\$PetscReal",
                                            "PetscInt"=>"\$PetscInt",
                                            "PetscComplex"=>"\$PetscComplex",
                                            ); 

function init_extract_parameters(typename::String, name::String, function_name::String, name_ccall::String, isarray::Bool, isoutput::Bool, stars::Int64) 
    init_arg    = ""  
    extract_arg = ""  
    @show function_name, typename, name, isarray, isoutput, stars
    if !isarray && isoutput && typename in CustomTypes
        # a custom type is being initialized
        typename_c  = replace(typename,"Petsc"=>"C", count=1)
        name_ccall  = "$(name)_"
        init_arg    = "$name_ccall = Ref{$typename_c}()"  
        extract_arg = "$name = $typename($(name_ccall)[], petsclib)"  

    elseif !isarray && !isoutput && stars==1 && typename in CustomTypes
        # a custom type is being deleted most likely
        typename_c  = replace(typename,"Petsc"=>"C", count=1)
        name_ccall  = "$(name)_"
        init_arg    = "$name_ccall = Ref($(name).ptr)"  
        extract_arg = "$(name).ptr = C_NULL"  
 
    elseif isarray && isoutput && stars==1
        # array that'll be an output
        name_ccall  = "$(name)_"
        init_arg    = "$name_ccall = Ref{Ptr{$typename}}()"  
        # Note: we make the implicit assumption here that the Vec is always called 'x'
        extract_arg = "$name = unsafe_wrap(Array, $name_ccall[], VecGetLocalSize(petsclib, x); own = false)"  

    elseif isarray && isoutput && stars==0
        # array that'll be an output
        name_ccall  = "$(name)"
        
        # Note: we make the implicit assumption here that the length is called "ni"; that may not hold always
        # I don't see a way to automatically correct this, so has to be done manually, which
        # is why a comment is left
        init_arg    = "$name_ccall = Vector{$typename}(undef, ni);  # CHECK SIZE!!"  
        extract_arg = ""  
    
    elseif isarray && !isoutput && stars==1
        # array that'll is an input
        name_ccall  = "$(name)_"
        init_arg    = "$name_ccall = Ref(pointer($name))"  
        extract_arg = ""  
    elseif !isarray && isoutput && contains(typename,"Type")
        #scalar
        name_ccall  = "$(name)_"
        init_arg    = "$name_ccall = Ref{$typename}()"  
        extract_arg = "$name = unsafe_string($(name_ccall)[])" 

    elseif !isarray && isoutput
        #scalar
        name_ccall  = "$(name)_"
        init_arg    = "$name_ccall = Ref{$typename}()"  
        extract_arg = "$name = $(name_ccall)[]" 

    end

    return init_arg, extract_arg, name_ccall
end

function is_this_an_output_argument(function_name::String, variablename::String, typename::String, stars::Int64, isarray::Bool, output_vars::Vector{String}, input_vars::Vector{String})
    output = false
    if stars==1
        output = true
    end
    if !isempty(output_vars)
        if variablename in output_vars
            output = true
        end
    end
    if !isempty(input_vars)
        if variablename in input_vars
            output = false
        end
    end

    # so the function type is NOT a simple thing such as PetscScalar, PetscBool, PetscReal 
    # and the function name does NOT contain "Create"
    if output
        if  (!contains(function_name, "Create") && 
             !contains(function_name, "Duplicate") && 
             !contains(function_name, "Type") && 
             !(typename in ["PetscScalar", "PetscBool", "PetscReal","PetscComplex","PetscInt"]))  ||
            contains(function_name, "RestoreArray") ||
            contains(function_name, "Copy") 
            
            output = false
        end
    end

    return output
end


# returns a julia struct 
function func_args(a::Py, function_name::String, input_vars=String[], output_vars=String[]) 
    stars      = pyconvert(Int64, a.stars)
    typename   = String(a.typename)
    typename   = replace_types(typename)    # some scrambling necessary
    name       = String(a.name)
    name       = replace_names(name)        # some scrambling necessary
    name_ccall = name;
    isconst    = Bool(a.const)
    isarray    = Bool(a.array)
    isoptional = Bool(a.optional)
    isfunction = Bool(a.isfunction)
    
    # default setting for output arguments
    @show input_vars, output_vars
    output = is_this_an_output_argument(function_name, name, typename, stars, isarray, output_vars, input_vars)

    # this is related to creating custom julia structs
    typename_ccall = typename
    if typename in CustomTypes
        typename_ccall = replace(typename_ccall,"Petsc"=>"C", count=1)
    end
    
    if stars==1
        # a single value, so can always be output
        # depending on what type of output, we need different strategies here
        #output      = true
        init_arg, extract_arg, name_ccall = init_extract_parameters(typename, name, function_name, name_ccall, isarray, output, stars) 
        ccall_str   = "Ptr{$typename_ccall}"
    else
        init_arg, extract_arg, name_ccall = init_extract_parameters(typename, name, function_name, name_ccall, isarray, output, stars) 

        #output      = false
        #init_arg    = ""
        #extract_arg = ""
        ccall_str   = "$typename_ccall"
    end
    if isarray
        # in case parameters are arrays
       # init_arg, extract_arg, name_ccall = init_extract_parameters(typename, name, name_ccall, isarray, output) 
       if typename=="Cchar"  
           typename     = "String" 
       else
           typename     = "Vector{$typename}" 
       end

        ccall_str    = "Ptr{$ccall_str}"
    end

    return f_args(name, name_ccall, typename, isarray, isconst, isoptional, stars, isfunction, Bool(a.stringlen), output, init_arg, extract_arg, ccall_str)
end

# creates a 
function julia_function_doc_header(args::Vector{f_args}, function_name::String)
    str_in   = ""; num_in=0
    str_out  = "";  str_out_doc  = ""; num_out=0
    for arg in args
        if arg.output # output argument
            str_out     *= num_out > 0 ? "," : ""
            str_out_doc *= num_out > 0 ? "," : ""
            str_out_doc *= "$(arg.name)::$(arg.typename)"
            str_out     *= "$(arg.name)"
            num_out += 1 
        else
            str_in *= num_in > 0 ? ", " : ""
            str_in *= "$(arg.name)::$(arg.typename)"
            num_in += 1
        end
    end

    if num_out>0
        if num_in>0
            str = "$str_out_doc = $function_name(petsclib::PetscLibType,$str_in)"
        else
            str = "$str_out_doc = $function_name(petsclib::PetscLibType)"
        end
    else
        if num_in>0
            str = "$function_name(petsclib::PetscLibType,$str_in)"
        else
            str = "$function_name(petsclib::PetscLibType)"
        end
    end

    return str, str_in, str_out, num_in, num_out
end

# returns 2 strings needed for the ccall routine
function julia_ccall_header(args::Vector{f_args})
    type_str, name_str = "(",""
    for (i,arg) in enumerate(args)
        if i>1
            type_str *= ", "
            name_str *= ", "
        end
        type_str *= arg.ccall_str
        name_str *= arg.name_ccall
    end
    if length(args)==1
        type_str *= ",)"
    else
        type_str *= ")"
    end
    return type_str, name_str, length(args)
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

# Retrieves input and output arguments of a petsc function
function process_function_arguments(arguments::Py, function_name::String, input_vars=String[], output_vars=String[])
    function_arguments  = f_args[]
    for arg in arguments
        push!(function_arguments,  func_args(arg, function_name, input_vars, output_vars))
    end
    return function_arguments
end

function write_initialize(arguments::Vector{f_args}, io = stdout)
    for arg in arguments
        if !isempty(arg.init_arg)
            println(io,"\t$(replace_dispatch_types(arg.init_arg))")
        end
    end
    println(io,"")
end

function write_extract(arguments::Vector{f_args}, io = stdout)
    println(io,"")
    for arg in arguments
        if !isempty(arg.extract_arg)
            println(io,"\t$(arg.extract_arg)")
        end
    end
end

"""
Writes a function to the specified IO stream
"""
function write_funcs(funcs_val::Py, io = stdout)
    if contains(funcs_val.name,"_") #|| isopaque(funcs_val)
        @info "Skipping function $(funcs_val.name)"
        return nothing
    end

    name          = String(funcs_val.name)

    # extract input/output arguments and doc string of function
    input_vars, output_vars, doc_str = extract_input_output_function(petsc_dir, name)
    
    # process all function arguments
    arguments     = process_function_arguments(funcs_val.arguments, name, input_vars, output_vars)
    
    #@info arguments
    # print function

    # Extract input/output arguments of function
    #@info name
    
    #julia_fct_str = julia_function_header(funcs_val)
    #julia_fct_str = replace_function_string(julia_fct_str)
    julia_doc_fct_str, str_in, str_out, num_in, num_out = julia_function_doc_header(arguments, name)
    ccall_type, ccall_name, n_arg = julia_ccall_header(arguments)
    
    # write help 
    println(io,"\"\"\"");
    println(io,"\t$julia_doc_fct_str ");
    if !isnothing(doc_str)
        for c in doc_str
            println(io,"$c")
        end
    end
    println(io,"");
    println(io,"# External Links");
    println(io,"\$(_doc_external(\"$(titlecase(String(funcs_val.mansec)))/$(funcs_val.name)\"))");
    println(io,"\"\"\"");
    if num_in>0
        println(io,"function $(funcs_val.name)(petsclib::PetscLibType, $str_in) end");
    else
        println(io,"function $(funcs_val.name)(petsclib::PetscLibType) end");
    end
    println(io,"");
    # write function itself (we will generate various dispatches using @for_petsc)
    if num_in>0
        println(io,"@for_petsc function $(funcs_val.name)(petsclib::\$UnionPetscLib, $(replace_dispatch_types(str_in)) )");
    else
        println(io,"@for_petsc function $(funcs_val.name)(petsclib::\$UnionPetscLib)");
    end
    # Write initial arguments 
    write_initialize(arguments, io)
    
    println(io,"    @chk ccall(");
    println(io,"               (:$(name), \$petsc_library),");
    println(io,"               PetscErrorCode,");
    println(io,"               $(replace_dispatch_types(ccall_type)),");
    if n_arg>0
    println(io,"               $ccall_name,");
    end
    println(io,"              )");

    # convert output args into correct values (if needed)
    write_extract(arguments, io)
    if num_out>0
        println(io,"\n\treturn $str_out");  # return them
    else
        println(io,"\n\treturn nothing");
    end
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



function write_typedefs(typedef_val::Py, io = stdout)
    type = String(typedef_val.value)
    type = replace_types(type)

    println(io,"const $(typedef_val.name) = $type")
    return nothing
end


function write_keys_to_file(filename::String, start_dir::String, object::Py, fn::Function; exclude=String[])
    open(joinpath(start_dir, filename), "w") do file
        # Call the write_enum function and pass the file as the io argument
        for val in object.keys()
            if !any(exclude .== String(val))
                fn(object[val], file)
            end
        end
    end
end

function write_skeys_to_file(filename::String, start_dir::String, object::Py; exclude=String[])
    open(joinpath(start_dir, filename), "w") do io
        # Call the write_enum function and pass the file as the io argument
        println(io, "# not quite sure yet how to deal with this")
        for val in object.keys()
            if !any(exclude .== String(val))
                println(io, "$(String(object[val].name))=Ptr{Cchar}")
            end
        end
    end
end

function write_functions_from_classes_to_file(filename::String, start_dir::String, classes::Py, class_name::String; exclude=String[])
    open(joinpath(start_dir, filename), "w") do file
        # write definitions to file (if needed)
        write_undefined_typearguments_from_classes_to_file(classes[class_name].functions, file=file)

        # Call the write_enum function and pass the file as the io argument
        for f in classes[class_name].functions
            name = String(f)
            if !any(exclude .== name)
                @info name
                write_funcs(classes[class_name].functions[name], file)
                #write_funcs(classes[function_name].functions[String(f)])

            end
        end
    end
end

function write_functions_to_file(filename::String, start_dir::String, funcs::Py; exclude=String[])
    open(joinpath(start_dir, filename), "w") do file
        # Call the write_enum function and pass the file as the io argument
        write_undefined_typearguments_from_classes_to_file(funcs, file=file)

        for f in funcs
            name = String(f)
            if !any(exclude .== name)
                #@info name
                write_funcs(funcs[name], file)
                #write_funcs(classes[function_name].functions[String(f)])
            end
        end
    end
end

function write_structs_to_file(filename::String, start_dir::String, structs::Py; exclude=String[])
    open(joinpath(start_dir, filename), "w") do file
        # Call the write_enum function and pass the file as the io argument
        for struct_val in structs
            if !any(exclude .== String(struct_val))
                write_struct(structs[String(struct_val)], file)
            end
        end
    end
end

function write_typedefs_to_file(filename::String, start_dir::String, typedefs::Py; exclude=String[])
    open(joinpath(start_dir, filename), "w") do file
        # Call the write_enum function and pass the file as the io argument
        for typedef_val in typedefs
            if !any(exclude .== String(typedef_val))
                write_typedefs(typedefs[String(typedef_val)], file)
            end
        end
    end 
end


function read_docstring_function(func::Py, petsc_dir::String)
    name = String(func.name)
    mansec = String(func.mansec)
end

# helps to quickly generate custom julia structs for PETSc types
function generate_custom_julia_struct(type_name::String)
    # Convert type name to appropriate formats
    ctype_name = "C$(type_name)"
    petsc_type_name = "Petsc$(type_name)"
    abstract_type_name = "AbstractPetsc$(type_name)"
    
    str= """
# ----- Custom Julia struct for PETSc $(type_name) -----
const $(ctype_name) = Ptr{Cvoid}
abstract type $(abstract_type_name){T} end
mutable struct $(petsc_type_name){PetscLib} <: $(abstract_type_name){PetscLib}
    ptr::$(ctype_name)
    age::Int
    
    # Constructor from pointer and age
    $(petsc_type_name){PetscLib}(ptr::$(ctype_name), age::Int = 0) where {PetscLib} = new{PetscLib}(ptr, age)
    
    # Constructor for empty $(type_name) (null pointer)
    $(petsc_type_name){PetscLib}() where {PetscLib} = new{PetscLib}(Ptr{Cvoid}(C_NULL), 0)
end

# Convenience constructor from petsclib instance
$(petsc_type_name)(lib::PetscLib) where {PetscLib} = $(petsc_type_name){PetscLib}()
$(petsc_type_name)(ptr::$(ctype_name), lib::PetscLib, age::Int = 0) where {PetscLib} = $(petsc_type_name){PetscLib}(ptr, age)
Base.convert(::Type{$(ctype_name)}, v::$(abstract_type_name)) = v.ptr
Base.unsafe_convert(::Type{$(ctype_name)}, v::$(abstract_type_name)) = v.ptr

# Custom display for REPL
function Base.show(io::IO, v::$(abstract_type_name){PetscLib}) where {PetscLib}
    if v.ptr == C_NULL
        print(io, "PETSc $(type_name) (null pointer)")
        return
    else
        print(io, "PETSc $(type_name) ")
    end
    return nothing
end
# ------------------------------------------------------
"""
    
    println(str)
    return nothing
end


"""
Helper function to simply finding a certrain PETSc function
"""
function find_functions(classes::Py, funcs::Py, function_name::String)

    found = false
    for class in classes.keys()
        # retrieve all functions
        funcs_list = String[]
        for f in classes[class].functions
            name = String(f)
            push!(funcs_list, name)
        end
        if any(funcs_list .== function_name)
            println("Found $function_name in class $(String(class))")
            found = true
        end
    end
    if !found
        println("Did not find $function_name in any class")
    end

    if !found
        funcs_list = String[]
        for f in funcs
            name = String(f)
            push!(funcs_list, name)
        end
         if any(funcs_list .== function_name)
            println("Found $function_name in funcs")
            found = true
        end
    end


    return nothing
end

function move_prologue(prologue_file="prologue.jl")
    # move prologue file to autowrapped directory
    src = joinpath(start_dir, prologue_file)
    dest = joinpath(start_dir, "../src/autowrapped/", "petsc_library.jl")
    cp(src, dest; force=true)
    return nothing
end

function extract_function_args_from_function(func::Py, function_name::String)
    arguments     = process_function_arguments(funcs_val.arguments, function_name, String[], String[])
    fn_args = String[]
    for arg in arguments
        if arg.isfunction
            push!(fn_args, arg.typename)
        end
    end
    return fn_args
end

# Extracts all type argiments from a class
function extract_function_typeargs_from_class(functions::Py)
    fn_args = Vector{String}()
    for f in functions
        name = String(f)
        arguments     = process_function_arguments(functions[name].arguments, name, String[], String[])
        for arg in arguments
            typename = arg.typename
            typename = replace(typename, "Vector{"=>"","}"=>"")

            push!(fn_args, typename)
        end
    end
    fn_args = unique(fn_args)
    for entry in CustomTypes
        fn_args = remove_entry(fn_args, entry)
    end

    return fn_args
end

function write_undefined_typearguments_from_classes_to_file(functions; file=stdout)
    fn_args = extract_function_typeargs_from_class(functions)

    println(file,"# autodefined type arguments for class ------")
    # Call the write_enum function and pass the file as the io argument
    for typearg in fn_args
       var_sym = Symbol(typearg)
        if !isdefined(Main, var_sym)
            # variable doesn't exist yet, so define it
            if contains(typearg,"Fn")
                println(file,"mutable struct $typearg end\n")
            else
                println(file,"mutable struct _n_$typearg end")
                println(file,"const $typearg = Ptr{_n_$typearg}\n")
            end
        end
    end
    println(file,"# -------------------------------------------------------")


end





#exclude=["KSPConvergedReason","PetscMemType"]
#write_keys_to_file("../src/autowrapped/enums_wrappers.jl",  start_dir,  enums, write_enum, exclude=exclude)  # Write enums to file
#write_skeys_to_file("../src/autowrapped/senums_wrappers.jl",start_dir, senums)              # Write string enums to file
#
#exclude = ["LandauCtx"]
# Adapt as needed
#write_structs_to_file("../src/autowrapped/struct_wrappers.jl", start_dir, structs,exclude=exclude)          # Write all structs to file
#exclude=["PetscGeom","PetscInt32"]
#write_typedefs_to_file("../src/autowrapped/typedefs_wrappers.jl", start_dir, typedefs,exclude=exclude)      # Write all typedefs to file

move_prologue("prologue.jl")


# Need to call this before writing any function files
# - it will error, but at that stage all necessary variables will be loaded, so you can comment this
#include("../src/autowrapped/petsc_library.jl")


#=
# Write KSP functions to file (this should be expanded to all other classes)
exclude=["KSPLSQRGetNorms"]
write_functions_from_classes_to_file("KSP_wrappers.jl",start_dir, classes, "KSP", exclude=exclude)     

# write general functions to file
exclude=["PetscHTTPSRequest",
"LandauKokkosJacobian","LandauKokkosDestroyMatMaps","LandauKokkosStaticDataSet","LandauKokkosStaticDataClear",
         "PetscSSLInitializeContext","PetscSSLDestroyContext","PetscHTTPSConnect",
         "PetscDataTypeToHDF5DataType","PetscHDF5DataTypeToPetscDataType","PetscHDF5IntCast",
         "PetscPostIrecvInt","PetscPostIrecvScalar",
         "PetscDTAltVInteriorPattern","LandauKokkosCreateMatMaps",
         "PetscCIntCast","PetscIntCast", "PetscBLASIntCast","PetscCuBLASIntCast","PetscHipBLASIntCast","PetscMPIIntCast",
         "PetscFixFilename",
         "PetscDTAltVApply","PetscDTAltVWedge","PetscDTSimplexQuadrature","PetscDTJacobiNorm",  # these have some comments that mess up julia docs        
         ]

=#
#exclude=["LandauKokkosJacobian","LandauKokkosDestroyMatMaps","LandauKokkosStaticDataSet","LandauKokkosStaticDataClear"]
#write_functions_to_file("Sys_wrappers.jl", start_dir, funcs, exclude=exclude)

# Write Vec functions to file (this should be expanded to all other classes)
exclude=[""]
#write_functions_from_classes_to_file("Vec_wrappers.jl",start_dir, classes, "Vec", exclude=exclude)     
#write_functions_from_classes_to_file("Vecs_wrappers.jl",start_dir, classes, "Vecs", exclude=exclude)     
#write_functions_from_classes_to_file("PetscOptions_wrappers.jl",start_dir, classes, "PetscOptions", exclude=exclude)     
write_functions_from_classes_to_file("PetscObject_wrappers.jl",start_dir, classes, "PetscObject", exclude=exclude)     

#exclude=["MatSolves","MatCreateVecs","MatCreateVecsFFTW"]
#exclude=[""]
#write_functions_from_classes_to_file("Mat_wrappers.jl",start_dir, classes, "Mat", exclude=exclude)     



# Todo: 
# - determine function parameters that have functions and handle those automatically
# 