# Once the huge library is created, we would like to create  
# more julia-friendly functions. The routines here read the header file,
# along with the PETSc source code and automatically create julia functions (with docstrings) 
# 
# Note that this works for a range of  functions, but not all. 
# Importantly, the key structures such as Mat, of DMStag still need to be declared manually,
# as we add additional functionality in the julia version. 

using Glob

hname = "headers/petscdmstag.h"
fname = "DMStagGetGlobalSizes"
fname = "DMStagGetDOF"
#fname = "DMStagGetEntriesPerElement"
#fname = "DMStagCreateCompatibleDMStag"
#fname = "DMStagCreateISFromStencils"
#fname = "DMStagCreate1d"
fname = "DMStagCreateISFromStencils"

# Works, but requires special attention:
# "DMStagCreate1d"                  # requires to specify options=true and add_petsclib=true
# "DMStagCreate2d"
# "DMStagCreate3d"
# "DMStagCreateCompatibleDMStag"    # works now, but we need to specify options=true
# "DMStagGetOwnershipRanges"        # needs some internal modification; as indicated with TODO comments
# DMStagVecSetValuesStencil         # simply check the length of input vectors and create a multiple dispatch
# DMStagVecGetValuesStencil         # needs some mangling with lengh of oputput vector
# DMStagMatGetValuesStencil         
# DMStagMatSetValuesStencil
# DMStagVecSplitToDMDA              # needs manual work
# DMStagStencilToIndexLocal         # needs some manual work as a vector is returned
# DMStagVecGetArray                 # needs manual work, to be improved
# DMStagVecRestoreArray             # 
# DMStagVecGetArrayRead
# DMStagVecRestoreArrayRead
# DMStagGetProductCoordinateArrays  # needs manual work (easy)
# DMStagRestoreProductCoordinateArraysRead

# Problematic ones are:

# DMStagCreateISFromStencils            # IS is not yet inmplemented, so this one doesn't work yet



"""
    download_petsc(petsc_version=v"3.22.0")

Automatically download the correct PETSc version and out it in the current directory
"""
function download_petsc(petsc_version=v"3.22.0")
    url = "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-$(petsc_version).tar.gz"
    Base.download(url,"petsc.tar.gz")
    run(`tar -xvf petsc.tar.gz`)
    run(`mv petsc-$(petsc_version) petsc`)
end

# This finds the C-file that contains the function we are looking for
function find_c_file(path::String, fct_name::String)
    
    files = glob("*.c", path)

    line_string = ""

    file_str = ""
    comment_block=""
    # Scan all files in the directory; find the line
    for file in files
        isfile = open(file, "r") do io
            contains(read(io, String), "PetscErrorCode $(fct_name)(")
        end
        
        if isfile
            file_str = file;
            line_string = read_c_function_args(file, fct_name)
            comment_block = read_c_function_docs(file, fct_name)
        end
    end
    
    return line_string, comment_block, file_str
end

# Finds the function in the correct PETSc c-file
function read_c_function_args(file, fct_name)
    line_string = ""
    open(file, "r") do f
        # read till end of file
        while ! eof(f)  
            line = readline(f)
            if contains(line, "PetscErrorCode $fct_name(")
                line_string = line
            end
        end
    end
    return line_string
end

function add_backticks(line::AbstractString)
    line_split = split(line,"- ")
    var_name = line_split[2]
    n = length(var_name)
    
    line_split[2] = rpad("`"*strip(var_name)*"`",n+2)
    
    line_new = join(line_split.*"- ")[1:end-2]
    return  line_new
end

# Find the docs of the c-string
function read_c_function_docs(file, fct_name)

    comment_block =  [""]
    open(file, "r") do f
        
        # read till end of file
        while ! eof(f)  
            line = readline(f)
            if contains(strip(line),"/*@")# start of comment block
                read_block = false
                line = readline(f)
              
                if split(strip(line))[1] == fct_name
                    read_block = true
                end
              
                while ! contains(line, "@*/") & read_block
                    line = strip(line)  # remove white space at beginningf
                    if !isempty(line)
                        if startswith(line,".") || startswith(line,"+")
                            line = "-"*line[2:end]
                        end
                    end
                 
                    # add backticks to input/output parameters
                    line = strip(line)
                    if startswith(line,'-')
                        if length(findall('-',line))>1
                            #  there are 2 dashes in lists with parameters
                            line = add_backticks(line)
                        end
                    end

                    push!(comment_block, line )

                    line = readline(f)
                end
            end
        end
    end

    comment_block = comment_block[2:end]
    if length(comment_block)>0
        comment_block[1] = strip(comment_block[1])
        comment_block[1] = strip(split(comment_block[1],"-")[2])
    end
    comment_block = replace.(comment_block, "\$" => "")
    comment_block = replace.(comment_block, "[](ch_stag)," => "")
    comment_block = replace.(comment_block, "[](ch_dmbase)," => "")
    comment_block = replace.(comment_block, "-seealso:  " => "See also: \n=== \n")
    comment_block = replace.(comment_block, "seealso:  " => "See also: \n=== \n")
    

    return comment_block
end



flip_arg(args::SubString) = split(args)[2]*"::"*split(args)[1]

wait_for_key(prompt) = (print(stdout, prompt); read(stdin, 1); nothing)

# This struct contains info about the petsc function, its arguments 
# its julia implementation
mutable struct PetscFct 
    name::String
    args_petsc::String
    deprecated::Bool
    docstr_petsc::Vector{String} # petsc docstring (C)
    doc_str::Vector{String}      # julia docstring
    body::Vector{String}         # julia body of fct
end

struct Petsc_Body 
    name::String
    body::Vector{String}         # julia body of fct
end


mutable struct PetscStruct
    name::String
    docstr_petsc::Vector{String}           # petsc docstring (C)
    doc_str::Vector{String}                # julia docstring
    body::Vector{String}                   # julia body of fct
end

mutable struct PetscEnum
    name::String
    docstr_petsc::Vector{String}           # petsc docstring (C)
    doc_str::Vector{String}                # julia docstring
    body::Vector{String}                   # julia body of fct
end

mutable struct PetscComment
    name::String
    docstr_petsc::Vector{String}        # petsc docstring (C)
    doc_str::Vector{String}             # julia docstring
end

# helper functions
add_to_namedtuple(tuple_list, name, fct) =  merge(tuple_list, NamedTuple{(Symbol(name),)}((fct,)))
add_to_namedtuple(petsc_fct_list, petsc_fct) = add_to_namedtuple(petsc_fct_list, petsc_fct.name, petsc_fct)

# extracts the function name and arguments from a line in the petsc header file
function extract_petsc_header_fct(s::String)
    ss          = split(s) # split into words
    depr        = contains(s,"PETSC_DEPRECATED_FUNCTION")  # deprecated or not?
     
    # Line that contains the function name:
    ind         = findall(contains.(ss,"PetscErrorCode"))[1]+1; 
    fct_name    = split(ss[ind],"(")[1] # function name
    id          = findfirst(fct_name,s)[end]+1    # index of function name
    fct_args    = split(s[id:end],")")[1]*")"     # arguments    

    petsc_fct =  PetscFct(fct_name, fct_args, depr, [""], [""],[""])
    return petsc_fct
end

# reads info about a struct of enum
function extract_petsc_header_struct(f,s)
    while ! contains(s,"}") # read till end of enum
        s = readline(f)
    end
    # name is given @ the end
    id = findfirst("}", s)[end];
    struct_name = strip(replace(s[id+1:end],";"=>""))
    return PetscStruct(struct_name,[""],[""],[""])
end

# reads info about a struct of enum
function extract_petsc_header_enum(f,s)
    while ! contains(s,"}") # read till end of enum
        s = readline(f)
    end
    # name is given @ the end
    id = findfirst("}", s)[end];
    enum_name = strip(replace(s[id+1:end],";"=>""))
    return PetscEnum(enum_name,[""],[""],[""])
end

# This reads a comment from a petsc header file
function extract_petsc_comment(f)
    s = readline(f)
    comment = String[s]

    name = split(s)[1]
    while ! contains(s,"E*/") # read till end of comment
        s = readline(f)
        push!(comment, s)
    end

    return PetscComment(name, comment,[""])
end

# Extracts a block
function extract_petsc_block(f,s)
    comment = String[s]
    while ! contains(s,"end") # read till end of comment
        s = readline(f)
        push!(comment, s)
    end
    return comment
end


function extract_petsc_fct(f,s)
    body = extract_petsc_block(f,s)
    name = strip(split(split(s,"@for_petsc function")[2],"(")[1]);
    return Petsc_Body(name,body)
end

function extract_petsc_struct(f,s)
    body = extract_petsc_block(f,s)
    name = strip(split(split(s,"mutable struct")[2],"(")[1]);
    name = strip(replace(name, "end" => ""))
    return Petsc_Body(name,body)
end

function extract_petsc_enum(f,s)
    body = extract_petsc_block(f,s)
    name = strip(split(split(s,"@enum")[2],"(")[1]);
    name = split(replace(name, "::" => " "))[1]
    return Petsc_Body(name,body)
end

# This analyzes a PETSc header file and extracts:
# - functions
# - enums
# - comments
# - structs 

function read_petsc_header(hname, fname)
    @assert isfile(hname)

    petsc_fct_list = NamedTuple()
    #struct_list =  NamedTuple()
    #enum_list =  NamedTuple()
    #comment_list =  NamedTuple()
    function_names = String[]
    open(hname) do f

        # do stuff with the open file instance 'f'
        line = 0  

        # read till end of file
        while ! eof(f)  

            s = readline(f)     # read a new / next line for every iteration         

            # Retrieve function names and arguments
            if contains(s, "PetscErrorCode") # petsc function
                petsc_fct = extract_petsc_header_fct(s);
                push!(function_names, petsc_fct.name)
                if petsc_fct.name==fname
                    petsc_fct_list = add_to_namedtuple(petsc_fct_list, petsc_fct)
                    break
                end
                #
            end

            #=
            # Retrieve struct/enum names
            if contains(s, "typedef struct") 
                petsc_struct = extract_petsc_header_struct(f,s)
                # Add to list of structures
                struct_list = add_to_namedtuple(struct_list, petsc_struct)
            end

            # Retrieve struct/enum names
            if contains(s, "typedef enum")  
                petsc_enum = extract_petsc_header_enum(f,s)
                enum_list = add_to_namedtuple(enum_list, petsc_enum)
            end

            # Retrieve comments
            if contains(s, "/*E") # comment
                petsc_comment = extract_petsc_comment(f)
                
                # Add to list of structures
                comment_list = add_to_namedtuple(comment_list, petsc_comment)
            end
            =#

            line += 1
        end

    end
    return petsc_fct_list, function_names 
end


function split_input_output(C_fct::AbstractString)
    @info C_fct
    # Get function and arguments
    C_fct = split(C_fct,"PetscErrorCode ")[2]
    args_string = split(C_fct,"(")[2][1:end-1]

    # Function arguments
    args_C = strip.(split(args_string,","))
    args_C_orig = copy(args_C)

    # determine if we have vectors. They will be input to the function if both const and *name are present
    if any(contains.(args_C,"[]")) || any(contains.(args_C,"*") .&& contains.(args_C,"const"))
        id_vec = findall(contains.(args_C,"[]") .|| 
                         (contains.(args_C,"*") .&& contains.(args_C,"const")))
        for id in id_vec
            # For now, I assume that vectors will always start with const
            str = args_C[id]
            if contains(str,"[]")
                str = strip(replace(str,"const "=>""))
                str_split = split(rsplit(str, "[]"; limit=2)[1])
                str_split[1] =  "Vector{"*str_split[1]*"}"
                str = strip(join(str_split.*" "))
            end
            if  (contains(str,"*") .&& contains(str,"const"))
                str = strip(replace(str,"const "=>""))
                str_split = split(rsplit(str, "[]"; limit=2)[1])
                str_split[1] =  "Vector{"*str_split[1]*"}"
                str = strip(join(str_split.*" "))

                #str = replace(str,"const "=>"Vector{")
                #str = replace(str," "=>"} ")
            end
            args_C[id] = str     
            #str = replace(str,"const "=>"Vector{")
            #str = replace(str," "=>"} ")
            args_C[id] = str     
        end   
    end

    # flip arguments and make it julia-like
    args = flip_arg.(args_C)

    
    # whatever has * is an output; rest input:
    # NOTE: this is wrong...; if we have vectors as input, they can be declared as, for example, const PetscScalar *val 

    output = String[]
    input  = String[]
    for (i,arg) in enumerate(args)
        if contains(args_C_orig[i],"*") && !contains(args_C_orig[i],"const")
            push!(output,arg)
        else
            push!(input,arg)
        end
    end

    # Replace *
    input  = replace.(input,"*"=>"")
    output = replace.(output,"*"=>"")

    # Replace input variables with julia types
    input_split = split.(input,"::")
    for id in eachindex(input_split)
        inp = input_split[id]
        if inp[2] == "DM"
            if dmtype=="DMStag"
                inp[2] = "AbstractDMStag{PetscLib}"
            elseif dmtype=="DM"
                inp[2] = "AbstractDM{PetscLib}"
            else
                println("inp[2]=$(inp[2]), dmtype=$dmtype")
                error("correct code")
            end
        elseif inp[2] == "Vec"
            inp[2] = "AbstractVector"
        elseif inp[2] == "Mat"
            inp[2] = "AbstractMatrix"
        elseif inp[2] == "DMBoundaryType"
            inp[2] = "DMBoundaryType"
        elseif inp[2] == "KSP"
            inp[2] = "AbstractKSP{PetscLib}"
        elseif inp[2] == "Vector{char}"
            inp[2] = "Vector{Char}"
        elseif inp[2] == "int"
            inp[2] = "Int"
        end
        input[id] = inp[1]*"::"*inp[2]

        # Replace PetscInt with Integer
        input[id] = replace(input[id],"PetscInt"=>"Int")
        input[id] = replace(input[id],"::PetscScalar"=>"::AbstractFloat")
        input[id] = replace(input[id],"::PetscReal"=>"::AbstractFloat")
        
        input[id] = replace(input[id],"Vector{DMStagStencil}"=>"Vector{DMStagStencil{Int}}")
        
    end

    return input, output
end


function print_petsc_function(input, output, fname, comment_block, type="DMStag"; options=false, add_petsclib=false, io=nothing, output_file="wrapped_functions.jl", have_tests=[""])
    # header of julia function
    if add_petsclib
        pushfirst!(input,"petsclib::PetscLib")
    end
    @show input

    function_header_main = "$(fname)("*join(input,",")*")"                       
    function_header = "function "*function_header_main                     
    function_header *= " where {PetscLib}"                       

    # initialize output variables
    input_split  = split.(input,"::")
    output_split = split.(output,"::")

    # Return string in doc_string
    return_doc = ""
    for outp in output_split
        return_doc *= "$(outp[1])"*","
    end
    return_doc = return_doc[1:end-1]
    if isempty(output)
        return_doc = "nothing"
    end
    
    if options
        # add options
        function_options =",dmsetfromoptions::Bool=true" 
        function_options *=",dmsetup::Bool=true" 
        function_options *=",options...)" 
        function_header = replace(function_header,")"=>function_options)
        function_header_main = replace(function_header_main,")"=>function_options)
    end

    # Print 
    if isnothing(io)
        io = open(output_file, "w")
    end
    println(io, "\"\"\"");
    if !any(contains.(have_tests, fname))
        println(io, "\t UNTESTED !!!"); # we don't have a test yet
    end

    if !isempty(output)
        println(io, "\t$return_doc = $function_header_main");
    else
        println(io, "\t $function_header_main");
    end
    println(io, "");
    
    # add lines of comment block
    for line in comment_block
        if !isempty(line)
            if split(line)[1] == "Output" 
                if options
                    println(io, "- dmsetfromoptions - call set from options")
                    println(io, "- dmsetup          - call setup")
                    println(io, "- options          - additional options")
                    println(io, "")
                end
            end
        end


        println(io, "$(line)")
        if !isempty(line)
            if split(line)[1] == "Input" 
                println(io, "===")
                if add_petsclib
                    println(io, "- petsclib     - the PETSc library")
                end
                
            elseif split(line)[1] == "Output" || split(line)[1] == "Options"
                println(io, "===")
            end
        end

    end # end of writing comment block
    println(io, "");
    println(io, "# External Links");
    println(io, "\$(_doc_external(\"$(dmtype)/$(fname)\"))");
    println(io, "\"\"\"");
    println(io, "$function_header");

    # If we have vectors as input, we likely also have the size of the vectors as input
    # It would make sense to add an assert statement here
    if any(contains.(input,"Vector"))
        println(io, "\t# TODO: you have vectors as input; make sure to test that the size of the vectors fits with something like: ");
        println(io, "\t# @assert length() == n ");
        println(io, "\t# You can likely also write a multiple dispatch version of this function where vector length is determined automatically ");
        @warn "test the size of input vectors to the generated function"
    end

    # Declare types
    if any(contains.(output,"PetscInt"))
        println(io, "\tPetscInt = PetscLib.PetscInt")
    elseif any(contains.(output,"PetscScalar")) 
        println(io, "\tPetscScalar = PetscLib.PetscScalar")
    end

    @show output_split
    for (id,outp) in enumerate(output_split)
        if outp[2]=="PetscInt" || outp[2]=="PetscScalar"
            println(io, "\t$(outp[1]) = [$(outp[2])(1)]")
        elseif outp[2]=="PetscReal"
                println(io, "\tPetscReal = PetscLib.PetscReal")
                println(io, "\t$(outp[1]) = [$(outp[2])(1)]")
        elseif outp[2]=="PetscBool"
            println(io, "\t$(outp[1]) = Ref{PetscBool}()")
        elseif outp[2]=="DMBoundaryType"
            println(io, "\t$(outp[1]) = Ref{DMBoundaryType}(DM_BOUNDARY_NONE)")
        elseif outp[2]=="DMStagStencilType"
            println(io, "\t$(outp[1]) = Ref{DMStagStencilType}()")
        elseif outp[2]=="ISColoringType"
            println(io, "\t$(outp[1]) = Ref{ISColoringType}()")
        elseif outp[2]=="DM"    
            # NOTE: I wrote this for DMStag; DMPlex and DMDA will probably use DM as well

            # array
            if !add_petsclib
                println(io,"\tpetsclib = getlib(PetscLib)")
            end
            
            # change this when options are an input
            if options
                println(io,"\topts = Options(petsclib; options...)")
            else
                println(io,"\topts = Options(petsclib)")
            end
            if dmtype=="DMStag"
                println(io, "\t$(outp[1]) = DMStag{PetscLib}(C_NULL, opts, petsclib.age)")
            elseif dmtype=="DM"
                println(io, "\t$(outp[1]) = DM{PetscLib}(C_NULL, opts, petsclib.age)")    
            else
                error("correct code")
            end

            
        elseif outp[2]=="void"
            # array
            @warn "outputting a vector - you likely have to manually change the routine"
          
            println(io, "\tPetscScalar = PetscLib.PetscScalar")
            println(io, "\t#TODO: your output is a vector; ensure that the size is correct!")
            println(io, "\t#It may involve: dims = DMStagGetGhostCorners(dm)[4:6]")
            println(io, "\t#dims = DMStagGetGhostCorners(dm)[4:6]    # dimensions including ghost values; set to 0 if not 2D/3D")
            println(io, "\t#dims = getindex(dims,findall(dims.>0))   # retrieve non-zero values")
            println(io, "\t#dmE  = DMStagGetEntriesPerElement(dm)    # dof per element" )
            
            println(io, "\tdims = (X,)")
            println(io, "\tr_$(outp[1]) = PETSc_RefPtr(dims, PetscScalar)")

        elseif outp[2]=="Vector{PetscInt}"
            println(io, "\tr_$(outp[1]) = Ref{Ptr{PetscInt}}(C_NULL)")
        
        elseif outp[2]=="DMType"
            println(io, "\tr_$(outp[1]) = Ref{PETSc.DMType}()")

        elseif outp[2]=="KSPType"
            println(io, "\tr_$(outp[1]) = Ref{PETSc.CKSPType}()")

        elseif outp[2]=="Vec"
            println(io, "\tr_$(outp[1]) = Ref{CVec}()")

        elseif outp[2]=="IS"
            # to be checked
            println(io, "\t$(outp[1]) = LibPETSc.IS()")

        elseif outp[2]=="PetscDS"
            # to be checked
            println(io, "\t$(outp[1]) = LibPETSc.PetscDS()")

        elseif outp[2]=="DMLabel"
            # to be checked
            println(io, "\t$(outp[1]) = LibPETSc.DMLabel()")

        elseif outp[2]=="PetscObject"
            # to be checked
            println(io, "\t$(outp[1]) = LibPETSc.PetscObject()")

        elseif outp[2]=="PetscSection"
            # to be checked
            println(io, "\t$(outp[1]) = LibPETSc.PetscSection()")

        elseif outp[2]=="MatOrderingType"
            # to be checked
            println(io, "\t$(outp[1]) = LibPETSc.MatOrderingType()")

        elseif outp[2]=="PetscBT"
            # to be checked
            println(io, "\t$(outp[1]) = LibPETSc.PetscBT()")

        elseif outp[2]=="ISLocalToGlobalMapping"
            # to be checked
            println(io, "\t$(outp[1]) = LibPETSc.ISLocalToGlobalMapping()")

        elseif outp[2]=="ISColoring"
            # to be checked
            println(io, "\t$(outp[1]) = LibPETSc.ISColoring()")
        
        elseif outp[2]=="DMBlockingType"
            # to be checked
            println(io, "\t$(outp[1]) = LibPETSc.DMBlockingType()")

        elseif outp[2]=="DMField"
            # to be checked
            println(io, "\t$(outp[1]) = LibPETSc.DMField()")
        
        elseif outp[2]=="PetscSF"
            # to be checked
            println(io, "\t$(outp[1]) = LibPETSc.PetscSF()")

        elseif outp[2]=="PetscSF"
            # to be checked
            println(io, "\t$(outp[1]) = LibPETSc.PetscSF()")

        elseif outp[2]=="VecType"
            println(io, "\t$(outp[1]) = Ref{VecType}()")

        elseif outp[2]=="MatType"
            println(io, "\t$(outp[1]) = Ref{MatType}()")

        else
            @warn "Don't know how to declare the type of: $(outp)"
        end
    end

    println(io, "");

    # Print call to lib
    space = "\t"
    if options
        println(io, "\twith(dm.opts) do")
        space *= "\t"
    end
            
    println(io, "$(space)LibPETSc.$fname(")
    if !add_petsclib
        println(io, "$(space)\tPetscLib,")
    end

    # input arguments
    for inp in input_split
        println(io, "$(space)\t$(inp[1]),")
    end

    # output (should have Ref if they are scalars)
    require_ptr_to_vec = false
    for outp in output_split
        if outp[2]=="PetscInt" || outp[2]=="PetscScalar"
            str = "Ref($(outp[1]),1)"
        elseif outp[2]=="PetscBool" || outp[2]=="DMBoundaryType"
            str = "$(outp[1])"
        elseif outp[2]=="DMStagStencilType" ||  outp[2]=="DM" || outp[2]=="Mat"
            str = "$(outp[1])"
        elseif outp[2]=="ISColoringType" 
            str = "$(outp[1])" 
        elseif outp[2]=="Vec" 
            str = "r_$(outp[1])"
            require_ptr_to_vec = true
        elseif outp[2]=="Vector{PetscInt}" || outp[2]=="void" || outp[2]=="DMType" || outp[2]=="VecType" || outp[2]=="MatType" || outp[2]=="KSPType"
            str = "r_$(outp[1])"
            require_ptr_to_vec = true
        elseif outp[2]=="PetscDS"  || outp[2]=="DMLabel" || outp[2]=="IS" || outp[2]=="ISLocalToGlobalMapping" 
            # to be checked
            str = "$(outp[1])"

        elseif outp[2]=="ISColoring"  || outp[2]=="DMBlockingType" || outp[2]=="PetscSection" || outp[2]=="PetscObject" 
              # to be checked
            str = "$(outp[1])"

        elseif outp[2]=="DMField"   || outp[2]=="PetscBT" || outp[2]=="PetscSF" || outp[2]=="PetscReal"
            # to be checked
          str = "$(outp[1])"

        else
            str = "$(outp[1])"
            @warn "check how to write this to the output of docstring : $(outp)"

#            error("stop here")
        end
        println(io, "$(space)\t$str,")
    end
    println(io, "$(space))")
    if options
        println(io, "\tend")
        
        println(io, "\tdmsetfromoptions && setfromoptions!($(output_split[end][1]))")
        println(io, "\tdmsetup && setup!($(output_split[end][1]))")
    end
    
    if require_ptr_to_vec
        println(io, "")
        for outp in output_split
            if outp[2]=="void"
                  str = "$(outp[1]) = PETSc_unsafe_wrap(r_$(outp[1]), dims; own=false)"
                  println(io, "\t$str")
            end

        end
    end
    
    println(io, "");

    # 
    return_string = "\treturn "
    for outp in output_split
        if outp[2]=="PetscInt" || outp[2]=="PetscScalar" || outp[2]=="PetscReal"
            return_string *= "$(outp[1])[1]"
        elseif outp[2]=="PetscBool"
            return_string *= "$(outp[1])[] == PETSC_TRUE"
        elseif outp[2]=="DM" 
            return_string *= "$(outp[1])"

        elseif outp[2]=="Vec"
            println(io, "\t$(outp[1]) = VecPtr(PetscLib, r_$(outp[1])[], false)");
            return_string *= "$(outp[1])"

        elseif outp[2]=="DMBoundaryType"
            return_string *= "$(outp[1])[]"
        elseif outp[2]=="DMStagStencilType"
            return_string *= "string($(outp[1])[])"
        elseif outp[2]=="Vector{PetscInt}"
            @warn "you likely need to set the size of the array >1"
            println(io, "\tn_$(outp[1]) = 1; # TODO: modify this piece of code" );
            println(io, "\t$(outp[1]) = unsafe_wrap(Array, r_$(outp[1])[], n_$(outp[1]); own = false)" );

            return_string *= "$(outp[1])"
        elseif outp[2]=="void"
            return_string *= "$(outp[1])"

        elseif outp[2]=="DMType" ||  outp[2]=="VecType"  || outp[2]=="KSPType"
            #   return unsafe_string(t_r[])
            println(io, "\t$(outp[1]) = unsafe_string(r_$(outp[1])[])" );
            return_string *= "$(outp[1])"

        elseif outp[2]=="IS" || outp[2]=="PetscDS" || outp[2]=="DMLabel" || outp[2]=="PetscObject"
            # to be checked
            return_string *= "$(outp[1])"

        elseif outp[2]=="PetscSection"  || outp[2]=="MatOrderingType"
            # to be checked
            return_string *= "$(outp[1])"

        else
            return_string *= "$(outp[1])"
            @warn "check how to write this to the return statement: $(outp); not writing it"
        end
        return_string *= ","
    end
    if isempty(output)
        return_string *= "nothing,"
    end

    return_string = return_string[1:end-1]

    println(io, "$return_string");
    println(io, "end")
    println(io, " ")
    println(io, " ")
    

   # close(io)

    function_header_str = "$fname("
    for inp in input_split
        function_header_str *= inp[1]*","
    end
    function_header_str = function_header_str[1:end-1]*")"

    
    # Write how the test could look like 
    println("\t\t out = PETSc.$function_header_str")   
    println("\t\t @test out == ")


    return io
end



"""
    wrap_petsc_function(headername, function_names; options_functionnames=[""], excluded=[""], have_tests=[""], output="wrapped_functions.jl", dmtype="DMStag")

Routine to create julia-like wrappers for the PETSc header file `headername`.
We automatically retrieve 

Input
===

- `headername`  - full name of the PETSc header file including directory that you want to wrap (ends with `*.h`)
- `function_names` - name of the function(s) you want to wrap. If you only want to wrap a single function, provide a string; for more, a vector of strings. If you want to process all, say `function_names=:all`.
- `path_within_petsc` - directory within the petsc repository that contains the relevant source files (e.g. "petsc/src/dm/impls/stag/" for DMSTAG)
- `excluded` - list with excliuded functions
- `options_functionnames` - list of routines where we add options
- `addpetsclib_functionnames` - list of routines that don't have `PETScLib` as type signature
- `dmtype` type of DM structure 

"""
function wrap_petsc_function(headername, function_names, path_within_petsc; options_functionnames=[""], addpetsclib_functionnames=[""], excluded=[""], have_tests=[""], output_file="wrapped_functions.jl", dmtype="DMStag")
                    
    if !isdir("petsc/")
        error("I don't find petsc in the current directory. please download the correct version with `download_petsc()`")
    end

    if function_names == :all
        _, function_names = read_petsc_header(headername, "bono")
        @info "found $(length(function_names)) PETSc functions to be wrapped"
    end
    if isa(function_names, String)
        function_names = [function_names]
    end
  
    io = nothing
    for (i,fname) in enumerate(function_names)
        @info i, fname 

        process_fct = true
        if any(contains.(excluded, fname))
            process_fct = false
        end
        if process_fct
            C_fct, comment_block, file_str = find_c_file(path_within_petsc, fname)
            if !isempty(C_fct)
                input, output = split_input_output(C_fct)       
                
                options = false;
                if any(contains.(options_functionnames, fname))
                    options = true
                end
                
                add_petsclib = false
                if any(contains.(addpetsclib_functionnames, fname))
                    add_petsclib = true
                end

            
                @show options add_petsclib
                io = print_petsc_function(input, output, fname, comment_block, options=options, add_petsclib=add_petsclib, io=io, output_file=output_file, have_tests=have_tests)
            
            else
                @info "deprecated function or cannot find it"
                
            end
        end

    end
    if !isnothing(io)
        close(io)
        @info "Wrote all functions to: $output_file"
    end
    
end



#=
# Wrap DMSTAG
headername                  =   "headers/petscdmstag.h"
#function_names              =   :all #["DMStagCreateISFromStencils"]

#function_names              =   ["DMStagSetUniformCoordinatesProduct"]
function_names              =   ["DMStagCreateISFromStencils"]
path_within_petsc           =   "petsc/src/dm/impls/stag/"
#output_file                 =   "../src/dmstag_wrapped1.jl"

output_file                 =   "wrapped_functions.jl"

options_functionnames       =   ["DMStagCreate1d","DMStagCreate2d","DMStagCreate3d","DMStagCreateCompatibleDMStag"]
addpetsclib_functionnames   =   ["DMStagCreate1d","DMStagCreate2d","DMStagCreate3d"]
#excluded                    =   ["DMStagCreateISFromStencils"]
excluded                    =   []
=#

#=
# Wrap DM
headername                  =   "headers/petscdm.h"
function_names              =   ["DMGetType"]
function_names = :all

path_within_petsc           =   "petsc/src/dm/interface/"
output_file                 =   "wrapped_functions.jl"

options_functionnames       =   [""]
addpetsclib_functionnames   =   ["DMCreate"]
dmtype                      =   "DM"

# lots of excluded files; most because they call another function
excluded                    =   ["DMInitializePackage","DMRegister","DMCoarsenHookAdd","DMCoarsenHookRemove","DMRefineHookAdd","DMRefineHookRemove","DMGenerateRegister","DMGenerateRegisterAll","DMGenerateRegisterDestroy","DMGlobalToLocalHookAdd","DMLocalToGlobalHookAdd",
"DMSubDomainHookAdd","DMSubDomainHookRemove","DMSetApplicationContextDestroy","DMSetVariableBounds","DMFinalizePackage","DMSetNullSpaceConstructor","DMGetNullSpaceConstructor","DMSetNearNullSpaceConstructor","DMGetNearNullSpaceConstructor","DMAddBoundary","DMProjectFunction",
"DMProjectFunctionLocal","DMProjectFunctionLabel","DMProjectFunctionLabelLocal",
"DMProjectFieldLocal","DMProjectFieldLabel","DMProjectFieldLabelLocal",
"DMProjectBdFieldLabelLocal","DMComputeL2Diff","DMComputeL2GradientDiff",
"DMComputeL2FieldDiff","DMMonitorSet","DMMonitorSetFromOptions","DMCreateFieldIS",
"DMCreateMatrix","DMCreateSuperDM","DMCreateSectionSuperDM","DMCreateFieldDecomposition","DMCreateDomainDecomposition","DMCreateDomainDecompositionScatters",
"DMReorderSectionGetDefault","DMReorderSectionGetType","DMGetOutputSequenceNumber","DMOutputSequenceLoad",
"DMCreateFEDefault","DMCompareLabels"
]
=#

# Wrap KSP routines
headername                  =   "headers/petscksp.h"
function_names              =   ["KSPGetDM"]
function_names = :all

path_within_petsc           =   "petsc/src/ksp/ksp/interface/"
output_file                 =   "wrapped_functions.jl"

options_functionnames       =   [""]
addpetsclib_functionnames   =   ["DMCreate"]
dmtype                      =   "DM"

# lots of excluded files; most because they call another function (which we fix manually), or because our wrapper cannot deal with it yet
excluded                    =   ["KSPInitializePackage","KSPFinalizePackage","KSPSetMatSolveBatchSize","KSPGetMatSolveBatchSize","KSPRegister",
                                "KSPMonitorRegister","KSPSetPreSolve","KSPSetPostSolve",
                                "KSPMonitorSet","KSPSetOptionsPrefix","KSPConvergedReasonViewSet",
                                "KSPSetConvergenceTest","KSPGetConvergenceTest","KSPGetAndClearConvergenceTest",
                                "KSPGuessRegister",
                                "KSPCheckSolve","KSPSetPCSide","KSPSetPC","KSPSetSupportedNorm","KSPComputeRitz",
                                "KSPMonitorLGCreate","KSPMonitorResidualDrawLGCreate","KSPMonitorTrueResidualDrawLGCreate",
                                "KSPMonitorErrorDrawLGCreate","KSPMonitorSolutionDrawLGCreate","KSPMonitorSingularValueCreate",
                                "KSPMonitorDynamicToleranceDestroy","KSPMonitorDynamicToleranceCreate",
                                "KSPMonitorDynamicToleranceSetCoefficient","KSPConvergedDefaultDestroy","KSPConvergedDefaultCreate",
                                "KSPGuessView","KSPGuessDestroy","KSPGuessCreate","KSPGuessSetType","KSPGuessGetType",
                                "KSPGuessSetTolerance","KSPGuessSetUp","KSPGuessUpdate","KSPGuessFormGuess"]

# functions to be checked check separately:
#  KSPComputeOperatorsFn, 

# functions that have tests already
have_tests = ["KSPGetIterationNumber","KSPGetResidualNorm","KSPGetType","KSPSetTolerances","KSPGetTolerances","KSPGetTotalIterations","KSPSetOperators","KSPGetSolution","KSPSetDM"]

wrap_petsc_function(headername, function_names, path_within_petsc; options_functionnames=options_functionnames, addpetsclib_functionnames=addpetsclib_functionnames, excluded=excluded, output_file=output_file, have_tests=have_tests)


