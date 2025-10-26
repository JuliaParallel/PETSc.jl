# this contains functions that search the PETSc source directory for function definitions
# and return the corresponding docstrings.
#
# We use that in the julia wrapper interface to:
#   1) determine which are input and which are output parameters of a function
#   2) extract the docstring for each function parameter, to generate a nicely formatted julia one

"""
    get_all_c_files(path::String)

Returns a vector of all *.c files found recursively in the given path.
"""
function get_all_c_files(petsc_dir::String)
    path = petsc_dir*"/src"
    c_files = String[]
    for (root, dirs, files) in walkdir(path)
        for file in files
            if endswith(file, ".c")
                push!(c_files, joinpath(root, file))
            end
        end
    end
    
    return c_files
end


function filter_files(files::Vector{String})
    ind = contains.(files,"benchmarks") .== false
    files = files[ind]
    ind = contains.(files,"tutorials") .== false
    files = files[ind]
    ind = contains.(files,"tests") .== false
    files = files[ind]
    ind = contains.(files,"petsc4py") .== false
    files = files[ind]
    return files
end

"""
    find_c_file(petsc_dir::String, fct_name::String)

Find the files in which the C function `fct_name` may be defined within the PETSc source directory `petsc_dir`. 
"""
function find_c_file(petsc_dir::String, fct_name::String)
    
    files = get_all_c_files(petsc_dir)
    files = filter_files(files)

    #line_string = ""

    file_str = String[]
    #comment_block=""
    # Scan all files in the directory; find the line
    for file in files
        isfile = open(file, "r") do io
            contains(read(io, String), "PetscErrorCode $(fct_name)(")
        end
        
        if isfile
            push!(file_str, file)
            #line_string = read_c_function_args(file, fct_name)
            #comment_block = read_c_function_docs(file, fct_name)
        end
    end
    
    return file_str #line_string, comment_block, 
end

# Finds the function in the correct PETSc c-file
function read_c_function_args(file::String, fct_name)
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


function  read_c_function_docs(files::Vector{String}, fct_name)
    for f in files
        comment_block = read_c_function_docs(f, fct_name)
        if length(comment_block)>0
            return comment_block
        end
    end
end

function read_c_function_docs(file::String, fct_name)

    comment_block =  [""]
    open(file, "r") do f
        
        # read till end of file
        while ! eof(f)  
            line = readline(f)
            if contains(strip(line),"/*@")# start of comment block
                read_block = false
                line = readline(f)
                
                if !isempty(line)
                    if split(strip(line))[1] == fct_name
                        read_block = true
                    end
                end
              
                while ! contains(line, "@*/") & read_block
                    line = strip(line)  # remove white space at beginning
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

function add_backticks(line::AbstractString)
    line_split = split(strip(line), "- ")
    if length(line_split)>1
        var_name = line_split[2]
        
        n = length(var_name)
        
        line_split[2] = rpad("`"*strip(var_name)*"`",n+2)
        
        line_new = join(line_split.*"- ")[1:end-2]
    else
        line_new = line
    end
    return  line_new
end


function get_docs_from_function(petscdir, fct_name)
    files = find_c_file(petscdir, fct_name)
    comment_block = read_c_function_docs(files, fct_name)
    return comment_block
end


function get_last_line(input_line_start, comment; last_str="")
    for l in input_line_start:length(comment)
        if isempty(strip(comment[l]))
            input_line_end = l-1
            return input_line_end
        end
    end
end

function get_last_line_keyword(comment, keyword="")
    l = findfirst(contains.(comment, keyword))
    if !isnothing(l)
        return l-1
    end
end

function extract_variable_from_string(str::String)
    if  contains(str,'`')
        var_name = split(str, '`')[2]
        return var_name
    else
        warning("problem with line $str; no variable found")
    end
end

function extract_variable_from_string(str::String)
    if  contains(str,'-')
        parts = split(str, '-')
        if length(parts) < 2
            return ""
        end
        var_part = strip(parts[2])
        var_part = replace(var_part, "`" => "")
        var_name = split(var_part, ' ')[1]
        return var_name
    else
        return ""
    end
end


function extract_input_output_vars(comment::Vector{String})
    hasinput  = any(startswith.(comment,"Input Parameter"))
    hasoutput = any(startswith.(comment,"Output Parameter"))
    
    input_vars  = String[]
    output_vars = String[]

    if hasinput
        input_line_start  = first(findall(startswith.(strip.(comment),"Input Parameter")))+1
        input_line_end    = get_last_line(input_line_start, comment)
        
        for l = input_line_start:input_line_end
            var_name = extract_variable_from_string(comment[l])
            if !isempty(var_name)
                push!(input_vars, var_name)
            end
        end
    end

    if hasoutput
        output_line_start = first(findall(startswith.(strip.(comment),"Output Parameter")))+1
        output_line_end   = get_last_line(output_line_start, comment)

        for l = output_line_start:output_line_end
            var_name = extract_variable_from_string(comment[l])
            if !isempty(var_name)
                push!(output_vars, var_name)
            end
        end
    end

    return input_vars, output_vars
end

"""
     input_vars, output_vars = extract_input_output_function(petsc_dir, fct_name)
This retrieves the input and output variables for a given function `fct_name` in the PETSc source directory `petsc_dir`.
"""
function extract_input_output_function(petsc_dir, fct_name)
    comment_block = get_docs_from_function(petsc_dir, fct_name)
    if isnothing(comment_block)
        input_vars, output_vars = String[], String[]
    else
        input_vars, output_vars = extract_input_output_vars(comment_block)
    end
    comment_block = remove_notes_from_comment(comment_block)
    return input_vars, output_vars, comment_block
end

# sometimes, PETSc docstrings contains math or other stuff in notes.
# this messes up the julia docstrings, so we remove them here
function remove_notes_from_comment(comment::Vector{String})
    comment = remove_notes_from_comment(comment, "Note:")
    comment = remove_notes_from_comment(comment, "Notes:")
    comment = remove_notes_from_comment(comment, "Developer Note:")
    comment = remove_notes_from_comment(comment, "Fortran Notes:")
    comment = remove_notes_from_comment(comment, "Example Usage:")
    comment = remove_notes_from_comment(comment, "-vb")
    comment = remove_notes_from_comment(comment, "Example Usage\\:")
    comment = remove_weird_signs_from_comment(comment)  # remove some weirdities in the PETSc docstrings that julia dislikes
    return comment
end

function remove_weird_signs_from_comment(comment::Vector{String})
    for i in eachindex(comment)
        comment[i] = replace(comment[i], "\\:" => ":","``"=>"`") 
    end
    return comment
end
remove_weird_signs_from_comment(comment::Nothing) = comment

function remove_notes_from_comment(comment::Vector{String}, keyword::String)
    note_start = findfirst(contains.(comment, keyword))
    if !isnothing(note_start)
        #note_end   = get_last_line(note_start, comment)
        note_end   = get_last_line_keyword(comment,"-seealso:")
        if !isnothing(note_end)
            deleteat!(comment, note_start:note_end)
        end
    end
    return comment
end
remove_notes_from_comment(comment::Nothing) = comment

