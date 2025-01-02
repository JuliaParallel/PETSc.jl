using Clang.Generators
using Clang.Generators.JLLEnvs
using CSTParser
using JuliaFormatter

## rewrite passes
struct Edit{T}
    loc::T
    text::String
end
mutable struct State
    offset::Int
    edits::Vector{Edit}
end

function pass(x, state, f = (x, state) -> nothing)
    res = f(x, state)
    if res && length(x) > 0
        for a in x
            pass(a, state, f)
        end
    else
        state.offset += x.fullspan
    end
    state
end

# insert `@chk` before each function with a `ccall` returning a checked type`
checked_types = ["PetscErrorCode"]
function add_chk(x, state)
    if x isa CSTParser.EXPR && x.head == :function
        _, def, body, _ = x
        @assert body isa CSTParser.EXPR && body.head == :block
        @assert length(body) == 1

        # Clang.jl-generated ccalls should be directly part of a function
        # definition
        call = body.args[1]
        if call isa CSTParser.EXPR &&
           call.head == :call &&
           call[1].val == "ccall"

            # get the ccall return type
            rv = call[5]

            if rv.val in checked_types
                push!(
                    state.edits,
                    Edit(
                        state.offset + x[1].span + 1 + x[2].span + 1 + 4,
                        "@chk ",
                    ),
                )
            end
        end
    end
    true
end

# make every function have $UnionPetscLib as the first argument
function add_UnionPetscLib(x, state)
    if x isa CSTParser.EXPR && x.head == :function
        _, def, body, _ = x
        @assert body isa CSTParser.EXPR && body.head == :block
        @assert length(body) == 1

        # Clang.jl-generated ccalls should be directly part of a function
        # definition
        call = body.args[1]
        if call isa CSTParser.EXPR &&
           call.head == :call &&
           call[1].val == "ccall"
            push!(
                state.edits,
                Edit(
                    state.offset + x[1].span + 1 + x[2][1].span + 1,
                    "::\$UnionPetscLib, ",
                ),
            )
        end
    end
    true
end

# make every function have $UnionPetscLib as the first argument
function add_for_petsc(x, state)
    if x isa CSTParser.EXPR && x.head == :function
        _, def, body, _ = x
        @assert body isa CSTParser.EXPR && body.head == :block
        @assert length(body) == 1

        # Clang.jl-generated ccalls should be directly part of a function
        # definition
        call = body.args[1]
        if call isa CSTParser.EXPR && (
            (call.head == :call && call[1].val == "ccall") ||
            (call.head == :macrocall && call[1].val == "@chk")
        )
            push!(state.edits, Edit(state.offset, "@for_petsc "))
        end
    end
    true
end

# interpolate the `val` in the functions with `@for_petsc`
function add_interpolations(x, state, vals)
    # skip structs
    if x isa CSTParser.EXPR && x.head == :struct
        return false
    end
    if x isa CSTParser.EXPR && x.head == :IDENTIFIER && x.val in vals
        push!(state.edits, Edit(state.offset, "\$"))
    end
    true
end

function add_struct_types(x, state, name, val)
    # skip structs
    if x isa CSTParser.EXPR && x.head == :struct
        return false
    end
    if x isa CSTParser.EXPR && x.head == :IDENTIFIER && x.val == val
        push!(state.edits, Edit(state.offset, "\$"))
    end
    true
end

# Make the `structs` have type values `vals`
function add_struct_type(x, state, structs, vals)
    # XXX: Maybe not the best way to do this...
    if x isa CSTParser.EXPR && x.head == :IDENTIFIER && x.val in structs
        push!(
            state.edits,
            Edit(
                state.offset + length(x.val),
                "{" * reduce((x, y) -> "" * x * "," * y, vals) * "}",
            ),
        )
    end
    true
end

# make petsc_library

function apply(text, edit::Edit{Int})
    string(text[1:(edit.loc)], edit.text, text[nextind(text, edit.loc):end])
end

function mangle_functions(output_file)
    let file = output_file
        text = read(file, String)

        for mod in (
            add_UnionPetscLib,
            add_chk,
            add_for_petsc,

            # Deal with structures that require type declarations:

            # Note some structs such as "_p_DMPlexPointQueue" needs some manual massaging 
            (x, s) -> add_struct_type(x, s, 
                    ("MatStencil","DMStagStencil","DMDALocalInfo","_p_DMPlexPointQueue"), 
                    ("PetscInt",)),
            (x, s) -> add_struct_type(x, s, 
                    ("DMDACoor2d","DMDACoor3d"), 
                    ("PetscScalar",)),
            (x, s) -> add_struct_type(x, s, 
                    ("PetscFVCellGeom",), 
                    ("PetscReal",)),
            (x, s) -> add_struct_type(x, s, 
                    ("PetscFVFaceGeom",), 
                    ("PetscReal","PetscScalar")),

            #(x, s) -> add_struct_type(x, s, 
            #        ("MatStencil","_n_PetscLayout","PetscViewerAndFormat","_PetscCDIntNd","_PetscCoarsenData","_PetscFormKey","_p_DMPlexPointQueue"), 
            #        ("PetscInt",)),
            #(x, s) -> add_struct_type(x, s, 
            #        ("petsc_mpiu_scalar_int",), 
            #        ("PetscScalar, PetscInt",)),
            #(x, s) -> add_struct_type(x, s, 
            #        ("petsc_mpiu_real_int",), 
            #        ("PetscReal, PetscInt",)),

            (x, s) -> add_interpolations(
                x,
                s,
                ("PetscReal", "PetscInt", "PetscComplex", "PetscScalar"),
            ),
        )
            state = State(0, Edit[])
            ast = CSTParser.parse(text, true)

            state.offset = 0
            pass(ast, state, mod)

            state.offset = 0
            sort!(
                state.edits,
                lt = (a, b) -> first(a.loc) < first(b.loc),
                rev = true,
            )

            for i in 1:length(state.edits)
                text = apply(text, state.edits[i])
            end
        end

        write(output_file, text)
    end
end

function final_mangling(output_file)
    # In a few cases, the procedure fails so deal with that manually  
    text = read(output_file, String)
    
    text_new = replace(text,"_p_DMPlexPointQueue{\$PetscInt}"=>"_p_DMPlexPointQueue")

    write(output_file, text_new)
end


function wrap(output_file, petsc_include_dir, mpi_include_dir)
    petsc_h = joinpath(petsc_include_dir, "petsc.h")
    @assert isfile(petsc_h)

    #petscis_h = joinpath(petsc_include_dir, "petscis.h")
    #@assert isfile(petscis_h)

    mpi_h = joinpath(mpi_include_dir, "mpi.h")
    @assert isfile(mpi_h)

    options = load_options(joinpath(@__DIR__, "generator.toml"))
    options["general"]["output_file_path"] = output_file

    args = get_default_args()

    push!(args, "-I$petsc_include_dir")
    push!(args, "-isystem$mpi_include_dir")

#    header_files = [petsc_h, petscis_h]
    header_files = [petsc_h]

    ctx = create_context(header_files, args, options)

    build!(ctx)
    @info "Finished building the context"
    mangle_functions(output_file)
    @info "mangled the functions"
    format(output_file)
    @info "formatted output file"

    # Deal with a few failed cases (so we don't have to do that by hand)
    final_mangling(output_file)
    @info "and performed final mangling"


    return output_file
end
