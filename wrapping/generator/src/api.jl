# Loading of the JSON API snapshot written by getapi_dump.py

struct Arg
    name::String
    typename::String       # raw C type name (before any Julia mapping)
    stars::Int
    array::Bool
    isconst::Bool
    optional::Bool
    isfunction::Bool
    stringlen::Bool
end

struct Fn
    name::String
    mansec::String
    class::String          # "" for a standalone function
    args::Vector{Arg}
end

struct API
    version::String
    layout::String
    functions::Dict{String,Fn}            # all functions by name (standalone + class methods)
    classes::Dict{String,Vector{String}}  # class name => sorted function names
    standalone::Vector{String}            # sorted names of functions not attached to a class
    enums::Dict{String,Vector{String}}    # enum name => values ("NAME" or "NAME = value")
    senums::Dict{String,Vector{String}}   # string-enum name => keys
    typedefs::Dict{String,String}         # typedef name => C type
    structs::Dict{String,Vector{String}}  # struct name => raw record strings
end

_str(x) = x === nothing ? "" : String(x)

function _arg(a)
    name = _str(a[:name])
    occursin(' ', strip(name)) && (name = String(last(split(name))))   # getAPI.py mis-parse such as "int cmp"
    Arg(name, _str(a[:typename]), Int(a[:stars]), Bool(a[:array]), Bool(a[:const]),
        Bool(a[:optional]), Bool(a[:isfunction]), Bool(a[:stringlen]))
end

_fn(f, class) = Fn(_str(f[:name]), _str(f[:mansec]), class, Arg[_arg(a) for a in f[:arguments]])

function load_api(path::AbstractString)
    js = JSON3.read(read(path, String))
    functions = Dict{String,Fn}()
    classes = Dict{String,Vector{String}}()
    for (cname, c) in pairs(js[:classes])
        names = String[]
        for (fname, f) in pairs(c[:functions])
            fn = _fn(f, String(cname))
            functions[fn.name] = fn
            push!(names, fn.name)
        end
        classes[String(cname)] = sort!(names)
    end
    standalone = String[]
    for (fname, f) in pairs(js[:funcs])
        fn = _fn(f, "")
        if !haskey(functions, fn.name)
            functions[fn.name] = fn
            push!(standalone, fn.name)
        end
    end
    sort!(standalone)
    enums = Dict(String(k) => String[_str(v) for v in e[:values]] for (k, e) in pairs(js[:enums]))
    senums = Dict(String(k) => sort!(String[String(kk) for kk in keys(e[:values])]) for (k, e) in pairs(js[:senums]))
    typedefs = Dict(String(k) => _str(t[:value]) for (k, t) in pairs(js[:typedefs]))
    structs = Dict(String(k) => String[_str(r[:type]) for r in s[:records]] for (k, s) in pairs(js[:structs]))
    API(_str(js[:petsc_version]), _str(js[:getapi_layout]), functions, classes, standalone,
        enums, senums, typedefs, structs)
end
