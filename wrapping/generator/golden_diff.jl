# Function-level comparison of two autowrapped directories.
#
#   julia wrapping/golden_diff.jl GOLDEN_DIR NEW_DIR [--show NAME] [--list] [--normalize]
#
# Files are split into blocks keyed by the wrapper/enum/typedef/struct/type they define, so
# the comparison is independent of the order of functions within a file and of which file a
# type declaration lives in. With --normalize, trailing whitespace and blank lines are ignored.

include(joinpath(@__DIR__, "src", "blocks.jl"))
load_dir(dir) = load_blocks(dir)

function main(args)
    length(args) >= 2 || error("usage: golden_diff.jl GOLDEN NEW [--show NAME] [--list] [--normalize]")
    g = load_dir(args[1]); n = load_dir(args[2])
    show = nothing; list = false; norm = false; cat = false
    i = 3
    while i <= length(args)
        if args[i] == "--show"; show = args[i+1]; i += 2
        elseif args[i] == "--list"; list = true; i += 1
        elseif args[i] == "--normalize"; norm = true; i += 1
        elseif args[i] == "--categorize"; cat = true; norm = true; i += 1
        else error("unknown $(args[i])") end
    end
    if show !== nothing
        key = haskey(g, show) ? show : "fn:$show"
        a = haskey(g, key) ? g[key].text : "<missing in golden>"
        b = haskey(n, key) ? n[key].text : "<missing in new>"
        ta = tempname(); tb = tempname()
        write(ta, a); write(tb, b)
        run(ignorestatus(`diff -u --label golden --label new $ta $tb`))
        return
    end
    cmp(s) = norm ? normalize(s) : s
    only_g = sort!(collect(setdiff(keys(g), keys(n))))
    only_n = sort!(collect(setdiff(keys(n), keys(g))))
    both = intersect(keys(g), keys(n))
    differ = sort!([k for k in both if cmp(g[k].text) != cmp(n[k].text)])
    same = length(both) - length(differ)
    println("golden blocks: $(length(g))   new blocks: $(length(n))")
    println("identical: $same   differing: $(length(differ))   only in golden: $(length(only_g))   only in new: $(length(only_n))")
    bycat(keys) = begin
        d = Dict{String,Int}()
        for k in keys; d[split(k, ':')[1]] = get(d, split(k, ':')[1], 0) + 1; end
        d
    end
    println("  differing by kind: ", bycat(differ))
    println("  only-golden by kind: ", bycat(only_g), "   only-new by kind: ", bycat(only_n))
    # differing functions per golden file
    perfile = Dict{String,Int}()
    for k in differ
        perfile[g[k].file] = get(perfile, g[k].file, 0) + 1
    end
    println("  differing per file: ", sort(collect(perfile), by = last, rev = true))
    if cat
        cats = Dict{String,Vector{String}}()
        for k in differ
            startswith(k, "fn:") || continue
            c = categorize(g[k].text, n[k].text)
            push!(get!(cats, c, String[]), k[4:end])
        end
        println("\n--- categories:")
        for (c, v) in sort(collect(cats), by = x -> -length(x[2]))
            println(rpad(c, 20), length(v), "   e.g. ", join(first(v, 6), " "))
        end
        outdir = get(ENV, "GOLDEN_DIFF_OUT", tempdir())
        for (c, v) in cats
            write(joinpath(outdir, "cat_$c.txt"), join(sort(v), "\n") * "\n")
        end
        println("lists written to $outdir/cat_*.txt")
    end
    if list
        println("\n--- differing:"); foreach(println, differ)
        println("\n--- only in golden:"); foreach(println, only_g)
        println("\n--- only in new:"); foreach(println, only_n)
    end
end

main(ARGS)
