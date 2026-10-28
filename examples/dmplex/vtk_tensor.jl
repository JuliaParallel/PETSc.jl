# EXCLUDE FROM TESTING
# Included by ex17.jl and ex62b.jl, not run on its own.

"""
    vtk_merge_tensor(fname, names...)

Post-process a VTK `.vtu` file written by PETSc's VTK viewer to merge 9 separate
scalar `DataArray`s named `name.0`…`name.8` (produced when a field has more than 3
components) into a single `NumberOfComponents="9"` DataArray, for each name given.
Also sets `Tensors="name1 name2 …"` on the `<PointData>` tag so ParaView treats the
arrays as 3×3 tensors. The file is rewritten in place.
"""
function vtk_merge_tensor(fname::AbstractString, names::AbstractString...)
    for name in names
        vtk_merge_one_tensor(fname, name)
    end
    return nothing
end

function vtk_merge_one_tensor(fname::AbstractString, name::AbstractString)
    isfile(fname) || return
    raw = read(fname)

    needle = b"<AppendedData"
    app_i  = let found = nothing
        for i in 1:length(raw)-length(needle)
            if @views raw[i:i+length(needle)-1] == needle
                found = i; break
            end
        end
        found
    end
    app_i === nothing && return
    gt  = findnext(isequal(UInt8('>')), raw, app_i + length(needle))
    gt  === nothing && return
    sep = findnext(isequal(UInt8('_')), raw, gt + 1)
    sep === nothing && return

    xml = String(raw[1:sep-1])
    bin = raw[sep+1:end]

    hdr_size = contains(xml, "header_type=\"UInt64\"") ? 8 : 4
    read_hdr(off) = hdr_size == 8 ?
        Int(read(IOBuffer(bin[off+1:off+8]), UInt64)) :
        Int(read(IOBuffer(bin[off+1:off+4]), UInt32))

    pat = Regex("""<DataArray\\b[^>]*\\bName=\"$(name)\\.(\\d+)\"[^>]*/?>""")
    ms  = collect(eachmatch(pat, xml))
    length(ms) == 9 || return

    entries = [(comp    = parse(Int, m[1]),
                xml_beg = m.offset,
                xml_end = m.offset + length(m.match) - 1,
                bin_off = parse(Int, match(r"offset=\"(\d+)\"", m.match)[1]))
               for m in ms]

    xml_first = minimum(e.xml_beg for e in entries)
    xml_last  = maximum(e.xml_end for e in entries)

    sorted   = sort(entries; by = e -> e.comp)
    offsets  = [e.bin_off for e in sorted]

    val_size = contains(ms[1].match, "Float64") ? 8 : 4
    RT       = val_size == 8 ? Float64 : Float32

    function read_scalar(off)
        nbytes = read_hdr(off)
        nvals  = nbytes ÷ val_size
        io     = IOBuffer(bin[off+hdr_size+1:off+hdr_size+nbytes])
        [read(io, RT) for _ in 1:nvals]
    end

    scalars = [read_scalar(off) for off in offsets]
    npoints = length(scalars[1])

    merged = Vector{RT}(undef, npoints * 9)
    for i in 0:npoints-1, k in 0:8
        merged[i*9+k+1] = scalars[k+1][i+1]
    end

    merged_nbytes = UInt64(npoints * 9 * val_size)
    merged_block  = let buf = IOBuffer()
        hdr_size == 8 ? write(buf, merged_nbytes) : write(buf, UInt32(merged_nbytes))
        for v in merged; write(buf, v); end
        take!(buf)
    end

    last_off    = offsets[9]
    old_bin_end = last_off + hdr_size + npoints * val_size
    Δ           = length(merged_block) - (old_bin_end - offsets[1])

    new_bin  = vcat(bin[1:offsets[1]], merged_block, bin[old_bin_end+1:end])

    type_str = RT == Float64 ? "Float64" : "Float32"
    new_da   = "<DataArray Name=\"$name\" NumberOfComponents=\"9\"" *
               " type=\"$type_str\" format=\"appended\" offset=\"$(offsets[1])\"/>"

    post = xml[xml_last+1:end]
    post_fixed = replace(post, r"offset=\"(\d+)\"" => function(s)
        m2  = match(r"\"(\d+)\"", s)
        off = parse(Int, m2[1])
        off > last_off ? "offset=\"$(off + Δ)\"" : s
    end)

    pre = xml[1:xml_first-1]
    # Julia's replace passes SubString (not RegexMatch) to the function, so we
    # must call match() explicitly to extract capture groups.
    pre = replace(pre, r"<PointData\b([^>]*)>" =>
        function(s)
            if contains(s, "Tensors=")
                m2 = match(r"Tensors=\"([^\"]*)\"", String(s))
                replace(String(s), m2.match => "Tensors=\"$(m2[1]) $name\"")
            else
                replace(String(s), "<PointData" => "<PointData Tensors=\"$name\"")
            end
        end)

    open(fname, "w") do io
        write(io, pre)
        write(io, new_da)
        write(io, post_fixed)
        write(io, UInt8('_'))
        write(io, new_bin)
    end
    return nothing
end
