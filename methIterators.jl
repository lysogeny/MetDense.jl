struct MethIterator
    mdf::MetDenseFile
    cpgInd::UnitRange
    cells::Vector{UInt32}
    chrom::String
    chrPositions::Vector{Int64}
    rowwise::Bool
end
struct MethIteratorCols
    mdf::MetDenseFile
    cpgInd::UnitRange
    chrom::String
    chrPositions::Vector{Int64}
end
struct MethIteratorRows
    mdf::MetDenseFile
    cells::Vector{UInt32}
end

function Base.getindex(x::MetDenseFile, gi::GenomicInterval)
    cpgInd, pos = get_interval(x, gi)
    return MethIteratorCols(x, cpgInd, gi.chrom, pos)
end
function Base.getindex(x::MetDenseFile, cells::Vector{Int})
    return MethIteratorRows()
end
function Base.getindex(x::MetDenseFile, gi::Union{Colon, GenomicInterval},
    cells::Union{Colon, Vector{Int}})
    #do something with colons
    return MethIterator(x, get_interval(x, gi), cells)
end

function Base.iterate(mi::MethIteratorCols)
    return MethIterator(
        mi.mdf,
        mi.cpgInd[1]:mi.cpgInd[1],
        1:length(mi.mdf.cell_names),
        mi.chrom,
        [mi.chrPositions[1]],
        true), 2
end
function Base.iterate(mi::MethIteratorCols, state::Int64)
    if length(mi.cpgInd) < state
        return nothing
    else
        return MethIterator(
            mi.mdf,
            mi.cpgInd[state]:mi.cpgInd[state],
            1:length(mi.mdf.cell_names),
            mi.chrom,
            mi.chrPositions[state],
            true), state + 1
    end
end

function Base.iterate(mi::MethIterator)
    #assume rowwise for now
    word = read_word(mi.mdf, mi.cpgInd[1], mi.cells[1])
    return MethCall(word & 0x03), (mi, word, 1, 1)
end
function Base.iterate(mi::MethIterator, state)
    mi, word, cpgInd, cellInd = state
    cellInd += 1
    if cellInd > length(mi.cells) #go to the next row
        cpgInd += 1
        cellInd = mi.cells[1]
    end
    if cpgInd > length(mi.cpgInd)
        return nothing
    end
    #now we need to figure out whether to use the current word or read the new one
    #just read the new word for test purposes
    word = read_word(mi.mdf, mi.cpgInd[cpgInd], mi.cells[cellInd])
    return MethCall(word & 0x03), (mi, word, cpgInd, cellInd)
end

function read_word(mdf::MetDenseFile, cpgInd, cellInd)
    newPosition = mdf.offset_data_block +
        (cpgInd - 1) * ceil(UInt64, length(mdf.cell_names)/16) * 4 + ceil(UInt64, cellInd/16)
    if position(mdf.f) != newPosition
        seek(mdf.f, newPosition)
    end
    return read(mdf.f, UInt32) >> (((cellInd - 1) % 16) * 2)
end

df = MetDenseFile("data/test.metdense")
pIterator = df[GenomicInterval("2", (3058898, 4050898)), :]

for p in df[GenomicInterval("1", (3823430, 3823431))] #p corresponds to a single position and all cells
    #println(p)
    for m in p #c is a value for a specified cell and position
        print("$m ")
    end
end

get_interval(df, GenomicInterval("1", (3062977, 3823431)))

for m in df[GenomicInterval("2", (3058898, 4050898)), :]
end

for c in df[[1, 2, 3]]
    for m in c[GenomicInterval(...)]
    end
end

for f in df.cell_names
    fh = GZip.open("/home/tyranchick/mnt/mnt/raid/sveta/dcm/data/covs/$f")
    println(f)
    println(readline(fh))
    close(fh)
end
