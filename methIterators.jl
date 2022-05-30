struct PositionMethIterator
    mdf::MetDenseFile
    cpgInds::UnitRange
    chrom::String
    positions::Vector{UInt32}
end
struct CellMethIterator
    mdf::MetDenseFile
    cells::Vector{UInt32}
end
struct FixedPositionMethIterator
    mdf::MetDenseFile
    position::GenomicPosition
    cells::Vector{UInt32}
end
struct FixedCellMethIterator
    mdf::MetDenseFile
    cell::UInt32
    chrom::String #I think, these three should form a separate struct
    positions::Vector{UInt32}
    cpgInds::UnitRange
end

function Base.getindex(x::MetDenseFile, gi::GenomicInterval)
    cpgInds, pos = get_interval(x, gi)
    return PositionMethIterator(x, cpgInds, gi.chrom, pos)
end
function Base.getindex(x::MetDenseFile, cells::Vector{Int})
    return CellMethIterator(x, cells)
end

function Base.iterate(mi::PositionMethIterator)
    return FixedPositionMethIterator(
        mi.mdf,
        GenomicPosition(mi.chrom, mi.positions[1], mi.cpgInds.start),
        1:length(mi.mdf.cell_names)
        ), 2
end
function Base.iterate(mi::PositionMethIterator, state::Int64)
    if length(mi.positions) < state
        return nothing
    else
        return FixedPositionMethIterator(
            mi.mdf,
            GenomicPosition(mi.chrom, mi.positions[state], mi.cpgInds[state]),
            1:length(mi.mdf.cell_names)
        ), state + 1
    end
end

function Base.iterate(mi::FixedPositionMethIterator)
    word = read_word(mi.mdf, mi.position, mi.cells[1])
    return (call = MethCall(word & 0x03), cell = mi.cells[1]), (word, 2)
end
function Base.iterate(mi::FixedPositionMethIterator, state::Tuple{UInt32, Int64})
    word, ind = state
    if length(mi.cells) < ind
        return nothing
    end
    if (mi.cells[ind] % 16  == 1) | (mi.cells[ind] - mi.cells[ind - 1] > 15)
        word = read_word(mi.mdf, mi.position, mi.cells[ind])
    else
        word = word >> ((mi.cells[ind] - mi.cells[ind - 1]) * 2)
    end
    return (call = MethCall(word & 0x03), cell = mi.cells[ind]), (word, ind + 1)
end

function read_word(mdf::MetDenseFile, gp::GenomicPosition, cellInd)
    newPosition = mdf.chroms_filepos[gp.chrom].data[gp.ind] + floor(UInt64, cellInd/4)
    if position(mdf.f) != newPosition
        seek(mdf.f, newPosition)
    end
    return read(mdf.f, UInt32) >> (((cellInd - 1) % 16) * 2)
end

df = MetDenseFile("data/test.metdense")
pIterator = df[GenomicInterval("2", (3058898, 4050898)), :]

for p in df[GenomicInterval("1", (3823430, 3823500))] #p corresponds to a single position and all cells
    println("Position $(p.position.pos)")
    for m in p #c is a value for a specified cell and position
        if m.call != nocall
            println("Cell: $(df.cell_names[m.cell]), call: $(m.call)")
        end
    end
    println("\n")
end

inds, pos = get_interval(df, GenomicInterval("1", (3823430, 3823500)))

MethCall(read_word(df, 2186, 1) & 0x03)


for m in df[GenomicInterval("2", (3058898, 4050898)), :]
end

for c in df[[1, 2, 3]]
    for m in c[GenomicInterval(...)]
    end
end
