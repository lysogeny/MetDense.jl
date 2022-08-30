struct PositionMethIterator
    mdf::MetDenseFile
    cpgInds::UnitRange
    chrom::String
    positions::Vector{UInt32}
    cells::Vector{UInt32}
    PositionMethIterator(
        mdf::MetDenseFile,
        cpgInds::UnitRange,
        chrom::String,
        positions::Vector{UInt32}
    ) = new(mdf, cpgInds, chrom, positions, 1:length(mdf.cell_names))
    PositionMethIterator(
        mi::PositionMethIterator, cells::Vector{Int}
    ) = new(mi.mdf, mi.cpgInds, mi.chrom, mi.positions, cells)
    PositionMethIterator(
        mi::PositionMethIterator, cells::BitVector
    ) = new(mi.mdf, mi.cpgInds, mi.chrom, mi.positions, collect(1:length(mi.mdf.cell_names))[cells])
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

function Base.getindex(x::MetDenseFile, gp::GenomicPosition)
    return FixedPositionMethIterator(x, get_position(x, gp), 1:length(x.cell_names))
end

function Base.getindex(x::MetDenseFile, cells::Vector{Int})
    return CellMethIterator(x, cells)
end

function Base.getindex(mi::PositionMethIterator, cells::Union{Vector{Int}, BitVector})
    return PositionMethIterator(mi, cells)
end

function Base.iterate(mi::PositionMethIterator)
    if length(mi.positions) == 0
        return nothing
    else 
        return FixedPositionMethIterator(
            mi.mdf,
            GenomicPosition(mi.chrom, mi.positions[1], mi.cpgInds.start),
            mi.cells
            ), 2
    end
end

function Base.iterate(mi::PositionMethIterator, state::Int64)
    if length(mi.positions) < state
        return nothing
    else
        return FixedPositionMethIterator(
            mi.mdf,
            GenomicPosition(mi.chrom, mi.positions[state], mi.cpgInds[state]),
            mi.cells
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
    if ((mi.cells[ind] - 1) ÷ 16 != (mi.cells[ind - 1] - 1) ÷ 16) || (mi.cells[ind] < mi.cells[ind - 1] )
        word = read_word(mi.mdf, mi.position, mi.cells[ind])
    else
        word = word >> ((mi.cells[ind] - mi.cells[ind - 1]) * 2)
    end
    return (call = MethCall(word & 0x03), cell = mi.cells[ind]), (word, ind + 1)
end

function Base.length(mi::FixedPositionMethIterator)
    return length(mi.cells)
end

function read_word(mdf::MetDenseFile, gp::GenomicPosition, cellInd)
    newPosition = mdf.chroms_filepos[gp.chrom].data[gp.ind] + floor(UInt64, cellInd/4)
    if position(mdf.f) != newPosition
        seek(mdf.f, newPosition)
    end
    return read(mdf.f, UInt32) >> (((cellInd - 1) % 16) * 2)
end
