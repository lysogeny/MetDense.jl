function read_value(mdf::MetDenseFile, cpgInd, cellInd)
    n = length(mdf.cell_names)
    seek(mdf.f, mdf.offset_data_block + (cpgInd - 1) * ceil(Int64, n/16) * 4 + ceil(Int64, cellInd/16))
    word = read(mdf.f, UInt32)
    return meth .* word .- unmeth .* word
end

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
function Base.getindex(x::MethDeseFile, gi::Union{Colon, GenomicInterval},
    cells::Uninon{Colon, Vector{Int}})
    #do something with colons
    return MethIterator(x, get_interval(x, gi), cells)
end

function Base.iterate(mi::MethIteratorCols)
    return MethIterator(
        mi.mdf,
        mi.cpgInd[1]:mi.cpgInd[1],
        1:length(mi.mdf.cell_names),
        mi.chrom,
        mi.chrPositions[1],
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

function Base.iterate(mi::MethIterator) = read_value(mi.mdf, mi.cpgInd[1], mi.cells[1])
function Base.iterate(mi::MethIterator, state)

end


df = MetDenseFile("data/test.metdense")
pIterator = df[GenomicInterval("2", (3058898, 4050898)), :]

for p in df[GenomicInterval("2", (3058898, 4050898))] #p corresponds to a single position and all cells
    for m in p[cells] #c is a value for a specified cell and position
    end
end

for m in df[GenomicInterval("2", (3058898, 4050898)), :]
end

for c in df[[1, 2, 3]]
    for m in c[GenomicInterval(...)]
    end
end
