struct GenomicPosition
    chrom :: String
    pos :: UInt32
    ind :: Int64
end
GenomicPosition(chrom::String, pos::Int) = GenomicPosition(chrom, UInt32(pos), -1)

function Base.isless( gp1 ::GenomicPosition, gp2 ::GenomicPosition )
    if gp1.chrom == gp2.chrom
        gp1.pos < gp2.pos
    else
        gp1.chrom < gp2.chrom
    end
end

struct EOFMarker
end

GenomicPositionOrEOF = Union{ GenomicPosition, EOFMarker }


Base.isless( gp1 ::EOFMarker, gp2 ::GenomicPositionOrEOF ) = false
Base.isless( gp1 ::GenomicPosition, gp2 ::EOFMarker ) = true

struct GenomicInterval
    chrom :: String
    iv :: Tuple{UInt32, UInt32}
end

function extend(gi::GenomicInterval; left = 0, right = 0)
    if Int64(gi.iv[2]) - Int64(gi.iv[1]) - left + right <= 0
        @warn "The resulting interval length is smaller than zero"
        return gi
    end

    if gi.iv[1] - left <= 0
        left = gi.iv[1] - 1
    end

    return GenomicInterval(gi.chrom, (gi.iv[1] - left, gi.iv[2] + right))
end

