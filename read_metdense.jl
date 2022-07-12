@enum MethCall ::UInt32 begin
    nocall = 0x00
    unmeth = 0x01
    meth = 0x02
    ambig = 0x03
end

struct GenomicPosition
    chrom :: String
    pos :: UInt32
    ind :: Int64
end
GenomicPosition(chrom::String, pos::Int) = GenomicPosition(chrom, UInt32(pos), -1)


struct GenomicInterval
    chrom :: String
    iv :: Tuple{UInt32, UInt32}
end

struct ChromFileposInfo
    pos ::StepRange{UInt64}  # File positions in Positions Block
    data ::StepRange{UInt64}  # File positions in Data Block
end

struct MetDenseFile
    f ::IOStream
    cell_names ::Vector{String}
    chroms_filepos ::Dict{ String, ChromFileposInfo }
    offset_data_block ::UInt64
end

function MetDenseFile( filename ::String )
    f = open( filename )

    # Read and check magic string
    s = Vector{UInt8}(undef,8)
    readbytes!( f, s )
    @assert s == b"MetDense"

    # Read version
    version = ( read( f, UInt32 ), read( f, UInt32 ) )
    @assert version == ( 0, 1 )

    # Get offsets
    offset_data_block = read( f, UInt64 )
    offset_chroms_block = read( f, UInt64 )

    # Read Cells block
    n_cells = read( f, UInt32 )
    names_cells = [ readline( f ) for i in 1:n_cells ]

    # Read Chromosomes block
    seek( f, offset_chroms_block )
    n_chroms = read( f, UInt32 )
    offsets_chroms = [ read( f, UInt64 ) for i in 1:n_chroms ]
    names_chroms = [ readline( f ) for i in 1:n_chroms ]
    chroms_filepos = Dict{ String, ChromFileposInfo }()
    start_in_data = offset_data_block
    rowlength = ceil( UInt64, n_cells / 16 ) * 4 # length of a site in the Data block
    
    for i in 1:n_chroms        

        filepos_in_positions_block = 
            range( offsets_chroms[i],
                i < n_chroms ? offsets_chroms[i + 1] - 1 : offset_chroms_block - 1;
                step = UInt64(4) )
        
        stop_in_data = start_in_data + rowlength * ( length( filepos_in_positions_block ) - 1 )
        filepos_in_data_block =
            range( start_in_data, stop_in_data; step = rowlength )
        start_in_data = stop_in_data + rowlength
        
        chroms_filepos[ names_chroms[i] ] = ChromFileposInfo( 
            filepos_in_positions_block, filepos_in_data_block )
    end            

    MetDenseFile( f, names_cells, chroms_filepos, offset_data_block )
end

function read_at_position( f, pos, type )
    seek( f, pos )
    read( f, type )
end

function get_interval( mdf::MetDenseFile, gi::GenomicInterval )
    start = searchsortedfirst(
        mdf.chroms_filepos[ gi.chrom ].pos, gi.iv[1];
        lt = (x, y) -> read_at_position( mdf.f, x, UInt32 ) < y )
    stop = searchsortedlast(
        mdf.chroms_filepos[ gi.chrom ].pos, gi.iv[2];
        lt = (x, y) -> x < read_at_position( mdf.f, y, UInt32 ) )
    if start > length(mdf.chroms_filepos[ gi.chrom ].pos)
        return start:stop, Vector{UInt32}(undef, 0)
    else
        seek( mdf.f, mdf.chroms_filepos[ gi.chrom ].pos[ start ])
        v = Vector{UInt32}( undef, (stop - start + 1) )
        read!( mdf.f, v )
        return start:stop, v
    end
end

# ok, this should be somehow reorganized not to keep file-related info and more general stuff together
function get_position(mdf::MetDenseFile, gp::GenomicPosition)
    if gp.ind !== -1
        return gp
    end
    ind, p = get_interval(mdf, GenomicInterval(gp.chrom, (gp.pos, gp.pos)))
    if length(p) == 0
        return nothing
    end
    @assert gp.pos == p[1]
    println(ind)
    return GenomicPosition(gp.chrom, gp.pos, ind[1])
end

function main_simon()
    mdf = MetDenseFile("data/gastrulation.metdense")
    iv, bpps = get_interval( mdf, GenomicInterval( "2", (4200000, 4700000) ))
    #println( iv )
    mypos = findfirst( x -> x == 4220582, bpps )
    println( "Now looking at position 2:", bpps[mypos] )
    seek( mdf.f, mdf.chroms_filepos["2"].data[ iv[mypos] ] )
    word = 99
    for cell in 1:length( mdf.cell_names )
        if (cell-1) % 16 == 0
            word = read( mdf.f, UInt32 )
        end
        call = MethCall( word & 0x03 )
        if call != nocall
            println( "  Cell $cell ($(mdf.cell_names[cell])): $call" )
        end
        word >>= 2
    end
end

#main_simon()
