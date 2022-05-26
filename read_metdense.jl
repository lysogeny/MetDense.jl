@enum MethCall ::UInt32 begin
    nocall = 0x00
    unmeth = 0x01
    meth = 0x02
    ambig = 0x03
end

struct GenomicPosition
    chrom :: String
    pos :: UInt32
end

struct GenomicInterval
    chrom :: String
    iv :: Tuple{UInt32, UInt32}
end

struct MetDenseFile
    f ::IOStream
    cell_names ::Vector{String}
    chroms_dict ::Dict{ String, StepRange{UInt64} }
    offset_data_block ::UInt64
end

function MetDenseFile( filename ::String )
    f = open( filename )

    # Read and check magic string
    s = Vector{UInt8}(undef,8)
    readbytes!( f, s )
    @assert s == Vector{UInt8}("MetDense")

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

    chroms_dict = Dict( names_chroms[i] =>
        range( offsets_chroms[i],
            i < n_chroms ? offsets_chroms[i + 1] - 1 : offset_chroms_block - 1;
            step = UInt64(4))
        for i = 1:n_chroms )

    MetDenseFile( f, names_cells, chroms_dict, offset_data_block )
end

function read_at_position( f, pos, type )
    seek( f, pos )
    read( f, type )
end

function get_interval( mdf::MetDenseFile, gi::GenomicInterval )
    start = searchsortedfirst(
        mdf.chroms_dict[ gi.chrom ], gi.iv[1];
        lt = (x, y) -> read_at_position( mdf.f, x, UInt32 ) < y )
    stop = searchsortedlast(
        mdf.chroms_dict[ gi.chrom ], gi.iv[2];
        lt = (x, y) -> x < read_at_position( mdf.f, y, UInt32 ) )
    seek( mdf.f, mdf.chroms_dict[ gi.chrom ][start])
    v = Vector{UInt32}( undef, (stop - start + 1) )
    read!( mdf.f, v )
    return start:stop, v
end


df = MetDenseFile("data/test.metdense")
int, pos = get_interval( df, GenomicInterval( "2", (3058898, 4050898) ))
