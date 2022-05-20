import GZip

@enum MethCall ::UInt32 begin
    nocall = 0x00
    unmeth = 0x01
    meth = 0x02
    ambig = 0x03
end

struct GenomicPosition
    chrom :: String
    pos :: Int32
end

struct EOFMarker
end

GenomicPositionOrEOF = Union{ GenomicPosition, EOFMarker }

function Base.isless( gp1 ::GenomicPosition, gp2 ::GenomicPosition )
    if gp1.chrom == gp2.chrom
        gp1.pos < gp2.pos
    else
        gp1.chrom < gp2.chrom
    end
end

Base.isless( gp1 ::EOFMarker, gp2 ::GenomicPositionOrEOF ) = false
Base.isless( gp1 ::GenomicPosition, gp2 ::EOFMarker ) = true

struct MethRecord
    gpos :: GenomicPositionOrEOF
    call :: MethCall   
end

function line_to_methrec( line )
    if line == ""
        return MethRecord( EOFMarker(), nocall )
    end
    fields = split( line, "\t" )
    gp = GenomicPosition( fields[1], parse( Int32, fields[2] ) )
    count_meth = parse( Int, fields[3] )
    count_unmeth = parse( Int, fields[4] )
    if count_meth == 0
        if count_unmeth == 0
            call = nocall
        else
            call = unmeth
        end
    else # count_meth > 0
        if count_unmeth == 0
            call = meth
        else
            call = ambig
        end
    end
    MethRecord( gp, call )
end

const data_block_pos_offset = 16

function write_header_block( fout )
    write( fout, "MetDense" )
    write( fout, UInt32(0) )  # major version
    write( fout, UInt32(0) )  # minor version
    @assert position(fout) == data_block_pos_offset
    write( fout, Int32(-1) )  # placeholder for position of Data block   
    write( fout, Int32(-1) )  # placeholder for position of Chromosomes block   
end

function write_cells_block( fout, cellnames )
    # Number of cells
    write( fout, UInt32( length(cellnames) ) )
    # Cell names
    for s in cellnames  
        write( fout, s, "\n" )
    end
    # Padding
    for i in 1:( position(fout) % 4 ) # Padding
        write( fout, UInt8(0) )
    end    
end

function handle_input_eof( inchannel:: Channel{MethRecord} ) ::Channel{MethRecord}
    f = function( ch::Channel )
        while isready( inch )
            put!( ch, take!(inchannel) )
        end
        put!( ch, MethRecord( GenomicPosition( nothing, 0 ), nocall ) )
    end
    Channel( f )
end

function write_data_block( fout, indata, tmp_filename )
    chrom_none_yet = "___none___just_starting___"
    word = UInt32(0)
    bitpos = 0
    current_recs :: Array{MethRecord} = take!.( indata )
    prev_chrom = chrom_none_yet
    fouttmp = open( tmp_filename, "w" )
    chroms = []
    while true

        # Get current position and write it out to temp file
        current_gpos = minimum( mr.gpos for mr in current_recs )

        # Are we starting a new chromosome?
        if current_gpos == EOFMarker() || current_gpos.chrom != prev_chrom
            if current_gpos == EOFMarker()
                break
            end
            @assert prev_chrom < current_gpos.chrom || prev_chrom == chrom_none_yet 
            if bitpos > 0
                # We have an unfinished word left over from the previous chromosome
                write( fout, UInt32(word) )
                word  = UInt32(0)
                bitpos = 0
            end
            print( "Processing chromosome $(current_gpos.chrom)\n")
            push!( chroms, ( name = current_gpos.chrom, filepos = position(fouttmp) ) )
            prev_chrom = current_gpos.chrom
        end

        # Get current position and write out current position to temp file
        write( fouttmp, current_gpos.pos )


        # Go through the cells and record calls for this position
        for i in 1:length(indata) 

            # Get call for current cell
            if current_recs[i].gpos != current_gpos
                @assert current_recs[i].gpos > current_gpos
                call = nocall
            else
                call = current_recs[i].call

                prev_pos = current_recs[i].gpos
                current_recs[i] = take!( indata[i] )
                if current_recs[i].gpos <= prev_pos
                    error( "File $i is not correctly sorted: " *
                       "$(prev_pos.chrom):$(prev_pos.pos) is followed by " *
                       "$(current_recs[i].gpos.chrom):$(current_recs[i].gpos.pos)." )
                end
            end

            # Add call to word
            word |= ( UInt32(call) << bitpos )
            bitpos += 2
    
            # Is the word full? If so, write it
            if bitpos > 30
                write( fout, UInt32(word) )
                word  = UInt32(0)
                bitpos = 0
            end

        end
    end    
    close( fouttmp )
    chroms
end

function copy_positions_block( fout, temp_filename )
    # We simply copy over the positions block from the tmpfile
    fin = open( temp_filename )
    buffer = Vector{UInt8}(undef, 8000)
    while( !eof(fin) )
        nb = readbytes!( fin, buffer, 8000 )
        write( fout, buffer[1:nb] )
    end
    close( fin )
end

function write_chromosomes_block( fout, chroms, start_positions_block )
    # Number of chromosomes
    write( fout, UInt32( length( chroms ) ) )
    for a in chroms
        write( fout, UInt32( a.filepos + start_positions_block ) )
    end
    for a in chroms
        write( fout, a.name, "\n" )
    end
end

function make_methrec_channel( fin )
    f = function(ch::Channel)
        while !eof( fin )
            put!( ch, line_to_methrec( readline( fin ) ) )
        end
    end
    Channel( f )
end

function make_metdense_file( outfilename, inputs, cellnames )
    fout = open( outfilename, "w" )
    temp_filename = outfilename * ".tmp"

    write_header_block( fout )   
    write_cells_block( fout, cellnames )    
    start_data_block = position( fout )
    chroms = write_data_block( fout, inputs, temp_filename )    
    start_positions_block = position( fout )
    copy_positions_block( fout, temp_filename )
    start_chromosomes_block = position( fout )
    write_chromosomes_block( fout, chroms, start_positions_block )
    seek( fout, data_block_pos_offset )
    write( fout, UInt32(start_data_block) )
    write( fout, UInt32(start_chromosomes_block) )

    close( fout )
end

function main()
    methcalls_dir = "/home/anders/w/metdense/gastrulation/raw_data"
    #methcalls_dir = "/home/anders/w/metdense/testshort/"
    methcalls_filenames = readdir( methcalls_dir )[1:40]
    cellnames = replace.( methcalls_filenames, ".tsv.gz"=>"")

    fins = GZip.open.( methcalls_dir * "/" .* methcalls_filenames )
    readline.( fins )  # Skip header
    inputs = make_methrec_channel.( fins )

    make_metdense_file( "test.metdense", inputs, cellnames )
end

main()