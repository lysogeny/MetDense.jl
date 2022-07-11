import GZip
import Pipe: @pipe

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

struct DataCols
    chrom::Int16
    pos::Int16
    meth::Int16
    unmeth::Int16
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

function line_to_methrec( line, data_cols, round )
    if line == ""
        return MethRecord( EOFMarker(), nocall )
    end
    fields = split( line, "\t" )
    gp = GenomicPosition( fields[data_cols.chrom], parse( UInt32, fields[data_cols.pos] ) )
    count_meth = parse( Int, fields[data_cols.meth] )
    count_unmeth = parse( Int, fields[data_cols.unmeth] )
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
            if round == "floor"
                call = unmeth
            elseif round == "ceil"
                call = meth
            elseif round == "round"
                if count_meth > count_unmeth 
                    call = meth
                elseif count_meth < count_unmeth
                    call = unmeth
                else
                    call = ambig
                end
            else
                call = ambig
            end
        end
    end
    MethRecord( gp, call )
end

const data_block_pos_offset = 16

function write_header_block( fout )
    write( fout, "MetDense" )
    write( fout, UInt32(0) )  # major version
    write( fout, UInt32(1) )  # minor version
    @assert position(fout) == data_block_pos_offset
    write( fout, Int64(-1) )  # placeholder for position of Data block
    write( fout, Int64(-1) )  # placeholder for position of Chromosomes block
end

function write_cells_block( fout, cellnames )
    # Number of cells
    write( fout, UInt32( length(cellnames) ) )
    # Cell names
    for s in cellnames
        write( fout, s, "\n" )
    end
    # Padding
    for i in 1:( 4 - position(fout) % 4 ) # Padding
        write( fout, UInt8(0) )
    end
end

function write_data_block( fout, indata, tmp_filename )
    chrom_none_yet = "___none___just_starting___"  # placeholder constant
    prev_chrom = chrom_none_yet
    fouttmp = open( tmp_filename, "w" )
    chroms = []
    numdots = 0
    word = UInt32(0)
    bitpos = 0
    current_recs :: Array{MethRecord} = take!.( indata )
    @assert position(fout) % 4 == 0

    while true

        # Get current position and write it out to temp file
        current_gpos = minimum( mr.gpos for mr in current_recs )

        # Are we starting a new chromosome?
        if current_gpos == EOFMarker() || current_gpos.chrom != prev_chrom
            if current_gpos == EOFMarker()
                break
            end
            @assert prev_chrom < current_gpos.chrom || prev_chrom == chrom_none_yet
            print( "\nProcessing chromosome $( rpad( current_gpos.chrom, 3) ) ")
            push!( chroms, ( name = current_gpos.chrom, filepos = position(fouttmp) ) )
            prev_chrom = current_gpos.chrom
            numdots = 0
        end

        # Progress indication
        if current_gpos.pos ÷ 10000000 > numdots
            newdots = current_gpos.pos ÷ 10000000 - numdots
            print( "." ^ newdots )
            numdots += newdots
        end

        # Write out current position to temp file, for later includion into Positions block
        write( fouttmp, current_gpos.pos )

        # Go through the cells and record calls for this position
        for i in 1:length(indata)

            # Get methylation call for current cell
            if current_recs[i].gpos != current_gpos
                @assert current_recs[i].gpos > current_gpos
                call = nocall
            else
                call = current_recs[i].call

                # Check for correct sorting
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
            if bitpos >= 32
                write( fout, UInt32(word) )
                word  = UInt32(0)
                bitpos = 0
            end

        end

        # Unless we have just written out the word, we still need to do that.
        if bitpos > 0
            write( fout, UInt32(word) )
            word  = UInt32(0)
            bitpos = 0
        end

    end
    close( fouttmp )
    print( "\n" )
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
    rm( temp_filename )
end

function write_chromosomes_block( fout, chroms, start_positions_block )
    # Number of chromosomes
    write( fout, UInt32( length( chroms ) ) )
    for a in chroms
        write( fout, UInt64( a.filepos + start_positions_block ) )
    end
    for a in chroms
        write( fout, a.name, "\n" )
    end
end

function make_methrec_channel( fin, data_cols, round )
    f = function(ch::Channel)
        while !eof( fin )
            put!( ch, line_to_methrec( readline( fin ), data_cols, round ) )
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
    write( fout, UInt64(start_data_block) )
    write( fout, UInt64(start_chromosomes_block) )

    close( fout )
end

function main(methcalls_dir::String, output::String, data_cols::DataCols; 
    file_end = "", skip_header = true, round = "none")
    # data_cols specifies the position of columns
    # 1 - chromosome
    # 2 - position
    # 3 - number of methylated reads
    # 4 - number of unmethylated reads
    # round can be "none", "floor", "ceil", "round"

    methcalls_filenames = @pipe readdir( methcalls_dir ) |>
            filter(x -> occursin(Regex(file_end * "\$"), x), _)
    cellnames = replace.( methcalls_filenames, file_end => "")
    
    fins = GZip.open.( methcalls_dir * "/" .* methcalls_filenames )
    if(skip_header)
        readline.( fins )  # Skip header
    end
    inputs = make_methrec_channel.( fins, Ref(data_cols), Ref(round) )

    make_metdense_file( output, inputs, cellnames )
end

main(
    "/home/tyranchick/mnt/mnt/raid/sveta/dcm/data/covs", 
    "test.metdense", 
    DataCols(1, 2, 5, 6), 
    file_end = "-DCM.cov.gz", 
    skip_header = false)
