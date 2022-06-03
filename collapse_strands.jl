using Pipe
using Profile

include("read_metdense.jl")
include("methIterators.jl")

genome_dir = "data/mouse_genome_m39/"
files = readdir(genome_dir)
genome_dict = Dict(map(x -> split(x, ".")[end - 1], files) .=> (genome_dir .* files))
genome_dict["__other"] = genome_dict["nonchromosomal"]

strandStep = UInt32(2)
input = "data/dcm.metdense"
output = "data/dcm_collapsed.metdense"
tempfile = "data/temp.metdense"

function find_chrom(name, genome_dict)
    if name in keys(genome_dict)
        genome = open(genome_dict[name])
        header = readline(genome)
    else
        genome = open(genome_dict["__other"])
        header = readline(genome)
        while !eof(genome) & (split(header[2:end], " ")[1] != name)
            header = readline(genome)
            while !eof(genome) & !startswith(header, ">")
                header = readline(genome)
            end
        end
    end
    @assert startswith(header, ">")
    @assert split(header[2:end], " ")[1] == name

    return genome
end

function write_line(pos::UInt32, data, fout, fout_tmp)
    write(fout, data)
    write(fout_tmp, pos)
end

function write_chrom(name, df, fout, fout_tmp, genome_dict, strandStep::UInt32, input)    
    genome = find_chrom(name, genome_dict)
    genome_starts = UInt32(position(genome))

    #data from df, positions from an extra IOStream
    posInput = open(input, "r")
    seek(posInput, df.chroms_filepos[name].pos[1])
    seek(df.f, df.chroms_filepos[name].data[1])
    collapsed = 0
    nonC = 0
    strayG = 0
    i = 0
    skip = false
    bpPos = -1
    C_data = Vector{UInt8}(undef, df.chroms_filepos[name].data.step)
    G_data = Vector{UInt8}(undef, df.chroms_filepos[name].data.step)

    while i < length(df.chroms_filepos[name].data)
    #while i < 1000
        if !skip
            bpPos = read(posInput, UInt32)
            read!(df.f, C_data)
        else
            skip = false
        end
        seek(genome, genome_starts + (bpPos - 1) + (bpPos - 1) ÷ 60)
        nct = read(genome, Char)

        if nct == 'C'
            nextBp = read(posInput, UInt32)
            read!(df.f, G_data)
            if nextBp - bpPos == strandStep
                C_data = C_data .| G_data
                write_line(bpPos, C_data, fout, fout_tmp)
                collapsed += 1
                i += 1
            else
                bpPos = nextBp
                C_data .= G_data
                skip = true
            end
        elseif nct == 'G'
            seek(genome,  genome_starts + (bpPos - strandStep - 1) + (bpPos - strandStep - 1) ÷ 60)
            pnct = read(genome, Char)
            if pnct != 'C'
                strayG += 1
            end
            write_line(bpPos - strandStep, C_data, fout, fout_tmp)
        else
            nonC += 1
        end
        i += 1
    end

    println("Chr. $name processed. $collapsed collapsed, $nonC non-C positions, $strayG unpaiered reverse-strand Cs")
    close(genome)
    close(posInput)

    return (name = name, filepos = position(fout_tmp))
end

df = MetDenseFile(input)

fout = open(output, "w")
fout_tmp = open(tempfile, "w")

write_header_block(fout) #from make_metdense
write_cells_block(fout, df.cell_names) #from make_metdense
start_data_block = position(fout)

chrom_names = df.chroms_filepos |>
    keys |> collect |> sort 

#temporary. Remove DCM manually. Later absent chromosomes should be copied without changes
chrom_names = chrom_names[chrom_names .!= "DCM"]

chroms = [write_chrom(chr, df, fout, fout_tmp, genome_dict, strandStep, input) for chr in chrom_names] #should return tuples (chromname, postion in temp)
close(fout_tmp)
start_positions_block = position(fout)
copy_positions_block(fout, tempfile) #from make_metdense
start_chromosomes_block = position(fout)
write_chromosomes_block(fout, chroms, start_positions_block) #from make_metdense
seek(fout, data_block_pos_offset)
write(fout, UInt64(start_data_block))
write(fout, UInt64(start_chromosomes_block))

close(fout)

# Play-around section (lot's of trying stuff)


function alternating(input, chrom_filepos)
    f1 = open(input, "r")
    f2 = open(input, "r")
    seek(f1, chrom_filepos.pos[1])
    seek(f2, chrom_filepos.data[1])
    for i in 1:length(chrom_filepos.pos)
        positon = read(f1, UInt32)
        data = read(f2, chrom_filepos.data.step)
    end
    close(f1)
    close(f2)
end

function doubleStream(input, chrom_filepos)
    f = open(input, "r")
    for i in 1:length(chrom_filepos.pos)
        seek(f, chrom_filepos.pos[i])
        position = read(f, UInt32)
        seek(f, chrom_filepos.data[i])
        data = read(f, chrom_filepos.data.step)
    end
    close(f)
end

using BenchmarkTools

@btime alternating($input, $(df.chroms_filepos["1"]))

@btime doubleStream($input, $(df.chroms_filepos["1"]))

g = find_chrom("GL456367.1", genome_dict)
g_start = position(g)
for i in 61:120
    seek(g, g_start + (i - 1) + (i - 1) ÷ 60 )
    print(read(g, Char))
end
println("")

function try_collapse()
    fout = open(output, "w")
    fout_tmp = open(tempfile, "w")
    
    write_header_block(fout) #from make_metdense
    write_cells_block(fout, df.cell_names) #from make_metdense
    start_data_block = position(fout)
    
    chrom_names = df.chroms_filepos |>
        keys |> collect |> sort 
    
    chroms = [write_chrom(chr, df, fout, fout_tmp, genome_dict, strandStep, input) for chr in chrom_names[1:2]] #should return tuples (chromname, postion in temp)
    chroms
end

using StatProfilerHTML
@profilehtml try_collapse()

tmp = open(tempfile, "r")

while !eof(tmp)
    println(convert(Int64, read(tmp, UInt32)))
end

p = get_interval(df, GenomicInterval("1", (194495341, 194495350)))

convert(Int64, p[2][1])

seek(df.f, df.chroms_filepos["1"].pos[end])
bpPos = convert(Int64, read(df.f, UInt32))