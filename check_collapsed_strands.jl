orig_df = MetDenseFile("/home/tyranchick/mnt/mnt/raid/sveta/dcm/data/dcm_cpg.metdense")
collapsed_df = MetDenseFile("/home/tyranchick/mnt/mnt/raid/sveta/dcm/data/dcm_cpg_collapsed.metdense")

#check that collapsed file has no neighbouring positions
nonzero = 0
for (name, chr) in collapsed_df.chroms_filepos
    println(name)
    if length(chr.pos) > 0
        seek(collapsed_df.f, chr.pos[1])
        prevPos = read(collapsed_df.f, UInt32)
        i = 1
        checked = 0
        neighbours = 0
        while i < length(chr.pos)
            i += 1
            curPos = read(collapsed_df.f, UInt32)
            if curPos - prevPos == UInt32(1)
                neighbours += 1
            end
            checked += 1
            prevPos = curPos
        end
        println("$checked positions checked, $neighbours neighbouring ones found")
        if neighbours > 0
            nonzero += 1
        end
    end
end
println("total nonzero: $nonzero")

#check that for every position in the original file there is one in the collapsed file
nonzero = 0
for (name, chr) in collapsed_df.chroms_filepos
    print("Chromosome $name ")
    notFound = 0
    if length(chr.pos) > 0
        seek(collapsed_df.f, chr.pos[1])
        seek(orig_df.f, orig_df.chroms_filepos[name].pos[1])
        i = 1
        colPos = read(collapsed_df.f, UInt32)
        while i <= length(orig_df.chroms_filepos[name].pos)
            origPos = read(orig_df.f, UInt32)
            if (origPos != colPos) & (origPos - UInt32(1) != colPos) & !eof(collapsed_df.f)
                colPos = read(collapsed_df.f, UInt32)
            end
            if (origPos != colPos) & (origPos - UInt32(1) != colPos)
                notFound += 1
            end
            i += 1
        end
    end
    if notFound > 0
        nonzero += 1
    end
    print("$notFound positions are missing\n")
end
println("Check for $nonzero chromosomes failed")

#check that for every position in the collapsed file there is a corresponding one in the original
nonzero = 0
for (name, chr) in collapsed_df.chroms_filepos
    print("Chromosome $name ")
    notFound = 0
    if length(chr.pos) > 0
        seek(collapsed_df.f, chr.pos[1])
        seek(orig_df.f, orig_df.chroms_filepos[name].pos[1])
        i = 1
        origPos = read(orig_df.f, UInt32)
        while i <= length(chr.pos)
            colPos = read(collapsed_df.f, UInt32)
            flag = false
            while ((origPos == colPos) | (origPos - UInt32(1) == colPos)) & !eof(orig_df.f)
                flag = true
                origPos = read(orig_df.f, UInt32)
            end
            if !flag
                notFound += 1
            end
            i += 1
        end
    end
    if notFound > 0
        nonzero += 1
    end
    print("$notFound positions are missing\n")    
end
println("Check for $nonzero chromosomes failed")

int_col = get_interval(collapsed_df, GenomicInterval("3", (1, 4000000)))
int_orit = get_interval(orig_df, GenomicInterval("3", (1, 4000000)))

gp = get_position(collapsed_df, GenomicPosition("3", 3000042))

for m in orig_df[GenomicPosition("3", 3000042)]
    print("$(m.call) ")
end

get_interval(collapsed_df, GenomicInterval("3", (3000042, 3000042)))