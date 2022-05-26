# This script checks a file (of the format used in the gastrulation data set )
# for correct sorting. Specify the filename as command line argument

import GZip

function check_sorted( filename )

    f = GZip.open( filename )
    readline( f )
    prev_chrom = ""
    while !eof(f)
        line = readline( f )
        if line == ""
            continue
        end
        firsttab = findfirst( '\t', line )
        if isnothing( firsttab )
            print( ">", line, "<" )
        end
        chrom = line[ 1 : firsttab-1 ]
        if chrom < prev_chrom
            return false
        end
        prev_chrom = chrom
    end
    true
end

sorted = check_sorted( ARGS[1] )
print( ARGS[1], " ", sorted, "\n" )