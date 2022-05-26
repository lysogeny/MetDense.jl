# Reads a methylation-calls file (of the format used in the gastrulations dataset)
# from stdin, sorts it, and prints it to stdout. Use only with small files.

function by( line )
    fields = split( line, "\t" )
    fields[1], parse( Int, fields[2] )
end

function sort()

    header = readline(stdin)
    a = readlines(stdin)
    sort!( a, by = by )
    println(header)
    println.(a)

end

sort()