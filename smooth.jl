# First attempt to get a smoothed heatmap

using Plots

include( "read_metdense.jl" )
include( "methIterators.jl" )

# Load metdense file
mdf = MetDenseFile( "data/gastrulation.metdense" )

# smoothing kernel
function tricube(x)
    if abs(x) < 1
        ( 1 - abs(x)^3 )^3
    else
        0
    end
end

# Load table assigning pseudotimes to cells
pseudotime = fill( NaN, length(mdf.cell_names) )
fpt = open( "data/DPT_to_mesoderm.csv" )
readline( fpt )
while !eof( fpt )
    fields = split( readline( fpt ), "," )
    if fields[2] != "NA"
        cellidx = findfirst( isequal(fields[1]), mdf.cell_names )
        if !isnothing( cellidx )
            pseudotime[ cellidx ] = parse( Float32, fields[2] )
        end
    end
end


# A grid of time points, and a bandwidth
pt_grid = range( minimum( filter( x -> !isnan(x), pseudotime ) ), maximum( filter( x -> !isnan(x), pseudotime ) ), length=100 )
bw_pt = 2.0

# Smoothing across cells

posIterator = mdf[GenomicInterval( "10", ( 61383530 - 50000, 61383530 + 50000 ) )]

numerators = fill( 0., ( length(posIterator.positions), length(pt_grid) ) )
denominators = fill( 0., ( length(posIterator.positions), length(pt_grid) ) )

for (i, pos) in enumerate(posIterator[.!isnan.(pseudotime)])
    for m in pos
        if m.call != nocall && m.call != ambig 
            weights = tricube.( ( pt_grid .- pseudotime[m.cell] ) ./ bw_pt )
            @assert any( weights .!= 0 )
            denominators[ i, : ] .+= weights
            if m.call == meth
                numerators[ i, : ] .+= weights
            end
        end
    end
end

heatmap( posIterator.positions, pt_grid, transpose( numerators ./ denominators ) )


# A grid of positions, and a smoothing bandwidth
bppos_grid = range( minimum(posIterator.positions), maximum(posIterator.positions), length=300 )
bw_pos = 800

# Now smoothing across position
numerators2 = fill( 0., ( length(bppos_grid), length(pt_grid) ) )
denominators2 = fill( 0., ( length(bppos_grid), length(pt_grid) ) )
for i in 1:length(posIterator.positions)
    weights = tricube.( ( bppos_grid .- Float64(posIterator.positions[i]) ) ./ bw_pos )
    for j in 1:length(bppos_grid)
        if weights[j] > 0
            numerators2[j,:] .+= numerators[ i, : ] * weights[j]
            denominators2[j,:] .+= denominators[ i, : ] * weights[j]
        end
    end
end

# New heatmap
heatmap( bppos_grid, pt_grid, transpose( numerators2 ./ denominators2 ) )

