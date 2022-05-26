# First attempt to get a smoothed heatmap

using Plots

include( "read_metdense.jl" )

mdf = MetDenseFile( "data/gastrulation.metdense" )

n_cells = length( mdf.cell_names )

iv, bpps = get_interval( mdf, GenomicInterval( "10", ( 61383530 - 50000, 61383530 + 50000 ) ) )


seek( mdf.f, mdf.chroms_filepos["10"].data[ iv.start ] )
n_calls = zeros( length(bpps) )
n_meth = zeros( length(bpps) )
for i in 1:length(bpps)
    word = 99
    @assert position(mdf.f) == mdf.chroms_filepos["10"].data[ iv.start+i-1 ]
    for cell in 1:length( mdf.cell_names )
        if (cell-1) % 16 == 0
            word = read( mdf.f, UInt32 )
        end
        call = MethCall( word & 0x03 )
        if call != nocall && call != ambig 
            n_calls[i] += 1
            if call == meth
                n_meth[i] += 1
            end
        end
        word >>= 2
    end
    n_calls / n_meth
end

scatter( bpps, n_meth ./ n_calls )


function tricube(x)
    if abs(x) < 1
        ( 1 - abs(x)^3 )^3
    else
        0
    end
end


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

pt_grid = range( minimum( filter( x -> !isnan(x), pseudotime ) ), maximum( filter( x -> !isnan(x), pseudotime ) ), length=100 )

numerators = fill( 0., ( length(bpps), length(pt_grid) ) )
denominators = fill( 0., ( length(bpps), length(pt_grid) ) )

bw_pt = 2.0
seek( mdf.f, mdf.chroms_filepos["10"].data[ iv.start ] )
for i in 1:length(bpps)
    word = 99
    @assert position(mdf.f) == mdf.chroms_filepos["10"].data[ iv.start+i-1 ]
    for cell in 1:length( mdf.cell_names )
        if (cell-1) % 16 == 0
            word = read( mdf.f, UInt32 )
        end
        call = MethCall( word & 0x03 )
        if call != nocall && call != ambig 
            if !isnan( pseudotime[ cell ] )
                weights = tricube.( ( pt_grid .- pseudotime[ cell ] ) ./ bw_pt )
                @assert any( weights .!= 0 )
                denominators[ i, : ] .+= weights
                if call == meth
                    numerators[ i, : ] .+= weights
                end
            end

        end
        word >>= 2
    end
end

heatmap( bpps, pt_grid, transpose( numerators ./ denominators ) )


bppos_grid = range( minimum(bpps), maximum(bpps), length=300 )
bw_pos = 800
numerators2 = fill( 0., ( length(bppos_grid), length(pt_grid) ) )
denominators2 = fill( 0., ( length(bppos_grid), length(pt_grid) ) )
for i in 1:length(bpps)
    weights = tricube.( ( bppos_grid .- Float64(bpps[i]) ) ./ bw_pos )
    for j in 1:length(bppos_grid)
        if weights[j] > 0
            numerators2[j,:] .+= numerators[ i, : ] * weights[j]
            denominators2[j,:] .+= denominators[ i, : ] * weights[j]
        end
    end
end

heatmap( bppos_grid, pt_grid, transpose( numerators2 ./ denominators2 ) )