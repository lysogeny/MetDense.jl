module MetDense

include("positions.jl")
include("read_metdense.jl")
include("iterators.jl")
include("make_metdense.jl")

export MetDenseFile, 
    get_interval, 
    MethCall, 
    make_metdense_file,
    extend,
    GenomicInterval

end # module
