module MetDense

include("read_metdense.jl")
include("iterators.jl")
include("make_metdense.jl")

export MetDenseFile, get_interval, seek, MethCall

end # module
