struct TestData
    chromosomes::Vector{String}
    positions::Vector{Vector{Int64}}
end

function TestData(chromosomes::Vector{String}, length::Int)
    chromosomes = sort(chromosomes)
    positions = [
        Vector{Int64}(cumsum(1 .+ rand(UInt8, Int(rand(UInt16) % length))))
        for _ in chromosomes
    ]
    TestData(chromosomes, positions)
end

function channel_function(td::TestData)
    function content(x::Channel)
        for (i, chromosome) in enumerate(td.chromosomes)
            positions = td.positions[i]
            for position in positions
                genomic_pos = MetDense.GenomicPosition(chromosome, position)
                put!(x, MetDense.MethRecord(genomic_pos, MetDense.meth))
            end
        end
        put!(x, MetDense.MethRecord(MetDense.EOFMarker(), MetDense.nocall))
        # Last thing to be sent needs to be EOFMarker?!
    end
    Channel(content)
end

@testset "Write single-celled file" begin
    data = TestData(["1", "11", "2", "3", "5", "MT", "X"], 100)
    fname = tempname()
    make_metdense_file(fname, [channel_function(data)], ["AS-hello-LR-world"])
    @test isfile(fname)
    file = MetDenseFile(fname)
    @test file.cell_names == ["AS-hello-LR-world"]
end
