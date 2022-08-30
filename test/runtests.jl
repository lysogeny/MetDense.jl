using MetDense
using Test

testfiles = [
    "create_file.jl"
]

for testfile in testfiles
    try 
        include(testfile)
        println("PASSED $testfile")
    catch e
        println("FAILED $testfile")
        rethrow(e)
    end
end
