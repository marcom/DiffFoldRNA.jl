using Test

# show which testset is currently running
showtestset() = println(" "^(2 * Test.get_testset_depth()), "testing ",
                        Test.get_testset().description)

@testset "DiffFoldRNA" verbose=true begin
    showtestset()
    include("aqua.jl")
    include("common.jl")
    include("brute_force.jl")
    include("energy.jl")
    include("nussinov.jl")
    include("sampling.jl")
    include("vienna.jl")
end
