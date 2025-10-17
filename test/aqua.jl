import Aqua
using DiffFoldRNA

@testset "Aqua.test_all" begin
    showtestset()
    Aqua.test_all(DiffFoldRNA)
end
