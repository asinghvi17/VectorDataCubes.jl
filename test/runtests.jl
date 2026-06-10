using VectorDataCubes
using Test

@testset "VectorDataCubes.jl" begin
    include("construction.jl")
    include("selectors.jl")
    include("basics.jl")
end
