using Test
using VectorDataCubes

using Rasters, DimensionalData
using Rasters.Lookups
import DimensionalData as DD
import GeometryOps as GO, GeoInterface as GI
using Extents

# Four hand-made squares with exactly known spatial relations:
# sq2 touches sq1 along the edge x = 1, sq3 lies strictly within sq1,
# and sq4 is disjoint from everything else.
_square(x1, y1, x2, y2) =
    GI.Polygon([GI.LinearRing([(x1, y1), (x2, y1), (x2, y2), (x1, y2), (x1, y1)])])

sq1 = _square(0.0, 0.0, 1.0, 1.0)
sq2 = _square(1.0, 0.0, 2.0, 1.0)
sq3 = _square(0.25, 0.25, 0.75, 0.75)
sq4 = _square(10.0, 10.0, 11.0, 11.0)
squares = [sq1, sq2, sq3, sq4]

# STRtree query order is unspecified, so multi-index results are sorted before comparison.
selinds(gl, sel) = sort(Lookups.selectindices(gl, sel))

@testset "selectors on hand-made squares (tree = $treedesc)" for (treedesc, treekw) in
    (("STRtree", (;)), ("nothing", (; tree = nothing)))

    gl = GeometryLookup(squares; treekw...)
    dv = rand(Geometry(gl))

    @testset "standard indices pass through" begin
        @test dv[Geometry = 2] == dv[2]
        @test dv[Geometry = 1:2] == dv[1:2]
        @test Lookups.selectindices(gl, 3) == 3
        @test Lookups.selectindices(gl, 1:2) == 1:2
    end

    @testset "Contains(point)" begin
        @test selinds(gl, Contains((0.5, 0.5))) == [1, 3]
        @test selinds(gl, Contains((1.5, 0.5))) == [2]
        @test selinds(gl, Contains((50.0, 50.0))) == Int[]
        @test dv[Geometry(Contains((1.5, 0.5)))] == dv[[2]]
        @test dv[Geometry = Contains((1.5, 0.5))] == dv[[2]]
    end

    @testset "At(geometry)" begin
        @test Lookups.selectindices(gl, At(sq2)) == 2
        # equality is geometric (GO.equals), not object identity
        @test Lookups.selectindices(gl, At(_square(1.0, 0.0, 2.0, 1.0))) == 2
        @test dv[Geometry(At(sq3))] == dv[3]
        @test_throws ArgumentError Lookups.selectindices(gl, At(_square(5.0, 5.0, 6.0, 6.0)))
    end

    @testset "Near(point)" begin
        @test Lookups.selectindices(gl, Near((10.4, 10.5))) == 4
        @test Lookups.selectindices(gl, Near((2.5, 0.5))) == 2
        @test dv[Geometry(Near((2.5, 0.5)))] == dv[2]
        # only point geometries are supported for Near so far
        @test_throws AssertionError Lookups.selectindices(gl, Near(sq1))
    end

    @testset "(X(At), Y(At)) point lookup" begin
        @test Lookups.selectindices(gl, (X(At(1.5)), Y(At(0.5)))) == 2
        @test dv[Geometry = (X(At(1.5)), Y(At(0.5)))] == dv[2]
        # At requires an exact match, so points in no geometry are errors,
        # whether inside the overall extent or not
        @test_throws ArgumentError Lookups.selectindices(gl, (X(At(5.0)), Y(At(5.0))))
        @test_throws ArgumentError Lookups.selectindices(gl, (X(At(50.0)), Y(At(50.0))))
    end

    @testset "(X(Contains), Y(Contains)) point lookup" begin
        @test selinds(gl, (X(Contains(0.5)), Y(Contains(0.5)))) == [1, 3]
        @test selinds(gl, (X(Contains(50.0)), Y(Contains(50.0)))) == Int[]
        @test length(dv[Geometry = (X(Contains(0.5)), Y(Contains(0.5)))]) == 2
    end

    @testset "(X(Touches), Y(Touches)) extent lookup" begin
        @test selinds(gl, (X(Touches(0.9, 1.1)), Y(Touches(0.4, 0.6)))) == [1, 2]
        @test selinds(gl, (X(Touches(9.0, 12.0)), Y(Touches(9.0, 12.0)))) == [4]
        @test selinds(gl, (X(Touches(50.0, 60.0)), Y(Touches(50.0, 60.0)))) == Int[]
    end

    @testset "(X(a .. b), Y(a .. b)) interval covers lookup" begin
        @test selinds(gl, (X(-0.1 .. 2.1), Y(-0.1 .. 1.1))) == [1, 2, 3]
        @test selinds(gl, (X(-0.1 .. 1.1), Y(-0.1 .. 1.1))) == [1, 3]
        # the interval only intersecting a geometry is not enough - it must cover it
        @test selinds(gl, (X(0.4 .. 0.6), Y(0.4 .. 0.6))) == Int[]
        @test dv[Geometry = (X(-0.1 .. 2.1), Y(-0.1 .. 1.1))] == dv[1:3]
    end

    @testset "Touches(extent)" begin
        @test selinds(gl, Touches(GI.extent(sq1))) == [1, 2, 3]
        @test selinds(gl, Touches(GI.extent(sq4))) == [4]
    end

    @testset "Where with GeometryOps predicates" begin
        @test selinds(gl, Where(GO.intersects(sq1))) == [1, 2, 3]
        # GO.equals has no curried form, so build the Fix2 directly
        @test selinds(gl, Where(Base.Fix2(GO.equals, sq2))) == [2]
        @test selinds(gl, Where(GO.contains(sq3))) == [1, 3]
        @test selinds(gl, Where(GO.within(sq1))) == [1, 3]
        @test selinds(gl, Where(GO.covers(sq3))) == [1, 3]
        @test selinds(gl, Where(GO.coveredby(sq1))) == [1, 3]
        @test selinds(gl, Where(GO.touches(sq1))) == [2]
        @test selinds(gl, Where(GO.disjoint(sq1))) == [4]
        @test dv[Geometry = Where(GO.disjoint(sq1))] == dv[[4]]
    end

    @testset "empty selections produce empty arrays" begin
        empty_dv = dv[Geometry(Contains((50.0, 50.0)))]
        @test isempty(empty_dv)
        @test empty_dv isa DD.AbstractDimVector
    end
end
