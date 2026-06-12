using Test

using Rasters, DimensionalData
using Rasters.Lookups
import DimensionalData as DD
import GeometryOps as GO, GeoInterface as GI
import SortTileRecursiveTree

_csquare(x1, y1, x2, y2) =
    GI.Polygon([GI.LinearRing([(x1, y1), (x2, y1), (x2, y2), (x1, y2), (x1, y1)])])

csq1 = _csquare(0.0, 0.0, 1.0, 1.0)
csq2 = _csquare(1.0, 0.0, 2.0, 1.0)
csq3 = _csquare(10.0, 10.0, 11.0, 11.0)
csquares = [csq1, csq2, csq3]

@testset "constructor variants" begin
    @testset "from a vector of geometries" begin
        gl = GeometryLookup(csquares)
        @test all(splat(GO.equals), zip(val(gl), csquares))
        @test gl.tree isa SortTileRecursiveTree.STRtree
        @test gl.manifold == GO.Planar()
    end

    @testset "from a table with geometrycolumn" begin
        tbl = (; geom = csquares, value = 1:3)
        gl = GeometryLookup(tbl; geometrycolumn = :geom)
        @test all(splat(GO.equals), zip(val(gl), csquares))
    end

    @testset "Union{Missing} eltype is narrowed" begin
        geoms = Union{Missing, eltype(csquares)}[csquares...]
        gl = GeometryLookup(geoms)
        @test !(Missing <: eltype(val(gl)))
        @test length(gl) == 3
    end

    @testset "error paths" begin
        @test_throws ArgumentError GeometryLookup([1, 2, 3])
        @test_throws ArgumentError GeometryLookup(csquares, (X(),))
        @test_throws ArgumentError GeometryLookup(csquares, (X(), Y(), Ti()))
        @test_throws ArgumentError GeometryLookup(csquares; tree = 42)
    end
end

@testset "tree keyword" begin
    @testset "tree = nothing disables the accelerator" begin
        gl = GeometryLookup(csquares; tree = nothing)
        @test gl.tree === nothing
    end

    @testset "tree as a type" begin
        gl = GeometryLookup(csquares; tree = SortTileRecursiveTree.STRtree)
        @test gl.tree isa SortTileRecursiveTree.STRtree
    end

    @testset "tree as a prebuilt instance" begin
        tree = SortTileRecursiveTree.STRtree(csquares)
        gl = GeometryLookup(csquares; tree)
        @test gl.tree === tree
    end
end

@testset "crs" begin
    @testset "no crs anywhere gives nothing" begin
        gl = GeometryLookup(csquares)
        @test GI.crs(gl) === nothing
        @test crs(gl) === nothing
    end

    @testset "explicit crs keyword" begin
        gl = GeometryLookup(csquares; crs = EPSG(4326))
        @test crs(gl) == EPSG(4326)
    end

    @testset "setcrs" begin
        gl = GeometryLookup(csquares; crs = EPSG(4326))
        gl2 = Rasters.setcrs(gl, EPSG(3857))
        @test crs(gl2) == EPSG(3857)
        # only the crs changed - data and tree are reused
        @test val(gl2) == val(gl)
        @test gl2.tree === gl.tree
    end
end

@testset "DimensionalData interface" begin
    gl = GeometryLookup(csquares)

    @testset "dims, order, parent" begin
        @test DD.name.(DD.dims(gl)) == (:X, :Y)
        @test DD.order(gl) == Lookups.Unordered()
        @test parent(gl) == val(gl)
        @test length(gl) == 3
    end

    @testset "Geometry dimension" begin
        @test Geometry <: DD.Dimension
        @test DD.name(Geometry) == :Geometry
        dv = rand(Geometry(gl))
        @test DD.dims(dv, Geometry) isa Geometry
        @test val(DD.lookup(dv, Geometry)) == val(gl)
    end

    @testset "rebuild" begin
        # same data reuses the tree
        rb_same = DD.rebuild(gl; data = gl.data)
        @test rb_same.tree === gl.tree

        # new data rebuilds the tree
        rb_new = DD.rebuild(gl; data = csquares[1:2])
        @test length(rb_new) == 2
        @test rb_new.tree isa SortTileRecursiveTree.STRtree
        @test rb_new.tree !== gl.tree

        # empty data has no tree
        rb_empty = DD.rebuild(gl; data = empty(csquares))
        @test isempty(rb_empty)
        @test rb_empty.tree === nothing

        # an explicit tree type is honored
        rb_tree = DD.rebuild(gl; data = csquares[1:2], tree = SortTileRecursiveTree.STRtree)
        @test rb_tree.tree isa SortTileRecursiveTree.STRtree
    end

    @testset "show" begin
        dv = rand(Geometry(gl))
        str = sprint(show, MIME"text/plain"(), dv)
        @test occursin("Geometry", str)
        @test occursin("GeometryLookup", sprint(show, MIME"text/plain"(), gl))
    end
end
