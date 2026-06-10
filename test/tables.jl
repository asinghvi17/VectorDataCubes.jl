using Test
using VectorDataCubes
using Rasters, DimensionalData
using Rasters.Lookups
import DimensionalData as DD
import GeometryOps as GO, GeoInterface as GI
import Tables, DataAPI

# A small axis-aligned unit square at (x, y), like the idiom in test/selectors.jl.
square(x, y; s=1.0) =
    GI.Polygon([GI.LinearRing([(x, y), (x + s, y), (x + s, y + s), (x, y + s), (x, y)])])

@testset "Tables.jl integration" begin
    geoms = [square(0, 0), square(2, 0), square(0, 2), square(2, 2)]
    gl = GeometryLookup(geoms; crs=EPSG(4326))

    @testset "1-D cube (Geometry only)" begin
        A = rand(Geometry(gl))
        for tbl in (DD.DimTable(A), VectorDataCubeTable(A))
            ct = Tables.columntable(tbl)
            @test :Geometry in keys(ct)
            @test :value in keys(ct)
            @test length(ct.Geometry) == length(geoms)
            # The Geometry column holds the *actual* geometry objects.
            @test all(GI.isgeometry, ct.Geometry)
            @test all(splat(GO.equals), zip(ct.Geometry, geoms))

            rt = Tables.rowtable(tbl)
            @test length(rt) == length(geoms)
            @test GO.equals(rt[1].Geometry, geoms[1])
        end
    end

    @testset "2-D cube (Geometry × Ti)" begin
        A = rand(Geometry(gl), Ti(1:3))
        for tbl in (DD.DimTable(A), VectorDataCubeTable(A))
            ct = Tables.columntable(tbl)
            @test Set(keys(ct)) == Set((:Geometry, :Ti, :value))
            @test length(ct.Geometry) == length(geoms) * 3
            @test all(GI.isgeometry, ct.Geometry)
            @test Set(ct.Ti) == Set(1:3)
            # Every (geometry, time) combination appears exactly once.
            @test length(unique(zip(map(GO.centroid, ct.Geometry), ct.Ti))) == length(geoms) * 3
        end
    end

    @testset "DimStack cube (multiple value columns)" begin
        st = DimStack((a=rand(Geometry(gl)), b=rand(Geometry(gl))))
        for tbl in (DD.DimTable(st), VectorDataCubeTable(st))
            ct = Tables.columntable(tbl)
            @test Set(keys(ct)) == Set((:Geometry, :a, :b))
            @test length(ct.Geometry) == length(geoms)
            @test all(GI.isgeometry, ct.Geometry)
            @test eltype(ct.a) == Float64
            @test eltype(ct.b) == Float64
        end
    end

    @testset "VectorDataCubeTable interface & metadata" begin
        A = rand(Geometry(gl), Ti(1:2))
        tbl = VectorDataCubeTable(A)
        @test vectordatacubetable(A) isa VectorDataCubeTable

        @test Tables.istable(typeof(tbl))
        @test Tables.columnaccess(typeof(tbl))
        @test Set(Tables.columnnames(tbl)) == Set((:Geometry, :Ti, :value))
        @test Tables.getcolumn(tbl, :Ti) == Tables.getcolumn(DD.DimTable(A), :Ti)

        # Schema matches the inner DimTable.
        sch = Tables.schema(tbl)
        @test :Geometry in sch.names

        # CRS exposed via DataAPI.metadata and GI.crs.
        @test GI.crs(tbl) == EPSG(4326)
        @test DataAPI.metadatasupport(typeof(tbl)) == (read=true, write=false)
        @test "crs" in DataAPI.metadatakeys(tbl)
        @test DataAPI.metadata(tbl, "crs") == EPSG(4326)
        @test DataAPI.metadata(tbl, "missingkey", :fallback) == :fallback
        v, style = DataAPI.metadata(tbl, "crs"; style=true)
        @test v == EPSG(4326)
        @test style == :default
    end

    @testset "no-crs cube exposes no crs metadata" begin
        gl2 = GeometryLookup(geoms)  # hand-made polygons have no crs
        A = rand(Geometry(gl2))
        tbl = VectorDataCubeTable(A)
        @test isnothing(GI.crs(tbl))
        @test DataAPI.metadatakeys(tbl) == ()
        # Round-trips through Tables regardless of crs.
        ct = Tables.columntable(tbl)
        @test all(GI.isgeometry, ct.Geometry)
    end

    @testset "errors on non-vector cube" begin
        ras = rand(X(1:3), Y(1:3))
        @test_throws ArgumentError VectorDataCubeTable(ras)
    end
end
