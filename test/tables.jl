using Test
using VectorDataCubes
using Rasters, DimensionalData
using Rasters.Lookups
import DimensionalData as DD
import GeometryOps as GO, GeoInterface as GI
import Tables

# A small axis-aligned unit square at (x, y), like the idiom in test/selectors.jl.
square(x, y; s=1.0) =
    GI.Polygon([GI.LinearRing([(x, y), (x + s, y), (x + s, y + s), (x, y + s), (x, y)])])

@testset "Tables.jl integration" begin
    geoms = [square(0, 0), square(2, 0), square(0, 2), square(2, 2)]
    gl = GeometryLookup(geoms; crs=EPSG(4326))

    @testset "1-D cube (Geometry only)" begin
        A = rand(Geometry(gl))
        for tbl in (DD.DimTable(A), vectordatacubetable(A))
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
        for tbl in (DD.DimTable(A), vectordatacubetable(A))
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
        for tbl in (DD.DimTable(st), vectordatacubetable(st))
            ct = Tables.columntable(tbl)
            @test Set(keys(ct)) == Set((:Geometry, :a, :b))
            @test length(ct.Geometry) == length(geoms)
            @test all(GI.isgeometry, ct.Geometry)
            @test eltype(ct.a) == Float64
            @test eltype(ct.b) == Float64
        end
    end

    @testset "crs is recorded in the table's metadata" begin
        A = rand(Geometry(gl), Ti(1:2))
        tbl = vectordatacubetable(A)
        @test tbl isa DD.DimTable
        @test DD.metadata(parent(tbl))[:crs] == EPSG(4326)
        # The rebuild is cheap: data and lookups are reused, not copied.
        @test parent(parent(tbl)) === parent(A)
        @test DD.lookup(parent(tbl), Geometry) === gl
        # Existing metadata is preserved alongside the crs.
        B = DD.rebuild(A; metadata=Dict{Symbol,Any}(:title => "t"))
        mdB = DD.metadata(parent(vectordatacubetable(B)))
        @test mdB[:crs] == EPSG(4326)
        @test mdB[:title] == "t"
    end

    @testset "no-crs cube gets no crs metadata" begin
        gl2 = GeometryLookup(geoms)  # hand-made polygons have no crs
        A = rand(Geometry(gl2))
        tbl = vectordatacubetable(A)
        @test DD.metadata(parent(tbl)) isa DD.Lookups.NoMetadata
        # Round-trips through Tables regardless of crs.
        ct = Tables.columntable(tbl)
        @test all(GI.isgeometry, ct.Geometry)
    end

    @testset "errors on non-vector cube" begin
        ras = rand(X(1:3), Y(1:3))
        @test_throws ArgumentError vectordatacubetable(ras)
    end
end

@testset "vectordatacube: table -> cube" begin
    geoms = [square(0, 0), square(2, 0), square(0, 2), square(2, 2)]
    tbl = (geometry=geoms, name=["a", "b", "c", "d"], pop=[10, 20, 30, 40])

    @testset "attribute columns become layers over Geometry" begin
        cube = vectordatacube(tbl; crs=EPSG(4326))
        @test cube isa DD.DimStack
        @test keys(cube) == (:name, :pop)
        @test DD.lookup(cube, Geometry) isa GeometryLookup
        @test crs(DD.lookup(cube, Geometry)) == EPSG(4326)
        @test all(splat(GO.equals), zip(val(DD.lookup(cube, Geometry)), geoms))
        @test parent(cube[:name]) == tbl.name
        @test parent(cube[:pop]) == tbl.pop
    end

    @testset "attributes stay aligned under subsetting" begin
        cube = vectordatacube(tbl)
        sub = cube[Geometry=2:3]
        @test parent(sub[:name]) == ["b", "c"]
        @test length(DD.lookup(sub, Geometry)) == 2
        # spatial selectors too
        hit = cube[Geometry(Contains((2.5, 0.5)))]
        @test parent(hit[:name]) == ["b"]
    end

    @testset "layers and geometrycolumn keywords" begin
        cube = vectordatacube(tbl; layers=(:pop,))
        @test keys(cube) == (:pop,)
        renamed = (geom=geoms, name=tbl.name)
        @test_throws ArgumentError vectordatacube(renamed)
        cube2 = vectordatacube(renamed; geometrycolumn=:geom)
        @test parent(cube2[:name]) == tbl.name
    end

    @testset "round trip through vectordatacubetable" begin
        cube = vectordatacube(tbl; crs=EPSG(4326))
        ct = Tables.columntable(vectordatacubetable(cube))
        @test Set(keys(ct)) == Set((:Geometry, :name, :pop))
        @test all(splat(GO.equals), zip(ct.Geometry, geoms))
        @test ct.name == tbl.name
        @test ct.pop == tbl.pop
    end

    @testset "errors" begin
        @test_throws ArgumentError vectordatacube(geoms)             # not a table
        @test_throws ArgumentError vectordatacube((geometry=geoms,)) # no attribute columns
    end
end
