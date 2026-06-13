using Test

using Rasters, DimensionalData
using Rasters.Lookups
# Explicitly bind VectorDataCubes' own `zonal` (Rasters exports one too).
using VectorDataCubes: zonal
import DimensionalData as DD
import GeometryOps as GO, GeoInterface as GI
using Statistics: mean
using Dates

_zsquare(x1, y1, x2, y2) =
    GI.Polygon([GI.LinearRing([(x1, y1), (x2, y1), (x2, y2), (x1, y2), (x1, y1)])])

# A 10x10 raster on [0, 10]^2 with cell centers at 0.5, 1.5, ..., 9.5,
# where each cell's value is its x index. Zone values are exactly known:
# zoneA covers x cells 1:2, zoneB covers x cells 5:7, zoneC is off-raster.
zx = X(0.5:1.0:9.5)
zy = Y(0.5:1.0:9.5)
zti = Ti(DateTime(2020, 1, 1):Month(1):DateTime(2020, 3, 1))
ras2d = Raster([Float64(xi) for xi in 1:10, yi in 1:10], (zx, zy); name=:vals)
# 3D: value = x index * time index
ras3d = Raster([Float64(xi * t) for xi in 1:10, yi in 1:10, t in 1:3], (zx, zy, zti); name=:cube)

zoneA = _zsquare(0.0, 0.0, 2.0, 2.0)    # 2x2 cells, x values {1, 2}
zoneB = _zsquare(4.0, 2.0, 7.0, 5.0)    # 3x3 cells, x values {5, 6, 7}
zoneC = _zsquare(20.0, 20.0, 22.0, 22.0) # entirely outside the raster
zgl = GeometryLookup([zoneA, zoneB, zoneC])

@testset "zonal with a GeometryLookup" begin
    @testset "2D raster -> vector over Geometry" begin
        res = zonal(sum, ras2d; of=zgl, progress=false)
        @test res isa Raster
        @test DD.dims(res, Geometry) isa Geometry
        @test val(DD.lookup(res, Geometry)) == zgl.data
        @test res[1] == (1 + 2) * 2        # {1,2} over 2 y-cells
        @test res[2] == (5 + 6 + 7) * 3    # {5,6,7} over 3 y-cells
        @test ismissing(res[3])            # off-raster geometry
        # the result is a real vector data cube: spatial selectors work on it
        @test res[Geometry(Contains((1.0, 1.0)))] == res[[1]]
    end

    @testset "of = Geometry(lookup) behaves the same" begin
        res_lookup = zonal(sum, ras2d; of=zgl, progress=false)
        res_dim = zonal(sum, ras2d; of=Geometry(zgl), progress=false)
        @test isequal(res_lookup, res_dim)
    end

    @testset "non-lookup `of` forwards to Rasters.zonal" begin
        res = zonal(sum, ras2d; of=[zoneA, zoneB], progress=false)
        @test res isa Vector
        @test res == [(1 + 2) * 2, (5 + 6 + 7) * 3]
    end

    @testset "3D raster -> cube over (Ti, Geometry)" begin
        res = zonal(mean, ras3d; of=zgl, progress=false)
        @test size(res) == (3, 3)
        @test DD.dims(res, Ti) == DD.dims(ras3d, Ti)
        @test DD.dims(res, Geometry) isa Geometry
        @test res[Ti=1, Geometry=1] ≈ 1.5    # mean({1,2} * 1)
        @test res[Ti=2, Geometry=1] ≈ 3.0    # mean({1,2} * 2)
        @test res[Ti=3, Geometry=2] ≈ 18.0   # mean({5,6,7} * 3)
        @test all(ismissing, res[Geometry=3])
        @test !any(ismissing, res[Geometry=1:2])
    end

    @testset "spatialslices = false reduces over all dims" begin
        res = zonal(mean, ras3d; of=zgl, spatialslices=false, progress=false)
        @test size(res) == (3,)
        @test res[1] ≈ mean([xi * t for xi in 1:2, _ in 1:2, t in 1:3]) # 3.0
        @test ismissing(res[3])
    end

    @testset "RasterStack -> stack of cubes" begin
        st = RasterStack((flat=ras2d, cube=ras3d))
        res = zonal(mean, st; of=zgl, progress=false)
        @test res isa RasterStack
        @test size(res[:flat]) == (3,)
        @test size(res[:cube]) == (3, 3)
        @test res[:flat][1] ≈ 1.5
        @test res[:cube][Ti=2, Geometry=2] ≈ 12.0  # mean({5,6,7} * 2)
        @test ismissing(res[:flat][3])
        @test all(ismissing, res[:cube][Geometry=3])
    end

    @testset "emptyval fills empty slices" begin
        # data missing at t = 2 in zoneD's cells: that slice is empty under
        # skipmissing, so it gets emptyval while other slices are computed
        data = Array{Union{Missing,Float64}}([Float64(xi * t) for xi in 1:10, yi in 1:10, t in 1:3])
        data[9:10, 9:10, 2] .= missing
        rasm = Raster(data, (zx, zy, zti); name=:gappy)
        zoneD = _zsquare(8.0, 8.0, 10.0, 10.0) # x cells 9:10, y cells 9:10
        res = zonal(mean, rasm; of=GeometryLookup([zoneD]), emptyval=NaN, progress=false)
        @test isnan(res[Ti=2, Geometry=1])
        @test res[Ti=1, Geometry=1] ≈ 9.5
        @test res[Ti=3, Geometry=1] ≈ 28.5
    end

    @testset "emptyval=missing with sub-cell and off-raster geometries" begin
        # zoneE sits between cell centers: crop is non-empty but the mask
        # removes every cell, so each slice is empty -> emptyval. Its slice
        # results have a different eltype than zoneA's, which must not break
        # assembling the cube.
        zoneE = _zsquare(0.6, 0.6, 0.9, 0.9)
        gl = GeometryLookup([zoneA, zoneE, zoneC])
        res = zonal(mean, ras3d; of=gl, emptyval=missing, progress=false)
        @test size(res) == (3, 3)
        @test !any(ismissing, res[Geometry=1])
        @test res[Ti=1, Geometry=1] ≈ 1.5
        @test all(ismissing, res[Geometry=2])  # sub-cell -> emptyval per slice
        @test all(ismissing, res[Geometry=3])  # off-raster -> missing
    end

    @testset "all geometries off-raster keeps the cube shape" begin
        zoneC2 = _zsquare(30.0, 30.0, 32.0, 32.0)
        off_gl = GeometryLookup([zoneC, zoneC2])
        res = zonal(mean, ras3d; of=off_gl, progress=false)
        @test size(res) == (3, 2)  # (Ti, Geometry), not collapsed to (2,)
        @test DD.dims(res, Ti) == DD.dims(ras3d, Ti)
        @test all(ismissing, res)
        res2d = zonal(mean, ras2d; of=off_gl, progress=false)
        @test size(res2d) == (2,)
        @test all(ismissing, res2d)
    end

    @testset "empty lookup throws" begin
        empty_gl = DD.rebuild(zgl; data=empty(zgl.data))
        @test_throws ArgumentError zonal(sum, ras2d; of=empty_gl, progress=false)
    end
end

@testset "reproject" begin
    gl = GeometryLookup([zoneA, zoneB]; crs=EPSG(4326))
    gl3857 = reproject(EPSG(3857), gl)
    @test crs(gl3857) == EPSG(3857)
    # Web Mercator coordinates are in meters, so points are far from the origin
    @test GI.x(GI.getpoint(GI.getexterior(gl3857[1]), 2)) ≈ 2.0 * 20037508.34 / 180 rtol = 1e-3
    # round trip back to EPSG:4326 recovers the original coordinates (up to float error)
    glback = reproject(EPSG(4326), gl3857)
    @test all(zip(val(glback), val(gl))) do (a, b)
        all(zip(GI.getpoint(a), GI.getpoint(b))) do (p, q)
            isapprox(GI.x(p), GI.x(q); atol=1e-6) && isapprox(GI.y(p), GI.y(q); atol=1e-6)
        end
    end
    # reprojecting a whole vector data cube reprojects the lookup
    dv = rand(Geometry(gl))
    dv3857 = reproject(EPSG(3857), dv)
    @test crs(DD.lookup(dv3857, Geometry)) == EPSG(3857)
    # no crs is an error
    @test_throws ArgumentError reproject(EPSG(3857), GeometryLookup([zoneA]))
end
