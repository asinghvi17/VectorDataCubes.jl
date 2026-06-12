#=
# Building a vector data cube: North Carolina SIDS

This is a port of the classic North Carolina sudden-infant-death-syndrome (SIDS)
example used to introduce vector data cubes in both the R `stars` package
(https://r-spatial.org/r/2022/09/12/vdc.html) and Python's `xvec`
(https://xvec.readthedocs.io/en/stable/intro.html).

The dataset is a set of 100 county polygons with birth and SIDS counts for two
periods (1974 and 1979). Instead of keeping these as columns of a flat table,
we arrange them as a *vector data cube*: a `DimStack` whose first dimension is
indexed by the county geometries themselves (via [`GeometryLookup`](@ref)), and
whose second dimension is the year.

The data ships with R's `sf` package; we download the shapefile straight from
its GitHub repository (~150 kB).
=#

using VectorDataCubes
using Rasters, DimensionalData
using Rasters.Lookups
import DimensionalData as DD
import GeometryOps as GO, GeoInterface as GI
import Shapefile
using DataFrames
import DataAPI, Tables
using Downloads: download

datadir = joinpath(@__DIR__, "data")
mkpath(datadir)
for ext in ("shp", "dbf", "shx", "prj")
    file = joinpath(datadir, "nc.$ext")
    isfile(file) || download(
        "https://raw.githubusercontent.com/r-spatial/sf/main/inst/shape/nc.$ext", file
    )
end

counties = Shapefile.Table(joinpath(datadir, "nc.shp"))

#=
## Constructing the cube

The geometry dimension is a `Geometry` dimension wrapping a `GeometryLookup` of
the county polygons. The lookup builds a spatial tree over the geometries, so
spatial selectors on the cube are fast, and it carries the CRS (the data is in
NAD27, EPSG:4267).
=#

geoms = GO.tuples(Shapefile.shapes(counties))
gl = GeometryLookup(geoms; crs=EPSG(4267))

years = Dim{:year}([1974, 1979])
cube = DimStack((;
    births = DimArray([counties.BIR74 counties.BIR79], (Geometry(gl), years)),
    sids = DimArray([counties.SID74 counties.SID79], (Geometry(gl), years)),
    nonwhite_births = DimArray([counties.NWBIR74 counties.NWBIR79], (Geometry(gl), years)),
))

#=
## Spatial selectors

The cube is indexed by real geometries, so we can ask spatial questions
directly. Which county contains Raleigh?
=#

raleigh = (-78.6382, 35.7796)
wake = cube[Geometry(Contains(raleigh))]
@assert size(wake[:births]) == (1, 2)

# The selector machinery also works on the lookup itself, which lets us recover
# the county's attributes from the original table:
wake_idx = only(Lookups.selectindices(gl, Contains(raleigh)))
println("Raleigh is in $(counties.NAME[wake_idx]) county")
println("births: ", parent(wake[:births]))

#=
## Derived layers

Cube arithmetic works as usual — dimensions (including the geometry lookup) are
carried through broadcasting. Here is the SIDS rate per 1000 births:
=#

rate = cube[:sids] ./ cube[:births] .* 1000

worst_rate, worst_idx = findmax(rate[year=At(1974)])
println("highest 1974 SIDS rate: $(counties.NAME[worst_idx]) ",
    "($(round(worst_rate; digits=2)) per 1000 births)")

#=
## To a table

`VectorDataCubeTable` flattens the cube into a Tables.jl table with one row per
(county, year) pair. The `:Geometry` column holds the actual polygons, and the
CRS is exposed as DataAPI metadata, which `DataFrame` picks up on construction.
=#

df = DataFrame(VectorDataCubeTable(rate))
@assert nrow(df) == 100 * 2
@assert DataFrames.metadata(df, "crs") == EPSG(4267)
println(first(sort(df, :value; rev=true), 5))
