#=
# Point extraction: city temperature time series

This is a port of xvec's "Extracting points from a geospatial raster" tutorial
(https://xvec.readthedocs.io/en/stable/extract_pts.html), and of the
station × time cube archetype from the R `stars` / "Spatial Data Science" NO₂
example: sample a raster time series at point geometries and arrange the
result as a vector data cube over `(Ti, Geometry)`, where the geometry
dimension holds the points.

Data is the same NCEP air temperature NetCDF as in `02_zonal_countries.jl`,
sampled at Natural Earth populated places (1:110m).
=#

using VectorDataCubes
using Rasters, DimensionalData
using Rasters.Lookups
import DimensionalData as DD
import GeometryOps as GO, GeoInterface as GI
import NCDatasets
using NaturalEarth
using DataFrames
using Statistics: mean
using Downloads: download

datadir = joinpath(@__DIR__, "data")
mkpath(datadir)
airfile = joinpath(datadir, "air_temperature.nc")
isfile(airfile) || download(
    "https://raw.githubusercontent.com/pydata/xarray-data/master/air_temperature.nc",
    airfile,
)

air = RasterStack(airfile)[:air]
air = set(air, X => val(dims(air, X)) .- 360)
air = air .- 273.15

#=
## The points

Natural Earth's populated places, filtered to the raster's spatial coverage.
=#

places = DataFrame(naturalearth("populated_places_simple", 110))
xb, yb = bounds(air, X), bounds(air, Y)
incoverage(g) = xb[1] <= GI.x(g) <= xb[2] && yb[1] <= GI.y(g) <= yb[2]
cities = subset(places, :geometry => ByRow(incoverage))
println(nrow(cities), " cities in coverage: ", join(cities.name, ", "))

#=
## Sampling -> vector data cube

We sample the nearest grid cell for each city, then assemble the series into
a cube whose geometry dimension is a `GeometryLookup` over the city *points*.
(The lookup works for any geometry type, not just polygons.)
=#

series = map(cities.geometry) do pt
    air[X=Near(GI.x(pt)), Y=Near(GI.y(pt))]
end
citygl = GeometryLookup(GO.tuples(cities.geometry); crs=EPSG(4326))
citycube = Raster(
    reduce(hcat, parent.(series)),
    (dims(air, Ti), Geometry(citygl));
    name=:air,
)
@assert size(citycube) == (length(dims(air, Ti)), nrow(cities))

#=
## Selectors on a point lookup

`Near` finds the city closest to a query point — here, the closest city in
the cube to lower Manhattan:
=#

nyc_series = citycube[Geometry(Near((-74.0, 40.7)))]
println("nearest city series has ", length(nyc_series), " time steps")

# Mean annual temperature per city, warmest first:
annual = vec(parent(mean(citycube; dims=Ti)))
for i in sortperm(annual; rev=true)
    println(rpad(cities.name[i], 16), round(annual[i]; digits=1), " °C")
end

#=
## To a table

One row per (time, city), with real point geometries in the `:Geometry`
column — ready for DataFrames or any other Tables.jl consumer.
=#

df = DataFrame(vectordatacubetable(citycube))
@assert nrow(df) == length(citycube)
println(first(df, 3))
