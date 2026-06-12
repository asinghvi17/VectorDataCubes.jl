# VectorDataCubes.jl examples

Ports of well-known vector data cube tutorials from the Python (`xvec`) and R
(`stars`) ecosystems. Each script is self-contained: it downloads its (small,
open) datasets into `examples/data/` on first run.

Run them with this directory as the project:

```sh
julia --project=examples examples/01_intro_nc_sids.jl
julia --project=examples examples/02_zonal_countries.jl
julia --project=examples examples/03_extract_points.jl
julia --project=examples examples/04_taxi_od_cube.jl
```

| Script | Ports | Data |
| --- | --- | --- |
| `01_intro_nc_sids.jl` | The North Carolina SIDS introduction from [R stars / "Vector Data Cubes"](https://r-spatial.org/r/2022/09/12/vdc.html) and the [xvec intro](https://xvec.readthedocs.io/en/stable/intro.html): build a cube directly over a geometry dimension (county × year), spatial selectors, derived layers, table conversion. | NC counties shapefile shipped with R's `sf` (~150 kB) |
| `02_zonal_countries.jl` | The [xvec "Zonal statistics" tutorial](https://xvec.readthedocs.io/en/stable/zonal_stats.html) and stars' [`aggregate(raster, by = polygons)`](https://r-spatial.github.io/stars/): aggregate a gridded temperature time series over country polygons into a `(Ti, Geometry)` cube, handle partial/empty coverage, reproject, convert to a table. | xarray's `air_temperature.nc` tutorial NetCDF (~7 MB) + Natural Earth 1:110m countries |
| `03_extract_points.jl` | The [xvec "Extracting points" tutorial](https://xvec.readthedocs.io/en/stable/extract_pts.html) and the station × time cube archetype from the [Spatial Data Science book](https://r-spatial.org/book/13-Geostatistics.html): sample a raster time series at point geometries into a `(Ti, Geometry)` cube over points. | same NetCDF + Natural Earth populated places |
| `04_taxi_od_cube.jl` | The [xvec "Indexing" tutorial](https://xvec.readthedocs.io/en/stable/indexing.html): a cube with *two* geometry-indexed dimensions — taxi trip counts over (origin zone, destination zone) — sliced with spatial selectors on either axis. | NYC TLC yellow-taxi trips, Jan 2022 (~37 MB Parquet) + taxi-zone shapefile |
