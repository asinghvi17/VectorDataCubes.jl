# Tables.jl integration in both directions: [`vectordatacube`](@ref) lifts a
# flat table (or feature collection) into a `DimStack` over a `Geometry`
# dimension, one layer per attribute column; [`vectordatacubetable`](@ref)
# flattens a cube back to a `DimTable`, one row per geometry × other-dim
# coordinate, with real geometry objects in the `Geometry` column and the
# lookup's crs recorded in the cube's metadata.

import Tables

"""
    vectordatacube(table; geometrycolumn=nothing, layers=nothing, crs=nokw)

Convert a table with a geometry column (a GeoJSON `FeatureCollection`, a
`Shapefile.Table`, a `DataFrame`, ...) to a vector data cube: a `DimStack`
over a `Geometry` dimension carrying a [`GeometryLookup`](@ref) of the
geometries, with one layer per remaining column.

Because the attributes are layers over the same `Geometry` dimension,
subsetting the cube (by index or spatial selector) keeps them aligned with the
geometries — there is no separate attribute table to keep in sync.

# Keywords

- `geometrycolumn`: the name of the geometry column. Defaults to the table's
  own metadata (`GeoInterface.geometrycolumns`), which is `:geometry` for most
  formats. Other geometry-typed columns are kept as ordinary layers.
- `layers`: an iterable of column names to keep as layers. Defaults to every
  column except the geometry column.
- `crs`: the coordinate reference system of the geometries. Defaults to the
  crs of the table or its geometries, if they carry one.

[`vectordatacubetable`](@ref) is the inverse: it flattens the cube back into a
table, with the geometries in a `:Geometry` column.

# Example

```julia
using VectorDataCubes, NaturalEarth
countries = vectordatacube(naturalearth("admin_0_countries", 110))
countries[:NAME]                          # a DimArray over Geometry
countries[Geometry(Contains((9.0, 50.0)))][:NAME]  # the country containing a point
```
"""
function vectordatacube(table; geometrycolumn=nothing, layers=nothing, crs=nokw)
    Tables.istable(table) || throw(ArgumentError("""
    `vectordatacube` requires a Tables.jl-compatible table with a geometry column,
    but `Tables.istable` is false for the input ($(typeof(table))).
    To build a cube from a plain geometry vector, use `Geometry(GeometryLookup(geoms))` directly.
    """))
    cols = Tables.columns(table)
    colnames = Tables.columnnames(cols)
    # The same geometry column detection as `GeometryOpsCore.get_geometries`.
    geomcol = if isnothing(geometrycolumn)
        geomcols = GI.geometrycolumns(table)
        isnothing(geomcols) || isempty(geomcols) ? :geometry : first(geomcols)
    else
        geometrycolumn
    end
    geomcol in colnames || throw(ArgumentError("""
    No geometry column :$geomcol found in the table (columns: $(join(colnames, ", "))).
    Pass the right column name with the `geometrycolumn` keyword.
    """))
    if isnokw(crs)
        table_crs = GI.crs(table)
        isnothing(table_crs) || (crs = table_crs)
    end
    gl = GeometryLookup(collect(Tables.getcolumn(cols, geomcol)); crs)
    layernames = Tuple(isnothing(layers) ? (n for n in colnames if n != geomcol) : layers)
    isempty(layernames) && throw(ArgumentError("""
    The table has no columns other than the geometry column :$geomcol, so there are
    no layers to build. Use `Geometry(GeometryLookup(geoms))` directly instead.
    """))
    gdim = Geometry(gl)
    return DD.DimStack(NamedTuple{layernames}(map(layernames) do name
        DD.DimArray(collect(Tables.getcolumn(cols, name)), gdim; name)
    end))
end

"""
    vectordatacubetable(cube)

Convert a vector data cube (a `DimArray`/`Raster`/`DimStack`/`RasterStack`
with a `Geometry([`GeometryLookup`](@ref))` dimension) to a
`DimensionalData.DimTable`, with one row per geometry × other-dimension
coordinate: a `Geometry` column holding the actual geometry objects, one
column per other dimension (e.g. `Ti`), and one value column per layer.

The geometry crs (if any) is recorded under `:crs` in the metadata of the
table's parent cube, retrievable as `DimensionalData.metadata(parent(tbl))[:crs]`.
"""
function vectordatacubetable(cube::Union{DD.AbstractDimArray,DD.AbstractDimStack})
    DD.hasdim(cube, Geometry) || throw(ArgumentError("""
    `vectordatacubetable` requires a vector data cube with a `Geometry` dimension,
    but the input has dimensions $(DD.basedims(cube)).
    Wrap your geometries in a `Geometry(GeometryLookup(geoms))` axis first.
    """))
    lookup = DD.lookup(cube, Geometry)
    crs = lookup isa GeometryLookup ? GI.crs(lookup) : nothing
    isnothing(crs) && return DD.DimTable(cube)
    # A cheap rebuild: the same data and lookups, plus metadata carrying the crs.
    md = DD.metadata(cube)
    newmd = md isa DD.Lookups.NoMetadata ? Dict{Symbol,Any}() :
        Dict{Symbol,Any}(k => md[k] for k in keys(md))
    newmd[:crs] = crs
    return DD.DimTable(DD.rebuild(cube; metadata=newmd))
end
