# Tables.jl integration: `DimensionalData.DimTable` already materializes a
# vector data cube as one row per geometry × other-dim coordinate, with real
# geometry objects in the `Geometry` column. [`vectordatacubetable`](@ref)
# additionally records the lookup's crs in the cube's metadata.

import Tables

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
