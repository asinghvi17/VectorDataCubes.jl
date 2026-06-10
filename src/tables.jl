#=
# Tables.jl integration for vector data cubes

A "vector data cube" is a `DimArray`/`Raster`/`DimStack`/`RasterStack` that has a
`Geometry(GeometryLookup(geoms))` axis (possibly alongside other dimensions such
as `Ti` or `Band`).

`DimensionalData`'s `DimTable` already does the right thing for these cubes: it
materialises one row per geometry × other-dim-coordinate combination, and the
`Geometry` column holds the *actual geometry objects* (because the `GeometryLookup`
parent is the vector of geometries, which `DimTable` iterates as the lookup values).

The only thing `DimTable` does *not* do is surface the CRS of the geometry column.
We therefore define a thin, package-owned wrapper, [`VectorDataCubeTable`](@ref),
that forwards the entire Tables.jl interface to an inner `DimTable` and additionally
implements `DataAPI.metadata` to expose the CRS.

All methods here dispatch on types this package owns (`VectorDataCubeTable`,
`GeometryLookup`, `Geometry`), so there is no type piracy.
=#

import Tables
import DataAPI

"""
    VectorDataCubeTable(cube)
    VectorDataCubeTable(dimtable)

A Tables.jl-compatible table view of a vector data cube.

A vector data cube is a `DimArray`/`Raster`/`DimStack`/`RasterStack` carrying a
`Geometry(GeometryLookup(...))` dimension. Converting it to a table yields one
row per geometry × other-dimension-coordinate combination, with:

- a `Geometry` column containing the actual geometry objects,
- one column per other dimension (e.g. `Ti`, `Band`),
- one value column per layer (named `:value` for a single array, or the layer
  names for a `DimStack`/`RasterStack`).

This wraps a `DimensionalData.DimTable` and forwards the full Tables.jl interface
to it, additionally exposing the geometry-column CRS through `DataAPI.metadata`
(key `"crs"`).

# Example

```julia
using VectorDataCubes, Rasters, Tables
cube = rand(Geometry(GeometryLookup(geoms)), Ti(1:3))
tbl = VectorDataCubeTable(cube)        # or `vectordatacubetable(cube)`
Tables.columntable(tbl)                # NamedTuple of columns, incl. `:Geometry`
DataAPI.metadata(tbl, "crs")           # the CRS of the geometry column
```
"""
struct VectorDataCubeTable{T<:DD.DimTable,C}
    table::T
    crs::C
end

# Build from a DimTable directly: pull the crs out of the Geometry dimension's lookup.
function VectorDataCubeTable(table::DD.DimTable)
    VectorDataCubeTable(table, _table_crs(table))
end

# Build from a cube (DimArray/Raster/DimStack/RasterStack). Validate the Geometry dim.
function VectorDataCubeTable(cube::Union{DD.AbstractDimArray,DD.AbstractDimStack})
    DD.hasdim(cube, Geometry) || throw(ArgumentError("""
    `VectorDataCubeTable` requires a vector data cube with a `Geometry` dimension,
    but the input has dimensions $(DD.basedims(cube)).
    Wrap your geometries in a `Geometry(GeometryLookup(geoms))` axis first.
    """))
    VectorDataCubeTable(DD.DimTable(cube))
end

"""
    vectordatacubetable(cube)

Construct a [`VectorDataCubeTable`](@ref) from a vector data cube. Convenience
alias for the `VectorDataCubeTable` constructor.
"""
vectordatacubetable(cube) = VectorDataCubeTable(cube)

# Extract the crs from the Geometry dimension's GeometryLookup, if present.
function _table_crs(table::DD.DimTable)
    A = DD.parent(table) # the underlying DimStack/DimArray of the DimTable
    if DD.hasdim(A, Geometry)
        lookup = DD.val(DD.dims(A, Geometry))
        if lookup isa GeometryLookup
            return GI.crs(lookup)
        end
    end
    return nothing
end

# Allow grabbing the crs straight from a GeometryLookup-backed dimension too.
_table_crs(lookup::GeometryLookup) = GI.crs(lookup)

DD.dims(t::VectorDataCubeTable) = DD.dims(t.table)
GI.crs(t::VectorDataCubeTable) = t.crs

#=
## Tables.jl interface

We forward everything to the inner `DimTable`. `VectorDataCubeTable` is a
column-access table.
=#

Tables.istable(::Type{<:VectorDataCubeTable}) = true
Tables.columnaccess(::Type{<:VectorDataCubeTable}) = true
Tables.columns(t::VectorDataCubeTable) = Tables.columns(t.table)
Tables.schema(t::VectorDataCubeTable) = Tables.schema(t.table)
Tables.columnnames(t::VectorDataCubeTable) = Tables.columnnames(t.table)
Tables.getcolumn(t::VectorDataCubeTable, i::Int) = Tables.getcolumn(t.table, i)
Tables.getcolumn(t::VectorDataCubeTable, nm::Symbol) = Tables.getcolumn(t.table, nm)

# Row access falls back through Tables' generic machinery on `columns`, but we
# forward it explicitly for clarity / efficiency.
Tables.rows(t::VectorDataCubeTable) = Tables.rows(t.table)

#=
## DataAPI metadata: expose the geometry-column CRS

`DimTable` itself does not support `DataAPI.metadata`, so exposing the CRS is a
natural, non-pirating use of `DataAPI` here. We expose a single table-level key,
`"crs"`, with `:default` style.
=#

DataAPI.metadatasupport(::Type{<:VectorDataCubeTable}) = (read=true, write=false)

function DataAPI.metadatakeys(t::VectorDataCubeTable)
    isnothing(t.crs) ? () : ("crs",)
end

function DataAPI.metadata(t::VectorDataCubeTable, key::AbstractString; style::Bool=false)
    if key == "crs"
        return style ? (t.crs, :default) : t.crs
    end
    throw(KeyError(key))
end

function DataAPI.metadata(t::VectorDataCubeTable, key::AbstractString, default; style::Bool=false)
    if key == "crs" && !isnothing(t.crs)
        return style ? (t.crs, :default) : t.crs
    end
    return style ? (default, :default) : default
end
