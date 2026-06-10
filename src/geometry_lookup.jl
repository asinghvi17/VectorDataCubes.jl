import DimensionalData as DD
import GeometryOps as GO, GeometryOpsCore as GOCore
import GeoInterface as GI
import Rasters as RA
import Extents
import Missings
import SortTileRecursiveTree
import Proj # to load GeometryOps' Proj extension for `reproject`

using Rasters: isnokw, nokw, Lookups, val, Dimensions, dims

"""
    Geometry(lookup::GeometryLookup)
    Geometry(sel)

A dimension for geometry lookups in vector data cubes.

Construct a `DimArray` with a [`GeometryLookup`](@ref) like
`DimArray(data, Geometry(GeometryLookup(geometries)))`, and
index into it with selectors like `Geometry(Contains(point))`.
"""
DD.@dim Geometry "Geometry"

"""
    GeometryLookup(data, dims = (X(), Y()); geometrycolumn = nothing)

A lookup type for geometry dimensions in vector data cubes.

`GeometryLookup` provides efficient spatial indexing and lookup for 
geometries using an STRtree (Sort-Tile-Recursive tree). 

It is used as the lookup type for geometry dimensions in vector 
data cubes, enabling fast spatial queries and operations.

It spans the dimensions given to it in `dims`, as well as the dimension
 it's wrapped in - you would construct a DimArray with a GeometryLookup
like `DimArray(data, Geometry(GeometryLookup(data, dims)))`.
Here, `Geometry` is a dimension - but selectors in X and Y will also work!

# Examples

```julia
using Rasters

using NaturalEarth
import GeometryOps as GO

# construct the polygon lookup
polygons = NaturalEarth.naturalearth("admin_0_countries", 110).geometry
polygon_lookup = GeometryLookup(polygons, (X(), Y()))

# create a DimArray with the polygon lookup
dv = rand(Geometry(polygon_lookup))

# select the polygon with the centroid of the 88th polygon
only(dv[Geometry(Contains(GO.centroid(polygons[88])))]) == dv[88] # true
```
"""
struct GeometryLookup{T,A<:AbstractVector{T},D,M<:GO.Manifold,Tree,CRS} <: DD.Dimensions.Lookup{T, 1}
    manifold::M
    data::A
    tree::Tree
    dims::D
    crs::CRS
end
function GeometryLookup(data, dims=(DD.X(), DD.Y()); geometrycolumn=nothing, crs=nokw, tree=nokw)
    # First, retrieve the geometries - from a table, vector of geometries, etc.
    geometries = GOCore.get_geometries(data; geometrycolumn)
    geometries = Missings.disallowmissing(geometries)

    if isnokw(crs)
        crs = GI.crs(data)
        if isnothing(crs)
            crs = GI.crs(first(geometries))
        end
    end
    
    # Check that the geometries are actually geometries
    if any(!GI.isgeometry, geometries)
        throw(ArgumentError("""
        The collection passed in to `GeometryLookup` has some elements that are not geometries 
        (`GI.isgeometry(x) == false` for some `x` in `data`).
        """))
    end
    # Make sure there are only two dimensions
    if length(dims) != 2
        throw(ArgumentError("""
        The `dims` argument to `GeometryLookup` must have two dimensions, but it has $(length(dims)) dimensions (`$(dims)`).
        Please make sure that it has only two dimensions, like `(X(), Y())`.
        """))
    end
    # Build the lookup accelerator tree
    tree = if isnokw(tree)
        SortTileRecursiveTree.STRtree(geometries)
    elseif GO.SpatialTreeInterface.isspatialtree(tree)
        if tree isa Type
            tree(geometries)
        else
            tree
        end
    elseif isnothing(tree)
        nothing
    else
        throw(ArgumentError("""
        Got an argument for `tree` which is not a valid spatial tree (according to `GeometryOps.SpatialTreeInterface`)
        nor `nokw` or `nothing`

        Type is $(typeof(tree))
        """))
    end
    # TODO: auto manifold detection and best tree type for that manifold
    GeometryLookup(GO.Planar(), geometries, tree, dims, crs)
end

GI.crs(l::GeometryLookup) = l.crs
RA.setcrs(l::GeometryLookup, crs) = DD.rebuild(l; crs)

# Rasters has a no-op fallback `reproject(target, l::Lookup) = l`, but a
# GeometryLookup holds actual coordinates, so it has to be reprojected for
# real. This also makes `reproject` work on Geometry dims and whole vector
# data cubes through Rasters' dimension/array methods.
function RA.reproject(target::RA.GeoFormat, l::GeometryLookup)
    isnothing(l.crs) && throw(ArgumentError("Cannot reproject a `GeometryLookup` with no crs. Set one first with `setcrs`."))
    new_data = GO.reproject(l.data; source_crs=l.crs, target_crs=target, always_xy=true)
    return DD.rebuild(l; data=new_data, crs=target)
end

#=

## DD methods for the lookup

Here we define DimensionalData's methods for the lookup.
This is broadly standard except for the `rebuild` method, which is used to update the tree accelerator when the data changes.
=#

DD.dims(l::GeometryLookup) = l.dims
# This has to return itself
# DD.dims(d::DD.Dimension{<:GeometryLookup}) = dims(val(d))
DD.order(::GeometryLookup) = Lookups.Unordered()
DD.parent(lookup::GeometryLookup) = lookup.data
# TODO: format for geometry lookup
DD.Dimensions.format(l::GeometryLookup, D::Type, values, axis::AbstractRange) = l

# Make sure that the tree is rebuilt if the data changes
function DD.rebuild(
        lookup::GeometryLookup; 
        data=lookup.data, tree=nokw, 
        dims=lookup.dims, crs=nokw, 
        manifold=nokw, metadata=nokw
    )
    # TODO: metadata support for geometry lookup
    new_tree = if isnokw(tree)
        if data == lookup.data
            lookup.tree
        elseif isempty(data)
            nothing
        else
            SortTileRecursiveTree.STRtree(data)
        end
    elseif GO.SpatialTreeInterface.isspatialtree(tree)
        if tree isa Type
            tree(data)
        else
            tree
        end
    else
        SortTileRecursiveTree.STRtree(data)
    end
    new_crs = if isnokw(crs)
        data_crs = GI.crs(data)
        if isnothing(data_crs)
            lookup.crs
        else
            data_crs
        end
    else
        crs
    end

    new_manifold = isnokw(manifold) ? lookup.manifold : manifold

    return GeometryLookup(new_manifold, Missings.disallowmissing(data), new_tree, dims, new_crs)
end

# # Bounds - get the bounds of the lookup
# function Lookups.bounds(lookup::GeometryLookup)
#     if isempty(lookup.data)
#         Extents.Extent(NamedTuple{DD.name.(lookup.dims)}(ntuple(2) do i; (nothing, nothing); end))
#     else
#         if isnothing(lookup.tree)
#             mapreduce(GI.extent, Extents.union, lookup.data)
#         else
#             GI.extent(lookup.tree)
#         end
#     end
# end

# Return an `Int` or Vector{Bool}
# Base case: got a standard index that can go into getindex on a base Array
Lookups.selectindices(lookup::GeometryLookup, sel::Lookups.StandardIndices) = sel
# other cases: 
# - decompose selectors
function Lookups.selectindices(lookup::GeometryLookup, sel::DD.DimTuple)
    DD.selectindices(lookup, map(_val_or_nothing, DD.sortdims(sel, DD.dims(lookup))))
end
function Lookups.selectindices(lookup::GeometryLookup, sel::NamedTuple{K}) where K
    dimsel = map(DD.rebuild, map(DD.name2dim, K), DD.values(sel))
    DD.selectindices(lookup, dimsel) 
end
function Lookups.selectindices(lookup::GeometryLookup, sel::Tuple)
    if (length(sel) == length(DD.dims(lookup))) && all(map(s -> s isa At, sel))
        i = findfirst(x -> all(map(DD.Dimensions._matches, sel, x)), lookup)
        isnothing(i) && _coord_not_found_error(sel)
        return i
    else
        return [DD.Dimensions._matches(sel, x) for x in lookup]
    end
end
function Lookups.selectindices(lookup::GeometryLookup, sel::DD.Lookups.Contains)
    sel_ext = GI.extent(val(sel))
    potential_candidates = _maybe_get_candidates(lookup, sel_ext)
    filter(potential_candidates) do candidate
        GO.contains(lookup.data[candidate], val(sel))
    end
end
function Lookups.selectindices(lookup::GeometryLookup, sel::DD.Lookups.At)
    geom = val(sel)
    @assert GI.isgeometry(geom)
    candidates = _maybe_get_candidates(lookup, GI.extent(geom))
    x = findfirst(candidates) do candidate
        GO.equals(geom, lookup.data[candidate])
    end
    if isnothing(x)
        throw(ArgumentError("$sel not found in lookup"))
    else
        return candidates[x]
    end
end
function Lookups.selectindices(lookup::GeometryLookup, sel::DD.Lookups.Near)
    geom = val(sel)
    @assert GI.isgeometry(geom)
    # TODO: temporary
    @assert GI.trait(geom) isa GI.PointTrait "Only point geometries are supported for the near lookup at this point!  We will add more geometry support in the future."

    # Get the nearest geometry
    # TODO: this sucks!  Use some branch and bound algorithm
    # on the spatial tree instead.
    # if pointtrait
    return findmin(x -> GO.distance(geom, x), lookup.data)[2]
    # else
    #     findmin(x -> GO.distance(GO.GEOS(), geom, x), lookup.data)[2]
    # end 
    # this depends on LibGEOS being installed.

end
function Lookups.selectindices(lookup::GeometryLookup, sel::DD.Lookups.Touches)
    sel_ext = GI.extent(val(sel))
    potential_candidates = _maybe_get_candidates(lookup, sel_ext)
    return filter(potential_candidates) do candidate
        GO.intersects(lookup.data[candidate], val(sel))
    end
end
function Lookups.selectindices(
    lookup::GeometryLookup, 
    (xs, ys)::Tuple{Union{<:DD.Lookups.Touches}, Union{<:DD.Lookups.Touches}}
)
    target_ext = Extents.Extent(X = (first(xs), last(xs)), Y = (first(ys), last(ys)))
    potential_candidates = _maybe_get_candidates(lookup, target_ext)
    return filter(potential_candidates) do candidate
        GO.intersects(lookup.data[candidate], target_ext)
    end
end
function Lookups.selectindices(
    lookup::GeometryLookup, 
    (xs, ys)::Tuple{Union{<:DD.IntervalSets.ClosedInterval},Union{<:DD.IntervalSets.ClosedInterval}}
)
    target_ext = Extents.Extent(X = extrema(xs), Y = extrema(ys))
    potential_candidates = _maybe_get_candidates(lookup, target_ext)
    filter(potential_candidates) do candidate
        GO.covers(target_ext, lookup.data[candidate])
    end
end
function Lookups.selectindices(
    lookup::GeometryLookup, 
    (x, y)::Tuple{Union{<:DD.Lookups.At,<:DD.Lookups.Contains}, Union{<:DD.Lookups.At,<:DD.Lookups.Contains}}
)
    xval, yval = val(x), val(y)
    # `At` requires an exact match, so a point matching no geometry is an error;
    # `Contains` just returns the (possibly empty) set of matches.
    is_at = x isa DD.Lookups.At && y isa DD.Lookups.At
    potential_candidates = if isnothing(lookup.tree)
        collect(eachindex(lookup.data))
    else
        # TODO: this should be
        # lookup_ext = Lookups.bounds(lookup)
        # but that relies on https://github.com/rafaqz/DimensionalData.jl/pull/991/
        # so for now we need GI.extent(lookup.tree)
        lookup_ext = GI.extent(lookup.tree)

        # This is a specialized implementation to save on lookups
        if lookup_ext.X[1] <= xval <= lookup_ext.X[2] && lookup_ext.Y[1] <= yval <= lookup_ext.Y[2]
            GO.SpatialTreeInterface.query(lookup.tree, (xval, yval))
        else
            Int[]
        end
    end
    if is_at
        # If both selectors are At(), return a single index (first match) for speed and clarity
        for candidate in potential_candidates
            if GO.contains(lookup.data[candidate], (xval, yval))
                return candidate
            end
        end
        throw(ArgumentError("Point ($xval, $yval) not found in lookup"))
    else
        # For Contains selectors, return all matching indices
        filter(potential_candidates) do candidate
            GO.contains(lookup.data[candidate], (xval, yval))
        end
    end
end

for fname in (:equals, :intersects, 
                :contains, :within, :covers, 
                :coveredby, :touches)
    @eval begin
        function Lookups.selectindices(lookup::GeometryLookup, sel::DD.Lookups.Where{Base.Fix2{typeof(GO.$fname)}})
            sel_ext = GI.extent(val(sel).x)
            potential_candidates = _maybe_get_candidates(lookup, sel_ext)
            f = val(sel)
            return filter(potential_candidates) do idx
                f(lookup.data[idx])
            end
        end
    end
end
# Disjoint needs a specialized implementation, which will look at intersects instead.
function Lookups.selectindices(lookup::GeometryLookup, sel::DD.Lookups.Where{Base.Fix2{typeof(GO.disjoint)}})
    sel_ext = GI.extent(val(sel).x)
    potential_candidates = _maybe_get_candidates(lookup, sel_ext, Extents.intersects)
    f = val(sel)
    actual_intersections = filter(potential_candidates) do idx
        f(lookup.data[idx])
    end
    return setdiff(1:length(lookup.data), actual_intersections)
end

@inline Lookups.reducelookup(l::GeometryLookup) = Lookups.NoLookup(Base.OneTo(1))

function Lookups.show_properties(io::IO, mime, lookup::GeometryLookup)
    print(io, " ")
    show(IOContext(io, :inset => "", :dimcolor => 244), mime, DD.basedims(lookup))
end

# Dimension methods

@inline _reducedims(lookup::GeometryLookup, dim::DD.Dimension) =
    DD.rebuild(dim, [map(x -> zero(x), dim.val[1])])

function DD.format(dim::DD.Dimension{<:GeometryLookup}, axis::AbstractRange)
    # checkaxis(dim, axis)
    return dim
end

# Local functions
_val_or_nothing(::Nothing) = nothing
_val_or_nothing(d::DD.Dimension) = val(d)

# Get the candidates for the selector extent.  
# If the selector extent is disjoint from the tree rootnode extent,
# you should raise an error.  We should have an error type that can be
# plotted etc. to allow debugging and understanding.
function _maybe_get_candidates(lookup::GeometryLookup, selector_extent, operation::O) where O
    tree = lookup.tree
    isnothing(tree) && return 1:length(lookup)
    Extents.disjoint(GI.extent(tree), selector_extent) && return Int[]
    potential_candidates = GO.SpatialTreeInterface.query(
        tree,
        Base.Fix1(operation, selector_extent)
    )
    isempty(potential_candidates) && return Int[]
    return potential_candidates
end

_maybe_get_candidates(lookup::GeometryLookup, selector_extent) = _maybe_get_candidates(lookup, selector_extent, Extents.intersects)