#=
# Zonal statistics returning vector data cubes

This file extends `Rasters.zonal` so that passing a [`GeometryLookup`](@ref)
(or a [`Geometry`](@ref) dimension wrapping one) as the `of` keyword returns a
vector data cube instead of a plain vector:

```julia
zonal(mean, raster; of = GeometryLookup(geometries))
```

returns a `Raster` over a `Geometry` dimension carrying the lookup, so the
result can be indexed with spatial selectors like `Geometry(Contains(point))`.

If the input raster has more dimensions than the lookup spans (e.g. `Ti` or
`Band` on top of `X` and `Y`), `f` is applied to each spatial slice
(conceptually `mapslices(f, masked; dims = (X, Y))`), and the result is a
data cube over the leftover dimensions plus `Geometry`. This is controlled
by the `spatialslices` keyword:

- `spatialslices = true` (the default): apply `f` per spatial slice, return
  a cube over `(otherdims..., Geometry)`.
- `spatialslices = false`: apply `f` to the whole masked raster per geometry,
  return a vector over `Geometry` only.
- `spatialslices = dims`: slice over the given dimensions instead of the
  lookup's spatial dims.

All cropping, masking, missing-value handling, threading, and progress
reporting is reused from Rasters' own zonal machinery - none of it is
reimplemented here.
=#

# Entry points: hook into Rasters' `_zonal(f, x, of; kw...)` dispatch.
RA._zonal(f, x::RA.RasterStackOrArray, of::Geometry; kw...) =
    RA._zonal(f, x, val(of); kw...)
# Stacks fan out by layer, so layers with different dimensions each produce
# a cube of the right shape. (This also sidesteps Rasters' stack-zonal path
# for geometries outside the raster, where `map(_ -> missing, st)` no longer
# produces a layer-wise NamedTuple on DimensionalData 0.30.)
function RA._zonal(f, st::RA.AbstractRasterStack, lookup::GeometryLookup; kw...)
    K = keys(st)
    layers = map(K) do k
        RA._zonal(f, st[k], lookup; kw...)
    end
    return RA.RasterStack(NamedTuple{K}(layers))
end
function RA._zonal(f, x::RA.AbstractRaster, lookup::GeometryLookup;
    spatialslices=true, emptyval=nokw, skipmissing=true, progress=true, threaded=true, kw...
)
    geoms = lookup.data
    isempty(geoms) && throw(ArgumentError("Cannot compute zonal statistics with an empty `GeometryLookup`."))
    slicedims = _zonal_slicedims(spatialslices, x, lookup)
    zs = if isnothing(slicedims)
        RA._zonal(f, x, nothing, geoms; skipmissing, emptyval, progress, threaded, kw...)
    else
        # When slicing, `emptyval` has to be applied per slice inside the
        # wrapper, not by Rasters per geometry, where it would produce a
        # scalar instead of a slice-shaped result. The per-geometry loop is
        # also run here rather than through Rasters' allocation path, which
        # types its result vector from the first geometry and so cannot hold
        # slice results whose eltype differs between geometries (e.g. an
        # all-`emptyval` result for a geometry smaller than a grid cell).
        inner = _SpatialSliceify(f, slicedims, emptyval)
        _zonal_eachgeom(inner, x, geoms; skipmissing, progress, threaded, kw...)
    end
    otherdims = isnothing(slicedims) ? () : DD.otherdims(x, slicedims)
    return _geometry_cube(x, zs, Geometry(lookup), otherdims)
end

# Like Rasters' `_zonal(f, x, ::Nothing, geoms)`, reusing its per-geometry
# crop/mask path and `_run` threading/progress, but collecting into an
# untyped vector that is narrowed afterwards.
function _zonal_eachgeom(f, x, geoms; skipmissing, progress, threaded, kw...)
    zs = Vector{Any}(undef, length(geoms))
    RA._run(eachindex(zs), threaded, progress, "Applying $f to each geometry...") do i
        zs[i] = RA._zonal(f, x, geoms[i]; skipmissing, emptyval=nokw, kw...)
    end
    return map(identity, zs)
end

_zonal_slicedims(spatialslices::Bool, x, lookup) =
    spatialslices ? DD.dims(x, DD.dims(lookup)) : nothing
_zonal_slicedims(spatialslices, x, lookup) = DD.dims(x, spatialslices)

#=
## Spatial slicing

`_SpatialSliceify` wraps `f` so that, instead of reducing the whole masked
raster to one value, it reduces each spatial slice and returns a `Raster`
over the remaining dimensions. Rasters' zonal pipeline calls the wrapped
function either with the masked raster directly (`skipmissing = false`) or
with `skipmissing(masked)` - the `Base.SkipMissing` method below unwraps the
latter and re-applies `skipmissing` per slice instead.
=#
struct _SpatialSliceify{F,D,E}
    f::F
    dims::D
    emptyval::E
end

(s::_SpatialSliceify)(x::DD.AbstractDimArray) =
    _mapspatialslices(_empty_aware(s.f, s.emptyval), x, s.dims)
(s::_SpatialSliceify)(sm::Base.SkipMissing) =
    _mapspatialslices(_empty_aware(s.f, s.emptyval) ∘ Base.skipmissing, sm.x, s.dims)

# If `emptyval` was passed, return it for empty (e.g. fully-masked) slices
# instead of calling `f` on an empty iterator.
function _empty_aware(f, emptyval)
    isnokw(emptyval) && return f
    return el -> isempty(el) ? emptyval : f(el)
end

function _mapspatialslices(g, x::DD.AbstractDimArray, slicedims)
    otherdims = DD.otherdims(x, slicedims)
    isempty(otherdims) && return g(x)
    slices = eachslice(x; dims=otherdims)
    return DD.rebuild(x; data=[g(slice) for slice in slices], dims=DD.dims(slices), refdims=())
end

#=
## Concatenation along the Geometry dimension

Per-geometry results come back from Rasters as a `Vector` holding scalars,
`Raster`s (when slicing), or `missing` (geometry entirely outside the
raster). These are assembled into a vector data cube along a `Geometry`
dimension that carries the lookup.
=#
function _geometry_cube(x::RA.AbstractRaster, zs::AbstractVector, geomdim::Geometry, otherdims::Tuple)
    i = findfirst(z -> z isa DD.AbstractDimArray, zs)
    if isnothing(i)
        # Scalar results: a vector over Geometry only...
        isempty(otherdims) && return RA.Raster(zs, (geomdim,); name=DD.name(x))
        # ...unless slicing over `otherdims` was requested and every geometry
        # was outside the raster - then keep the cube shape, so the output
        # dimensionality doesn't depend on data coverage.
        data = Base.stack(map(z -> fill(z, length.(otherdims)), zs))
        return RA.Raster(data, (otherdims..., geomdim); name=DD.name(x))
    end
    # Sliced results: concatenate along a new last dimension.
    # Geometries entirely outside the raster came back as `missing` and are
    # expanded to missing-filled slices.
    template = zs[i]
    arrays = map(zs) do z
        z isa DD.AbstractDimArray ? parent(z) : fill(z, size(template))
    end
    data = Base.stack(arrays)
    return RA.Raster(data, (DD.dims(template)..., geomdim); name=DD.name(x))
end
