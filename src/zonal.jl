# Zonal statistics that return vector data cubes. `VectorDataCubes.zonal` is a
# package-owned function (not a method of `Rasters.zonal`, which can't dispatch
# on its `of` keyword): geometry-lookup `of`s are handled here, everything else
# forwards to `Rasters.zonal`.

"""
    VectorDataCubes.zonal(f, x; of, kw...)

Calculate zonal statistics like `Rasters.zonal`, but return a vector data cube
when `of` is a [`GeometryLookup`](@ref) or a dimension wrapping one: a `Raster`
over a `Geometry` dimension carrying the lookup, so the result can be indexed
with spatial selectors like `Geometry(Contains(point))`. Any other `of` is
forwarded to `Rasters.zonal` unchanged.

If `x` has more dimensions than the lookup spans (e.g. `Ti` or `Band` on top
of `X` and `Y`), `f` is applied to each spatial slice (conceptually
`mapslices(f, masked; dims = (X, Y))`) and the result is a cube over the
leftover dimensions plus `Geometry`. The `spatialslices` keyword controls
this: `true` (the default) slices over the lookup's spatial dims, `false`
applies `f` to the whole masked raster per geometry, and a tuple of dims
slices over those dims instead.

All other keywords (`skipmissing`, `emptyval`, `progress`, `threaded`, ...)
are passed through to Rasters' zonal machinery, which does all the cropping,
masking, and missing-value handling.

This function is deliberately not exported, since Rasters exports `zonal` too;
call it qualified, or bind it explicitly with `using VectorDataCubes: zonal`.
"""
zonal(f, x; of, kw...) = _zonal(f, x, of; kw...)

# Any `of` without a geometry lookup is Rasters' business.
_zonal(f, x, of; kw...) = RA.zonal(f, x; of, kw...)
_zonal(f, x, of::DD.Dimension{<:GeometryLookup}; kw...) = _zonal(f, x, val(of); kw...)
# Stacks fan out by layer, so layers with different dimensions each produce
# a cube of the right shape.
function _zonal(f, st::RA.AbstractRasterStack, lookup::GeometryLookup; kw...)
    K = keys(st)
    layers = map(K) do k
        _zonal(f, st[k], lookup; kw...)
    end
    return RA.RasterStack(NamedTuple{K}(layers))
end
function _zonal(f, x::RA.AbstractRaster, lookup::GeometryLookup;
    spatialslices=true, skipmissing=true, emptyval=nokw, kw...
)
    geoms = lookup.data
    isempty(geoms) && throw(ArgumentError("Cannot compute zonal statistics with an empty `GeometryLookup`."))
    # The same `open`/`_prepare_for_burning` preamble as `Rasters.zonal`.
    return Base.open(x) do o
        xp = RA._prepare_for_burning(o)
        slicedims = _zonal_slicedims(spatialslices, xp, lookup)
        # When slicing, `emptyval` has to be applied per slice inside the
        # wrapper, not by Rasters per geometry, where it would produce a
        # scalar instead of a slice-shaped result.
        inner, inner_emptyval = if isnothing(slicedims)
            f, emptyval
        else
            _SpatialSliceify(f, slicedims, emptyval), nokw
        end
        zs = RA._zonal(inner, xp, nothing, geoms; skipmissing, emptyval=inner_emptyval, kw...)
        otherdims = isnothing(slicedims) ? () : DD.otherdims(xp, slicedims)
        _geometry_cube(xp, zs, Geometry(lookup), otherdims)
    end
end

_zonal_slicedims(spatialslices::Bool, x, lookup) =
    spatialslices ? DD.dims(x, DD.dims(lookup)) : nothing
_zonal_slicedims(spatialslices, x, lookup) = DD.dims(x, spatialslices)

# `_SpatialSliceify` wraps `f` to reduce each spatial slice instead of the
# whole masked raster, returning a `Raster` over the remaining dims. Rasters
# passes the wrapped function `skipmissing(masked)` when `skipmissing=true`;
# that is unwrapped and `skipmissing` re-applied per slice.
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

# Assemble the per-geometry results (scalars, `Raster`s when slicing, or
# `missing` for geometries entirely outside the raster) into a vector data
# cube along a `Geometry` dimension carrying the lookup.
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
    # Geometries entirely outside the raster came back as `missing` and are
    # expanded to missing-filled slices.
    template = zs[i]
    arrays = map(zs) do z
        z isa DD.AbstractDimArray ? parent(z) : fill(z, size(template))
    end
    data = Base.stack(arrays)
    return RA.Raster(data, (DD.dims(template)..., geomdim); name=DD.name(x))
end
