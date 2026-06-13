# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this package is

VectorDataCubes.jl makes a **vector data cube**: a `DimensionalData`/`Rasters` array or
stack whose `Geometry` dimension is backed by a `GeometryLookup` (a geometry vector + an
STRtree). The payoff is spatial indexing — `cube[Geometry(Contains(point))]`,
`cube[X(a..b), Y(c..d)]` — on top of everything DimensionalData/Rasters already gives you.

Three entry points, one per `src/` file:

- `geometry_lookup.jl` — the `Geometry` dimension and `GeometryLookup` (the spatial-indexing core).
- `zonal.jl` — `VectorDataCubes.zonal`, aggregating a raster over geometries into a cube.
- `tables.jl` — `vectordatacube` / `vectordatacubetable`, round-tripping a table ↔ cube.

## Commands

Julia workspace package: root `Project.toml` has `[workspace] projects = ["docs", "test"]`; `test/` and `docs/` carry their own `Project.toml`.  When running code that wants to use this package, use the `docs/` environment, since that contains the ecosystem packages too.

```sh
# Full suite
julia --project=. -e 'using Pkg; Pkg.test()'
# Faster iteration (no sandbox build)
julia --project=test test/runtests.jl
# A single test file — each is self-contained (does its own imports), so just include it.
julia --project=test -e 'include("test/zonal.jl")'
# An example (each self-downloads its data to examples/data/)
julia --project=docs examples/02_zonal_countries.jl
```

`runtests.jl` includes each test file via `SafeTestsets.@safetestset`, so every file runs
in its own module — each must do its own `using VectorDataCubes` (and other imports);
nothing leaks in from `runtests.jl`. A new test file must follow suit.

`test/basics.jl` and the examples need network access. CI runs Julia 1.10 and `1`.

For iterative work, prefer the persistent Julia REPL via the `mcp__julia__*` tools — it
avoids re-paying Julia's per-process compile latency on every run.

## Architecture

- **`Lookups.selectindices` (in `geometry_lookup.jl`) is the heart.** Every spatial
  selector resolves to indices there: narrow with an STRtree extent query, then refine
  with an exact GeometryOps predicate. A new selector means a new method here.
- **The lookup spans `(X(), Y())` *and* the `Geometry` dim wrapping it** — that's why both
  `Geometry(...)` and `X()/Y()` selectors work on one axis. `DD.rebuild` rebuilds the
  STRtree when the data changes.
- **`zonal` is package-owned, not a `Rasters.zonal` method** (Rasters can't dispatch on
  `of`), and not exported (call it qualified). A `GeometryLookup` `of` yields a cube;
  anything else forwards to `Rasters.zonal`.

## Conventions

- **Import aliases**, consistent everywhere: `DD`, `GO`, `GOCore`, `GI`, `RA`, plus
  `Extents`, `Missings`, `SortTileRecursiveTree`.
- **`nokw` / `isnokw`** (from Rasters) is the "keyword not supplied" sentinel, distinct
  from a meaningful `nothing` (e.g. no CRS / no tree).
