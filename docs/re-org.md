# bgcArgoDMQC

Document to plan out new organizational structure for version 0.3.0. Some
constituent files have become too long to reasonably work on.

## core

Main functions that can be accessed directory from the parent package. These
should be the "front-facing" functions, those that serve the principal
purpose of the package.

## `configure`

Functions to set up permanent directories, metadata, etc. to improve workflow.

## `io`

Functions responsible for file input/output, as well as file fetching and
downloading. Includes netCDF file manipulation.

- fetch
- nc

## `interp`

Interpolation functions, mainly taking reference data and mapping it to a
float trajectory.

## `fplt`

Plotting and data visualization.

note: re-name to `plot`.

## `unit`

For converting between various units.

## `util`

Ad-hoc functions for indexing, etc.

## `lut`

Tempertaure/boundary layer thickness response time lookup table.
