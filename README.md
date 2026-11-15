# PAD Python: Precipitation Attribution Distance

Precipitation Attribution Distance (PAD) measures spatial displacement between
precipitation fields on a planar grid. Its results describe how much
precipitation is attributed between locations, the distances between those
locations, and the amounts left unattributed in each field.

This Python library provides a NumPy interface to the C++ implementation.
Attribution uses randomized nearest-neighbor searches with k-d trees.
Normalization, overlap removal and distance limits are controlled from Python.

The package is distributed and imported as `pad_python`.

GitHub repository: [skokg/PAD](https://github.com/skokg/PAD).

## Installation

Python 3.10 or newer is required. Install the package with:

```sh
pip install pad_python
```

To install from source:

```sh
pip install --no-binary=pad_python pad_python
```

Compiling from source requires a C++ compiler: GCC or Clang on Linux, MSVC Build
Tools on Windows, or Xcode Command Line Tools on macOS.

The package can also be installed from conda-forge:

```sh
conda install -c conda-forge pad_python
```

Alternatively, download the source from the
[GitHub repository](https://github.com/skokg/PAD) and run `make` in the
project's root directory on Linux. This requires GNU Make and a C++ compiler.
Running `make` creates `libpad_python.so` in the `src` folder.
Copy `pad_python.py` and `libpad_python.so` from `src` into the folder
containing your Python script. NumPy must also be installed.
No compiled library is included in this source checkout.

## Usage

Import the library with:

```python
import pad_python as pad
```

Inputs are precipitation amounts in two equally shaped, nonempty 2D NumPy
arrays. A field has shape `(y, x)` and represents a rectangular grid with equal
spacing along both axes. Values must be finite and nonnegative, with at least
one positive value in each field. Masked arrays are not supported.

The two main calculation functions are:

- `pad.calculate_PAD_attributions()`: attribute precipitation between fields
  and return the attribution records and the two residual fields.
- `pad.calculate_PAD_distance_from_attributions()`: calculate the amount-weighted
  mean attribution distance from those records.

By default, `normalize=True` independently scales each field to sum to one,
`remove_overlap=True` attributes pointwise overlap at zero distance first,
and `distance_cutoff=None` imposes no distance limit.
A fixed `random_seed` makes repeated calculations
reproducible for identical inputs and the same native-library build.
The selected seed is printed when nearest-neighbor work is needed.

Distances and cutoffs are in grid cells. Multiply a distance by the common
grid spacing to express it in physical units. With `normalize=False`, amounts
and residuals retain the input units. For equal-area cells, rainfall depths
give the same weighting as water volumes; the library does not calculate
cell areas.

PAD averages only attributed amounts. A cutoff can leave precipitation
unmatched; if there are no attributions, the mean distance is undefined and
the distance function raises `ValueError`. Detailed parameters and returned
columns are described in **Function descriptions** below.

This interface supports equally shaped planar fields. For latitude/longitude
coordinates and spherical distances, see
[PAD on Sphere](https://github.com/skokg/PAD_on_Sphere).

### Idealized planar field example

This example uses a grid with 120 rows and 200 columns. A rectangular region
has 10 mm of observed rainfall; the forecast contains the same region shifted
60 grid cells to the right, leaving 20 dry columns between the rectangles.
An attribution cutoff of 50 grid cells leaves the outer edges unmatched.

The grid cells are assumed to have equal areas. Setting `normalize=False`
retains rainfall depths in mm, so the residual fields can be plotted directly.
Weighting by depth is equivalent to weighting by volume on this grid.

```python
import numpy as np
import pad_python as pad

# Rainfall on a planar grid with equal spacing in both directions.
cutoff_cells = 50
observed = np.zeros((120, 200))
observed[40:80, 30:70] = 10.0
forecast = np.roll(observed, shift=60, axis=1)  # No rainfall wraps at the boundary.

# Keep original amounts: attribution weights and residuals are in mm.
attributions, remaining_forecast, remaining_observed = pad.calculate_PAD_attributions(
    forecast, observed, normalize=False, distance_cutoff=cutoff_cells, random_seed=42,
)

distance_cells = pad.calculate_PAD_distance_from_attributions(attributions)
print(f"PAD: {distance_cells:.2f} grid cells")
```

To display the fields, attributions and unattributed precipitation, append:

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)
axes = axes.ravel()
ax_pdf = axes[2]

for ax, fields, title in zip(axes[:2],
        [(observed, forecast), (remaining_observed, remaining_forecast)],
        ["Rainfall and sampled attributions", "Unattributed precipitation"]):
    for rain, cmap, label in zip(fields, ["Blues", "Reds"], ["Observed", "Forecast"]):
        image = ax.imshow(np.ma.masked_equal(rain, 0), cmap=cmap,
                          vmin=0, vmax=10, alpha=0.6, origin="lower",
                          interpolation="nearest")
        fig.colorbar(image, ax=ax, location="right", label=f"{label} (mm)", shrink=0.7)
    ax.set(xlabel="Column (grid cells)", ylabel="Row (grid cells)",
           title=f"{title} (cutoff: {cutoff_cells} grid cells)")

# Plot at most 100 attributions; all are still used to calculate PAD.
indices = np.random.default_rng(42).choice(len(attributions),
                                          min(100, len(attributions)), replace=False)
for x1, y1, x2, y2 in attributions[indices, 2:]:
    axes[0].plot([x1, x2], [y1, y2], "k.-", linewidth=0.8, markersize=3)

# Amount-weighted PDF using all attributions (also volume-weighted here).
ax_pdf.hist(attributions[:, 0], bins=25,
             weights=attributions[:, 1], density=True, color="steelblue")
ax_pdf.axvline(distance_cells, color="black", linestyle="--",
               label=f"PAD distance: {distance_cells:.1f} grid cells")
ax_pdf.legend()
ax_pdf.set(xlabel="Attribution distance (grid cells)",
           ylabel="Probability density (per grid cell)",
           title=f"Attribution PDF (cutoff: {cutoff_cells} grid cells)")

# Amount-weighted two-dimensional displacement PDF.
displacement = attributions[:, 4:6] - attributions[:, 2:4]
limit = max(1., np.abs(displacement).max())
ax_displacement = axes[3]
hist = ax_displacement.hist2d(*displacement.T, bins=25, range=[(-limit, limit)] * 2,
                             weights=attributions[:, 1], density=True, cmap="Blues")
fig.colorbar(hist[3], ax=ax_displacement, label="Probability density (per grid cell²)")
ax_displacement.axhline(0, color="grey", linewidth=0.5)
ax_displacement.axvline(0, color="grey", linewidth=0.5)
ax_displacement.set(xlabel="Δx (grid cells)", ylabel="Δy (grid cells)",
                    title="2D attribution PDF", aspect="equal")
plt.show()
```

Observed rainfall is blue and forecast rainfall is red. Dry cells are
transparent. The top-left panel shows a reproducible random sample of up to
100 attributions in black; zero-distance matches would appear as points.
The top-right panel shows observed and forecast residuals, both as
nonnegative rainfall depths in mm. Each field panel has two colourbars on
the right. The bottom-left panel is an amount-weighted PDF histogram of all
attribution distances; the bottom-right panel shows their two-dimensional
displacement PDF. Displacements are `Δx = x2 - x1` and `Δy = y2 - y1`, from
the first input field to the second. Here the forecast is the first input,
so attribution lines point toward the observed rectangle on its left.
Unattributed precipitation is excluded from both histograms.
This planar example needs NumPy and Matplotlib, without geographic map data.

PAD is the amount-weighted mean attribution distance, equivalent to the
volume-weighted mean for these equal-area cells. It is not simply the
60-cell shift of the rectangle. Observed columns 30–39 and forecast columns
120–129 are more than 50 grid cells from any counterpart and remain
unattributed. Other residuals depend on the attribution. Returned residuals
are already in mm because normalization is disabled; with the default
`normalize=True`, they would instead be in normalized amount units.

## Further examples

The source distribution and GitHub repository include two scripts and their
sample NetCDF data in `examples/`. Examples and data are not installed into
`site-packages` by a wheel. Download and extract the source to use them.
The examples additionally require Matplotlib, netCDF4 and Cartopy.

Then run the examples from the downloaded source root:

```sh
python examples/example_PAD_01.py
python examples/example_PAD_02.py
```

| Script | Input | Output |
| --- | --- | --- |
| `example_PAD_01.py` | Two synthetic fields, shape `(300, 500)`, with displaced rainy rectangles. | Printed attributions and PAD; a 2×2 figure with fields and sampled attribution lines, residuals, distance PDF and displacement PDF. |
| `example_PAD_02.py` | `sample_data/sample_field_A.nc` and `sample_data/sample_field_B.nc`, relative to the script. | Printed attributions and PAD; the same 2×2 layout with maps in the top row. |

Both scripts display figures rather than save output files. They can also be
run from another working directory because the second resolves its data paths
relative to its own file. Longitude/latitude in the second example are used
only for map display; the PAD calculation still uses grid-cell distances.
The second script disables normalization and retains the original amount units;
it applies no distance cutoff by default.
Its precipitation colours use a square-root stretch, and the 2D attribution
PDF uses logarithmic colours to reveal weaker features.
Cartopy may download coastline data the first time the map example is used.

## Function descriptions

The following describes the two calculation functions used above. Attribution
records and residual fields are NumPy arrays. Input fields are never modified.

### `pad.calculate_PAD_attributions(fa, fb, remove_overlap=True, distance_cutoff=None, random_seed=None, normalize=True)`

**Inputs:** `fa` and `fb` are real numeric NumPy arrays shaped `(y, x)`, with
identical shapes and nonempty dimensions. Each field must contain only finite,
nonnegative values and at least one positive value. Masked arrays are rejected.
The fields represent the same rectangular grid, with equal spacing along x
and y; no separate coordinate arrays are supplied. Values are converted to
`float64` for calculation and must remain finite, with a positive amount in
each field after conversion.

| Parameter | Meaning |
| --- | --- |
| `remove_overlap` | Boolean, default `True`. Attribute the pointwise minimum of the two prepared fields at zero distance before nearest-neighbor attribution. `False` skips this preprocessing, but zero-distance matches can still occur. |
| `distance_cutoff` | `None` for no limit, or a finite nonnegative Euclidean distance in grid cells. Matches exactly on the cutoff are permitted. Zero allows only zero-distance matches. |
| `random_seed` | `None` or `-1` chooses a random seed. An integer from `0` to `4294967295` selects a seed explicitly; Boolean values are rejected. A fixed seed reproduces results for identical inputs and the same native-library build. |
| `normalize` | Boolean, default `True`. Independently scale each field to sum to one before overlap removal and attribution. `False` retains the original amounts. |

**Output:** A list `[attributions, remaining1, remaining2]`:

| Array | Contents |
| --- | --- |
| `attributions` | A `float64` array shaped `(N, 6)`, with one row per attributed amount. It has shape `(0, 6)` when no matches are possible under the cutoff. |
| `remaining1` | A `float64` array with the shape of `fa`, containing its remaining unattributed amounts. |
| `remaining2` | A `float64` array with the shape of `fb`, containing its remaining unattributed amounts. |

The six attribution columns are:

| Column | Contents |
| --- | --- |
| `0`: distance | Euclidean distance between the endpoints, in grid cells. |
| `1`: amount | Amount attributed between the endpoints. |
| `2`: x1 | Zero-based column of the endpoint in `fa`. |
| `3`: y1 | Zero-based row of the endpoint in `fa`. |
| `4`: x2 | Zero-based column of the endpoint in `fb`. |
| `5`: y2 | Zero-based row of the endpoint in `fb`. |

Endpoint coordinates represent integer grid indices stored as `float64` values.
For array indexing, convert them to integers and use `[y, x]` order.
Attribution amounts and both residual fields use normalized amount units when
`normalize=True`, and the original input units when `normalize=False`.
Returned arrays do not share storage with the inputs and remain valid after
later calculation calls.

**Functionality:** Validates the fields and options, optionally normalizes each
field, and optionally attributes their pointwise overlap. It then alternates
between fields, selecting points in randomized order and attributing amounts
to the nearest remaining precipitation in the opposite field. Each match
transfers the smaller available amount, consuming at least one endpoint.
Points without an acceptable neighbor under the cutoff retain their amounts
in the residual output. Attribution stops when one field has no points left
for further searches.

The returned rows describe the attribution, rather than the scalar PAD value;
pass them to `calculate_PAD_distance_from_attributions` for that value.
Nearest-neighbor work prints the selected random seed, including when a seed
was supplied explicitly. No seed is printed when overlap removal leaves no
nearest-neighbor work. Invalid inputs or options raise `TypeError` or
`ValueError`; native calculation failures are reported as `RuntimeError`.

### `pad.calculate_PAD_distance_from_attributions(PAD_attributions)`

**Input:** A nonempty six-column attribution array, normally returned by
`calculate_PAD_attributions`. Array-like input convertible to `float64` is
also accepted. Column 0 contains distances and column 1 contains attribution
amounts. Both columns must be finite and nonnegative, with at least one
strictly positive amount. The four endpoint columns are not used in the
distance calculation.

**Output:** A Python `float` containing the amount-weighted mean attribution
distance. For attributions produced by this library, the unit is grid cells.
The result is `0.0` when all positively weighted distances are zero.

**Functionality:** Computes `sum(distance * amount) / sum(amount)` using scaled
values to reduce overflow risk. Zero-weight rows contribute nothing and are
excluded before choosing the distance scale; their distances must still pass
input validation. The input array is not modified.

Only attributed amounts enter the mean. Residual precipitation is not included,
and the cutoff used to construct the attributions therefore affects which
amounts contribute. An empty array, incorrect number of columns, invalid
distances or amounts, or entirely zero amounts raises `ValueError`; PAD is
undefined without a positive attributed amount.

## Author and license

Gregor Skok, Faculty of Mathematics and Physics, University of Ljubljana,
Slovenia.

Email: [Gregor.Skok@fmf.uni-lj.si](mailto:Gregor.Skok@fmf.uni-lj.si).
Please contact the author if you have questions or trouble using the package.

Licensed under the MIT License. See LICENSE, which also includes
the modified third-party k-d tree notice.

## References

The paper below describes precipitation attribution distance; please cite it
when publishing work using PAD:

- Skok, G. (2023). Precipitation attribution distance. Atmospheric Research, 295.
  [https://doi.org/10.1016/j.atmosres.2023.106998](https://doi.org/10.1016/j.atmosres.2023.106998)
