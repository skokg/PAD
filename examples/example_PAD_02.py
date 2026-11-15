"""NetCDF fields, residuals and distance/displacement PDFs in a 2-by-2 figure.

Install pad_python, NumPy, Matplotlib, netCDF4 and Cartopy before running.
Cartopy may download coastline data when displaying the maps.
"""
from pathlib import Path

import cartopy.crs as ccrs
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, LogNorm, PowerNorm
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np

from pad_python import calculate_PAD_attributions, calculate_PAD_distance_from_attributions


def read_field(filename):
    path = Path(__file__).resolve().parent / "sample_data" / filename
    with Dataset(path, "r") as dataset:
        arrays = [dataset.variables[name][:] for name in ("lon", "lat", "precip")]
        if any(np.ma.is_masked(array) for array in arrays):
            raise ValueError("Example inputs must not contain masked values.")
        return [np.asarray(array) for array in arrays]


lon, lat, fa = read_field("sample_field_A.nc")
lon2, lat2, fb = read_field("sample_field_B.nc")
if (lon.ndim != 1 or lat.ndim != 1 or fa.shape != (lat.size, lon.size)
        or fb.shape != fa.shape or not np.array_equal(lon, lon2)
        or not np.array_equal(lat, lat2) or not np.all(np.isfinite(lon))
        or not np.all(np.isfinite(lat))):
    raise ValueError("The example requires matching fields on the same lon/lat axes.")

# Planar PAD still uses grid-cell distances; lon/lat are only for map display.
# Set a finite cutoff in grid cells to leave unmatched precipitation.
distance_cutoff = None
random_seed = 5489
attributions, remaining1, remaining2 = calculate_PAD_attributions(
    fa, fb, remove_overlap=True, normalize=False,
    distance_cutoff=distance_cutoff, random_seed=random_seed,
)
print(attributions)
PAD_distance = calculate_PAD_distance_from_attributions(attributions)
print(f"PAD distance: {PAD_distance:.2f} grid cells")

# Inputs, attribution amounts and residuals retain the original amount units.
fields = (fa, fb)
vmax = max(field.max() for field in fields)
field_norm = PowerNorm(gamma=0.5, vmin=0, vmax=vmax)
field_cmaps = [
    LinearSegmentedColormap.from_list("field_a", ["white", (0.25, 0.35, 1.)]),
    LinearSegmentedColormap.from_list("field_b", ["white", (1., 0.25, 0.25)]),
]
# The supplied lon/lat axes are regular; align image edges with the cell centres.
dlon, dlat = lon[1] - lon[0], lat[1] - lat[0]
image_extent = (lon[0] - dlon / 2, lon[-1] + dlon / 2,
                lat[0] - dlat / 2, lat[-1] + dlat / 2)
cutoff_label = "no cutoff" if distance_cutoff is None else f"cutoff: {distance_cutoff:g} grid cells"

fig = plt.figure(figsize=(14, 10), constrained_layout=True)
grid = fig.add_gridspec(2, 2)
projection = ccrs.PlateCarree()
axes = [fig.add_subplot(grid[0, column], projection=projection) for column in range(2)]
ax_pdf = fig.add_subplot(grid[1, 0])
ax_displacement = fig.add_subplot(grid[1, 1])
for ax, pair, title in zip(axes, [fields, (remaining1, remaining2)],
                           ["Fields and sampled attributions", "Unattributed precipitation"]):
    # Stretch small amounts and multiply the colours so neither field hides the other.
    rgb = np.ones((*fa.shape, 3))
    for field, cmap, label in zip(pair, field_cmaps, ["Field A", "Field B"]):
        rgb *= cmap(field_norm(field))[..., :3]
        fig.colorbar(ScalarMappable(norm=field_norm, cmap=cmap), ax=ax,
                     label=f"{label} (original amount units)", shrink=0.7)
    ax.imshow(rgb, origin="lower", extent=image_extent, transform=projection,
              interpolation="nearest")
    ax.set_extent((lon.min(), lon.max(), lat.min(), lat.max()), crs=projection)
    ax.coastlines(resolution="50m", color="grey", linewidth=0.6)
    ax.gridlines(draw_labels={"bottom": "x", "left": "y"}, linewidth=0.5, alpha=0.4)
    ax.set_title(f"{title} ({cutoff_label})")

# Sample 300 attribution lines, weighted by amount; convert grid endpoints to lon/lat.
probabilities = attributions[:, 1] / attributions[:, 1].sum()
indices = np.random.default_rng(random_seed).choice(len(attributions), size=300, p=probabilities)
x1, y1, x2, y2 = attributions[indices, 2:6].astype(np.intp).T
axes[0].plot(np.stack((lon[x1], lon[x2])), np.stack((lat[y1], lat[y2])),
             "k.-", transform=projection, alpha=0.3, linewidth=0.8, markersize=2)

# All distances and the weighted mean remain in grid cells, not kilometres.
ax_pdf.hist(attributions[:, 0], bins=100, weights=attributions[:, 1],
            density=True, color="steelblue")
ax_pdf.axvline(PAD_distance, color="black", linestyle="--",
              label=f"PAD distance: {PAD_distance:.1f} grid cells")
ax_pdf.set(xlabel="Attribution distance (grid cells)",
           ylabel="Probability density (per grid cell)",
           title=f"Attribution PDF ({cutoff_label})", xlim=(0, None))
ax_pdf.legend()

# Displacements use grid endpoints, not differences in latitude/longitude.
displacement = attributions[:, 4:6] - attributions[:, 2:4]
limit = max(1., np.abs(displacement).max())
density, x_edges, y_edges = np.histogram2d(
    *displacement.T, bins=49, range=[(-limit, limit)] * 2,
    weights=attributions[:, 1], density=True,
)
# As in the original example, show densities down to one-thousandth of the peak.
hist = ax_displacement.pcolormesh(
    x_edges, y_edges, np.ma.masked_equal(density.T, 0), cmap="Blues",
    norm=LogNorm(vmin=density.max() * 0.001, vmax=density.max()), shading="auto",
)
fig.colorbar(hist, ax=ax_displacement, label="Probability density (per grid cell²)")
ax_displacement.axhline(0, color="grey", linewidth=0.5)
ax_displacement.axvline(0, color="grey", linewidth=0.5)
ax_displacement.set(xlabel="Δx = x2 − x1 (grid cells)", ylabel="Δy = y2 − y1 (grid cells)",
                    title="2D attribution PDF", aspect="equal")
plt.show()
plt.close(fig)
