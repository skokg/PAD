"""Synthetic fields, residuals and distance/displacement PDFs in a 2-by-2 figure.

Install pad_python, NumPy and Matplotlib before running this example.
"""
import matplotlib.pyplot as plt
import numpy as np

from pad_python import calculate_PAD_attributions, calculate_PAD_distance_from_attributions

# Two 50-by-50 rainy squares, displaced by 100 columns on a 300-by-500 grid.
fa = np.zeros((300, 500))
fb = np.zeros_like(fa)
fa[100:150, 100:150] = 1
fb[100:150, 200:250] = 1

# Set a finite cutoff in grid cells to leave unmatched precipitation.
distance_cutoff = None
random_seed = 5489
attributions, remaining1, remaining2 = calculate_PAD_attributions(
    fa, fb, remove_overlap=True, normalize=True,
    distance_cutoff=distance_cutoff, random_seed=random_seed,
)
print(attributions)
PAD_distance = calculate_PAD_distance_from_attributions(attributions)
print(f"PAD distance: {PAD_distance:.2f} grid cells")

# Display inputs and residuals in the same normalized amount units.
fields = []
for field in (fa, fb):
    scaled = field / field.max()
    fields.append(scaled / scaled.sum())
vmax = max(field.max() for field in fields)
cutoff_label = "no cutoff" if distance_cutoff is None else f"cutoff: {distance_cutoff:g} grid cells"

fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)
axes = axes.ravel()
for ax, pair, title in zip(axes[:2], [fields, (remaining1, remaining2)],
                           ["Fields and sampled attributions", "Unattributed precipitation"]):
    for field, cmap, label in zip(pair, ["Blues", "Reds"], ["Field A", "Field B"]):
        image = ax.imshow(np.ma.masked_equal(field, 0), cmap=cmap, vmin=0, vmax=vmax,
                          alpha=0.6, origin="lower", interpolation="nearest")
        fig.colorbar(image, ax=ax, label=f"{label} (normalized amount)", shrink=0.7)
    ax.set(xlabel="Column (grid cells)", ylabel="Row (grid cells)",
           title=f"{title} ({cutoff_label})")

# Sample 50 attribution lines with probability proportional to their amounts.
probabilities = attributions[:, 1] / attributions[:, 1].sum()
indices = np.random.default_rng(random_seed).choice(len(attributions), size=50, p=probabilities)
for x1, y1, x2, y2 in attributions[indices, 2:]:
    axes[0].plot([x1, x2], [y1, y2], "k.-", alpha=0.3, linewidth=0.8, markersize=2)

# All attributions contribute to the weighted histogram and PAD mean.
ax_pdf = axes[2]
ax_pdf.hist(attributions[:, 0], bins=100, weights=attributions[:, 1],
            density=True, color="steelblue")
ax_pdf.axvline(PAD_distance, color="black", linestyle="--",
              label=f"PAD distance: {PAD_distance:.1f} grid cells")
ax_pdf.set(xlabel="Attribution distance (grid cells)",
           ylabel="Probability density (per grid cell)",
           title=f"Attribution PDF ({cutoff_label})", xlim=(0, None))
ax_pdf.legend()

# Displacement is from the first field to the second: (x2 - x1, y2 - y1).
displacement = attributions[:, 4:6] - attributions[:, 2:4]
limit = max(1., np.abs(displacement).max())
ax_displacement = axes[3]
hist = ax_displacement.hist2d(*displacement.T, bins=49, range=[(-limit, limit)] * 2,
                             weights=attributions[:, 1], density=True, cmap="Blues")
fig.colorbar(hist[3], ax=ax_displacement, label="Probability density (per grid cell²)")
ax_displacement.axhline(0, color="grey", linewidth=0.5)
ax_displacement.axvline(0, color="grey", linewidth=0.5)
ax_displacement.set(xlabel="Δx (grid cells)", ylabel="Δy (grid cells)",
                    title="2D attribution PDF", aspect="equal")
plt.show()
plt.close(fig)
