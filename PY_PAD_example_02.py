# ------------------------------
# Import the needed libraries
# ------------------------------

# Import the PAD python library
# The PAD python library file (PY_PAD_library.py) as well as the PAD C++ shared library file (PAD_Cxx_shared_library.so) 
# need to be in the same folder as the script using them
import numpy as np
from pathlib import Path
from PY_PAD_library import calculate_PAD_attributions, calculate_PAD_distance_from_attributions

# Import matplotlib 
import matplotlib as matplotlib
import matplotlib.pyplot as plt

# Import netcdf library 
from netCDF4 import Dataset

# Import cartopy library 
import cartopy.crs as ccrs

# ------------------------------
# Read sample netcdf fields
# ------------------------------

def read_field(filename):
    with Dataset(Path(__file__).resolve().parent / filename, 'r') as dataset:
        arrays = [dataset.variables[name][:] for name in ("lon", "lat", "precip")]
        if any(np.ma.is_masked(array) for array in arrays):
            raise ValueError("Example inputs must not contain masked coordinates or precipitation.")
        return [np.asarray(array) for array in arrays]


lon, lat, fa = read_field("PY_PAD_example_02_field_A.nc")
lon2, lat2, fb = read_field("PY_PAD_example_02_field_B.nc")
if (lon.ndim != 1 or lat.ndim != 1 or fa.shape != (lat.size, lon.size)
        or fb.shape != fa.shape or not np.array_equal(lon, lon2) or not np.array_equal(lat, lat2)
        or not np.all(np.isfinite(lon)) or not np.all(np.isfinite(lat))):
    raise ValueError("The example requires matching 2D fields on the same 1D lon/lat axes.")

# ------------------------------------------------------
# Calculate the PAD value and the PAD attribution PDFs
# ------------------------------------------------------

# Planar distances are in grid cells; lon/lat are used only for the map below.
# True removes same-grid overlap first; False skips overlap preprocessing.
PAD_attributions, remaining1, remaining2 = calculate_PAD_attributions(
    fa, fb, remove_overlap=True, normalize=True, distance_cutoff=None, random_seed=5489,
)
print(PAD_attributions)
# The six columns are distance, normalized amount, x1, y1, x2, y2.
x1, y1, x2, y2 = PAD_attributions[:, 2:6].astype(np.intp).T

# Calculate the PAD attribution PDF (Probability Density Function) 
PAD_distance = calculate_PAD_distance_from_attributions(PAD_attributions)
print("PAD distance: " + str(PAD_distance))

# ------------------------------------------------------------------------------------------------------------------------
# Draw the attribution PDF 
# ------------------------------------------------------------------------------------------------------------------------

# Import matplotlib library 
import matplotlib
import matplotlib.pyplot as plt

# ---  PDF
hist = np.histogram(PAD_attributions[:,0], bins=100, range=(0,max(1.0,np.max(PAD_attributions[:,0]))), density=True, weights=PAD_attributions[:,1])
PAD_PDF = np.asarray([hist[1][:-1], hist[0][:]]).transpose(1,0)
fig, ax = plt.subplots()
plt.fill_between(PAD_PDF[:,0], PAD_PDF[:,1], 0, linestyle='-')
plt.axvline(PAD_distance, color="navy", label="PAD_distance = "+str(PAD_distance), linestyle='--')
plt.ylim(bottom = 0, top = 1.2*np.max(PAD_PDF[:,1]))
plt.xlim(left = 0)
plt.xlabel("Attribution distance (grid cells)")
plt.ylabel("PDF (per grid cell)")
leg=plt.legend(loc = "upper right")
plt.show()
plt.close()


# Display the two-dimensional PDF 
dx = x2 - x1
dy = y2 - y1
maxdistance = np.max(fa.shape)
hist2d = np.histogram2d(dx,dy,weights=PAD_attributions[:,1], bins=[50-1,50-1], range=[(-maxdistance,maxdistance),(-maxdistance,maxdistance)], density=True)
fig = plt.figure(figsize=(5, 5), linewidth = 3)
ax = fig.add_subplot(1, 1, 1)
cmap_b = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',["white","blue"])
norm_b = matplotlib.colors.LogNorm(vmin = np.max(hist2d[0]*0.001))
img_extent = (-maxdistance, maxdistance, -maxdistance, maxdistance)
img = ax.imshow(np.transpose(hist2d[0],(1,0)), interpolation='nearest', origin='lower', cmap = cmap_b, extent=img_extent, norm=norm_b)
ax.grid(which='major', color='grey', alpha=0.5, linestyle=':', linewidth=1)
ax.axvline(0, color='grey', alpha=0.8, linestyle='-', linewidth=1)
ax.axhline(0, color='grey', alpha=0.8, linestyle='-', linewidth=1)
ax.set_xlabel(r"$\Delta x$" , fontsize = 13)
ax.set_ylabel(r"$\Delta y$" , fontsize = 13)
ax.tick_params(axis='both', labelsize = 13)
plt.show()
plt.close()



# ----------------------------------------------------------
# Visualize PAD attributions
# ----------------------------------------------------------

# number of attribution lines shown in the figure
number_of_shown_attributions = 300
# cumulative distribution
cumulative = np.cumsum(PAD_attributions[:, 1]/np.sum(PAD_attributions[:, 1]))
cumulative[-1] = 1.0  # Protect searchsorted from cumulative-sum roundoff.
# randomly select attributions - the probability of selection is affected by the attribution value
rand = np.random.default_rng(5489).random(number_of_shown_attributions)
ind = np.searchsorted(cumulative,rand)

cmap_b = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',["white",(0.3,0.3,1.0)],512)
cmap_r = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',["white",(1.0,0.3,0.3)],512)
norm_b = matplotlib.colors.Normalize()
norm_r = matplotlib.colors.Normalize()
# convert to rgb values using colormaps
fax = cmap_r(norm_r(fa))
fbx = cmap_b(norm_b(fb))
# get the combined colors for the image using the multiply effect
fx = fax*fbx
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
img_extent = (np.min(lon), np.max(lon), np.min(lat), np.max(lat))
img = plt.imshow(fx, transform=ccrs.PlateCarree(), interpolation='nearest', origin='lower', extent=img_extent)
# Convert endpoint grid coordinates to map coordinates and display lines.
ax.plot(np.stack((lon[x1[ind]], lon[x2[ind]])), np.stack((lat[y1[ind]], lat[y2[ind]])), '-ok', alpha=0.3, markersize=2, mfc='black', mec='black', transform=ccrs.PlateCarree())
ax.coastlines(resolution='50m', color='grey', linestyle='-', alpha=1)
plt.show()
plt.close()
