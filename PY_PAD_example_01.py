# ------------------------------
# Import the needed libraries
# ------------------------------

# Import the PAD python library
# The PAD python library file (PY_PAD_library.py) as well as the PAD C++ shared library file (PAD_Cxx_shared_library.so) 
# need to be in the same folder as the script using them
from PY_PAD_library import *

# Import matplotlib 
import matplotlib as matplotlib
import matplotlib.pyplot as plt

# ------------------------------
# Construct sample fields
# ------------------------------

# Set domain size
dimx=500
dimy=300

# Construct two fields filled with zeroes
fa=np.zeros((dimy, dimx))
fb=fa.copy()

# In the fields, define two identically sized (50 x 50 grid points) non-zero regions with values equal to one, which are displaced by 100 grid points. 
fa[100:150,100:150]=1
fb[100:150,200:250]=1

# ------------------------------------------------------
# Calculate the PAD value and the PAD attribution PDFs
# ------------------------------------------------------

# Calculate and print the PAD value (the output for the sample fields will be approximately 100 grid points)
PAD_attributions = calculate_PAD_attributions(fa,fb)
print(PAD_attributions)

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
hist = np.histogram(PAD_attributions[:,0], bins=100, range=(0,np.max(PAD_attributions[:,0])), density=True, weights=PAD_attributions[:,1])
PAD_PDF = np.asarray([hist[1][:-1], hist[0][:]]).transpose(1,0)
fig, ax = plt.subplots()
plt.fill_between(PAD_PDF[:,0], PAD_PDF[:,1], 0, linestyle='-')
plt.axvline(PAD_distance, color="navy", label="PAD_distance = "+str(PAD_distance), linestyle='--')
plt.ylim(bottom = 0, top = 1.2*np.max(PAD_PDF[10:,1]))
plt.xlim(left = 0)
plt.xlabel("Attribution distance")
plt.ylabel("PDF")
leg=plt.legend(loc = "upper right")
plt.show()
plt.close()


# Display the two-dimensional PDF 
dx = PAD_attributions[:,4] - PAD_attributions[:,2] 
dy = PAD_attributions[:,5] - PAD_attributions[:,3] 
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
ax.set_xlabel("$\Delta x$" , fontsize = 13)
ax.set_ylabel("$\Delta y$" , fontsize = 13)
ax.tick_params(axis='both', labelsize = 13)
plt.show()
plt.close()

# ----------------------------------------------------------
# Visualize PAD attributions
# ----------------------------------------------------------

# number of attribution lines shown in the figure
number_of_shown_attributions = 50
# cumulative distribution
cumulative = np.cumsum(PAD_attributions[:, 1]/np.sum(PAD_attributions[:, 1]))
# randomly select attributions - the probability of selection is affected by the attribution value
rand = np.random.rand(number_of_shown_attributions)	
ind = np.searchsorted(cumulative,rand)
PAD_attributions_selected=PAD_attributions[ind]

cmap_b = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',["white",(0.3,0.3,1.0)],512)
cmap_r = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',["white",(1.0,0.3,0.3)],512)
norm_b = matplotlib.colors.Normalize()
norm_r = matplotlib.colors.Normalize()
# convert to rgb values using colormaps
fax = cmap_b(norm_b(fa))
fbx = cmap_r(norm_r(fb))
# get the combined colors for the image using the multiply effect
fx = fax*fbx
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1)
img = plt.imshow(fx, interpolation='nearest', origin='lower')
ax.plot(PAD_attributions_selected[:,[2,4]].transpose(), PAD_attributions_selected[:,[3,5]].transpose() ,'-ok',  alpha=0.3, markersize = 2, mfc='black', mec='black')
#ax.coastlines(resolution='110m', color='grey', linestyle='-', alpha=1)
plt.show()
plt.close()


