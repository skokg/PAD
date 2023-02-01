# ------------------------------
# Import the needed libraries
# ------------------------------

# Import the PAD python library
# The PAD python library file (PY_PAD_library.py) as well as the PAD C++ shared library file (PAD_Cxx_shared_library.so) 
# need to be in the same folder as the script using them
from PY_PAD_library import *

# Import matplotlib 
import matplotlib.pyplot as plt

# ------------------------------
# Construct sample fields
# ------------------------------

# Set domain size
dimx=300

# Construct two fields filled with zeroes
fa=np.zeros((dimx, dimx))
fb=fa.copy()

# In the fields, define two identically sized (50 x 50 grid points) non-zero regions with values equal to one, which are displaced by 100 grid points. 
fa[100:150,100:150]=1
fb[100:150,200:250]=1

# ------------------------------------------------------
# Calculate the PAD value and the PAD attribution PDF
# ------------------------------------------------------

# Calculate and print the PAD value (the output for the sample fields will be approximately 100 grid points)
PAD_value = calculate_PAD(fa,fb)
print(PAD_value)

# Calculate the PAD attribution PDF (Probability Density Function) 
PAD_PDF = calculate_PAD_attribution_PDF(fa, fb)

# Display the PAD attribution PDF on a graph
fig, ax = plt.subplots()
plt.fill_between(PAD_PDF[:,0], PAD_PDF[:,1], 0, linestyle='-')
mean_value=np.sum(PAD_PDF[:,0]*PAD_PDF[:,1])/np.sum(PAD_PDF[:,1])
plt.axvline(mean_value, color="navy", label="mean", linestyle='--')
plt.ylim(bottom = 0)
plt.xlabel("attribution distance")
plt.ylabel("PDF")
leg=plt.legend(loc = "upper right")
plt.show()
plt.close()

# ----------------------------------------------------------
# Calculate the list of PAD attributions and visualize them 
# ----------------------------------------------------------

# Calculate and print the PAD attribution list
PAD_attributions = calculate_PAD_attributions(fa,fb)
print(PAD_attributions)

# Make a simple visualization of the attributions 
# Only a randomly chosen subset of attributions is shown since showing all attributions 
# would make the visualization unintelligible
number_of_shown_attributions = 10 

np.random.shuffle(PAD_attributions)
PAD_attributions_subset = PAD_attributions[:number_of_shown_attributions]
fig, ax = plt.subplots()
fx = np.transpose(np.asarray([1-fa/np.max(fa),(1-fa/np.max(fa))*(1-fb/np.max(fb)), 1-fb/np.max(fb)]),(1, 2, 0))
plt.imshow(fx,origin='lower')
plt.plot(PAD_attributions_subset[:,[2,4]].transpose(), PAD_attributions_subset[:,[3,5]].transpose() ,'-ok',  alpha=0.3, markersize = 2, mfc='C1', mec='C1')
plt.show()
plt.close()

