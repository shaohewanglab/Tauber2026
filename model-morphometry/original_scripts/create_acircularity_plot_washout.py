import numpy as np

import matplotlib.pyplot as plt
import matplotlib

import PlotLibrary as plotlib

from matplotlib.ticker import (MultipleLocator)

cmap = matplotlib.colormaps['Greys']

# Plot Figure settings
figSpecs = plotlib.FigureSettings()
figSpecs.set_journal('PNAS')
figSpecs.set_figureHeight(35)

# What is the size of this figure relative to the journal specs
xFraction = 0.5
yFraction = 1.0

width = xFraction*figSpecs.singleColumn
height = yFraction*figSpecs.figureHeight

# Figure Settings
fig = plt.figure(figsize=(width,height)) 
ax = plt.subplot()

######################################################################
# Functions
######################################################################

def circularity(P, A):
	return (2*np.sqrt(np.pi*A))/ P

def acircularity(P, A):
	return 1/circularity(P, A) - 1

def compute_area(R):
	# Extract x and y coordinates
	x = R[:, 0]
	y = R[:, 1]

	# Append the starting node to close the loop with the chord
	x = np.append(x, x[0])
	y = np.append(y, y[0])

	# Compute area using the shoelace formula
	Area = 0.5 * np.abs(np.sum(x[:-1] * y[1:]) - np.sum(y[:-1] * x[1:]))

	return Area

def compute_perimeter(R):
	# Compute distances between consecutive nodes
	distances = np.sqrt(np.sum((R[1:] - R[:-1])**2, axis=1))

	# Sum up the distances to get the total length
	Perimeter = np.sum(distances)

	return Perimeter

######################################################################
# Script for plotting data
######################################################################

# Washout simulation
z_list = np.array([0.015])
for cnt, z in enumerate(z_list):
    
    # Load data from npz file          
    data = np.load(f"output/washout_paper/shape_array_r_5_z_{z}_kP_0.01_kn_5e-06_tg_1e+06_tgi_2e+06_M_1e-07_dt_0.1_expV0.1_washout.npz")
    shapes = data["data"]

    Nframes = np.shape(shapes)[2]

    area_array = np.zeros(Nframes)
    perimeter_array = np.zeros(Nframes)
    time_array = np.arange(Nframes)*0.01

    for i in range(Nframes):
        shape = shapes[:,:,i]
        area_array[i] = compute_area(shape)
        perimeter_array[i] = compute_perimeter(shape)
	
    chi = acircularity(perimeter_array,area_array) 

    ax.plot(time_array, 100*(chi-chi[0]), c="firebrick",label=f"Washout",lw=1,zorder=3)

# Reference simulations
z_list = np.array([0.001,0.015])
c_list = np.array([0.44,0.70])
for cnt, z in enumerate(z_list):
    
    # Load data from npz file          
    data = np.load(f"output/thickness_paper/shape_array_r_5_z_{z}_kP_0.01_kn_5e-06_tg_1e+06_tgi_2e+06_M_1e-07_dt_0.1_expV0.1.npz")
    shapes = data["data"]

    Nframes = np.shape(shapes)[2]

    area_array = np.zeros(Nframes)
    perimeter_array = np.zeros(Nframes)
    time_array = np.arange(Nframes)*0.01

    for i in range(Nframes):
        shape = shapes[:,:,i]
        area_array[i] = compute_area(shape)
        perimeter_array[i] = compute_perimeter(shape)
	
    chi = acircularity(perimeter_array,area_array) 

    print(c_list[cnt])
    ax.plot(time_array, 100*(chi-chi[0]), c=cmap(c_list[cnt]),label=f"$z={z_list[cnt]:.2g}$",lw=1)

######################################################################
# Final layout settings
######################################################################

# Final adjustments to the figure
plotlib.set_box(ax, alt_style=True)

ax.set_xlabel(r'$ \tau / \tau_g$')
ax.set_ylabel(r'$\Delta \chi_{2D} [\%]$')

ax.set_xlim([0,1.5])
ax.set_ylim([0,150])

leg = plotlib.set_legend(ax,pos=0)
ax.xaxis.set_major_locator(MultipleLocator(0.5))
ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.yaxis.set_major_locator(MultipleLocator(50))
ax.yaxis.set_minor_locator(MultipleLocator(10))

plotlib.set_position(ax,x=0.21,y=0.20,width=0.71,height=0.70)

# Save as PDF
plt.savefig(f"results/acircularity_washout_paper_exp.pdf", format="pdf", bbox_inches="tight",transparent=True)

plt.show()