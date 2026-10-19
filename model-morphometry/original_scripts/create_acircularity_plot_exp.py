import numpy as np

import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import (MultipleLocator)

import PlotLibrary as plotlib

# Plot Figure settings
figSpecs = plotlib.FigureSettings()
figSpecs.set_journal('PNAS')
figSpecs.set_figureHeight(30)

# What is the size of this figure relative to the journal specs
xFraction = 60/178
yFraction = 1.0

width = xFraction*figSpecs.doubleColumn
height = yFraction*figSpecs.figureHeight

# Figure Settings
fig = plt.figure(figsize=(width,height)) 
ax = plt.subplot()

cmap = matplotlib.colormaps['rainbow']

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

z_list = np.array([0.0005, 0.001, 0.003, 0.005, 0.008, 0.01])
N = len(z_list)

for cnt, z in enumerate(z_list):
	
	# Load data from npz file          
	if cnt < 1:
		data = np.load(f"output/thickness_paper/shape_array_r_5_z_{z}_kP_0.01_kn_5e-06_tg_1e+06_tgi_2e+06_M_1e-07_dt_0.01_expV0.1.npz")
	else:
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

	ax.plot(time_array, 100*(chi-chi[0]), c=cmap(cnt/(N-1)),label=f"${z_list[cnt]:.1g}$",lw=1)

# ballooning transition
r0 = 5
h0 = 1
tau_g = 1e6
tau_g_i = 2*tau_g

L0 = 2*np.pi*(r0+h0)
A0 = np.pi*(r0+h0)**2
A0I = np.pi*(r0)**2
A0O = A0 - A0I
P0 = L0

Agrowth = A0I*np.exp(time_array*tau_g/tau_g_i) + A0O*np.exp(time_array)
Pgrowth = P0*np.exp(time_array)

Pmin_growth = 2*np.pi*np.sqrt(Agrowth/np.pi)
Pstretch = (Pmin_growth - P0)/P0

chi0 = P0/(2*np.sqrt(np.pi*A0))  # = 1 if initially circular

chi_upper_model = chi0*(np.exp(time_array)/np.sqrt(Agrowth/A0) - 1)

ax.fill_between(time_array,100*chi_upper_model,150*np.ones(len(time_array)),color="lightgray",zorder=0)

######################################################################
# Final layout settings
######################################################################

# Final adjustments to the figure
plotlib.set_box(ax, alt_style=True)

ax.set_xlabel(r'$ \tau / \tau_g$')
ax.set_ylabel(r'$\Delta \chi_{2D} [\%]$')

ax.set_xlim([0,1.5])
ax.set_ylim([0,150])

ax.xaxis.set_major_locator(MultipleLocator(0.5))
ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.yaxis.set_major_locator(MultipleLocator(50))
ax.yaxis.set_minor_locator(MultipleLocator(10))

plotlib.set_position(ax,x=0.17,y=0.27,width=0.78,height=0.67)

# Save as PDF
plt.savefig(f"results/acircularity_area_evolution_exp.pdf", format="pdf", bbox_inches="tight",transparent=True)

plt.show()