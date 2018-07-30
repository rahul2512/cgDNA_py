import numpy as np
import scipy
from scipy import io
import sys
import os
import matplotlib
import matplotlib.pyplot as plt
path = os.getcwd()
sys.path.insert(0,os.getcwd())
from vector2shapes import vec2shapes
from nondim2cur import nondim2cur
import matplotlib.pylab as pylab


"""
--------------------------------------------------------------------------
 Output[2 plots one for each intra and inter bp coordinates] = cgDNA2dplot('seq1','seq2',current_dir,'degrees')
--------------------------------------------------------------------------
 This function plots the ground-state coordinates to the screen.

 Input:
1. Label of the input sequence(s).
2. Name of working dircetory.
3. Unit of Y-axis i.e degrees or radians

 Output: 2 X 6 panels    plot of each intra- and interbasepair coordinate along the molecule
"""


def cgDNA2dplot(*arg):

	n_arg = len(arg)

	#### This is just the axes label ... tex code for Angstrom symbol

	t_unit = ' ($\AA$)' 
	if str(arg[n_arg - 1]) == 'degrees':
		r_unit = ' ($^\circ$)'
	if str(arg[n_arg - 1]) == 'radians':
		r_unit = ' (rad/5)'

	#### Loop to plot all the input sequences

	for k in range(n_arg - 2):

		ss = scipy.io.mmread(path + '/../' + arg[n_arg-2] + '/shape_' + arg[k] + '.mtx')
		nbp, intra_r ,intra_t ,inter_r ,inter_t = vec2shapes(ss)

		#### The if statements to get the plot as per input units ####

		if str(arg[n_arg - 1]) == 'degrees':
			intra_r, intra_t, inter_r, inter_t = (36/np.pi)*intra_r, intra_t, (36/np.pi)*inter_r, inter_t      

		#### transpose to make it compatible for plotting ### 

		intra_r , intra_t = np.transpose(intra_r), np.transpose(intra_t)
		inter_r , inter_t = np.transpose(inter_r), np.transpose(inter_t)

		### X-axis #### x_intra+1 as indexing starts from 0 in python ###

		x_intra, x_inter = 1 + np.array(range(nbp)), 1 + np.array(range(nbp-1))

		#### Subplots ####

		axes = ["ax1", "ax2", "ax3", "ax4", "ax5", "ax6"]
		position = [321, 322, 323, 324, 325, 326]		
		coords_intra = [intra_r[0], intra_t[0], intra_r[1], intra_t[1], intra_r[2], intra_t[2]]
		ylabel_fig1 = ["buckle" + r_unit, "shear" + t_unit, "propeller" + r_unit,
				"stretch" + t_unit, "opening" + r_unit, "stagger" + t_unit]
		coords_inter = [inter_r[0], inter_t[0], inter_r[1], inter_t[1], inter_r[2], inter_t[2]]
		ylabel_fig2 = ["tilt" + r_unit, "shift" + t_unit, "roll" + r_unit,
				"slide" + t_unit, "twist" + r_unit, "rise" + t_unit]

		#### Two different plots for inter and intra coordinates ####

		fig1 = plt.figure(1)

		count = 0
		for ax, pos, y_intra, ylabel in zip(axes, position, coords_intra, ylabel_fig1):
			count = count + 1
			ax = plt.subplot(pos)
			plt.plot(x_intra, y_intra ,label = arg[k])
			ax.margins(x=0)
			plt.legend(loc = 'best', fontsize = 'x-small',frameon = False)
			if count == 5 or count == 6:
				ax.set(xlabel = 'basepair',ylabel = ylabel)
			else:
				ax.set(ylabel = ylabel)
			if count == 1:
				plt.title('Intra-rotational')
			if count == 2:
				plt.title('Intra-translational')

		fig2=plt.figure(2)

		count = 0
		for ax, pos, y_intra, ylabel in zip(axes, position, coords_inter, ylabel_fig2):
			count = count + 1
			ax = plt.subplot(pos)
			plt.plot(x_inter, y_intra, label = arg[k])
			ax.margins(x = 0)
			plt.legend(loc = 'best', fontsize = 'x-small',frameon = False)
			if count == 5 or count == 6:
				ax.set(xlabel = 'basepair',ylabel = ylabel)
			else:
				ax.set(ylabel = ylabel)
			if count == 1:
				plt.title('Inter-rotational')
			if count == 2:
				plt.title('Inter-translational')
			
	### differntly optimized for python major version 2 and 3 #########
	
	if sys.version_info[0] < 3: 
		fig1.tight_layout(h_pad = 0, w_pad = 0)
		fig1.subplots_adjust(left = 0.08, right = 0.99, top = 0.95, bottom = 0.06, wspace = 0.22)
		fig2.tight_layout(h_pad = 0,w_pad = 0)
		fig2.subplots_adjust(left = 0.08, right = 0.99, top = 0.95, bottom = 0.06, wspace = 0.22)

	if sys.version_info[0] > 2:
		fig1.tight_layout(h_pad = 0, w_pad = 0)
		fig1.subplots_adjust(left = 0.11, right = 0.98, top = 0.94, bottom = 0.1, wspace = 0.3)		
		fig2.tight_layout(h_pad = 0, w_pad = 0)
		fig2.subplots_adjust(left = 0.11, right = 0.98, top = 0.94, bottom = 0.1, wspace = 0.3)

	params = {
		'legend.fontsize': 'large', 
		'figure.figsize': (15, 5),
		'axes.labelsize': 'large',
		'axes.titlesize':'large',
		'xtick.labelsize':'large',
		'ytick.labelsize':'large'
		}

	pylab.rcParams.update(params)
	plt.show()

	#### For publication quality plot, uncomment following ####

	#fig1.savefig(path+'/../'+arg[l-2]+'/intra.png',dpi = 1200)
	#fig2.savefig(path+'/../'+arg[l-2]+'/inter.png',dpi = 1200)

	return
