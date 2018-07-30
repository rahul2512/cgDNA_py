import numpy as np
import scipy
from scipy import io
import os
import sys
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
path = os.getcwd()
sys.path.insert(0,os.getcwd())
from frames import frames


"""
--------------------------------------------------------------------------
  Output[plot] = plot3D('label1','label2'....,'labelN',working_dir,'Axes options') # if multiple sequences
	       = or  plot3D('label1',working_dir)  # if one input sequence
--------------------------------------------------------------------------
 This function plots the 3D reconstruction of the ground states as rigid
 bodies colored according to the sequence and the Crick Watson pairing.

 Input:  "Label"        structure array as output by the function cons_shape_stiff
	 Other Input options: Name of the Working directory, Axes on and off option

 Output:   Figure 1    plot of rigid bodies displayed as patches. 
	   Color scheme: red A, blue T, green G, and yellow C. 

   Note 1: If MyData contains M sequences, then M structures will be
           displayed
"""


def getRigidBodyConfig(*arg):

	##### sample input: getRigidBodyConfig(rrw, RRw, base) for original strand
	##### sample input: getRigidBodyConfig(rrc, RRc, baseC, 'C') for complementary strand, the last argument has
	##### no meaning, it is just to differentiate between the original and complementary strand.  

	base=arg[2]

	dz, dx, Y = 0.36, 4.20, -1.2   # this function generates the 3D frame for the bases

	if base == 'A' or base == 'G':
		dy = 6
	if base == 'T' or base == 'C':
		dy = 4.2	
	if len(arg) == 4:		
		dy, Y = -dy, -Y

	verteces_RB = [ [-0.5*dx, -0.5*Y      , -0.5*dz],
			[ 0.5*dx, -0.5*Y      , -0.5*dz],
			[ 0.5*dx,  dy - 0.5*Y , -0.5*dz],
			[-0.5*dx,  dy - 0.5*Y , -0.5*dz],
			[-0.5*dx, -0.5*Y      ,  0.5*dz],
			[ 0.5*dx, -0.5*Y      ,  0.5*dz],
			[ 0.5*dx,  dy - 0.5*Y ,  0.5*dz],
			[-0.5*dx,  dy - 0.5*Y ,  0.5*dz]]

	if base == 'A':
		color = "r"
	if base == 'T':
		color = "b"
	if base == 'C':
		color = "y"
	if base == 'G':
		color = "g"

	RB_conf = np.add(np.transpose(np.matmul(arg[1], np.transpose(verteces_RB))), np.tile(np.transpose(arg[0]), [8, 1]))

	return color, RB_conf


def getRigidBodies(label,working_dir):   

	sseq = scipy.io.mmread(path + '/../' + working_dir + '/shape_' + label + '.mtx')  
	nbp, R, r, Rc, rc, Rw, rw = frames(sseq)
	color, colorC, RB, RBC = {}, {}, {}, {} 

	for i in range(int(nbp)):
		rrw, rrc, RRw, RRc = rw[i], rc[i], Rw[i], Rc[i]
		fuf = open(path + '/../' + working_dir + '/complete_' + label + '.txt', 'r')
		sus = fuf.read()
		fuf.close()
		base = sus[i]

		if base == 'A':
			baseC = 'T'
		if base == 'T':
			baseC = 'A'
		if base == 'C':
			baseC = 'G'
		if base == 'G':
			baseC = 'C'

		#### The variable with C is complementary strand,c i.e. color and colorC ####
	
		color[i], RB[i] = getRigidBodyConfig(rrw, RRw, base) # call the above function for both the watson 
		colorC[i], RBC[i] = getRigidBodyConfig(rrc, RRc, baseC, 'C') # and crick strands

	return nbp, RBC, colorC, RB, color


def plot3D(*arg):  
	
	#this function plots the 3D structure of the input sequences. 
	#sample input: plot3D('seq1','seq2', current_dir, 'off')

	n_arg = len(arg)
	working_dir = arg[n_arg - 2]
	axes_opt = arg[n_arg - 1] # Show axis on or off, provided in the input
	fig = plt.figure()
	ax = Axes3D(fig)
	nbp = np.zeros(n_arg - 2)

	for ii in range(n_arg-2):
		label = arg[ii]
		nbp[ii], RBC, colorC, RB, color = getRigidBodies(label, working_dir)


		#### The following loop will create rectangular patches which will represent the faces for base ####
		#### The two loops are for watson and crick strands ####

		for i in range(len(colorC)):
			k, c, kc, cc = RB[i], color[i], RBC[i], colorC[i]

			x1, y1, z1 = [k[0][0], k[1][0], k[2][0], k[3][0]], [k[0][1], k[1][1], k[2][1], k[3][1]], [k[0][2], k[1][2], k[2][2], k[3][2]]
			x2, y2, z2 = [k[4][0], k[5][0], k[6][0], k[7][0]], [k[4][1], k[5][1], k[6][1], k[7][1]], [k[4][2], k[5][2], k[6][2], k[7][2]]
			f1, f2 = [list(zip(x1, y1, z1))], [list(zip(x2, y2, z2))]

			for i in f1, f2:
				alpha = 0.5
				fc = c
				pc = Poly3DCollection(i, alpha = alpha, facecolors = fc, linewidths = 1)
				ax.add_collection3d(pc)

			x1, y1, z1 = [kc[0][0], kc[1][0], kc[2][0], kc[3][0]], [kc[0][1], kc[1][1], kc[2][1], kc[3][1]], [kc[0][2], kc[1][2], kc[2][2], kc[3][2]]
			x2, y2, z2 = [kc[4][0], kc[5][0], kc[6][0], kc[7][0]], [kc[4][1], kc[5][1], kc[6][1], kc[7][1]], [kc[4][2], kc[5][2], kc[6][2], kc[7][2]]

			f1, f2 = [list(zip(x1, y1, z1))], [list(zip(x2, y2, z2))]
			for i in f1,f2:
				alpha = 0.5
				fc = cc
				pc = Poly3DCollection(i, alpha = alpha, facecolors = fc, linewidths = 1)
				ax.add_collection3d(pc)
	bb = np.amax(nbp)
	ax.set_xlim(-1.75*bb, 1.75*bb)  # it optimizes the 3D plot view. 
	ax.set_ylim(-1.75*bb, 1.75*bb)
	ax.set_zlim(0, 3.5*bb)
	plt.axis(axes_opt) # ON or OFF
#	plt.savefig(path +'/'+ '3D.png', format = 'png', dpi = 1200)
	plt.show()
