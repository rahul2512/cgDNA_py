import numpy as np
import numpy
import numpy.linalg
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import scipy
from scipy import io
import scipy.linalg
from scipy.sparse.linalg import inv
from scipy.stats import multivariate_normal
import os
import os.path
from scipy import io
import sys
sys.path.insert(0,os.getcwd())
from frames import frames
import time
path = os.getcwd()


"""
Output[l_p,l_d] = cgDNA_MonteCarlo('label',Number of samples,working_dir)

This function draws a given number of configurations from the cgDNA probability  density 
function for the sequence(s) in the given Data structure.

Apparent and dynamic (see ref below for definitions) persistence lengths are computed from the 
sampling for each sequence, and are added as fields in the structure Data.

This is a python version of the much more efficient cgDNAmc code 
http://lcvmwww.epfl.ch/software/cgDNAmc/doc/index.html
Jonathan S. Mitchell, Jaroslaw Glowacki, Alexandre E. Grandchamp, Robert S. Manning, and John H. Maddocks, 
Sequence-dependent persistence lengths of DNA J. Chem. Theory Comput.

If only one sequence is provided as input, this function also plots all the base-pairs 
positions of all sampled configuration in blue and the ground-state in red.

Input : - label           Data structure as output by the cgDNA.m function
        - nbr_samples    Number of wanted sampled configuration. We advise to not use above 500 samples 
			 or 500 bp long sequence. 
	- Working directory  The location of shape, sequence and stiffness matrix
Output : - Data          apparent and dynamic persistence length for all the input sequences
			 plot of all all the base-pairs positions of all sampled configuration 
			 in blue and the ground-state in red (only if one sequence is provided)
"""


def Compute_pl(nbp,dot_d3,dot_d3_intr):  #this function computes the persistence length

	x = np.zeros((nbp,1))  

	for i in range(nbp):
		x[i][0] = i+1   

	y = np.log(dot_d3)

	######## least square fit of the log(tangents) 
	tmp_pl = np.divide(numpy.matmul(np.transpose(x),y),numpy.matmul(np.transpose(x),x)) 
	pl = -1/tmp_pl
	y = np.subtract(np.log(dot_d3), np.log(dot_d3_intr))
	tmp_pl = np.divide(numpy.matmul(np.transpose(x),y),numpy.matmul(np.transpose(x),x))
	pl_fact = -1/tmp_pl
	return pl, pl_fact, tmp_pl

def cgDNA_MonteCarlo(*arg):

	#####Sample input: cgDNA_MonteCarlo(seq_label(s),nbr_samples,working_dir)

	n_arg = len(arg)
	nbr_samples = int(arg[n_arg-2])
	l_p, l_d = {}, {}

	for label in range(n_arg-2):

		fuf = open(path+'/../'+arg[n_arg-1]+'/complete_'+arg[label]+'.txt', 'r')
		sus = fuf.read()
		fuf.close()

		shape=scipy.io.mmread(path+'/../'+arg[n_arg-1]+'/shape_'+arg[label]+'.mtx')
		stiff=scipy.io.mmread(path+'/../'+arg[n_arg-1]+'/stiff_'+arg[label]+'.mtx')

		nbp=len(sus.strip())
		dim=12*nbp-6

		MC = np.random.multivariate_normal(np.zeros(dim),np.eye(dim,dtype=int),nbr_samples)
		ttc = np.zeros((nbp,1))

		L = scipy.linalg.cholesky(stiff.todense())

		#### if only one seq is provided, it will plot the different MC configurations####
		if n_arg == 3:
			fig = plt.figure()
			ax = plt.axes(projection='3d')

		for jj in range(nbr_samples):
			x = np.linalg.solve(L,np.transpose(MC[jj])) + shape
			nbp, R, r, Rc, rc, Rw, rw = frames(x)
			XYZ = np.zeros((nbp,3))

			for k in range(nbp): 
				XYZ[k] = r[k]
				ttc[k][0] =  np.divide(R[k][2][2],nbr_samples) + ttc[k][0] 

			if n_arg == 3:	
				XYZ = np.transpose(XYZ)
				ax.plot3D(XYZ[0], XYZ[1], XYZ[2], 'blue')
				XYZ = np.transpose(XYZ)

		nbp, R, r, Rc, rc, Rw, rw = frames(shape)
		ttc_int = np.zeros((nbp,1))

		for k in range(nbp):
			XYZ[k] = r[k]
			ttc_int[k][0] = R[k][2][2]

		if n_arg == 3:
			XYZ = np.transpose(XYZ)
			ax.plot3D(XYZ[0], XYZ[1], XYZ[2], 'red')

		l_p[label],l_d[label],not_imp = Compute_pl(nbp,ttc,ttc_int) 

		print("For ",arg[label], " l_p = ", l_p[label]," and l_d = ", l_d[label])

	if n_arg == 3:
		plt.show()
		plt.savefig(path+'Mon.png', format='png', dpi=1200)
