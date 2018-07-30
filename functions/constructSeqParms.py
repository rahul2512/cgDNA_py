import numpy as np
import scipy
import scipy.io
from scipy import sparse
from scipy.sparse import csc_matrix, lil_matrix
from scipy.sparse.linalg import spsolve
import os
import os.path
path = os.getcwd()
from seq_edit import seq_edit


"""
####Important note: indexing in python starts from 0   #### 

###   INPUT   ####

This function takes input 'sequence label',parameterset number,working directory name 
If you want to do multiple sequence: label1,label2,....,labelN,parameterset number,working directory name

###  OUTPUT  #####

 1. Expanded sequence = saved in the working directory with name complete_label.txt
 2. Ground state = saved in the working directory with name shape_label.txt (N X 1 matrix)
 3. Stiffness matrix = saved in the working directory with name stiff_label.txt (N X N matrix)
   where N =  12*nbp - 6 and nbp is the length of the sequence seq (number of basepairs)

###  Parameterset format ####  Parameterset4 is located at cgDNA_py/Parametersets/ps4 ###
  - stiff_end3: 18x18 stiffness 3' end blocks,
  - stiff_end5: 18x18 stiffness 5' end blocks,
  - stiff_int: 18x18 stiffness interior blocks,
  - sigma_end3: 18x1 weighted shape 3' end vector blocks,
  - sigma_end5: 18x1 weighted shape 5' end vector blocks,
  - sigma_int: 18x1 weighted shape interior vector blocks,
  All the above vector are avaiable for 36 dinucleotides 
  Example : sigma_end3.AT is a 18x1 vector

#### The entries in shapes and stiff are consistent with the following ordering of the structural coordinates

      y_1, z_1, ..., y_{nbp-1}, z_{nbp-1}, y_{nbp}
      where for each a=1,2,3,... we have 
      y_a = (Buckle,Propeller,Opening,Shear,Stretch,Stagger)_a
      z_a = (Tilt,Roll,Twist,Shift,Slide,Rise)_a.

For example: 

 shapes((a-1)*12+1) = Buckle at basepair a
 shapes((a-1)*12+2) = Propeller at basepair a
  ...
 shapes((a-1)*12+6) = Stagger at basepair a
 shapes((a-1)*12+7) = Tilt at junction a, a+1
 shapes((a-1)*12+8) = Roll at junction a, a+1
  ...
 shapes((a-1)*12+12) = Rise at junction a, a+1.

Correspondingly, we have

 stiff(i,j) = stiffness coefficient for the pair of
              coordinates shapes(i) and shapes(j).

If you find this code useful, please cite:
"""

#### *arg means a function with any number of arguments  ####
####Important note: indexing in python starts from 0 unlike MATLAB which starts from 1 ####


def cons_shape_stiff(*arg):

	n_arg = len(arg)
	ps = scipy.io.loadmat(path + '/..' + '/Parametersets/cgDNAps' + str(arg[n_arg-2]) + '.mat')

	#### Following loop take every input sequence and construct shape and stiff matrix ###

	for label in range(n_arg-2):

		###Read the input sequence ###
		f = open(arg[label] + '.txt', 'r') 
		seq_inp = f.read()
		f.close()

		#### seq_edit function expand the sequence from shorter format like A_5 to AAAAA ####
		s_seq = seq_edit(seq_inp)
		nbp = len(s_seq.strip())
		N = 12*nbp - 6

		#### Initialise stiff and sigma matrix ###		

		k =  lil_matrix((N, N),shape = None, dtype = 'd', copy = False)
		s = np.zeros((N,1))

		#### Error report if sequence provided is less than 2 bp #### 

		if nbp < 2:
			print("sequence length must be greater than or equal to 2")
			sys.exit() 

		#### treated nbp=2 as separate case #######
		## Note first two zeroes while calling ps are redundant, because of the way python treats .mat file ###

		if nbp == 2:
			k[0:18, 0:18] = k[0:18, 0:18] +  ps['stiff_end5'][s_seq[0:2]][0][0][0:18, 0:18]
			s[0:18] = s[0:18] + ps['sigma_end5'][s_seq[0:2]][0][0][0:18]

		#### the loops below will create the stiffness matrix k, weighted shape (s)  ########

		if nbp > 2 :
			
			#### 5' end ####
			k[0:18, 0:18] = k[0:18, 0:18] + ps['stiff_end5'][s_seq[0:2]][0][0][0:18,0:18]
			s[0:18] = s[0:18] + ps['sigma_end5'][s_seq[0:2]][0][0][0:18]

			#### interior blocks  ###
			for i in range(1, nbp-2, 1):
				k[12*i:12*i+18, 12*i:12*i+18] = k[12*i:12*i+18, 12*i:12*i+18] + ps['stiff_int'][s_seq[i:i+2]][0][0][0:18, 0:18]
				s[12*i:12*i+18] = s[12*i:12*i+18] + ps['sigma_int'][s_seq[i:i+2]][0][0][0:18]
			
			#### 3' end ####
			k[N-18:N, N-18:N] = k[N-18:N, N-18:N] + ps['stiff_end3'][s_seq[nbp-2:nbp]][0][0][0:18, 0:18]
			s[N-18:N] = s[N-18:N] + ps['sigma_end3'][s_seq[nbp-2:nbp]][0][0][0:18]

		#### Groudstate calculation ####
		shape = spsolve(csc_matrix(k), s) 

		#### The following command write the expanded sequence, stiffness, shape matrix ####

		h = open(path+'/../'+arg[n_arg - 1]+'/complete_'+arg[label] + '_ps' + str(arg[n_arg - 2])+'.txt', "w")
		h.write(s_seq)
		h.close()

		scipy.io.mmwrite(path + '/../' + arg[n_arg - 1] + '/shape_' + arg[label]+'_ps' + str(arg[n_arg - 2]), 
				np.asmatrix(shape))
		scipy.io.mmwrite(path + '/../' + arg[n_arg - 1] + '/stiff_' + arg[label] + '_ps' + str(arg[n_arg - 2]),
				 k)

		#### If you want sigma matrix, please uncomment the following command ####

#		scipy.io.mmwrite(path + '/../' + arg[n_arg - 1] + '/sigma_' + arg[label] + '_ps' + str(arg[n_arg - 2]), np.asmatrix(s))
		
		#### If you need shape, stiff or sigma in MATLAB/Octave format, uncomment follwoing ##### 

#		scipy.io.savemat(path+'/../'+arg[n_arg-1]+'/out_'+arg[label]+'_ps' + str(arg[n_arg - 2]), 
#				mdict={
#					'sigma': np.transpose(np.matrix(s)),
#					'shape': np.transpose(np.matrix(shape)),
#					'stiff': k,
#					'nbp': nbp,
#					'ps': arg[n_arg-2],
#					'sequence': seq_inp 
#					}, 
#				appendmat=True)		

	return
