import numpy as np
import os
import sys
import scipy.io
from scipy.linalg import sqrtm, inv
path = os.getcwd()
sys.path.insert(0,os.getcwd())
from vector2shapes import vec2shapes



"""
--------------------------------------------------------------------------
 cgDNA function: Output [nbp, R, r, Rc, rc, Rw, rw] = frames(shapes)
--------------------------------------------------------------------------
 Given a ground-state coordinate vector in non-dimensional Curves+ form, this function constructs a 
reference point and frame for each base on each strand of the DNA according to the Tsukuba convention, and 
stores the result in the basepair structure (see Note 1). The reference point and frame vectors for each 
base are expressed relative to a fixed lab frame, which is assumed to coincide with the first basepair frame.

 Input: 
   shapes   ground-state coordinate vector 
            in non-dimensional Curves+ form (see Note 2)
            [size N x 1].
 Output:
   nbp, R, r, Rc, rc, Rw, rw       structure with reference point and frame
                                   for each base on each strand (see Note 1).

 Note 1: All the following outputs are (1 x nbp) struct arrays:
   - 'R' : the frame of the basepair;
   - 'r' : the coordinates of the basepair;
   - 'Rw' : the frame of the base on the reading strand;
   - 'rw' : the coordinates of the base on the r. s.;
   - 'Rc': the frame of the base on the complementary strand;
   - 'rc': the coordinates of the base on the c. s.;

   Reference point coordinates are 3x1 vectors, while frames 
   are 3x3 matrices, with the frame coordinate vectors stored
   as columns.  'nbp' is the length of the sequence.

Note 2: The entries of the input variable shapes must be consistent with the following ordering of the 
        structural coordinates

    y_1, z_1, ..., y_{nbp-1}, z_{nbp-1}, y_{nbp}
   where for each a=1,2,3,... we have
    y_a = (Buckle,Propeller,Opening,Shear,Stretch,Stagger)_a
    z_a = (Tilt,Roll,Twist,Shift,Slide,Rise)_a.

   For example

    shapes((a-1)*12+1) = Buckle at basepair a
    shapes((a-1)*12+2) = Propeller at basepair a
     ...
    shapes((a-1)*12+6) = Stagger at basepair a
    shapes((a-1)*12+7) = Tilt at junction a, a+1
    shapes((a-1)*12+8) = Roll at junction a, a+1
     ...
    shapes((a-1)*12+12) = Rise at junction a, a+1.

If you find this code useful, please cite:
"""


def cay(k1):

	I=np.identity(3)
	k = np.divide(k1,10)
	X = [[0,-k[2],k[1]],[k[2],0,-k[0]],[-k[1],k[0],0]]
	Q = np.matmul(inv(np.subtract(I,X)),np.add(I,X))
	return Q

def frames(s):
	
	#relative coordinates of the oligomer
	(nbp,intra_r,intra_t,inter_r,inter_t) = vec2shapes(s)
	
	# absolute coordinates for the first bp
	G = np.identity(3)   
	q = np.transpose(np.zeros(3))

	R, r, Rc, rc, Rw, rw = {}, {}, {}, {}, {}, {}

	for i in range(nbp):
		R[i] = G #merging everything at origin
		r[i] = q #merging everything at origin
		L = cay(intra_r[i]) #basepair
		Gw = np.matmul(G,(np.transpose(intra_t[i]))) #basepair
		Rc[i] = np.matmul(G,np.transpose(np.real(sqrtm(L))))   #complementray strand
		rc[i] = np.subtract(q,np.multiply(0.5,Gw))   #complementray strand
		Rw[i] = np.matmul(Rc[i],L) 	#original strand
		rw[i] = np.add(rc[i],Gw)        #original strand

		if i < nbp-1:
			ru = cay(inter_r[i])
			H = np.matmul(G,np.real(sqrtm(ru)))	

			################## Next base pair #################
			G = np.matmul(G,ru)
			q = np.add(q,np.matmul(H,np.transpose(inter_t[i])))

	return nbp, R, r, Rc, rc, Rw, rw
