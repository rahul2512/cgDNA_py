import numpy as np
import os
import sys
import math
path = os.getcwd()
sys.path.insert(0,os.getcwd())
from vector2shapes import vec2shapes
from shapes2vec import shapes2vec


"""
-------------------------------------------------------
 cgDNA function: [curshapes] = nondim2cur(shapes)
-------------------------------------------------------
 This function transforms the ground-state coordinates from non-dimensional 
Curves+ form to the standard (dimensional) Curves+ form.

 Input: 

   shapes  ground-state coordinate vector 
           in non-dimensional Curves+ form
           [size N x 1].

 Output:

   curshapes  ground-state coordinate vector 
              in standard Curves+ form
              [size N x 1]

   where N = 12*nbp - 6 and nbp is the length 
   of the DNA sequence (number of basepairs).

If you find this code useful, please cite:

D. Petkeviciute, M. Pasi, O. Gonzalez and J.H. Maddocks. 
 cgDNA: a software package for the prediction of 
 sequence-dependent coarse-grain free energies of B-form 
 DNA. Nucleic Acids Research 2014; doi: 10.1093/nar/gku825.
"""



def nondim2cur(shape):

	(nbp,intra_r,intra_t,inter_r,inter_t) = vec2shapes(shape)
	cay = intra_r/10
	ncay = np.transpose(np.tile(np.sqrt(np.square(cay).sum(axis=1)),(3,1)))
	angle = 2*np.arctan(ncay)/np.pi*180
	etac = np.multiply(angle,np.divide(cay,ncay))	
	cay = inter_r/10
	ncay = np.transpose(np.tile(np.sqrt(np.square(cay).sum(axis=1)),(3,1)))
	angle = 2*np.arctan(ncay)/np.pi*180
	uc = np.multiply(angle,np.divide(cay,ncay))
	curshapes = shapes2vec(etac, intra_t, uc, inter_t)

	return curshapes
