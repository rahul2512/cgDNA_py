import numpy as np
import os
import sys
import scipy
from scipy import io
path = os.getcwd()
sys.path.insert(0,os.getcwd())
from vector2shapes import vec2shapes
from frames import frames


"""
--------------------------------------------------------
 cgDNA function: output[label.pdb] = makePDB(seq_label, working_dir)
--------------------------------------------------------
 Given the reference point and frame of each base,
 this function constructs the ideal coordinates of 
 the non-hydrogen atoms of each base according to the
 Tsukuba definition, and writes the output to a PDB 
 file (backbone atoms are not included).  The atomic
 coordinates are expressed relative to a fixed lab
 frame which coincides with the first basepair frame.

 Input: 

  Data             shape as output by constructSeqParms.py: It will extract the shape if label is provided. 

 Auxiliary input: 

   ideal_bases.txt  text file with the ideal coordinates (in base frame) of 
                    the non-hydrogen atoms of the four bases T, A, C, G.
Output: label.pdb will be stored in working_dir. 

Note 1:

   - 'Rw' : the frame of the base on the reading strand (Watson);
   - 'rw' : the coordinates of the base on the r. s.;
   - 'Rc': the frame of the base on the complementary strand (Crick);
   - 'rc': the coordinates of the base on the c. s.;

   Reference point coordinates are 3x1 vectors, while frames 
   are 3x3 matrices, with the frame coordinate vectors stored
   as columns.  'nbp' is the length of the sequence.

If you find this code useful, please cite:

D. Petkeviciute, M. Pasi, O. Gonzalez and J.H. Maddocks. 
 cgDNA: a software package for the prediction of 
 sequence-dependent coarse-grain free energies of B-form 
 DNA. Nucleic Acids Research 2014; doi: 10.1093/nar/gku825.
"""


def comp(base):

        if base == 'A':
                baseC = 'T'

        if base == 'G':
                baseC = 'C'

        if base == 'C':
                baseC = 'G'

        if base == 'T':
                baseC = 'A'

        return baseC


def makepdb(*arg):
	n_arg = len(arg)
	for ii in range(n_arg - 1):
		fuf = open(path + '/../' + arg[n_arg - 1] + '/complete_' + arg[ii] + '.txt', 'r')
		sus = fuf.read()
		fuf.close()
		out = open(arg[ii] + '.pdb', 'w')
		shape = scipy.io.mmread(path + '/../' + arg[n_arg - 1] + '/shape_' + arg[ii] + '.mtx')
		nbp, R, r, Rc, rc, Rw, rw = frames(shape)

		ideal = np.loadtxt(path + '/../functions/' + 'ideal_bases.txt', 
					skiprows=1,
					delimiter=',',
					usecols = (1,2,3))

		atoms = np.loadtxt(path + '/../functions/' + 'ideal_bases.txt',
					skiprows=1,
					delimiter=',',
					usecols=(0),
					dtype=str)

		kk=0
		for i in range(nbp): 
			base = sus[i]	
			cord = {}
			if base == 'A':
				for j in range(11):
					kk = kk + 1
					cord[j] = rw[i] + np.matmul(Rw[i],np.transpose(ideal[j]))
					cord[j] = rw[i] + np.matmul(Rw[i],ideal[j])
					out.write("ATOM  "+str(kk)+"  "+str(atoms[j])+"  "+base+"  "+str(i+1)+"  "+'   '.join(map(str, cord[j]))+"\n")
			if base == 'T':
				for j in range(10): 
					kk = kk+1
					cord[j] = rw[i] + np.matmul(Rw[i],np.transpose(ideal[j+11]))
					out.write("ATOM  "+str(kk)+"  "+str(atoms[j+11])+"  "+base+"  "+str(i+1)+"  "+'   '.join(map(str, cord[j]))+"\n")
			if base == 'G':
				for j in range(12):
					kk = kk+1
					cord[j] = rw[i] + np.matmul(Rw[i],np.transpose(ideal[j+11+10]))
					out.write("ATOM  "+str(kk)+"  "+str(atoms[j+11+10])+"  "+base+"  "+str(i+1)+"  "+'   '.join(map(str, cord[j]))+"\n")
			if base == 'C':
				for j in range(9):
					kk = kk+1
					cord[j] = rw[i] + np.matmul(Rw[i],np.transpose(ideal[j+11+10+12]))
					out.write("ATOM  "+str(kk)+"  "+str(atoms[j+11+10+12])+"  "+base+"  "+str(i+1)+"  "+'   '.join(map(str, cord[j]))+"\n")


		out.write('TER  \n')
		F = [[1,0,0],[0,-1,0],[0,0,-1]]
		for i in range(nbp-1,-1,-1):
			base=comp(sus[i])
			cord = {}
			if base == 'A':
				for j in range(11):
					kk = kk+1
					cord[j] = rc[i] + np.matmul(Rc[i],np.matmul(F,np.transpose(ideal[j])))
					out.write("ATOM  "+str(kk)+"  "+str(atoms[j])+"  "+base+"  "+str(2*nbp-i)+"  "+'   '.join(map(str, cord[j]))+"\n")
			if base == 'T':
				for j in range(10):
					kk = kk+1
					cord[j] = rc[i] + np.matmul(Rc[i],np.matmul(F,np.transpose(ideal[j+11])))
					out.write("ATOM  "+str(kk)+"  "+str(atoms[j+11])+"  "+base+"  "+str(2*nbp-i)+"  "+'   '.join(map(str, cord[j]))+"\n")
			if base == 'G':
				for j in range(12):
					kk = kk+1
					cord[j] = rc[i] + np.matmul(Rc[i],np.matmul(F,np.transpose(ideal[j+11+10])))
					out.write("ATOM  "+str(kk)+"  "+str(atoms[j+11+10])+"  "+base+"  "+str(2*nbp-i)+"  "+'   '.join(map(str, cord[j]))+"\n")
			if base == 'C':
				for j in range(9):
					kk = kk+1
					cord[j] = rc[i] + np.matmul(Rc[i],np.matmul(F,np.transpose(ideal[j+11+10+12])))
					out.write("ATOM  "+str(kk)+"  "+str(atoms[j+11+10+12])+"  "+base+"  "+str(2*nbp-i)+"  "+'   '.join(map(str, cord[j]))+"\n")
		out.close()
	return
