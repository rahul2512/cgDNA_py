import sys
import os
os.getcwd()
sys.path.insert(0,os.getcwd()+ '/../functions')
from constructSeqParms import cons_shape_stiff
from cgDNA2dplot import cgDNA2dplot

from cgDNA_MonteCarlo import cgDNA_MonteCarlo
from cgDNA3dplot import plot3D
from makepdb import makepdb
#----------------------------------------------------------------------------------------------------------
#Note1: the Parameterset directory, functions directory and the current working directory should be at the 
#same level i.e. in the same directory: for example it's in cgDNA_py
#----------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------
#Note2: The sequence file should be named as label.txt and should be located in the same directory as input.py
#if you want to repeat any oligomer --- type [sequence]_N : it will repeat that oligomer N times
#example: A_5 = AAAAA, A_2G_3= AAGGG, A_2[GG]_2=AAGGGG, A_2[T_2]_2G_2=AATTTTGG  
#Note that the code first reads the brackets and then underscore i.e  A_2[T_2]_2G_2 = A_2T_2T_2G_2 = AATTTTGG
#For more details, see the original code,seq_edit.py in functions dircetory.  
#-----------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------
#Following changes need to be done in this input file

#  1. Need to change it everytime you work in different directory name: for example the current directory 

#  2. cons_shape_stiff ('label of the sequence file',paramaters set number,'name of the current directory')
#  The function cons_shape_stiff will predict the groundstate and stiffness matrix and save in the working 
#  directory with the name shape_label.mtx and stiff_label.mtx respectively. 
#  You can compute multiple sequences at the sam time. Just provide the names of label. 

#### Now in all the functions below, one must specify the parameterset used, this would also faciliate comparing different parametersets. ###

#  3. cgDNA2dplot('label1'---'labelN',current_dir,'degrees or radians')  
#  It also allows to plot the graph in both degrees and radians (cgDNA coordinates)

#  4. plot3D('label1',---,'labelN',current_dir,'Axes_option: ON or OFF')
#  This function will plot the groundstate of the provided sequence in 3D plot. 
#  Axes option will allow to turn ON/OFF the axes on the 3D graph. 

#  5. cgDNA_MonteCarlo('label1',---,'labelN',nbr_samples,current_dir)  
#  Here, we recommend to use a nbr_samples = 500 for a sequence of around 300 bp. For higher number of samples, 
#  please use much more efficient C++ code: https://lcvmwww.epfl.ch/research/cgDNA/downloads.php 

#  6. makepdb('label1'---'labelN',current_dir) 
#  This function will save a pdb file in the current_dir named labelN.pdb




current_dir = 'Examples'                              #1 .. change to be made if requires



cons_shape_stiff ('seq1','seq2','seq3',1,current_dir)               #2 .. changes to be made if requires
cons_shape_stiff ('seq1','seq2','seq3',2,current_dir)
cons_shape_stiff ('seq1','seq2','seq3',3,current_dir)
cons_shape_stiff ('seq1','seq2','seq3',4,current_dir)

cgDNA2dplot('seq1_ps3','seq1_ps2',current_dir,'degrees')             #3 .. degrees or radians
cgDNA2dplot('seq1_ps4','seq2_ps4',current_dir,'degrees')
cgDNA2dplot('seq1_ps4','seq2_ps4','seq3_ps4',current_dir,'degrees')
cgDNA2dplot('seq3_ps4',current_dir,'radians')            

plot3D('seq1_ps4','seq1_ps2',current_dir,'off')       #4 axis option = on/off
plot3D('seq3_ps4','seq3_ps2',current_dir,'off')
plot3D('seq1_ps4','seq2_ps4','seq3_ps4',current_dir,'off')

cgDNA_MonteCarlo('seq1_ps1','seq1_ps4',100,current_dir)        #5 .. For long sequence, recommended nbr_samples = 500
makepdb('seq1_ps2','seq2_ps4',current_dir)                           #6 .. It will create seq1.pdb
