import numpy as np
import sys


"""
This function expands the compact sequence given as input. 
For example:    A_5G_2 = AAAAAGG
	        [A_2]_3 = AAAAAA
		[A_2G]_3[T_2]_2 = AAGAAGAAGTTTT
Note1: The code first reads the [] and then non-bracket. In case of multiple brackets, it will read the 
       outermost bracket first ..... to .... innnermost bracket in the last.  

For example:  [A_2G]_3[T_2]_2 = A_2GA_2GA_2GT_2T_2 = AAGAAGAAGTTTT
              [[[AT]_2g]_2C]_2 = [[AT]_2g]_2C[[AT]_2g]_2C = [AT]_2g[AT]_2gC[AT]_2g[AT]_2gC = 
				= ATATgATATgCATATgATATgC
Note2: Two or more _ can't be put consequently. Instead, You can use the brackets. 
	i.e. A_2_2 can be written as [A_2]_2
"""


def finder(seq):
	istart = []  
	end = {}
	start = []
	for i, c in enumerate(seq):
		if c == '[':
			istart.append(i)
			start.append(i)
		if c == ']':
			try:
				end[istart.pop()] = i
			except IndexError:
				print('Too many closing parentheses')
	if istart:  # check if stack is empty afterwards
		print('Too many opening parentheses')
	return end, start

def mult(seq):
	i =seq.rfind('_') 
	if seq[i+1].isdigit():
		a = seq[i+1]
		if seq[i+2].isdigit():
			a = a + seq[i+2]
			if seq[i+3].isdigit():
				a = a + seq[i+3]
				if seq[i+4].isdigit():
					a = a + seq[i+4]
					if seq[i+5].isdigit():
						a = a + seq[i+5]
	return a
def seq_edit(seq):
	s = seq.upper()
	while s.rfind('_')>0:
		if s[s.rfind('_')-1].isdigit():
			print("Please write the input sequence correctly. Two or more _ can't be put consequently. You can use the brackets. i.e. A_2_2 can be written as [A_2]_2")
			exit()
		if s[s.rfind('_')-1] != ']':
			a = int(mult(s))
			s = s[:s.rfind('_')-1]+ s[s.rfind('_')-1]*a +  s[s.rfind('_')+1+len(str((a))):]
		if s[s.rfind('_')-1] == ']':
			end,start = finder(s)
			ka=(2,len(start))
			h=np.zeros(ka)
			for i in range(len(start)):
				h[0][i] = start[i]
				h[1][i] = end[start[i]]	
			ss=  int(max(h[1]))
			ee=  int(h[0][np.argmax(h[1])])
			a = int(mult(s))
			s =  s[0:ee] + s[ee+1:ss]*a + s[ss+2+len(str((a))):] 
	return s	
