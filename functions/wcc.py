import sys
"""
This function take the sequence input and output the complementary sequence in 5' to 3' fashion.
output[complementary] = wcc[input_sequence]
Example: 
wcc('TAGCTACTCAGACGACCAGTAGATAAAAG') will give "CTTTTATCTACTGGTCGTCTGAGTAGCTA"
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
def wcc(s):
	c = []
	for i in range(len(s.strip())):
		c.append(comp(s[i]))
	sys.stdout.write(" Note that the output sequence is from 5' to 3' \n")
	c = c[::-1]   #invert the sequence i.e. make it 5' to 3'
	for element in c:
		sys.stdout.write(str(element))	
	sys.stdout.write("\n")
