import scipy.io
import numpy as np
ps = scipy.io.loadmat('cgDNAps4.mat',mat_dtype=1)
k = ps['stiff_end5']['AA']
print(np.ma.shape(k))
f = [1,2]
a = [4, 3, 5, 7, 6, 8]
print(k[0][0][0][0])
