import numpy as np
from scipy import io as spio

size_mat = np.random.randint(5,10,size=1)[0]
A = np.random.normal(20, 30, size=(size_mat, size_mat))
b = np.random.normal(20, 30, size=(size_mat,1))
x = np.linalg.solve(A,b).reshape(size_mat,1)

spio.mmwrite("A_test",A)
spio.mmwrite("b_test",b)
spio.mmwrite("x_test",x)

det_orig = np.linalg.det(A)

print("%.16e" % det_orig)
