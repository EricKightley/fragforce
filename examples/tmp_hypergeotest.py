import numpy as np
import fragforce as frag

# compute a hypergeometric integral
aV = np.array([3., 2., 1.])
a1,a2,a3 = aV
bV = np.array([3/2., 3/2., 1/2.])
b1,b2,b3 = bV
pV = np.ascontiguousarray(np.array([a1**2,a2**2,a3**2, b1,b2,b3]))

import time

N = 150

s1 = time.time()
for i in range(N):
    hyperint = frag.py_integrate_hypergeo3(pV)
e1 = time.time()

s2 = time.time()
for i in range(N):
    chi = frag.py_set_chi(np.array([a1,a2,a3]))
    X1,X2,X3 = chi
    result = (X2-X1) / (a1**2-a2**2)
e2 = time.time()

#print(hyperint, result, np.abs(hyperint - result) / hyperint)

sc = 10000 / 60.
print(sc*(e1-s1), sc*(e2-s2))



