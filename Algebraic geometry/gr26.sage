from itertools import combinations
from math import factorial
import numpy as np

# build up the ring and group the variables nicely
R = PolynomialRing(QQ, [
    "x" + str(i) + str(j) 
        for i in range(1, 7) for j in range(i + 1, 7)
])
x12, x13, x14, x15, x16, x23, x24, x25, x26, x34, x35, x36, x45, x46, x56 = R.gens()
_x = [[x12, x13, x14, x15, x16], [x23, x24, x25, x26], [x34, x35, x36], [x45, x46], [x56]]
x = lambda i, j: _x[i - 1][j - i - 1]

# construct the ideal describing Gr(2, 6)
polys = []
for seq in combinations([1, 2, 3, 4, 5, 6], 4):
    (i, j, u, v) = sorted(seq)
    p = x(i,j) * x(u,v) + x(i,v) * x(j,u) - x(i,u) * x(j,v)
    polys.append(p)

I = R.ideal(polys)
print(I.hilbert_polynomial()) 
# 1/2880*t^8 + 1/120*t^7 + 41/480*t^6 + 39/80*t^5 + 541/320*t^4 + 291/80*t^3 + 3401/720*t^2 + 101/30*t + 1