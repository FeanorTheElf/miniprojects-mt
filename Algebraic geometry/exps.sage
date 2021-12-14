from itertools import combinations
from functools import reduce
from math import factorial
import numpy as np

n = 8

codes = [str(i + 1) for i in range(n)]

# build up the ring and group the variables nicely
R = PolynomialRing(QQ, [
    "x" + codes[i] + codes[j]
        for i in range(n) for j in range(i + 1, n)
])

_x = []
vars = list(R.gens())
for i in range(n - 1, 0, -1):
    _x.append(vars[:i])
    vars = vars[i:]
x = lambda i, j: _x[i - 1][j - i - 1]

def gen_f(a, b, c, d):
    (i, j, u, v) = sorted((a, b, c, d))
    return x(i,j) * x(u,v) + x(i,v) * x(j,u) - x(i,u) * x(j,v)

# construct the ideal describing Gr(2, 6)
polys = []
for seq in combinations(range(1, n + 1), 4):
    polys.append(gen_f(*seq))

I = R.ideal(polys)
for f in I.groebner_basis():
    print(f)