from itertools import combinations
from functools import reduce
from math import factorial
import numpy as np

n = 6

codes = ["a", "j", "u", "b", "c", "d"]

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
for seq in combinations([1, 2, 3, 4, 5, 6], 4):
    polys.append(gen_f(*seq))

I = R.ideal(polys)

a, j, u, b, c, d = 1, 2, 3, 4, 5, 6
i, v = a, d
S = x(b, c) * gen_f(i, j, u, v) - x(j, u) * gen_f(a, b, c, d)
print(S)
Sp = S + x(1, 3) * gen_f(2, 4, 5, 6)
print(Sp)
Spp = Sp - x(a, j) * gen_f(u, d, b, c)
print(Spp)
Sppp = Spp - x(b, d) * gen_f(a, c, j, u)
print(Sppp)
Spppp = Sppp + x(c, d) * gen_f(a, b, j, u)
print(Spppp)