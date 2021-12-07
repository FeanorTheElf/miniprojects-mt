from itertools import combinations
from math import factorial
import numpy as np

# build up the ring and group the variables nicely
R = PolynomialRing(QQ, [
    "x" + str(i) + str(j) 
        for i in range(1, 6) for j in range(i + 1, 6)
])
x12, x13, x14, x15, x23, x24, x25, x34, x35, x45 = R.gens()
_x = [[x12, x13, x14, x15], [x23, x24, x25], [x34, x35], [x45]]
x = lambda i, j: _x[i - 1][j - i - 1]

# construct the ideal describing Gr(2, 5)
polys = []
for seq in combinations([1, 2, 3, 4, 5], 4):
    (i, j, u, v) = sorted(seq)
    p = x(i,j) * x(u,v) + x(i,v) * x(j,u) - x(i,u) * x(j,v)
    polys.append(p)

I = R.ideal(polys)
dimension = I.dimension() - 1
assert dimension == 6
assert I.is_prime()

hyperplane_vectors = [
    np.random.randint(-4, 4, R.ngens(), int) 
        for i in range(dimension + 1)
]
for vecs in combinations(hyperplane_vectors, dimension):
    eqs = [
        sum(map(lambda t: t[0] * t[1], zip(vec, R.gens())))
            for vec in vecs
    ]
    J = I + R.ideal(eqs)
    # the number of intersection points is clearly equal to the
    # dimension of S(X)_d for large enough d
    hp = J.hilbert_polynomial()
    degree = hp.leading_coefficient()
    print(degree) # usually prints 5
