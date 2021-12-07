from itertools import combinations
from math import factorial

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
codimension = 9 - dimension
assert codimension == 3
assert I.is_prime()

unit_vector = lambda i: [1 if j == i else 0 for j in range(10)]
hyperplane_vectors = [unit_vector(i) for i in range(10)]
for vecs in combinations(hyperplane_vectors, codimension):
    eqs = [
        sum(map(lambda t: t[0] * t[1], zip(vec, R.gens())))
            for vec in vecs
    ]
    J = I + R.ideal(eqs)
    print(J.hilbert_numerator())
    hp = J.hilbert_polynomial()
    print(hp)
    degree = hp.leading_coefficient()
    print(degree)
