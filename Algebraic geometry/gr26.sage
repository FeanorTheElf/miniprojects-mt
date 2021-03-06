from itertools import combinations

codes = [str(i) for i in range(1, 7)]
# build up the ring and group the variables nicely
R = PolynomialRing(QQ, [
    "x" + codes[i] + codes[j]
        for i in range(6) for j in range(i + 1, 6)
])

_x = []
vars = list(R.gens())
for i in range(5, 0, -1):
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
print(I.hilbert_polynomial())