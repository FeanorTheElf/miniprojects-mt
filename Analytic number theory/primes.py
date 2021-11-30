import itertools
import math

def primes():
    yield 2
    found_primes = [2]
    for n in itertools.count(3):
        for p in found_primes:
            if n % p == 0:
                break
            elif p >= math.sqrt(n):
                yield n
                found_primes.append(n)
                break

def primes_leq(n):
    return itertools.takewhile(lambda p: p <= n, primes())

for i in range(1, 8):
    print("Consider interval [1, 10**" + str(i) + "]")
    print("  Number of primes = 1 mod 4 is " + str(
        sum(1 for p in primes_leq(10**i) if (p - 1) % 4 == 0)
    ))
    print("  Number of primes = 3 mod 4 is " + str(
        sum(1 for p in primes_leq(10**i) if (p - 3) % 4 == 0)
    ))
    print()