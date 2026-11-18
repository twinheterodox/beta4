"""Check one fixed seed at a chain of prime exponents using the Chinese remainder theorem."""

import json
import math
from kernel_coefficients import N, ORDERS, coefficients_by_products


def is_prime(n):
    return n >= 2 and all(n % k for k in range(2, math.isqrt(n) + 1))


def exponent_chain():
    chain = [193]
    while 10 * chain[-1] < 40_000_000:
        candidate = 27 * chain[-1] // 25
        while not is_prime(candidate):
            candidate -= 1
        assert candidate > chain[-1]
        chain.append(candidate)
    return chain


def normalized_seed(modulus):
    # Only the local polar factors and harmonic denominators are inverted.
    # Every one is a unit for these primes, even when it divides a binomial.
    derivative = 20
    top = 9 * max(ORDERS)
    bins = [0] + [math.comb(i + derivative - 1, derivative) for i in range(1, top + 1)]
    harmonic = [0] * (top + derivative + 1)
    seed = 0
    for a in range(N + 1):
        row = coefficients_by_products(a, 9, modulus)
        for i, coeff in enumerate(row, 1):
            seed -= (-1)**a * bins[i] * coeff * harmonic[i + derivative]
        seed %= modulus
        if a == N:
            break
        inverse = pow(2 * a + 1, -1, modulus)
        term = pow(inverse, derivative, modulus)
        for s in range(derivative + 1, top + derivative + 1):
            term = term * inverse % modulus
            harmonic[s] = (harmonic[s] + (-1)**a * term) % modulus
    return seed


def main():
    primes = exponent_chain()
    controls = [10007, 10009]
    modulus = math.prod(set(primes + controls))
    seed = normalized_seed(modulus)
    assert [seed % p for p in controls] == [1491, 7947]
    valid = [p for p in primes if seed % p]
    assert valid and valid[0] == primes[0]
    intervals = [((299 * p + 32) // 33, 10 * p) for p in valid]
    assert all(left[1] + 1 >= right[0] for left, right in zip(intervals, intervals[1:]))
    assert intervals[-1][1] >= 40_000_000
    print(json.dumps({
        "power": 9, "derivative": 20, "first_q": intervals[0][0],
        "last_q": intervals[-1][1], "tested_primes": len(primes),
        "prime_exponents_and_normalized_residues": [[p, seed % p] for p in primes],
        "control_residues": [[p, seed % p] for p in controls],
    }, indent=2))


if __name__ == "__main__":
    main()
