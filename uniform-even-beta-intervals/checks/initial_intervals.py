"""Certify the finite initial interval cover using exact Laurent jets modulo 65537."""

from array import array
import json
import math

from flint import nmod_poly


N = 96
PRIME = 65537
FIRST_PRIME_BRIDGE = 1749
ORDERS = [1 + min(a, N - a, 30) for a in range(N + 1)]


def cover():
    """Each tuple is (power, derivative, first original q, last original q)."""
    cases = []
    q = 1
    while q < FIRST_PRIME_BRIDGE:
        j = (33 * q) // 299
        if j >= 1 and q <= 10 * j:
            power, derivative, end = 9 * j, 20 * j, 10 * j
        else:
            power, derivative, end = q, 2 * q, q
        assert derivative + 2 >= 2 * end + 2
        assert 31 * power + derivative <= 33 * q
        cases.append((power, derivative, q, end))
        q = end + 1
    return cases


def harmonic_table(maximum, p):
    """Table of sigma_a(s), including a=0. Each row uses four bytes per residue."""
    result = [array("I", [0]) * (maximum + 1)]
    for a in range(N):
        row = array("I", result[-1])
        inverse = pow(2 * a + 1, -1, p)
        term = (-1)**a
        for s in range(1, maximum + 1):
            term = term * inverse % p
            row[s] = (row[s] + term) % p
        result.append(row)
    return result


def local_fraction(a, p):
    """Unit numerator and denominator of the pole-removed seed; central shift separately."""
    x = N - 2 * a
    numerator = nmod_poly([1], p)
    for v in range(-239, 240, 2):
        numerator *= nmod_poly([x - v, 1], p)
    numerator **= 3
    if x:
        numerator *= nmod_poly([x, 1], p)
    denominator = nmod_poly([1], p)
    for b, order in enumerate(ORDERS):
        if a != b:
            denominator *= nmod_poly([2 * (b - a), 1], p)**order
    return numerator, denominator, int(x == 0)


def certificate():
    p = PRIME
    assert all(p % k for k in range(2, math.isqrt(p) + 1))
    cases = cover()
    maximum = max(31 * power + d for power, d, _, _ in cases)
    assert p > max(480, maximum)
    harmonic = harmonic_table(maximum, p)
    inverses = [0] + [pow(i, -1, p) for i in range(1, maximum + 1)]
    bins = {}
    for power, d, _, _ in cases:
        weights = array("I", [0, 1])
        for i in range(1, 31 * power):
            weights.append(weights[-1] * (i + d) * inverses[i] % p)
        bins[power, d] = weights

    constants = {(power, d): 0 for power, d, _, _ in cases}
    verified = 0
    for a in range(N // 2 + 1):
        numerator, denominator, shift = local_fraction(a, p)
        largest = max(power for power, _, _, _ in cases) * (ORDERS[a] - shift)
        local = (numerator * denominator.inverse_series_trunc(largest)).truncate(largest)
        log_den = numerator * denominator
        log_num = numerator.derivative() * denominator - denominator.derivative() * numerator
        for power, d, _, _ in cases:
            length = power * (ORDERS[a] - shift)
            jet = local.truncate(length).pow_trunc(power, length)
            # Y=(P/Q)^power. The constant and this differential equation
            # uniquely certify every coefficient used, since Q(0)P(0) is a
            # unit and length < p. The central w^power was divided out.
            assert int(jet[0]) == pow(int(numerator[0]) * pow(int(denominator[0]), -1, p), power, p)
            assert not (log_den * jet.derivative() - power * log_num * jet).truncate(length - 1)
            verified += length
            weights = bins[power, d]
            total = 0
            for h, elem in enumerate(jet.coeffs()):
                i = length - h
                value = harmonic[a][i + d]
                if a != N - a:
                    # Reflection c_(N-a,i)=(-1)^i c_(a,i), with N even.
                    reflected = harmonic[N - a][i + d]
                    value += reflected if i % 2 == 0 else -reflected
                total += int(elem) * weights[i] * value
            constants[power, d] = (constants[power, d] - (-1)**a * total) % p
    rows = []
    for power, d, start, end in cases:
        residue = constants[power, d]
        assert residue, (power, d, start, end)
        rows.append({"power": power, "derivative": d, "first_q": start,
                     "last_q": end, "normalized_seed_residue": residue})
    return {"prime": p, "first_q": rows[0]["first_q"], "last_q": rows[-1]["last_q"],
            "largest_index_in_binomials": maximum, "number_of_certified_seeds": len(rows),
            "verified_unit_jet_coefficients_before_reflection": verified, "seeds": rows}


if __name__ == "__main__":
    print(json.dumps(certificate(), indent=2))
