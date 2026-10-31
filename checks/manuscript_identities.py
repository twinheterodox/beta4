"""Exact checks of the Rodrigues norm calculation added to the manuscript.

These finite checks use polynomial coefficients and rational integration,
separately from the archived determinant implementations. They are not
an independent mathematical audit or formal verification.
"""
from fractions import Fraction as F
from math import comb, factorial


def determinant(rows):
    a = [list(map(F, row)) for row in rows]
    result = F(1)
    for i in range(len(a)):
        pivot = next((j for j in range(i, len(a)) if a[j][i]), None)
        if pivot is None:
            return F(0)
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            result = -result
        v = a[i][i]
        result *= v
        for j in range(i + 1, len(a)):
            c = a[j][i] / v
            for k in range(i + 1, len(a)):
                a[j][k] -= c * a[i][k]
            a[j][i] = F(0)
    return result


def run():
    norms, orthogonality = [], 0
    for j in range(13):
        # Differentiate y^j(1-y)^(j+1) j times.
        derivative = [F((-1)**k * comb(j+1, k) * factorial(j+k), factorial(k))
                      for k in range(j+2)]
        quotient = []
        acc = F(0)
        for c in derivative[:-1]:
            acc += c
            quotient.append(acc)
        assert derivative[-1] == -quotient[-1]
        scale = F((-1)**j * factorial(j+1), factorial(2*j+1))
        p = [c * scale for c in quotient]
        assert len(p) == j + 1 and p[-1] == 1
        for k in range(j):
            moment = sum((c / ((i+k+1)*(i+k+2)) for i, c in enumerate(p)), F(0))
            assert moment == 0
            orthogonality += 1
        norm = sum((a*b / ((i+k+1)*(i+k+2))
                    for i, a in enumerate(p) for k, b in enumerate(p)), F(0))
        expected = F(1, (2*j+2)*comb(2*j+1, j)**2)
        assert norm == expected
        norms.append(norm)
    products = []
    product = F(1)
    for n in range(1, 10):
        product *= norms[n-1]
        direct = determinant([[F(1, (i+j+1)*(i+j+2)) for j in range(n)]
                              for i in range(n)])
        assert direct == product
        products.append({'n': n, 'determinant': str(direct)})
    # Primitive content example in Section 2, including a noninteger multiplier.
    assert [F(5, 6)*F(a, 5) for a in (12, 6)] == [2, 1]
    return {'monic_polynomials': 13, 'orthogonality_equalities': orthogonality,
            'exact_norms': [str(v) for v in norms],
            'gram_product_checks': products, 'primitive_example': 'passed'}


if __name__ == '__main__':
    import json
    print(json.dumps(run(), indent=2))
