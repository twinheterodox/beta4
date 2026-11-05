"""Integer polynomial basis for the uniform odd-index seed construction."""
from fractions import Fraction as F
from math import comb, gcd
from base import (Parameters, partials_newton, partials_product, constant_term,
                  constant_term_direct, determinant, evaluate, modq, chi)


def mul_linear(poly, constant, slope):
    out = [0]*(len(poly)+1)
    for i, c in enumerate(poly):
        out[i] += constant*c
        out[i+1] += slope*c
    return out


def divide_linear(poly, constant, slope):
    out = [0]*(len(poly)-1)
    remainder = poly[:]
    for i in range(len(out)-1, -1, -1):
        assert remainder[i+1] % slope == 0
        out[i] = remainder[i+1]//slope
        remainder[i+1] -= slope*out[i]
        remainder[i] -= constant*out[i]
    assert not any(remainder)
    return out


def gram_basis(B):
    nodes = sorted({F(2*j-1, b) for b in range(1, B+1)
                    for j in range(1, b+1)})
    Q = [1]
    for a in nodes:
        Q = mul_linear(Q, a.numerator**2, a.denominator**2)
    quotients = {a: divide_linear(Q, a.numerator**2, a.denominator**2)
                 for a in nodes}
    hs = []
    for b in range(1, B+1):
        h = [0]*len(nodes)
        for j in range(1, b+1):
            a = F(2*j-1, b)
            coefficient = -2*(-1)**(b-j)*a.numerator*a.denominator
            for k, c in enumerate(quotients[a]):
                h[k] += coefficient*c
        hs.append(h)
    return nodes, Q, hs


def poly_eval(poly, x):
    value = 0
    for a in reversed(poly):
        value = value*x+a
    return value


def local_multiplier(poly, J, n, h):
    # W_n(t)=sum_k poly[k](-1)^k n^(2J-2k)(t+2n)^(2k+1).
    c = 2*n-4*h
    out = [0]*(2*J+2)
    for k, a in enumerate(poly):
        a *= (-1)**k*n**(2*J-2*k)
        for v in range(2*k+2):
            out[v] += a*comb(2*k+1, v)*c**(2*k+1-v)
    return out


def modified_coefficients(P, a, poly, J, p, n=1):
    out = {}
    for h in range(n+1):
        w = [c % p for c in local_multiplier(poly, J, n, h)]
        for i in range(1, P.S):
            out[i, h] = sum(c*a.get((i+k, h), 0) for k, c in enumerate(w)) % p
    return out


def H_centered(b, x):
    return sum(-2*(-1)**(b-j)*F(2*j-1, b)/(x+F(2*j-1, b)**2)
               for j in range(1, b+1))


def poly_mul(a, b):
    out = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def moment_matrix(P):
    """Exact rational Gram formula for the unnormalized seed matrix.

    Returns J and a rational scalar such that M_raw=scalar*J.
    J uses moments scaled by 2^(4S-2). Their integrality is checked,
    not needed as an imported theorem in the general proof.
    """
    nodes, Q, hs = gram_basis(P.B)
    alpha = sorted([F(2)-c for c in P.constants(1) if c < 2])
    assert len(alpha) == P.D*P.L//2 and len(set(alpha)) == len(alpha)
    assert set(nodes) <= set(alpha)
    scale, T = 1, [1]
    for a in alpha:
        scale *= a.denominator**2
        if a not in nodes:
            T = mul_linear(T, a.numerator**2, a.denominator**2)
    assert scale == __import__('math').prod(theta.denominator**P.L for theta in P.theta)
    assert P.K % scale == 0
    Jdegree = len(nodes)-1
    degree = 1+len(T)-1+2*Jdegree
    assert degree < P.S
    moments = [comb(2*P.S-2, P.S-1)]
    for k in range(degree):
        numerator = 4*(2*k+1)*moments[-1]
        denominator = 2*P.S-2*k-3
        assert numerator % denominator == 0
        moments.append(numerator//denominator)
    J = []
    for hb in hs:
        row = []
        for hc in hs:
            poly = [0]+poly_mul(T, poly_mul(hb, hc))
            row.append(sum(a*m for a, m in zip(poly, moments)))
        J.append(row)
    sign = (-1)**(P.D*P.L//2+1)
    scalar = F(sign, scale*2**(4*P.S-2))
    return J, scalar, {'paired_root_count': len(alpha),
                       'kernel_poles_cancel_in_numerator': True,
                       'remaining_positive_polynomial_degree': len(T)-1,
                       'largest_moment': degree,
                       'moment_count': len(moments),
                       'common_sign': sign}


def seed_matrix(P, p=1000003, two_methods=True):
    nodes, Q, hs = gram_basis(P.B)
    J = len(nodes)-1
    assert len(nodes) <= 2*P.D and P.s > 2*J+1
    for b in range(1, P.B+1):
        for x in (F(0), F(1, 2), F(2), F(7)):
            assert poly_eval(hs[b-1], x)/F(poly_eval(Q, x)) == H_centered(b, x)
    a = partials_newton(P, p, n=1)
    if two_methods:
        assert a == partials_product(P, p, n=1)
    M = [[0]*P.B for _ in range(P.B)]
    for j, poly in enumerate(hs):
        aj = modified_coefficients(P, a, poly, J, p)
        assert all(aj[i, 1] == (-1)**i*aj[i, 0] % p for i in range(1, P.S))
        assert all(sum(aj[i, h] for h in range(2)) % p == 0 for i in range(1, P.S, 2))
        for t in (F(P.r+1), F(2*P.r+3, 2)):
            W = (t+2)*poly_eval(poly, -(t+2)**2)
            expected = evaluate(P, t, p, n=1)*modq(W, p) % p
            actual = sum(value*pow(modq(t+4*h, p), -i, p)
                         for (i, h), value in aj.items()) % p
            assert actual == expected
        for b in range(1, P.B+1):
            M[b-1][j] = constant_term(P, aj, b, p, n=1)
            assert M[b-1][j] == constant_term_direct(P, aj, b, p, n=1)
            # Reflection and zeros give q_b=Res_(t=0) T(t) H_b(t).
            residue = sum(chi(k)*sum((-1)**(i-1)*aj[i, 0]
                          *pow(b*pow(k, -1, p) % p, i, p)
                          for i in range(1, P.S))
                          for k in range(1, 4*b, 2)) % p
            assert residue == M[b-1][j]
    assert all(M[i][j] == M[j][i] for i in range(P.B) for j in range(P.B))
    det = determinant(M, p)
    return {'parameters': {'B': P.B, 's': P.s, 'r': P.r},
            'prime': p, 'seed_index': 1, 'node_count': len(nodes),
            'multiplier_degree': 2*J+1, 'dilation_exponent': P.s-2*J-1,
            'unnormalized_matrix': M, 'symmetric': True,
            'determinant_residue': det,
            'two_partial_fraction_algorithms_agree': two_methods,
            'two_constant_algorithms_agree': True,
            'residue_identity_checks': P.B**2,
            'rational_kernel_identity_checks': 4*P.B,
            'reflection_and_odd_cancellation': True,
            'modified_reconstruction_checks': 2*P.B}
