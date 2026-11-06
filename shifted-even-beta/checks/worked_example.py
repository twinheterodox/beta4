"""Exact self-check of the worked B=2, m=1 Gram example."""
from fractions import Fraction as F
from pathlib import Path
import json
from checks import (
    parameters, kernels, mul, moment, forms, determinant, lcm
)


def main():
    assert __debug__
    B, m, r, a = 2, 1, 2, 70
    _, D, L, K, _ = parameters(B, r)
    C = 2**a * K
    assert (D, L) == (3, 5) and C == 2**85 * 3**16
    hs, distances = kernels(B, m)
    assert hs == ((-4, -9), (4, -9)) and distances == (F(2, 3),)

    # w(y)=C/(18*pi) times this polynomial in y^2 over (1+y^2)^70.
    weight = [F(1)]
    for alpha in [F(4, 3), F(2), F(8, 3), F(10, 3), F(4), F(14, 3)]:
        weight = mul(weight, [alpha**2, 1])
    mu = [F(C, 9)*sum(c*moment(a, j+ell) for j, c in enumerate(weight))
          for ell in range(3)]
    u, v, w = mu
    N = [[16*u+72*v+81*w, -16*u+81*w],
         [-16*u+81*w, 16*u-72*v+81*w]]
    columns = [forms(B, m, r, a, 1, h)[0] for h in hs]
    M = [list(row) for row in zip(*columns)]
    assert M == [[-x for x in row] for row in N]
    det = determinant(N)
    assert det == 5184*(u*w-v*v) > 0

    coefficients = [N[0][1], -N[0][0]]
    assert [sum(M[i][j]*coefficients[j] for j in range(2))
            for i in range(2)] == [0, det]
    scale = 1
    for value in coefficients:
        scale = lcm(scale, value.denominator)
    integer_h = tuple(int(scale*sum(coefficients[j]*hs[j][i] for j in range(2)))
                      for i in range(2))
    seed = forms(B, m, r, a, 1, integer_h)[0]
    assert seed == (0, scale*det) and seed[1] != 0

    record = {
        "scope": "Exact self-check of the post-audit exposition example; not an independent audit.",
        "parameters": {"B": B, "m": m, "r": r, "a": a},
        "h1": hs[0], "h3": hs[1],
        "normalization": "C_1=2^85*3^16; weight=C_1/(18*pi)*product of six factors/(1+y^2)^70.",
        "gram_matrix_matches_seed_partial_fractions": True,
        "determinant_identity": True,
        "rational_combination_has_seed": [0, "det(N)>0"],
        "integer_combination_has_nonzero_last_seed_only": True
    }
    Path(__file__).with_suffix(".json").write_text(json.dumps(record, indent=2)+"\n")
    print("Worked Gram example: all exact checks passed.")


if __name__ == "__main__":
    main()
