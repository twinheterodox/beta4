"""Block profiles and balanced products used by the proof checkers.

These two functions are retained verbatim from the original research
module; unrelated determinant experiments are omitted.
"""
from __future__ import annotations
from flint import fmpz_poly


def truncated_product(factors, precision):
    """Balanced multiplication avoids multiplying huge early prefixes."""
    factors=list(factors)
    if not factors:
        return fmpz_poly([1])
    while len(factors)>1:
        following=[factors[i].mul_low(factors[i+1],precision) for i in range(0,len(factors)-1,2)]
        if len(factors)%2:
            following.append(factors[-1])
        factors=following
    return factors[0]


def profile(n: int, pole_order: int, kind: str = "uniform") -> list[int]:
    """Symmetric odd-total block profiles, with explicit finite rounding."""
    if kind == "uniform":
        return [pole_order] * (n + 1)
    if kind == "krisman35":
        if pole_order != 35:
            raise ValueError("krisman35 requires pole_order=35")
        offsets = [0] * 8 + list(range(1, 9)) + list(range(10, 47, 2))
    elif kind == "reduced33":
        if pole_order != 33:
            raise ValueError("reduced33 requires pole_order=33")
        offsets = [0] * 8 + list(range(1, 9)) + list(range(10, 43, 2))
    else:
        raise ValueError(kind)
    # ceil rounding produces a finite analogue; only n divisible by 160 has
    # exactly the published/asymptotic rational offsets.
    aa = [(a * n + 159) // 160 for a in offsets]
    return [sum(a <= k <= n - a for a in aa) for k in range(n + 1)]
