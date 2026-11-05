"""Exact rational specialization of Sorokin's displayed finite estimate.

This checks the numerical comparison in Appendix A, not the source proof
or optimality of any endpoint.
"""
from fractions import Fraction as F
from pathlib import Path
import json


def main():
    if not __debug__:
        raise RuntimeError("Run without -O: assertions are part of this check.")
    e_upper, kernel_upper = F(68, 25), F(32361, 10000)
    assert F(163, 60) + F(7, 4320) < e_upper
    assert kernel_upper > 1 and (kernel_upper - 1)**2 > 5
    a, b, r = 628, 5, 11
    assert a % 2 == 0 and b % 2 == 1 and 2*b*r <= a
    c = a//2 + b
    upper = (e_upper**(2*(a+b-1)) * 2**a * F(r)**(r*c)
             * (2*kernel_upper*(r+1))**(2*b*(r+1))
             / F(r+2)**((r+2)*c))
    assert upper < 1
    result = {
        "a": a, "b": b, "r": r, "endpoint": b+a-1,
        "kernel_upper": str(kernel_upper), "e_upper": str(e_upper),
        "rate_upper_less_than_one": True,
        "rate_upper_illustrative_decimal": float(upper),
        "arithmetic": "Exact rational; decimal is illustrative only.",
        "scope": "Appendix A specialization of the source estimate; not a proof audit or optimality claim.",
    }
    Path(__file__).with_name("comparison.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
