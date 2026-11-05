"""Exact finite comparisons for the rounded and unrounded count bounds.

Only the standard library is used. The infinite-range inequalities are
proved in Section 9 of the manuscript; their constants are checked here.
"""
from fractions import Fraction as Q
from hashlib import sha256
import json
from math import gcd
from pathlib import Path

HERE = Path(__file__).resolve().parent
TERMS = 18
C_SQUARED_UPPER = Q(968, 6615)
SHARPER = {2: 160, 3: 364, 4: 588, 5: 1030, 6: 1296,
           7: 1992, 8: 2502, 9: 3260, 10: 3828}


def log_series(y):
    y = Q(y)
    assert 1 <= y <= 2
    z = (y - 1) / (y + 1)
    z_squared = z*z
    power, lower = z, Q(0)
    for j in range(TERMS):
        lower += 2*power/(2*j+1)
        power *= z_squared
    remainder = 2*power/((2*TERMS+1)*(1-z_squared))
    return lower, lower+remainder


LOG_TWO = log_series(2)


def log_interval(x):
    y, k = Q(x), 0
    assert y >= 1
    while y >= 2:
        y /= 2
        k += 1
    lower, upper = log_series(y)
    return lower+k*LOG_TWO[0], upper+k*LOG_TWO[1]


def ceiling(x):
    return -((-x.numerator)//x.denominator)


def endpoint(k, d):
    lower_log, upper_log = log_interval(k)
    lower = 2*ceiling(d*(33+Q(45, 2)*lower_log))
    upper = 2*ceiling(d*(33+Q(45, 2)*upper_log))
    assert lower == upper, ("endpoint not isolated", k)
    return lower


def arctangent_interval(z):
    assert 0 < z < 1
    lower = sum(((-1)**j*z**(2*j+1)/(2*j+1)
                 for j in range(8)), Q(0))
    return lower, lower+z**17/17


def verify():
    sums = [0]
    for k in range(1, 939):
        sums.append(sums[-1]+sum(gcd(j, k) == 1 for j in range(1, k+1)))
    assert endpoint(25, sums[25]) == 42170

    unrounded, minimum_unrounded = [], None
    for k in range(26, 256):
        s = endpoint(k, sums[k])
        t = s-2
        lower_log, _ = log_interval(t)
        margin = (k-1)**2*lower_log/t-C_SQUARED_UPPER
        assert margin > Q(24529, 100000000)
        if minimum_unrounded is None or margin < minimum_unrounded[0]:
            minimum_unrounded = (margin, k)
        unrounded.append({"k": k, "D_k": sums[k], "s_k": s, "passed": True})
    assert minimum_unrounded[1] == 29

    rounded, minimum_rounded = [], None
    previous = SHARPER[2]
    for k in range(3, 26):
        s = SHARPER[k] if k in SHARPER else endpoint(k, sums[k])
        assert previous < s
        previous = s
        t = s-2
        lower_log, _ = log_interval(t)
        margin = k*k*lower_log/t-C_SQUARED_UPPER
        assert margin > Q(14325, 100000000)
        if minimum_rounded is None or margin < minimum_rounded[0]:
            minimum_rounded = (margin, k)
        rounded.append({"k": k, "t_k": s, "passed": True})
    assert minimum_rounded[1] == 3
    assert previous == 42170

    # Constants in the proof for every k >= 256.
    lower_log, upper_log = log_interval(256)
    assert Q(11, 2) < lower_log <= upper_log < Q(111, 20)
    tail_upper = (2*(Q(111, 20)+Q(22, 15))
                  *((Q(111, 20)+2)*Q(484, 147)+2)*Q(256, 255**2))
    assert tail_upper < Q(3, 2)
    assert 45*Q(147, 484)*(Q(11, 2)+Q(22, 15)) > 95
    assert Q(49, 18) < Q(11, 4)
    assert Q(11, 4)**9 < 95**2
    assert Q(9, 2)-Q(44, 15) == Q(47, 30) > Q(3, 2)
    assert Q(2, 11)+Q(2, 15) < 1

    # Machin's formula and alternating-series remainders.
    atan5 = arctangent_interval(Q(1, 5))
    atan239 = arctangent_interval(Q(1, 239))
    pi_lower = 16*atan5[0]-4*atan239[1]
    pi_upper = 16*atan5[1]-4*atan239[0]
    assert Q(157, 50) < pi_lower < pi_upper < Q(22, 7)
    c_squared_lower = 2*Q(157, 50)**2/135
    small_floors = []
    for s, count in ((160, 2), (364, 3)):
        lo, hi = log_interval(s)
        assert c_squared_lower*s > count**2*hi
        assert C_SQUARED_UPPER*s < (count+1)**2*lo
        small_floors.append({"s": s, "floor": count})
    assert c_squared_lower > Q(19, 50)**2
    assert log_interval(160)[1]/160 < Q(9, 50)**2
    assert Q(19, 50)-Q(9, 50) == Q(1, 5)

    # The introduction's footnote: a smaller coefficient covers s >= 64.
    alternative_c = Q(179, 500)
    assert alternative_c**2 < c_squared_lower
    alternative_floors = []
    for s, count in ((64, 1), (158, 1), (364, 2)):
        lo, hi = log_interval(s)
        assert alternative_c**2*s > count**2*hi
        assert alternative_c**2*s < (count+1)**2*lo
        alternative_floors.append({"s": s, "floor": count})

    inverses = []
    for s, b in ((42170, 25), (100000, 37), (1000000, 108), (100000000, 937)):
        _, upper_b = log_interval(b)
        lower_next, _ = log_interval(b+1)
        assert sums[b]*(66+45*upper_b) <= s
        assert sums[b+1]*(66+45*lower_next) > s
        inverses.append({"s": s, "finite_inverse_count": b})

    return {
        "status": "all exact endpoint comparisons passed",
        "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "logarithm_series_terms": TERMS,
        "unrounded_cutoff_even_s": 42170,
        "rounded_cutoff_even_s": 160,
        "unrounded_checks": unrounded,
        "rounded_checks": rounded,
        "unrounded_minimum_margin_index": minimum_unrounded[1],
        "unrounded_margin_strict_lower_bound": "24529/100000000",
        "rounded_minimum_margin_index": minimum_rounded[1],
        "rounded_margin_strict_lower_bound": "14325/100000000",
        "infinite_range_constants": {
            "all_passed": True, "first_k": 256,
            "T_256_rational_upper": str(tail_upper),
            "G_256_strict_lower": "47/30",
        },
        "small_floors": small_floors,
        "alternative_coefficient_footnote": {
            "coefficient": str(alternative_c),
            "rounded_cutoff_even_s": 64,
            "smaller_than_main_coefficient": True,
            "endpoint_floors": alternative_floors,
            "range_argument": "Monotonicity and I(64)>=1 cover 64<=s<=158; the main theorem covers every even s>=160.",
        },
        "coefficient_after_rounding_error_greater_than": "1/5",
        "finite_inverse_examples": inverses,
    }


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run without -O: this checker uses assertions.")
    result = verify()
    (HERE / "endpoint-verification.json").write_text(
        json.dumps(result, indent=2)+"\n", encoding="utf-8"
    )
    print(result["status"])
    print("230 unrounded intervals, 23 rounded intervals, all tail constants")
    print("The rounded bound holds for every even s >= 160.")
    print("The 0.358 footnote checks pass at 64, 158 and 364.")
