"""Exact comparisons for Appendix A, at the even-beta endpoint 364.

Run with Python 3.8+ and without -O. No third-party packages or private
files are required. These checks concern the estimates and parameter
families specified in the appendix, not refinements of the constructions
or complete audits of the cited proofs. Floating-point extrema are
illustrative only; every certificate decision uses integers or Fractions.
"""

from fractions import Fraction as F
from math import factorial, log
from pathlib import Path
import json


def elementary_bounds():
    # The positive exponential series gives e > 8/3. For k >= 6,
    # consecutive terms have ratio at most 1/7, so the remaining tail
    # is at most (1/6!)/(1-1/7).
    assert sum((F(1, factorial(k)) for k in range(4)), F(0)) == F(8, 3)
    e_upper = sum((F(1, factorial(k)) for k in range(6)), F(0)) + F(7, 4320)
    assert e_upper < F(68, 25)
    assert e_upper < 4  # log(2) > 1/2, hence log(4) > 1.
    assert e_upper**2 < 8  # log(8) > 2.
    t = F(7, 10)
    assert sum((t**k / factorial(k) for k in range(5)), F(0)) > 2
    return {"e_lower": "8/3", "e_upper": "68/25",
            "log_2_bounds": ["1/2", "7/10"], "log_8_gt_2": True}


def rivoal_zudilin():
    q = 366
    displayed, profile_upper = [], []
    for r in range(2, (q - 1) // 2 + 1, 2):
        # Appendix (A.1): exponentiate U(q,r)<2, then replace e by
        # the smaller 8/3. A strict integer inequality suffices.
        numerator = (8**(4*(q-1)) * 2**(q-2)
                     * (2*r+1)**(2*(2*r+1)) * r**(q*r))
        denominator = 3**(4*(q-1)) * (r+1)**(q*(r+1))
        assert numerator > denominator, (q, r)
        height = (2+log(2))*(q-1) + (2*r+1)*log(2*r+1)
        displayed.append(((q*log(r)+(q-2*r-1)*log(2))/height, r))
        profile_upper.append((q*(log(2)+(r+1)*log(r+1)-r*log(r))/height, r))
    assert len(displayed) == 91
    return {
        "source": "Rivoal-Zudilin (2003), Section 6, (26) and (29)",
        "q": q, "beta_endpoint": q-2, "even_r_cases": len(displayed),
        "all_profile_upper_bounds_lt_2": True,
        "displayed_max_illustrative": max(displayed),
        "profile_upper_max_illustrative": max(profile_upper),
        "extension_to_smaller_endpoints": "U(q,r) increases with q for fixed r>=2.",
        "scope": "Uniform family (26), using the published height estimate."
    }


def nishimoto():
    # With a<=364, 4|d, r>=2 and (a+d)/d>=3r, we have a>=5d.
    # Thus r+1 <= 364/12+4/3 <32 and Delta <10*log(2)+2 <9.
    assert F(364, 12) + F(4, 3) < 32
    assert 10*F(7, 10)+2 == 9
    # For d=4, A<36 and B>9*a>=180. For d>=8,
    # A<(8-2*d)*a+9*d <= d*(49-10*d)<0.
    assert F(1) + F(36, 180) == F(6, 5)
    assert 49-10*8 < 0
    return {
        "source": "Nishimoto (2011), Theorem 5.5",
        "max_beta_endpoint": 364, "periods": "positive multiples of 4",
        "period_4_dimension_estimate_lt": "6/5",
        "larger_periods_gain_numerator_negative": True,
        "scope": "Explicit expression in Theorem 5.5; algebraic reduction in Appendix A.2."
    }


def fischler_2018():
    a, rows = 365, []
    # Ordinary representation: 4|N. Twisted representation: N=2 mod 4.
    for N in range(2, a//2+1, 2):
        for r in range(1, (a-1)//(2*N)+1):
            # Raise alpha to N and replace e by 8/3.
            power = -a + N*(4*r+2)
            numerator = (8**(a*N) * 4**(a-N*r)
                         * (N+1)**(N*(2*r+2)) * r**max(power, 0))
            denominator = 3**(a*N) * r**max(-power, 0)
            assert numerator > denominator, (N, r)
            # At a=0, (N+1)^(2r+2)>4^r and r^(4r+2)>=1.
            assert (N+1)**(2*r+2) > 4**r
            value = (a+(a/N-r)*log(4)+(2*r+2)*log(N+1)
                     +(-a/N+4*r+2)*log(r))
            rows.append((value, N, r))
    assert len(rows) == 429
    return {
        "source": "Fischler (2018), Section 4.2 and Proposition 4.4",
        "a": a, "beta_endpoint": 364, "periods": "all positive even N",
        "parameter_pairs": len(rows), "all_alpha_gt_1": True,
        "minimum_log_alpha_illustrative": min(rows),
        "extension_to_smaller_a": "log(alpha) is affine in a and positive at both a=0 and a=365."
    }


def fischler_2020():
    a, rows = 365, []
    for N in range(4, a//3+1, 4):
        for r in range(1, (a-1)//(3*N)+1):
            # Raise alpha to N and replace e by 8/3.
            power = -(a+1) + 4*N*(r+1)
            numerator = (32**(a+1) * (2*N)**(N*(2*r+2))
                         * r**max(power, 0))
            denominator = 3**(a+1) * r**max(-power, 0)
            assert numerator > denominator, (N, r)
            # At a=0, all factors are >=1 and the first two are >1.
            assert 4*(r+1)-F(1, N) > 0
            value = (((a+1)/N)*(log(4)+1)+(2*r+2)*log(2*N)
                     +(-(a+1)/N+4*(r+1))*log(r))
            rows.append((value, N, r))
    assert len(rows) == 111
    return {
        "source": "Fischler (2020), Proposition 1, (1.1), author revision 10 October 2019",
        "a": a, "beta_endpoint": 364, "periods": "all positive multiples of 4",
        "parameter_pairs": len(rows), "all_alpha_gt_1": True,
        "minimum_log_alpha_illustrative": min(rows),
        "extension_to_smaller_a": "log(alpha) is affine in a and positive at both a=0 and a=365.",
        "scope": "Includes the enlarged periods used in the irrational-value count."
    }


def square_root_bound():
    # Appendix A.4 derives incompatible inequalities if R>=2:
    # (Omega*log(r)-3*X)/kappa <365/32 and >6+3*H_5.
    left_upper = F(365, 32)
    right_lower = 6 + 3*sum((F(1, k) for k in range(1, 6)), F(0))
    assert right_lower == F(257, 20)
    assert left_upper < right_lower
    return {
        "source": "twinheterodox (2026), (2.2), (5.2), (5.3) and Section 7",
        "max_a_plus_h": 365, "all_admissible_R_lt": 2,
        "left_upper": str(left_upper), "right_lower": str(right_lower),
        "scope": "Exact constants for the all-parameter argument in Appendix A.4."
    }


def main():
    if not __debug__:
        raise SystemExit("Run without -O: the certificate checks use assertions.")
    result = {
        "purpose": "Finite comparisons for the specified estimates in Appendix A.",
        "certificate_arithmetic": "integers and exact rational numbers",
        "decimal_values": "illustrative only",
        "elementary_bounds": elementary_bounds(),
        "rivoal_zudilin_2003": rivoal_zudilin(),
        "nishimoto_2011": nishimoto(),
        "fischler_2018": fischler_2018(),
        "fischler_2020": fischler_2020(),
        "square_root_bound": square_root_bound(),
    }
    output = Path(__file__).with_name("asymptotic-comparisons.json")
    output.write_text(json.dumps(result, indent=2)+"\n")
    print("All exact comparisons passed: 91 Rivoal-Zudilin cases, "
          "429 Fischler 2018 cases, 111 Fischler 2020 cases, "
          "and the rational constants for Nishimoto and the square-root bound.")


if __name__ == "__main__":
    main()
