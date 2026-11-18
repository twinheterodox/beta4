"""Check the seed-height and prime-interval inequalities used in the manuscript."""

from fractions import Fraction
import json
import math


def main():
    orders = [min(a, 96 - a, 30) + 1 for a in range(97)]
    common_gap = math.lcm(*range(1, 193))
    at_pole = [
        math.prod(abs(2 * (a - b)) ** (9 * order)
                  for b, order in enumerate(orders) if a != b)
        for a in range(97)
    ]
    denominator = math.lcm(*at_pole) * common_gap**299
    terms = 9 * sum(orders)
    height = denominator * 97**9 * 336**6480 * 2**299 * terms * 96
    assert terms == 18693
    assert height.bit_length() == 256201
    assert height.bit_length() * Fraction(7, 10) < 180000
    bands = []
    for start, low_log, width in [(4_000_000, 14, Fraction(1, 10)),
                                  (200_000_000, 18, Fraction(3, 1000))]:
        assert start > 3_594_641 and Fraction(11, 4)**low_log < start
        err = Fraction(1, 5 * low_log**2)
        mass = start * ((1 + width) * (1 - err) - (1 + err))
        assert mass > 180000
        bands.append({"start": start, "relative_width": str(width),
                      "theta_difference_exceeds": str(mass)})
    assert 299 * Fraction(11, 10) / 10 < 33
    assert 299 * Fraction(1003, 1000) / 20 < 15
    print(json.dumps({"height_upper_bound_bits": height.bit_length(),
                      "common_denominator_bits": denominator.bit_length(),
                      "bands": bands}, indent=2))


if __name__ == "__main__":
    main()
