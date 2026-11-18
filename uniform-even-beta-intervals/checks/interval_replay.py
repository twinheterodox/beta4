"""Replay finite nonvanishing by a second modulus and reverse harmonic summation."""

import json
import math
from pathlib import Path

from flint import nmod_poly

from initial_intervals import cover, local_fraction, N, ORDERS
from prime_bridge import is_prime


def binomial_weights(power, derivative, p):
    top = max(ORDERS) * power
    if top < p:
        weights = [0, 1]
        for i in range(1, top):
            weights.append(weights[-1] * (i + derivative) * pow(i, -1, p) % p)
        return weights
    return [0] + [math.comb(i + derivative - 1, derivative) % p for i in range(1, top + 1)]


def reverse_residues(p, cases):
    assert p > 192 and is_prime(p)
    totals = {case: 0 for case in cases}
    suffix = {case: nmod_poly([], p) for case in cases}
    bins = {case: binomial_weights(*case, p) for case in cases}
    largest_power = max(power for power, _ in cases)
    for a in range(N, 0, -1):
        # Directly construct the polynomial at every address, including the
        # right half. This route never infers an address by reflection.
        numerator, denominator, shift = local_fraction(a, p)
        degree = (ORDERS[a] - shift) * largest_power
        local = (numerator * denominator.inverse_series_trunc(degree)).truncate(degree)
        inv = pow(2 * a - 1, -1, p)
        for case in cases:
            power, derivative = case
            length = (ORDERS[a] - shift) * power
            sol = local.truncate(length).pow_trunc(power, length).coeffs()
            poly = [0] * (length + 1)
            for h, coeff in enumerate(sol):
                i = length - h
                poly[i] = int(coeff) * bins[case][i] % p
            suffix[case] += (-1)**a * nmod_poly(poly, p)
            # h=a-1. The suffix already contains every pole address >h.
            totals[case] -= (-1)**(a - 1) * pow(inv, derivative, p) * int(suffix[case](inv))
    return {case: residue % p for case, residue in totals.items()}


def main():
    root = Path(__file__).resolve().parent
    first_receipt = json.loads((root / "initial-intervals.json").read_text())
    bridge_receipt = json.loads((root / "prime-bridge.json").read_text())
    cases = [(power, d) for power, d, _, _ in cover()]
    assert len(cases) == len(set(cases))
    second_prime = 100003
    residues = reverse_residues(second_prime, cases)
    for row in first_receipt["seeds"]:
        case = row["power"], row["derivative"]
        assert row["normalized_seed_residue"] % first_receipt["prime"] and residues[case]
        for q in range(row["first_q"], row["last_q"] + 1):
            assert case[1] + 2 >= 2 * q + 2
            assert (31 * case[0] + case[1]) // 2 <= (33 * q) // 2

    bridge = bridge_receipt["prime_exponents_and_normalized_residues"]
    for p, reported in bridge + bridge_receipt["control_residues"]:
        assert reverse_residues(p, [(9, 20)])[(9, 20)] == reported != 0
    print(json.dumps({
        "second_prime": second_prime,
        "unreflected_reverse_seed_residues": [[*case, residues[case]] for case in cases],
        "individual_prime_bridge_checks": len(bridge),
        "additional_control_moduli_checked": 2,
        "all_interval_inclusions_checked": True,
    }, indent=2))


if __name__ == "__main__":
    main()
