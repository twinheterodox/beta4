"""Check the three-node, prime-seven illustration using exact fractions."""
from fractions import Fraction as F
from pathlib import Path
import json


def tail_sum(order, node):
    return sum((F(1, (4*j+1)**order) - F(1, (4*j+3)**order)
                for j in range(node)), F(0))


def valuation(value, prime):
    assert value
    numerator, denominator = value.numerator, value.denominator
    result = 0
    while numerator % prime == 0:
        numerator //= prime
        result += 1
    while denominator % prime == 0:
        denominator //= prime
        result -= 1
    return result


def residue(value, prime):
    assert value.denominator % prime
    return value.numerator * pow(value.denominator, -1, prime) % prime


def run():
    nodes = [1, 2, 3]
    matrix = [[-tail_sum(4, a) + tail_sum(3, a)/(12*a) if a == b
               else (b*tail_sum(3, b)-a*tail_sum(3, a))/(6*(b*b-a*a))
               for b in nodes] for a in nodes]
    scaled = [[residue(7**4*x, 7) for x in row] for row in matrix]
    assert scaled == [[0, 0, 0], [0, 1, 0], [0, 0, 1]]
    minors = [matrix[i][i]*matrix[j][j]-matrix[i][j]**2
              for i, j in [(0, 1), (0, 2), (1, 2)]]
    minor_residues = [residue(7**8*x, 7) for x in minors]
    assert minor_residues == [0, 0, 1]
    coefficient = sum(minors, F(0))
    assert valuation(coefficient, 7) == -8
    assert residue(7**8*coefficient, 7) == 1
    return {
        "status": "passed", "nodes": nodes, "prime": 7,
        "scaled_matrix_mod_7": scaled,
        "scaled_principal_minor_residues": minor_residues,
        "coefficient_of_X": str(coefficient),
        "coefficient_valuation": -8,
        "primitive_multiplier_valuation_lower_bound": 8,
        "scope": "One exact explanatory example, not a new asymptotic check or independent audit."
    }


if __name__ == "__main__":
    result = run()
    target = Path(__file__).with_name("expository-example.json")
    target.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
