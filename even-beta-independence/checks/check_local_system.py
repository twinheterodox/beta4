"""Exact checks supporting the local-system and rank derivations.

This is a same-assistant computational audit, not an independent mathematical
review or formal verification. No numerical value of beta enters these checks.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb, factorial, lcm
from pathlib import Path
import json
import sympy as sp

HERE = Path(__file__).resolve().parent
report = {}


def delta(d, k):
    return __import__("functools").reduce(
        lambda x, j: x * lcm(*range(1, k//j + 1)), range(1, d+1), 1)


# Direct definition of the multiplier versus the prime-valuation formula.
count = 0
for k in range(13):
    for d in range(1, 5):
        direct = 1
        for size in range(min(k, d)+1):
            for subset in combinations(range(1, k+1), size):
                product = 1
                for t in subset:
                    product *= t
                direct = lcm(direct, product)
        assert direct == delta(d, k)
        count += 1
report["distinct_products_lcm_checks"] = count

# Route 1: dense Taylor recurrence for the original differential system.
d, K = 7, 36
U = [[[F(int(i == j)) for j in range(d+1)] for i in range(d+1)]]
for k in range(K):
    nxt = [[F(0) for j in range(d+1)] for i in range(d+1)]
    for q in range(k+1):
        a10 = -sum((F(1, 2**j) for j in range(q+1)), F(0))
        prev = U[k-q]
        for j in range(d+1):
            nxt[1][j] += a10*prev[0][j] + prev[1][j]/2
            for i in range(2, d+1):
                nxt[i][j] += prev[i-1][j] + prev[i][j]/2
    for i in range(d+1):
        for j in range(d+1):
            nxt[i][j] /= k+1
    U.append(nxt)

# Route 2: solve the gauged triangular system by successive integrations,
# then multiply by (1-u)^(-1/2). This includes every initial-value column.
b = [F(comb(2*k, k), 4**k) for k in range(K+1)]
g = [-sum((b[r]/2**(k-r) for r in range(k+1)), F(0))
     for k in range(K)]
comparisons = 0
nonzero_coefficients = 0
for j in range(d+1):
    V = [[F(0) for k in range(K+1)] for i in range(d+1)]
    for i in range(d+1):
        V[i][0] = F(int(i == j))
    for k in range(1, K+1):
        V[1][k] = g[k-1]*V[0][0]/k
    for i in range(2, d+1):
        for k in range(1, K+1):
            V[i][k] = sum(V[i-1][:k], F(0))/k
    for i in range(d+1):
        for k in range(1, K+1):
            value = V[i][k] if i == 0 else sum(
                (b[r]*V[i][k-r] for r in range(k+1)), F(0))
            assert value == U[k][i][j], (i, j, k)
            assert (value * 4**k * delta(d, k)).denominator == 1
            comparisons += 1
            nonzero_coefficients += bool(value)
report["fundamental_matrix_route_comparisons"] = comparisons
report["nonzero_coefficients_compared"] = nonzero_coefficients
report["taylor_depth"] = d
report["max_taylor_order"] = K
print("Exact denominator and two-route Taylor checks passed", flush=True)

# Enlarged system: symbolic differentiation of the full binomial-log sum.
z, w = sp.symbols("z w")
identities = 0
for h in range(5):
    bb = 3
    gg = sp.symbols("g0:"+str(bb+h+1))
    def derivative(expr):
        out = sp.diff(expr, z)
        for j in range(1, bb+h+1):
            gp = (2*gg[0]/(z*(1-z)) if j == 1 else -gg[j-1]/z)-gg[j]/(2*z)
            out += sp.diff(expr, gg[j])*gp
        return sp.simplify(out)
    for q in range(h+1):
        y = []
        for i in range(bb+h+1):
            if i < h-q:
                y.append(sp.Integer(0))
            elif i <= h:
                y.append(sp.factorial(q)*(-sp.log(z))**(i+q-h)/sp.factorial(i+q-h))
            else:
                y.append(sum(sp.binomial(q,p)*(-sp.log(z))**(q-p)
                             *(-1)**p*sp.rf(i-h,p)*gg[i-h+p]
                             for p in range(q+1)))
        for i in range(bb+h+1):
            if i == 0:
                rhs = 0
            elif i <= h:
                rhs = -y[i-1]/z
            elif i == h+1:
                rhs = 2*y[i-1]/(z*(1-z))-y[i]/(2*z)
            else:
                rhs = -y[i-1]/z-y[i]/(2*z)
            assert sp.simplify((derivative(y[i])-rhs).subs(gg[0], 1)) == 0
            identities += 1
report["enlarged_system_symbolic_identities"] = identities
assert sp.simplify(w * 2/(w*w*(1-w*w)) * 2*w - 4/(1-w*w)) == 0
assert sp.residue(-2/w, w, 0) == -2
assert sp.residue(4/(1-w*w), w, 1) == -2
report["pullback_and_residue_checks"] = 3

# Direct differentiation versus the residue coefficient in the chain lemma.
residue_checks = 0
for a, point in [(-2/w, 0), (4/(1-w*w), 1)]:
    for D in range(6):
        top = sum((j+1)*w**j for j in range(D+1))
        lower = w**(D+2) + 3*w + 2
        current_top, current_lower = top, lower
        for _ in range(D+1):
            current_top, current_lower = (
                sp.diff(current_top, w),
                sp.diff(current_lower, w) + a*current_top)
        formula = sp.diff(lower,w,D+1) + sum(
            sp.diff(a*sp.diff(top,w,j),w,D-j) for j in range(D+1))
        assert sp.cancel(current_top) == 0
        assert sp.cancel(current_lower-formula) == 0
        assert sp.residue(formula,w,point) == sp.diff(top,w,D)*sp.residue(a,w,point)
        residue_checks += 1
report["descending_residue_identities"] = residue_checks

# Exact rational certificate for the intentionally coarse constant 0.21.
# log(39/10)=2*atanh(29/49); the omitted series is positive.
x = F(29,49)
log_lower = 2*sum((x**(2*j+1)/(2*j+1) for j in range(16)),F(0))
log_upper = log_lower + 2*x**33/(33*(1-x*x))
assert F(136,100) < log_lower < log_upper < F(1361,1000)
sqrt_upper = F(1167,1000)
assert sqrt_upper**2 > F(34,25)
r, kappa = F(39,10), F(529,50)
constant_lower = r*log_lower/((r*r*log_lower/2+kappa)*sqrt_upper)
assert constant_lower > F(21,100)
assert F(9,25)*(kappa-2*r) > 1
report["constant_certificate"] = {
    "log_lower": str(log_lower), "log_upper": str(log_upper),
    "sqrt_upper": str(sqrt_upper), "constant_lower": str(constant_lower),
    "constant_lower_decimal_diagnostic": float(constant_lower),
    "greater_than_0_21_exact": True, "rank_slope_greater_than_1_exact": True
}
report["status"] = "Exact finite supporting checks; same assistant; not formal or independent mathematical verification."
(HERE/"local-system-results.json").write_text(json.dumps(report,indent=2)+"\n")
print(json.dumps({k:v for k,v in report.items() if k != "constant_certificate"},indent=2))

