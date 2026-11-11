from fractions import Fraction
from pathlib import Path
import json
from flint import arb, ctx, fmpq, fmpq_poly

ROOT = Path(__file__).resolve().parent
ctx.prec = 512
z = fmpq_poly([0,1])
numerator = (z+fmpq(5,2))**3
denominator = z**3
for j in range(31):
    numerator *= z+fmpq(72+j,96)
    denominator *= z+fmpq(168-j,96)
polynomial = numerator-denominator
assert polynomial.degree() == 33
assert all(polynomial[j]>0 for j in range(3))
assert all(polynomial[j]<0 for j in range(3,34))
left = fmpq(378667493,500000000)
right = fmpq(757334987,1000000000)
r = fmpq(3,4)
assert polynomial(left-r)>0 and polynomial(right-r)<0
advertised_bracket = [str(left),str(right)]
for _ in range(80):
    middle = (left+right)/2
    if polynomial(middle-r)>0:
        left = middle
    else:
        right = middle
point = arb((left+right)/2)
point += arb(0, arb((right-left)/2))
def H(x):
    return x*x.log()
phase = 3*(H(point+1+arb(r))-H(point-arb(r)))
for j in range(31):
    delta = arb(fmpq(j,96))
    phase -= H(point+1-delta)-H(point+delta)
maximum_bound = fmpq(-18835112315041,10**12)
assert phase < arb(maximum_bound)
gamma = fmpq(341,16)
tail = 3*fmpq(5,2)*3/(3-r)-3*gamma/4
assert tail == fmpq(-383,64)
cleared = arb(fmpq(41,5))+15*arb(2).log()-arb(fmpq(1883,100))
assert cleared < arb(fmpq(-23,100))
arithmetic = json.loads((ROOT/'independent-arithmetic.json').read_text())
total = Fraction(arithmetic['total_upper'])
sharper = arb(fmpq(total.numerator,total.denominator))+15*arb(2).log()+phase
report = dict(polynomial_degree=33, positive_coefficients=[0,1,2],
              negative_coefficients=list(range(3,34)),
              maximum_point_bracket=advertised_bracket,
              refined_bracket=[str(left),str(right)],
              phase_interval=phase.str(50), maximum_upper=str(maximum_bound),
              simplified_cleared_exponent=str(cleared),
              full_certificate_cleared_exponent=str(sharper), tail_slope=str(tail))
(ROOT/'independent-phase.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
