from fractions import Fraction as Rational
from itertools import combinations
from pathlib import Path
import hashlib
import json
import subprocess

ROOT = Path(__file__).resolve().parent
offsets = [(j,0) for j in range(31)] + [(96-j,0) for j in range(31)]
offsets += [(-72,-48),(168,-48)]
# Reduce all affine crossings to arithmetic progressions in their numerators.
progressions = set()
for (s,t),(u,v) in combinations(offsets,2):
    difference = s-u
    if difference:
        sign = 1 if difference>0 else -1
        progressions.add((abs(difference), ((v-t)*sign)%96))
events = set(map(Rational, range(1,401)))
for denominator, residue in progressions:
    first = denominator + (residue-denominator)%96
    events.update(Rational(numerator,denominator)
                  for numerator in range(first, 400*denominator+1, 96))
events = sorted(events)
assert len(events) == 38372
midpoints = [(a+b)/2 for a,b in zip(events,events[1:])]
data = ''.join(f'{x.numerator} {x.denominator}\n' for x in midpoints)
result = subprocess.run([str(ROOT/'independent_floor')], input=data,
                        text=True, capture_output=True, check=True)
values = [int(value) for value in result.stdout.split()]
assert len(values) == len(midpoints)
scale = 2**96
total = 0
for a,b,value in zip(events,events[1:],values):
    term = Rational(value,16)*(1/a-1/b)*scale
    total += -((-term.numerator)//term.denominator)
body = Rational(total,scale)
tail = Rational(36,400)+Rational(221,16*2*399)
assert body+tail < Rational(41,5)
# Compare the independently enumerated partition and maxima only after computing them.
reference = ROOT/'replay/checks/fifteen-arithmetic-midpoints.tsv'
for x,value,line in zip(midpoints,values,reference.read_text().splitlines()):
    a,b,c,d = map(int,line.split())
    assert x == Rational(a,b)
    assert Rational(value,16) == Rational(c,d)
report = dict(event_points=len(events), intervals=len(values),
              body_upper=str(body), total_upper=str(body+tail),
              total_decimal=float(body+tail), rounding_bits=96,
              negative_cells=sum(v<0 for v in values),
              event_hash=hashlib.sha256('\n'.join(map(str,events)).encode()).hexdigest())
(ROOT/'independent-arithmetic.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
