from pathlib import Path
from math import comb, isqrt
import json
from flint import nmod_poly as Polynomial

ROOT = Path(__file__).resolve().parent
reports = []
for prime in (10007,10009):
    assert all(prime%d for d in range(2,isqrt(prime)+1))
    x = Polynomial([0,1],prime)
    orders = [min(a+1,97-a,31) for a in range(97)]
    numerator = x
    for zero in range(1,240,2):
        numerator *= (x*x-zero*zero)**3
    denominator = Polynomial([1],prime)
    for a,order in enumerate(orders):
        denominator *= (x-(96-2*a))**order
    rows = []
    for a,order in enumerate(orders):
        centre = 96-2*a
        pole = (x-centre)**order
        quotient,remainder = divmod(denominator,pole)
        assert not remainder
        gcd,inverse,_ = quotient.xgcd(pole)
        assert gcd == Polynomial([1],prime)
        local_polynomial = (numerator*inverse)%pole
        # Translate the polynomial residue using the binomial theorem.
        taylor = [sum(int(local_polynomial[j])*comb(j,k)*pow(centre,j-k,prime)
                      for j in range(k,order))%prime for k in range(order)]
        rows.append([0]+[taylor[order-i] for i in range(1,order+1)])
    constant = 0
    beta = [0]*34
    for a,row in enumerate(rows):
        for i,coefficient in enumerate(row[1:],1):
            weight = (-1)**a*i*(i+1)*coefficient
            beta[i+2] = (beta[i+2]+weight)%prime
            harmonic = sum((-1)**h*pow(2*h+1,-i-2,prime) for h in range(a))
            constant = (constant-weight*harmonic)%prime
    assert all(beta[s]==0 for s in range(3,34,2))
    reference = json.loads((ROOT/f'replay/checks/fifteen-seed-rows-p{prime}.json').read_text())
    assert rows == reference
    assert constant == {10007:3738,10009:5333}[prime]
    reports.append(dict(prime=prime, coefficients=sum(orders), residue=constant))
(ROOT/'independent-seed.json').write_text(json.dumps(reports,indent=2)+'\n')
print(json.dumps(reports,indent=2))
