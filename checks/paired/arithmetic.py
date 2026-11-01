"""An upper ceiling for the source-style termwise q-divisor estimates.

This is NOT an upper bound on actual common content. Ignoring all P costs,
it bounds the saving the specified forward-term/dual-row minimum estimates
can certify, using selected terms to upper-bound the full minima.
"""
from pathlib import Path
from fractions import Fraction as F
from math import floor
import json

HERE=Path(__file__).resolve().parent
INTEGER_EXP={8:1,9:2,11:2,12:1,14:-1,16:-1}
SHIFT_EXP={11:1,12:2,13:2,14:1,16:-1,17:-2,18:-2,19:-1}

def fterm(x,r,phase):
    v=sum(floor((r+A)*x-phase)-floor((r+A-10)*x-phase)-floor(10*x) for A in (12,13,14))
    v+=floor(15*x)-floor((1+r)*x)-floor((14-r)*x)
    v+=floor(15*x)-floor(r*x)-floor((15-r)*x)
    return v

def dlead(x,r,phase):
    numerator=floor((16+r)*x-phase)-floor((2+r)*x-phase)+floor((17+r)*x-phase)-floor((1+r)*x-phase)
    denominator=sum(floor(A*x) for A in (r,11-r,r-1,12-r,r-2,10-r))
    return numerator-denominator

def qceiling(x,phase):
    f=min(fterm(x,r,phase) for r in range(15))
    lead=min(dlead(x,r,phase) for r in range(3,10))
    lam=sum(e*floor(A*x) for A,e in INTEGER_EXP.items())+sum(e*floor(A*x-phase) for A,e in SHIFT_EXP.items())
    dual=lam+2*lead-1
    return max(f,dual),f,dual

if __name__=='__main__':
    # All functions use integer slopes <=28 and shifts 0,1/4,3/4.
    grid=sorted({F(0),F(1)}|{F(4*m+l,4*d) for d in range(1,29)
                   for l in (0,1,3) for m in range(d+1) if 0<=4*m+l<=4*d})
    cells=[]
    for a,b in zip(grid,grid[1:]):
        mid=(a+b)/2
        vals=[qceiling(mid,phase) for phase in (F(1,4),F(3,4))]
        assert vals==[qceiling((2*a+b)/3,phase) for phase in (F(1,4),F(3,4))]
        assert vals==[qceiling(mid+1,phase) for phase in (F(1,4),F(3,4))]
        assert all(v[0]>=0 for v in vals)
        cells.append((a,b,vals))
    assert cells[0][2][0][0]==cells[0][2][1][0]==0
    cutoff=40
    integrals=[F(0),F(0)]
    for shift in range(cutoff):
        for a,b,vals in cells:
            weight=1/(a+shift)-1/(b+shift) if a+shift else F(0)
            for i in range(2):integrals[i]+=vals[i][0]*weight
    maxima=[max(vals[i][0] for _,_,vals in cells) for i in range(2)]
    upper=[integrals[i]+F(maxima[i],cutoff) for i in range(2)]
    out={'scope':'Upper ceiling for the specified termwise certificate only, not for primitive content.',
         'period_cells':len(cells),'cutoff':cutoff,'integrals':list(map(str,integrals)),
         'period_maxima':maxima,'upper_bounds':list(map(str,upper)),
         'decimal_upper_bounds':[float(v) for v in upper],
         'combined_upper':str(sum(upper)),'combined_upper_decimal':float(sum(upper))}
    (HERE/'arithmetic-result.json').write_text(json.dumps(out,indent=2)+'\n')
    print('q-only divisor ceiling:',out['decimal_upper_bounds'],'sum',out['combined_upper_decimal'],flush=True)
    print('period cells',len(cells),'maxima',maxima,flush=True)
