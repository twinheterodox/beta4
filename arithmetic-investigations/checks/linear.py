"""Recompute the fixed linear-form certificate from its definitions.
Adapted from the separate integer-pattern/Taylor-recurrence implementation.
No frozen workspace, private path, saved large coefficient row, or generator import.
"""
from pathlib import Path
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations
import sys,json,math
from flint import fmpq,arb,ctx
sys.set_int_max_str_digits(0)
ctx.prec=256
ROOT=Path(__file__).resolve().parents[1]/'evidence/linear'

def events(cutoff):
    ds=[Q(j,116) for j in range(29)]
    offsets=set([(v,Q(0)) for v in ds]+[(1-v,Q(0)) for v in ds]
                +[(Q(-1),Q(-1,2)),(Q(2),Q(-1,2))])
    out={Q(k) for k in range(1,cutoff+1)}
    for (s,t),(ss,tt) in combinations(sorted(offsets),2):
        slope,shift=ss-s,t-tt
        if not slope: continue
        for k in range(math.floor(slope-shift)-1,
                       math.ceil(cutoff*slope-shift)+2):
            x=(shift+k)/slope
            if 1<x<cutoff: out.add(x)
    return sorted(out)

def envelope(x):
    a,b=x.numerator,x.denominator
    den=116*b; base=0; jumps={}
    def add(num,jump):
        residue=num%den
        if residue: jumps[residue]=jumps.get(residue,0)+jump
    for j in range(29):
        lo,hi=j*a,(116-j)*a
        base+=(-lo)//den-(-hi)//den
        add(lo,1);add(hi,-1)
    base-=2*((116*a+58*b)//den-(-232*a+58*b)//den)
    add(-116*a-58*b,-2);add(232*a-58*b,2)
    cap=116*(a-b) if a<2*b else den
    current=maximum=base
    for position,jump in sorted(jumps.items()):
        if position>cap: break
        current+=jump;maximum=max(maximum,current)
    return max(0,maximum+1)-16*(a//b)

def periodic_check():
    data=json.loads((ROOT/'periodic-arithmetic.json').read_text())
    cells=data['initial_records']+data['period_records']
    points=events(118)
    assert len(cells)==len(points)-1
    column_saving=Q(1)  # integral of 2/x^2 from 2 to infinity
    for (l,r,e),a,b in zip(cells,points,points[1:]):
        assert Q(l)==a and Q(r)==b
        x=(a+b)/2
        assert e==envelope(x)
        if b<=2:
            # Recover max(0,M+1); on this range it is positive.
            M=e+16*(x.numerator//x.denominator)-1
            gain=max(0,M+1)-max(0,M-1)
            column_saving+=gain*(1/a-1/b)
    assert column_saving==Q(223,116)
    total=0;scale=2**data['rounding_bits']
    for k in range(data['cycles']):
        recs=data['period_records']+(data['initial_records'] if k==0 else [])
        for l,r,e in recs:
            l,r=Q(l),Q(r)
            shift=116*k if l>=2 else 0
            term=e*(1/(l+shift)-1/(r+shift))*scale
            total+=ceildiv(term.numerator,term.denominator)
    upper=Q(total,scale)+Q(32,data['cutoff'])+Q(8,data['cutoff']-1)
    assert Q(total,scale)==Q(data['body_upper'])
    assert upper==Q(data['integral_upper'])
    profile=json.loads((ROOT/'profile-saving-6.json').read_text())
    assert [(Q(v['left']),Q(v['right'])) for v in profile['records']]==list(zip(events(6),events(6)[1:]))
    return dict(envelope_cells=len(cells),profile_cells=len(profile['records']),
                arithmetic_upper=str(upper),column_saving=str(column_saving))

def phase_check():
    ds=[Q(j,116) for j in range(29)]
    curvature=Q(2,3)-20+sum((1/(1+d)-1/(Q(21,10)-d) for d in ds),Q(0))
    assert curvature<0
    def convert(x):return arb(fmpq(x.numerator,x.denominator))
    def phase(x):
        y=convert(x)
        return 2*((y+2)*(y+2).log()-(y-1)*(y-1).log())-sum(
            ((y+1-convert(d))*(y+1-convert(d)).log()
             -(y+convert(d))*(y+convert(d)).log() for d in ds),arb(0))
    def slope(x):
        y=convert(x)
        return 2*((y+2)/(y-1)).log()-sum(
            (((y+1-convert(d))/(y+convert(d))).log() for d in ds),arb(0))
    left,right=Q(10001,10000),Q(11,10)
    assert slope(left)>0 and slope(right)<0
    for _ in range(160):
        mid=(left+right)/2
        if slope(mid)>0:left=mid
        else:
            assert slope(mid)<0
            right=mid
    x=(left+right)/2
    upper=phase(x)+abs(slope(x))*convert(max(x-1,Q(11,10)-x))
    saved=Q(json.loads((ROOT/'quantitative-bounds.json').read_text())['inherited_tail_phase_upper'])
    assert upper<convert(saved)
    g=convert(Q(22,1)/(Q(11,10)+1))-2*convert((Q(11,10)+2)/(Q(11,10)-1)).log()
    assert g>0
    assert Q(24,3)-Q(88,5)<-9
    return dict(critical_lower=str(left),critical_upper=str(right),
                tangent_upper_enclosure=str(upper),saved_upper=str(saved),
                curvature_upper=str(curvature))
def ceildiv(a,b):return -((-a)//b)

def integer_patterns(x):
    a,b=x.numerator,x.denominator;D=116*b;cap=min(D,116*(a-b))
    cuts={0,cap}
    residues=[j*a for j in range(29)]+[(116-j)*a for j in range(29)]
    residues.extend([232*a+58*b,-116*a+58*b])
    cuts.update(v%D for v in residues if 0<v%D<cap)
    cuts=sorted(cuts)
    samples=[2*v for v in cuts]+[l+r for l,r in zip(cuts,cuts[1:])]
    result=set()
    for v in samples:
        if v>=2*D:continue
        orders=tuple(sum(2*j*a<=v+2*D*k<=2*(116-j)*a for j in range(29))
                     for k in range((232*a-v)//(2*D)+1))
        if sum(bool(z) for z in orders)<2:continue
        roots=tuple(range(ceildiv(v-464*a-D,2*D),(v+232*a-D)//(2*D)+1))
        result.add((orders,roots))
    return result

@lru_cache(None)
def principal_parts(t):
    orders,roots=t;out=[]
    for a,B in enumerate(orders):
        c0=fmpq(1)
        for k in roots:c0*=fmpq(-2*a-2*k-1,2)**2
        for j,b in enumerate(orders):
            if j!=a:c0/=fmpq(j-a)**b
        coeff=[c0]
        power_sums=[fmpq(0)]
        for k in range(1,B):
            z=sum((2/fmpq(-2*a-2*r-1,2)**k for r in roots),fmpq(0))
            z-=sum((fmpq(b)/fmpq(j-a)**k for j,b in enumerate(orders) if j!=a),fmpq(0))
            power_sums.append((-1)**(k-1)*z)
            coeff.append(sum((power_sums[j]*coeff[k-j] for j in range(1,k+1)),fmpq(0))/k)
        out.append(tuple(coeff))
    return tuple(out)

@lru_cache(None)
def functionals(t,h):
    assert h in (0,1)
    orders,roots=t;series=principal_parts(t);cs=[]
    for a,B in enumerate(orders):
        c=[fmpq(0)]
        for i in range(1,B+1):
            k=B-i
            z=series[a][k] if h==0 else -a*series[a][k]+(series[a][k-1] if k else 0)
            c.append(z)
        cs.append(c)
    bb=sum(((-1)**a*c[1] for a,c in enumerate(cs)),fmpq(0))
    qq=fmpq(0)
    for v in range(len(orders)-1):
        for j in range(v+1,len(orders)):
            for i in range(1,len(cs[j])):
                qq+=(-1)**v*(-i)*cs[j][i]/fmpq(2*(j-v)-1,2)**(i+1)
    return bb,qq

def profile_check():
    d=json.loads((ROOT/'profile-saving-6.json').read_text());saving=Q(0);checks=0;gains={0:0,1:0,2:0}
    left=Q(d['window'][0])
    for rec in d['records']:
        l,r=Q(rec['left']),Q(rec['right']);assert l==left and l<r;left=r
        patterns=integer_patterns((l+r)/2)
        M=max(sum(t[0])-2*len(t[1]) for t in patterns)
        assert M==rec['M']
        near={t for t in patterns if sum(t[0])-2*len(t[1])>=M-1}
        saved={(tuple(t['orders']),tuple(t['half_root_k'])) for t in rec['templates']}
        assert near==saved
        g=rec['saving'];gains[g]+=1
        if not g:continue
        assert M>=g+1;k=fmpq(rec['kappa'])
        for t in patterns:
            F=sum(t[0])-2*len(t[1])
            for h in range(g-(M-F)):
                bb,qq=functionals(t,h);assert qq==k*bb;checks+=1
        saving+=g*(1/l-1/r)
    assert left==Q(d['window'][1]) and saving==Q(d['saving'])
    # All saved individual functional values also agree with the fresh recurrence.
    for e in d['functionals']:
        t=tuple(e['orders']),tuple(e['half_root_k'])
        bb,qq=functionals(t,e['h']);assert bb==fmpq(e['B']) and qq==fmpq(e['Q'])
    return dict(intervals=len(d['records']),template_equalities=checks,
                distinct_functionals=len(d['functionals']),gain_cells=gains,saving=str(saving))

def bounds_check():
    d=json.loads((ROOT/'quantitative-bounds.json').read_text())
    A=fmpq(d['arithmetic_upper']);S=fmpq(d['local_saving'])
    psi=6*arb(fmpq(3,2)).log()
    for j in range(29):
        x=fmpq(1,2)-fmpq(j,116);psi-=2*arb(x)*arb(x).log()
    H=arb(A)-arb(fmpq(223,116))-10*arb(2).log()+psi
    C=-arb(A)-12*arb(2).log()-arb(fmpq(d['inherited_tail_phase_upper']))+arb(S)
    assert H<arb(fmpq(24522,1000)) and H/5<arb(fmpq(491,100))
    assert C>arb(fmpq(60384,10000)) and C-arb(fmpq(491,100))>arb(fmpq(112,100))
    return dict(H_enclosure=str(H),C_enclosure=str(C),eleven_height_rate='4.91',cleared_decay_rate='1.12')

if __name__=='__main__':
    from argparse import ArgumentParser
    parser=ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result={'status':'PASS','periodic':periodic_check()}
    print('envelope and period certificate passed',flush=True)
    result['profile']=profile_check()
    print('local cancellation templates passed',flush=True)
    result['phase']=phase_check()
    result['bounds']=bounds_check()
    result['scope']='Recomputed finite envelope, exact local functionals, phase enclosure and constants. General prime sums, tail estimates and lattice theorems require mathematical proof.'
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='phase'},indent=2))
