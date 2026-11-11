"""Fresh exact-Q checks using global Hasse derivatives and long division.

Does not import or call manuscript checkers. No rational-series library is used.
"""
from pathlib import Path
import json,math,time
from flint import fmpz_poly,fmpq

HERE=Path(__file__).resolve().parent

def product(items):
    while len(items)>1:
        items=[items[i]*items[i+1] if i+1<len(items) else items[i]
               for i in range(0,len(items),2)]
    return items[0]

def hasse_values(poly,points,count):
    out=[[0]*count for _ in points]
    for k in range(count):
        for j,x in enumerate(points):out[j][k]=int(poly(x))
        coefficients=[int(c) for c in poly.derivative().coeffs()]
        assert all(c%(k+1)==0 for c in coefficients)
        poly=fmpz_poly([c//(k+1) for c in coefficients])
    return out

def convolution(a,b,cap):
    out=[fmpq(0) for _ in range(cap)]
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b[:cap-i]):
                if y:out[i+j]+=x*y
    return out

def vp(value,p):
    if not value:return math.inf
    a,b=abs(int(value.numerator)),int(value.denominator);v=0
    while a%p==0:a//=p;v+=1
    while b%p==0:b//=p;v-=1
    return v

def mod(value,p):
    assert int(value.denominator)%p
    return int(value.numerator)%p*pow(int(value.denominator),-1,p)%p

def decompose(n,h,blocks,trim,maxq=3):
    orders=[sum(m for a,b,m in blocks if a<=j<=b) for j in range(n+1)]
    poles=[n-2*j for j in range(n+1)]
    edge=n+2*h+1-2*trim
    zeros=list(range(-edge,edge+1,2))
    num=fmpz_poly([0,1])*product([fmpz_poly([-z,1]) for z in zeros])**3
    den=product([fmpz_poly([-a,1])**r for a,r in zip(poles,orders)])
    assert num.degree()<den.degree()
    maxr=max(orders)
    nj=hasse_values(num,poles,maxq*maxr)
    dj=hasse_values(den,poles,(maxq+1)*maxr)
    states={q:dict(beta=[fmpq(0)]*((maxr+2)*q+1),constant=fmpq(0),partials=[])
            for q in range(1,maxq+1)}
    tail=[fmpq(0)]*((maxr+2)*maxq+1)
    for j,(a,r) in enumerate(zip(poles,orders)):
        cap=maxq*r
        assert all(dj[j][k]==0 for k in range(r)) and dj[j][r]
        d=dj[j][r:r+cap];v=nj[j][:cap]
        base=[]
        for k in range(cap):
            value=fmpq(v[k])-sum((d[t]*base[k-t] for t in range(1,k+1)),fmpq(0))
            base.append(value/d[0])
        power=[fmpq(1)]+[fmpq(0)]*(cap-1)
        for q,state in states.items():
            power=convolution(power,base,cap)
            size=q*r
            aa=[fmpq(0)]+[power[size-i] for i in range(1,size+1)]
            state['partials'].append(aa)
            for i in range(1,size+1):
                s=i+2*q
                weight=(-1)**j*math.prod(range(i,s))*aa[i]
                state['beta'][s]+=weight
                state['constant']-=weight*tail[s]
        for s in range(3,len(tail)):
            tail[s]+=fmpq((-1)**j,(2*j+1)**s)
    point_checks=0;reflections=0
    for q,state in states.items():
        assert all(not state['beta'][s] for s in range(2*q+1,len(state['beta']),2))
        for j,aa in enumerate(state['partials']):
            assert all(aa[i]==(-1)**i*state['partials'][n-j][i] for i in range(1,len(aa)))
            reflections+=len(aa)-1
        for x in (fmpq(1,3),fmpq(2,7),fmpq(17,11)):
            actual=(num(x)/den(x))**q
            total=sum((aa[i]/(x-a)**i for a,aa in zip(poles,state['partials'])
                       for i in range(1,len(aa))),fmpq(0))
            assert total==actual
            point_checks+=1
    return dict(n=n,h=h,blocks=blocks,orders=orders,poles=poles,zeros=zeros,
                states=states,point_checks=point_checks,reflections=reflections)

def primes_through(n):
    return [p for p in range(2,n+1) if all(p%d for d in range(2,math.isqrt(p)+1))]

def exponent_bound(n,p):
    k=0;value=1
    while value*p<=n:value*=p;k+=1
    return k

def clearing(data):
    n,h=data['n'],data['h'];orders=data['orders'];S=max(orders)
    zeros=range(-n-h-1,h+1);M=4*n+2*h+3;Z=3*len(zeros)
    multiplier=fmpq(1);exp={}
    for p in primes_through(n):
        lam=exponent_bound(n if p==2 else M,p)
        costs=[]
        for j,B in enumerate(orders):
            if p==2:
                nu=1-Z-sum(orders[k]*vp(fmpq(k-j),2) for k in range(n+1) if k!=j)
                cost=-nu+S*lam
            elif p*p>M:
                eps=3*sum((2*j+1+2*l)%p==0 for l in zeros)+(int((n-2*j)%p==0))
                eps-=sum(orders[k] for k in range(n+1) if k!=j and (k-j)%p==0)
                cost=B+2-eps
            else:
                nu=3*sum(vp(fmpq(2*j+1+2*l),p) for l in zeros)
                nu-=sum(orders[k]*vp(fmpq(k-j),p) for k in range(n+1) if k!=j)
                cost=-nu+(S+2)*lam
            costs.append(cost)
        exp[p]=int(max(costs));multiplier*=fmpq(p)**exp[p]
    # F(x)=2^{-C} R(t), x=2t+n+1. The t derivative adds 2^{2q}.
    C=sum(orders)-3*len(zeros)
    checked=0;groups=0
    for q,state in data['states'].items():
        scalar=multiplier**q*fmpq(2)**(q*C+2*q)
        for x in [state['constant'],*state['beta']]:
            assert (scalar*x).denominator==1
            checked+=1
        for ell in range(1,n+1):
            point=n+1-2*ell
            total=sum((math.prod(range(i,i+2*q))*aa[i]/fmpq(point-a)**(i+2*q)
                       for a,aa in zip(data['poles'],state['partials'])
                       for i in range(1,len(aa))),fmpq(0))
            assert total==0;groups+=1
    return dict(n=n,h=h,exponents=exp,cleared_coefficients_checked=checked,
                grouped_zero_identities_checked=groups)

def dilation(seed,big,p):
    n,h=seed['n'],seed['h'];blocks=seed['blocks'];L=n+2*h
    W=sum((b-a)*m for a,b,m in blocks);S=sum(m for a,b,m in blocks)
    A=W-3*L+S-1;GN=W-3*L
    assert p>n+L+1
    output=[]
    for q,state in big['states'].items():
        tau=q*(A+2);selected=other=0
        for j,aa in enumerate(state['partials']):
            for i in range(1,len(aa)):
                if j%p==0:
                    lhs=aa[i]*fmpq(p)**(q*A-i)
                    rhs=(-1)**q*seed['states'][q]['partials'][j//p][i]
                    assert vp(lhs,p)>=0 and mod(lhs-rhs,p)==0
                    selected+=1
                else:
                    assert vp(aa[i],p)>=-q*GN+i
                    other+=1
        assert all(vp(b*fmpq(p)**tau,p)>=1 for b in state['beta'])
        lhs=state['constant']*fmpq(p)**tau
        rhs=(-1)**(q+(p-1)//2)*seed['states'][q]['constant']
        assert vp(lhs,p)>=0 and mod(lhs-rhs,p)==0
        output.append(dict(seed_n=n,seed_h=h,prime=p,power=q,selected=selected,
            other_pole_bounds=other,constant_residue=mod(lhs,p),
            seed_constant_residue=mod(seed['states'][q]['constant'],p),
            beta_coefficients_checked=len(state['beta'])))
    return output

if __name__=='__main__':
    start=time.monotonic();report=dict(clearing=[],dilation=[],point_checks=0,reflection_checks=0)
    geometries=[(2,2,[(0,2,17)],[11,13]),(4,4,[(0,4,9),(1,3,8)],[19]),
                (6,7,[(0,6,9),(1,5,8)],[29])]
    for n,h,blocks,primes in geometries:
        trimmed=decompose(n,h,blocks,True)
        full=decompose(n,h,blocks,False)
        report['clearing'].append(clearing(full))
        for obj in (trimmed,full):
            report['point_checks']+=obj['point_checks'];report['reflection_checks']+=obj['reflections']
        print('Completed seed geometry',n,h,flush=True)
        for p in primes:
            large=decompose(n*p,h*p,[(a*p,b*p,m) for a,b,m in blocks],False)
            report['point_checks']+=large['point_checks'];report['reflection_checks']+=large['reflections']
            report['dilation']+=dilation(trimmed,large,p)
            print('Completed dilation',n,'prime',p,'all powers 1,2,3',flush=True)
            (HERE/'fresh-rational.json').write_text(json.dumps(report,indent=2)+'\n')
    report['seconds']=time.monotonic()-start
    (HERE/'fresh-rational.json').write_text(json.dumps(report,indent=2)+'\n')
    print('All rational checks passed in',report['seconds'],'seconds',flush=True)
