"""Independent exact reconstructions of formulas newly printed in the supplement.

These bounded checks supplement the written proofs; they are not asymptotic
proofs by extrapolation. They do not modify the frozen research evidence.
"""
from fractions import Fraction as F
from math import factorial, gcd, lcm
from pathlib import Path
import json,datetime,sys
import sympy as s
from flint import arb,ctx,fmpq,fmpq_mat
x,t,X=s.symbols('x t X');counts={}
def count(name,n=1):counts[name]=counts.get(name,0)+n
def chi(j):return 0 if j%2==0 else (-1)**((j-1)//2)
def vp(q,p):
    q=F(q)
    if not q:return 10**9
    a,b=abs(q.numerator),q.denominator;v=0
    while a%p==0:v+=1;a//=p
    while b%p==0:v-=1;b//=p
    return v

def mod(q,p):
    q=F(q);return q.numerator*pow(q.denominator,-1,p)%p

def detmod(a,p):
    a=[row[:] for row in a];out=1
    for i in range(len(a)):
        j=next(j for j in range(i,len(a)) if a[j][i]%p)
        if j!=i:a[i],a[j]=a[j],a[i];out=-out
        v=a[i][i]%p;out=out*v%p
        for j in range(i+1,len(a)):
            z=a[j][i]*pow(v,-1,p)%p
            for k in range(i+1,len(a)):a[j][k]=(a[j][k]-z*a[i][k])%p
    return out%p

def primitive(poly):
    cs=s.Poly(poly,X).all_coeffs();d=s.ilcm(*(s.denom(c) for c in cs))
    zz=[int(d*c) for c in cs];g=0
    for c in zz:g=gcd(g,abs(c))
    return s.expand(poly*d/g)

# The resolvent is recovered by differentiating the rational weight, rather
# than trusting the recorded pole formula's normalization.
a=s.symbols('a',positive=True)
assert s.factor(s.diff(x**4/(a*a+x*x),x,3)/3-8*a**4*x*(a*a-x*x)/(a*a+x*x)**4)==0
count('rational_weight_derivative')
for K,h in ((3,1),(5,2),(7,3)):
    nodes=[s.Rational(2*j+1,4) for j in range(K)]
    D=s.prod(t+z*z for z in nodes);interp=s.interpolate([(-z*z,1/z) for z in nodes],t)
    for ell in range(2*h-1):
        quotient,remainder=s.div(interp*t**ell,D,t)
        residues=[s.cancel(((-z*z)**ell)/(z*s.diff(D,t).subs(t,-z*z))) for z in nodes]
        assert s.cancel(interp*t**ell/D-quotient-sum(c/(t+z*z) for c,z in zip(residues,nodes)))==0
        assert sum(c*z**3 for c,z in zip(residues,nodes))==0
        count('squared_partial_fraction_and_period_cancellation')
    cK=s.Rational(4**K,(2*K-1)*factorial(K-1)**2)
    assert s.LC(s.Poly(interp,t))==cK
    count('squared_leading_coefficient')

def unsquared(M,h,K):
    nodes=[F(2*j+1,4) for j in range(K)]
    tail=[sum((F(256,(4*k+1+2*(j%2))**4) for k in range(j//2)),F()) for j in range(K)]
    rows=[]
    for ell in range(2*h-1):
        A=B=F()
        for j,z in enumerate(nodes):
            der=F((-1)**j*factorial(j)*factorial(K-1-j),2**(K-1))
            V=F(1)
            for m in range(1,M+1):V*=(-z-m)**4
            r=V*(-z)**ell/der
            A+=128*(-1)**j*r*z**4;B-=r*z**4*tail[j]
        rows.append((A,B))
    return rows
first=unsquared(0,1,6)[0]
assert primitive(s.Rational(first[0])*X+s.Rational(first[1]))==118682030250*X-117369943777
count('unsquared_primitive_fixture')
for n in (1,3,4,6):
    p=8*n+5;assert s.isprime(p)
    rows=unsquared(n,2*n,8*n+4)
    for theta in (F(0),F(2,3),F(7,2)):
        scaled=[p**4*(a*theta+b) for a,b in rows]
        assert all(vp(v,p)>=0 for v in scaled)
        matrix=[[mod(scaled[i+j],p) for j in range(2*n)] for i in range(2*n)]
        assert detmod(matrix,p)
        count('unsquared_prime_nonvanishing')

# Exact common-base coefficients, using pairs of rational numbers for complex
# multiplication; compare both forms of the determinant and the dyadic bound.
def mul(z,w):return (z[0]*w[0]-z[1]*w[1],z[0]*w[1]+z[1]*w[0])
def conv(a,b,N):return [sum((a[j]*b[k-j] for j in range(k+1)),F()) for k in range(N+1)]
N=24;hh=[F()]*(N+1);ll=[F()]*(N+1)
for j in range(1,N//2+1):ll[2*j]=F(1,j*2**j)
l2=conv(ll,ll,N)
for m,c in ((2,F(225,136)),(3,F(-25,34)),(5,F(9,68))):
    q=(F(1,2),F(-1,2));z=(F(1),F())
    for _ in range(m):z=mul(z,q)
    z=(-z[0],-z[1]);power=(F(1),F());blocks={r:[F()]*(N+1) for r in (2,3,4)}
    for j in range(1,N//m+1):
        power=mul(power,z)
        for r in blocks:blocks[r][m*j]=power[1]/j**r
    l3=conv(ll,blocks[3],N);l22=conv(l2,blocks[2],N)
    for j in range(N+1):hh[j]+=c*(blocks[4][j]+F(m,2)*l3[j]+F(m*m,12)*l22[j])
qrat=lambda z:fmpq(z.numerator,z.denominator)
def determinant(a):return F(str(fmpq_mat([[qrat(z) for z in r] for r in a]).det()))
for n in range(1,13):
    values=[]
    for target in (F(0),F(1)):
        border=[[target]+[-hh[j] for j in range(1,n+1)]]
        border += [[hh[i]]+[hh[i+j]-hh[i+j-1] for j in range(1,n+1)] for i in range(1,n+1)]
        cc=[target-sum(hh[1:j+1],F()) for j in range(2*n+1)]
        H=[[cc[i+j] for j in range(n+1)] for i in range(n+1)]
        v=determinant(border);assert determinant(H)==(-1)**n*v
        values.append(v);count('border_hankel_identity')
    B=n*(n+1)//2+(6+2*(n.bit_length()-1))*(n+1)
    assert min(vp(values[0],2),vp(values[1]-values[0],2))>=-B
    count('common_base_dyadic_bound')
p2=t**4+4;p3=t**6-4*t**3+8;p5=t**10-8*t**5+32
D=s.div(p2*p3*p5,(t*t+2*t+2)**2,t)[0]
assert s.rem(t**120-2**60,(2-t*t)*D,t)==0
count('common_base_radial_star')

# Check the endpoint system against independently imposing vanishing of the
# original rational source, including b<4 and singular-kernel possibilities.
for N in (1,2,3):
 for aa in range(1,2*N):
  b=2*N-aa
  for base in (1+x*x,1-x+x*x):
    dn=s.Poly(base**N,x);poly=s.prod(t-j for j in range(aa+1,2*N+1))
    K=s.Matrix(4,4,lambda r,k:sum(d*s.diff(t**r*poly,t,k).subs(t,ell[0])/factorial(k) for ell,d in dn.terms()))
    Flow=sum(s.Rational(chi(j),j**4)*x**j for j in range(1,aa+1))
    low=-s.Poly(dn.as_expr()*Flow,x)
    low=sum(c*x**j[0] for j,c in low.terms() if j[0]<=aa)
    # exp(ell*y) times V_b: direct derivatives at zero by binomial convolution.
    T=[s.diff(1/(2*s.cosh(t)),t,j).subs(t,0) for j in range(b)]
    PP=[s.Poly(t**r*poly,t) for r in range(4)]
    vv=[]
    for P in PP:
        val=sum(c*P.eval(j[0]) for j,c in s.Poly(low,x).terms())
        for ell,d in dn.terms():
            for degree,coef in P.terms():
                z=degree[0]
                for m in range(4,min(z,b+3)+1):val+=d*coef*s.binomial(z,m)*ell[0]**(z-m)*T[m-4]
        vv.append(val)
    for ker in K.row_join(-s.Matrix(vv)).nullspace():
        gs=list(ker[:4]);A=ker[4]
        unknown=s.symbols('r0:'+str(b));R=A*low+sum(unknown[j]*x**(aa+1+j) for j in range(b))
        G=R/base**N
        # Solve b initial conditions by jets of R(e^y), using the desired
        # G jets including V_b for orders above three.
        eq=[]
        for order in range(b):
            desired=0
            for ell,d in dn.terms():
                for k in range(order+1):
                    gv=gs[k] if k<4 else -A*T[k-4]
                    desired+=d*s.binomial(order,k)*ell[0]**(order-k)*gv
            eq.append(sum(c*j[0]**order for j,c in s.Poly(R,x).terms())-desired)
        sol=s.solve(eq,unknown);R=s.expand(R.subs(sol));G=R/base**N
        e4=G
        for _ in range(4):e4=s.cancel(x*s.diff(e4,x))
        num,den=s.fraction(s.cancel(A/(1+x*x)+e4/x));P=s.Poly(num,x)
        assert s.rem(P.as_expr(),x**aa,x)==0 and s.rem(P.as_expr(),(x-1)**b,x)==0
        assert s.cancel(G.subs(x,1)-gs[0])==0
        count('hermite_reconstruction')

# Modular coefficient congruences computed by integral q->t composition are
# already invariant at the first nonintegral coefficient. Check that leading
# coefficient directly from the independent Lambert divisor sums.
for p in (3,5,7,11,13,17,19,23,29,31):
    LF=sum(F(chi(d)*(-1)**(p//d-1),d**4) for d in s.divisors(p) if d%2)
    LG=sum(F((-1)**((p//d-1)//2),d**4) for d in s.divisors(p) if d%2 and (p//d)%2)
    b0=401*16*LF+1192*LG;b1=298*16*LF+802*LG
    assert mod(p**4*b0,p)==(6416*chi(p)+1192)%p
    assert mod(p**4*b1,p)==(4768*chi(p)+802)%p
    count('modular_first_prime_layer',2)

# Signed sparse determinant collapse is checked directly over polynomial X.
for M in (3,5,7):
    degree=(M-1)//2; rr=sum(s.Rational(j+1,j+2)*t**j for j in range(1,degree+1));c=s.Rational(-2,7)
    hs=[s.expand(rr+c*t**M).coeff(t,j) for j in range(2*M-1)]
    n=M-1;B=s.Matrix([[X]+[-hs[j] for j in range(1,n+1)]]+[[hs[i]]+[hs[i+j]-hs[i+j-1] for j in range(1,n+1)] for i in range(1,n+1)])
    got=B.det(method='domain-ge');rv=rr.subs(t,1);expected=(-1)**(n*(n-1)//2+1)*c**(M-2)*((rv-c)*X-rv**2)
    assert s.expand(got-expected)==0;count('sparse_determinant_identity')

# Reconstruct the five-period rational function and aggregate its pole tails.
for n in (0,1):
    lam=4*n+1
    aj=[]
    for j in range(3*n+1):
        den=s.prod(x+q for q in range(n+j+1,4*n+1))*s.prod(x+q for q in range(4*n+2,4*n+j+2))
        aj.append(s.Rational(factorial(4*n),factorial(n))*(-1)**j*s.binomial(3*n,j)/den)
    I0=0
    for j in range(3*n+1):
        for l in range(3*n+1):
            L=2*(x+3*n+1+j);M=2*(x+2*n+1+l);T=2*(x+2*n+1+j+l);V=2*(x+1+j+l)
            I0+=2*aj[j]*aj[l]*(1/(L*T)+1/(L*M)+1/(M*V))
    C=s.prod(x+2*r for r in range(1,4*n+1))*s.Rational(1,2**(4*n)*factorial(4*n))
    R=s.cancel(4*factorial(n)**4*C*I0/(factorial(4*n)*(x+lam)**2))
    assert s.cancel(R.subs(x,-2*lam-x)-R)==0
    periods={('beta',i):s.Rational(0) for i in range(1,5)}
    periods.update({('eta',i):s.Rational(0) for i in range(1,5)})
    rat=s.Rational(0)
    for q in range(1,8*n+2):
        reg=s.cancel((x+q)**4*R)
        for order in range(1,5):
            coeff=s.diff(reg,x,4-order).subs(x,-q)/factorial(4-order)
            if q%2:
                h=(q-1)//2;sgn=(-1)**h;periods['beta',order]+=sgn*coeff
                rat-=coeff*sgn*sum(s.Rational((-1)**j,(2*j+1)**order) for j in range(h))
            else:
                h=(q-2)//2;sgn=s.Rational((-1)**h,2**order);periods['eta',order]+=sgn*coeff
                rat-=coeff*sgn*sum(s.Rational((-1)**(j-1),j**order) for j in range(1,h+1))
    assert all(periods[key]==0 for key in [('beta',1),('beta',3),('eta',2),('eta',4)])
    result=[rat,periods['beta',2],periods['beta',4],periods['eta',1],periods['eta',3]*s.Rational(3,4)]
    expected=[0,0,6,0,0] if n==0 else list(map(s.Rational,['301/384','-503/64','243/32','-73/24','27/32']))
    assert result==expected;count('coupled_period_reduction')

# Tail-multiplier identity and the dominant real saddle in the simultaneous family.
for z,P in [(s.I/2,t*t+4),((1+s.I)/4,t*t-4*t+8),((1-s.I)/8,t*t-8*t+32)]:
    assert s.simplify(P.subs(t,1/(z*x))-(1-x)*(1-z/s.conjugate(z)*x)/(z*z*x*x))==0
    count('compatible_pole_tail_identity')
fprime=s.diff(x*x*(1-x*x)**2,x)/(x*x*(1-x*x)**2)
assert s.factor(fprime-2*(1-3*x*x)/(x*(1-x*x)))==0
assert s.diff(fprime,x).subs(x,1/s.sqrt(3))==-18
count('compatible_multiplier_saddle')
for n in (1,2,3):
 for p in (5,7,11):
    degree=p+n-1;F0=sum(s.Rational(chi(j),j**4)*x**j for j in range(1,degree+1))
    # Truncate after each multiplication; ordinary derivative includes degree!.
    base=s.Poly(3*F0*(2-3*F0),x);base=s.Poly(sum(c*x**j[0] for j,c in base.terms() if j[0]<=degree),x)
    acc=s.Poly(1,x)
    for _ in range(n):
        acc=acc*base;acc=s.Poly(sum(c*x**j[0] for j,c in acc.terms() if j[0]<=degree),x)
    v=acc.nth(degree)*factorial(degree)/factorial(n)
    assert vp(v,p)==-3;count('niven_origin_prime_layer')

# The exact interleaved energy comparison is enclosed with Arb, not a bare float.
ctx.dps=60
L=lambda z:z*z*z.log() if z else arb(0)
def I(a,b,equal=False):return (L(a+b)+(arb(0) if equal else L(abs(a-b)))-2*L(a)-2*L(b))/2
c=arb(1)/5;energy=(I(arb(1),arb(1),True)+2*I(arb(1),1+c)+I(1+c,1+c,True))/4
assert energy>arb('1.6523') and energy<arb('1.6524') and energy<arb(44)/25
count('interleaved_energy_enclosure')
out={'status':'passed','run_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'counts':counts,'total_exact_or_interval_checks':sum(counts.values()),'interleaved_energy':energy.str(25),'scope':'Bounded exact reconstructions and one Arb inequality supporting the newly exposed formulas; no formal or independent-model verification.'}
Path(__file__).with_name('case-studies-result.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
