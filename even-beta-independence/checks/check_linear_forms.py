"""Exact two-route coefficient comparison. Same-assistant check, not an independent model review. Adapted unchanged from the initial investigation."""
import json
from math import factorial
from pathlib import Path
import sympy as sp

t,z=sp.symbols('t z')
HERE=Path(__file__).resolve().parent


def pole_coefficients(expr,shifts,depth):
    data={}
    for shift in shifts:
        analytic=sp.cancel((t+shift)**depth*expr)
        for s in range(1,depth+1):
            data[shift,s]=sp.diff(analytic,t,depth-s).subs(t,-shift)/factorial(depth-s)
    return data


def fundamental(depth,K):
    U=[sp.eye(depth+1)]
    for k in range(K):
        out=sp.zeros(depth+1)
        for q in range(k+1):
            A=sp.zeros(depth+1)
            A[1,0]=-sum(sp.Rational(1,2**h) for h in range(q+1))
            for i in range(1,depth+1):
                A[i,i]=sp.Rational(1,2)
                if i>1:
                    A[i,i-1]=1
            out+=A*U[k-q]
        U.append(out/(k+1))
    return U


n,a,M=1,4,2
F=1/(t**a*(t+1)**a)
parts=pole_coefficients(F,[sp.Integer(0),sp.Integer(1)],a)
P={i:sum(parts[j,i]*z**j for j in range(n+1)) for i in range(1,a+1)}
results=[]
for p in [0,1]:
    d=a+p
    U=fundamental(d,6)
    C=[sp.Integer(0)]*(d+1)
    for i in range(1,a+1):
        C[i+p]=(-1)**p*sp.rf(i,p)*z**M*P[i]
    Fp=sp.diff(F,t,p)
    for ell in [4,5,6]:
        qrow=sp.zeros(1,d+1)
        for v in range(ell+1):
            crow=sp.Matrix([[sp.diff(c,z,v).subs(z,-1)/factorial(v) for c in C]])
            qrow+=crow*U[ell-v]
        ode=[qrow[0,0]]+[(-2**(i+1)*qrow[0,i] if i%2==0 else 0) for i in range(1,d+1)]
        pos=sp.prod(M-t-1-j for j in range(ell))/factorial(ell)
        neg=sp.prod(M+t-j for j in range(ell))/factorial(ell)
        raw=sp.cancel(sp.Integer(-1)**(M-1-ell)*(pos*Fp.subs(t,t+sp.Rational(1,2))+neg*Fp.subs(t,-t-sp.Rational(1,2))))
        shifts=[sp.Rational(1,2)+j for j in range(-n,n+1)]
        coeff=pole_coefficients(raw,shifts,d)
        assert sp.cancel(raw-sum(v/(t+c)**s for (c,s),v in coeff.items()))==0
        direct=[sp.Integer(0)]*(d+1)
        for (c,s),v in coeff.items():
            j=int(c-sp.Rational(1,2))
            sign=sp.Integer(-1)**j
            direct[s]+=v*sign*2**s
            direct[0]-=v*sign*sum(sp.Rational((-1)**h,1)/(sp.Rational(1,2)+h)**s for h in range(M+j))
        assert all(sp.simplify(x-y)==0 for x,y in zip(ode,direct)),(p,ell,ode,direct)
        results.append({'p':p,'derivative_order':ell,'beta_coefficients':[str(x) for x in direct]})
        print('Passed',p,ell,flush=True)
(HERE/'half-shift-forms-checks.json').write_text(json.dumps({'fixtures':results,'status':'Exact two-route coefficient comparison; same-assistant implementation.'},indent=2)+'\n')
