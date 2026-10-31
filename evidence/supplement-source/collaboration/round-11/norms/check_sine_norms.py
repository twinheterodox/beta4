#!/usr/bin/env python3
"""Prescribed exact cyclotomic norm gate, with separate diagnostics."""
from math import gcd
from pathlib import Path
import json
import sympy as sp
import mpmath as mp

HERE=Path(__file__).resolve().parent
w,X=sp.symbols('w X')
mp.mp.dps=100
theta=mp.im(mp.polylog(4,1j))
delta=2-mp.pi**4/96
C=mp.pi**4/96

def exact(p):
    modulus=sp.Poly(sp.cyclotomic_poly(4*p,w),w,domain=sp.QQ)
    def red(z):return sp.rem(sp.Poly(z,w,domain=sp.QQ),modulus)
    def wp(k):return red(w**(k%(4*p)))
    base=wp(1)-wp(-1)
    U=sp.Poly(0,w,domain=sp.QQ)
    for n in range(1,p,2):
        denom=wp(n)-wp(-n)
        ratio=red(base*sp.invert(denom,modulus))
        assert all(c.q==1 for c in ratio.all_coeffs())
        U=red(U+(-1)**((n-1)//2)*red(ratio**4)*red((wp(n)+wp(-n))**3))
    assert all(c.q==1 for c in U.all_coeffs())
    def automorphism(k):
        return red(sum(c*w**((k*j[0])%(4*p)) for j,c in U.terms()))
    assert automorphism(-1)==U
    assert automorphism(2*p+1)==-U
    norm=sp.Poly(sp.resultant(modulus.as_expr(),X-U.as_expr(),w),X)
    unit,factors=sp.factor_list(norm)
    assert unit==1 and len(factors)==1 and factors[0][1]==2
    minimal=factors[0][0]
    assert minimal.degree()==p-1 and minimal.LC()==1
    poly=sp.Poly(minimal.as_expr().subs(X,8*X),X,domain=sp.ZZ)
    content,primitive=poly.primitive()
    h=mp.pi/(2*p)
    conjugates=[mp.sin(k*h)**4*mp.fsum(
        (-1)**((n-1)//2)*mp.cos(k*n*h)**3/mp.sin(k*n*h)**4
        for n in range(1,p,2)) for k in range(1,2*p,2) if k!=p]
    direct=mp.fprod(theta-aa for aa in conjugates)*mp.mpf(8)**(p-1)/int(content)
    evaluated=mp.polyval([mp.mpf(int(c)) for c in primitive.all_coeffs()],theta)
    assert abs(direct-evaluated)<max(1,abs(evaluated))*mp.mpf('1e-75')
    return {'p':p,'real_degree':p-1,'eight_alpha_polynomial_in_w':str(U.as_expr()),
            'minimal_polynomial_eight_alpha':str(minimal.as_expr()),
            'primitive_norm_polynomial':str(primitive.as_expr()),
            'norm_content':int(content),'constant_term_odd':bool(int(minimal.TC())%2),
            'log_abs_primitive_at_theta':mp.nstr(mp.log(abs(evaluated)),35)}

rows=[exact(p) for p in (3,5,7,11)]
diagnostics=[]
for p in (11,31,101):
    h=mp.pi/(2*p)
    vals=[]
    for k in range(1,2*p,2):
        if k==p:continue
        aa=mp.sin(k*h)**4*mp.fsum(
            (-1)**((n-1)//2)*mp.cos(k*n*h)**3/mp.sin(k*n*h)**4
            for n in range(1,p,2))
        scale=mp.cos(h)**3*(mp.sin(k*h)/mp.sin(h))**4
        assert delta*scale <= abs(aa) <= C*scale
        vals.append((k,aa))
    error=theta-vals[0][1]
    assert error>0
    logprod=mp.fsum(mp.log(abs(theta-aa)) for _,aa in vals)
    diagnostics.append({'p':p,'principal_error_times_p2':mp.nstr(error*p*p,35),
        'smallest_other_abs_conjugate':mp.nstr(min(abs(aa) for k,aa in vals if k not in(1,2*p-1)),35),
        'log_monic_norm_at_theta':mp.nstr(logprod,35),
        'log_norm_over_d_log_p':mp.nstr(logprod/((p-1)*mp.log(p)),35)})

out={'scope':'Four prescribed exact cases; three all-conjugate diagnostic checks. No parameter scan.',
     'exact':rows,'diagnostics':diagnostics,
     'diagnostic_status':'mpmath100-digit values; not certified intervals.'}
(HERE/'checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'exact':[{'p':r['p'],'content':r['norm_content'],
                          'log_primitive':r['log_abs_primitive_at_theta']} for r in rows],
                  'diagnostics':diagnostics},indent=2))
