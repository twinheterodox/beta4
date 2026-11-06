"""Exact endpoint certificates; floating point only locates rational trials."""
from fractions import Fraction as F
from math import ceil, exp, factorial, log
from pathlib import Path
import json
from checks import parameters


ELOW=F(163,60)
EUP=sum((F(1,factorial(j)) for j in range(21)),F(0))+F(1,20*factorial(20))


def stationary(B,m,r,a):
    _,D,L,_,_=parameters(B,r)
    def fun(x):return m*D*log(1+2*L/x)-a*log(1+2/(x+2*r))
    low,high=1e-100,1.
    while fun(high)>0:high*=2
    for _ in range(180):
        mid=(low+high)/2
        if fun(mid)>0:low=mid
        else:high=mid
    return high


def certificate(B,m,r,a,u):
    _,D,L,K,_=parameters(B,r)
    assert a%2==0 and a>m*D*L
    f=((u+2*L)/u)**(m*D)*((u+2*r)/(u+2*r+2))**a
    U=(EUP**(2*(a+m-1))*K**(2*m)*(u+2*L)**(2*m*D*L)
       *2**(2*a)*(F(2*r)**(2*r)/F(2*r+2)**(2*r+2))**a)
    assert f<1 and U<1
    return {'B':B,'m':m,'r':r,'a':a,'first_even_index':m+1,
            'last_even_index':a+m-1,'u':str(u),'f_lt_1':True,'U_lt_1':True,
            'f_approx':exp(log(f.numerator)-log(f.denominator)),
            'U_approx':exp(log(U.numerator)-log(U.denominator))}


def exp_bounds(x,N=80):
    assert 0<=x<N+2
    term=F(1);lower=term
    for j in range(1,N+1):term*=x/j;lower+=term
    nextterm=term*x/(N+1)
    return lower,lower+nextterm/(1-x/F(N+2))


def main():
    assert __debug__
    choices=[(1,1,18),(1,3,52),(1,5,88),(1,7,122),(1,19,334),(1,99,1740),
             (2,1,70),(2,3,210),(2,19,1326),
             (3,1,186),(3,3,556),(3,19,3520),
             (4,1,360),(4,3,1084),(4,19,6860),
             (5,1,554),(5,3,1662),(5,19,10528),
             (10,1,2794),(10,3,8382),(10,19,53092)]
    certificates=[]
    for B,m,a in choices:
        u=F(ceil(stationary(B,m,2,a)*10**12)+1,10**12)
        row=certificate(B,m,2,a,u)
        certificates.append(row)
        print('Endpoint passed',B,m,row['last_even_index'],row['U_approx'],flush=True)

    # Uniform B=1 endpoint a=2 ceil((879m-50)/100).
    assert 1001**3*F(401,601)**52<1
    assert EUP**929*2**250*F(1001,100)**250*F(4,27)**879<1
    assert F(27,4)<ELOW**2
    # Globally exclude every other integer r from beating a<=18m.
    exclusions=[]
    for r in [1,3,4,5,6,7,8]:
        L=2*r+1
        witness=ELOW**18*2**(2*L)*L**L*(F(r)**r/F(r+1)**(r+1))**18
        assert witness>1
        exclusions.append({'r':r,'profile_at_zero_for_a_18m_already_too_large':True})
    # r>=9 is excluded by the degree condition a>m(2r+1)>=19m.

    e0=F(68,25)
    assert F(1001,50)<ELOW**3
    assert F(4)**10*e0**19<F(27)**10
    assert 1001*F(401,601)**18<1
    assert e0<sum((F(11,10)**j/F(factorial(j)) for j in range(6)),F(0))
    uniform=[]
    for B in [1,2,3,4,5,10]:
        _,D,_,_,_=parameters(B,2);M=2*B-1
        if B==1:N=9
        else:
            N=ceil(D*(9+40/9*log(M)))
            low=F(9*(N-1-9*D),40*D)
            high=F(9*(N-9*D),40*D)
            assert exp_bounds(low)[1]<M<exp_bounds(high)[0]
        uniform.append({'B':B,'D_odd':D,'a_B':2*N,'endpoint_slope':2*N+1})
    record={'status':'All exact endpoint certificates passed; not an independent audit.',
            'e_upper':str(EUP),'e_upper_proof':'Sum through 20 plus 1/(20*20!); geometric tail bound.',
            'certificates':certificates,'uniform_general_B':uniform,
            'uniform_single_value':'a(m)=2*ceil((879*m-50)/100), endpoint a(m)+m-1; every positive odd m.',
            'single_value_uniform_comparisons_passed':True,
            'integer_zero_range_exclusions':exclusions,
            'r_at_least_9':'Degree condition forces a>19m, so cannot improve a<=18m.',
            'optimization_scope':'Same factorial denominator normalization and absolute-value exponential profile; integer zero-range r. No global obstruction for other constructions is asserted.'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(record,indent=2)+'\n')


if __name__=='__main__':main()
