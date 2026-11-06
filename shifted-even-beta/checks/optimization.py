"""Exploratory double-precision limits. Exact claims are in endpoints.py."""
from math import log, log1p
from pathlib import Path
import json
from checks import parameters


def stationary(mu,D,r):
    L=2*r+1
    def f(y):return D*log1p(L/y)-mu*log1p(1/(y+r))
    low,high=1e-100,1.
    while f(high)>0:high*=2
    for _ in range(180):
        mid=(low+high)/2
        if f(mid)>0:low=mid
        else:high=mid
    return high


def profile(mu,D,K,r):
    L=2*r+1;y=stationary(mu,D,r)
    return (mu+1+log(K)+D*L*log(2)+D*L*log(y+L)
            +mu*(r*log(y+r)-(r+1)*log(y+r+1)))


def main():
    records=[]
    for B in [1,2,3,4,5,10]:
        rows=[]
        for r in range(1,16):
            _,D,L,K,_=parameters(B,r)
            low=D*L+1e-6;high=2*low
            while profile(high,D,K,r)>0:high*=2
            for _ in range(65):
                mid=(low+high)/2
                if profile(mid,D,K,r)>0:low=mid
                else:high=mid
            rows.append({'r':r,'limiting_a_over_m':high,
                         'limiting_endpoint_slope':high+1,
                         'stationary_y':stationary(high,D,r)})
        best=min(rows,key=lambda row:row['limiting_a_over_m'])
        records.append({'B':B,'best_in_scan':best,'scan':rows})
    result={'status':'Numerical exploration, not an exact certificate.',
            'scope':'r=1,...,15; fixed factorial normalization and absolute-value profile.',
            'exact_single_value_scope':'endpoints.py proves every other integer r cannot improve a<=18m.',
            'results':records}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    for row in records:print(row['B'],row['best_in_scan'])


if __name__=='__main__':main()
