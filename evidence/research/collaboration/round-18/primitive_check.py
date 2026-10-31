"""Full primitive normalization of the explorer's existing 24 exact forms.
No new determinant construction, matrix size or parameter is introduced.
"""
from pathlib import Path
import json,math,sys
from flint import fmpq,arb,ctx
sys.set_int_max_str_digits(0)
here=Path(__file__).resolve().parent
data=json.loads((here/'explorer/exact-model-coefficients.json').read_text())
out=[]
for d in data:
    for name,ak,bk in [('direct','a4','b4'),('first_variation','wa','wb')]:
        a=fmpq(d[ak]);b=fmpq(d[bk]);D=math.lcm(int(a.denominator),int(b.denominator))
        Ai=int(a*D);Bi=int(b*D);g=math.gcd(Ai,Bi);A=Ai//g;B=Bi//g
        assert g>0 and math.gcd(A,B)==1
        assert a*D==A*g and b*D==B*g
        ctx.prec=max(abs(A).bit_length(),abs(B).bit_length())+512
        theta=(arb(4).zeta(arb(1)/4)-arb(4).zeta(arb(3)/4))/256
        value=A*theta+B
        assert not value.contains(0)
        logvalue=abs(value).log();n=d['n']
        row={'n':n,'model':name,'D':str(D),'g':str(g),
             'A':str(A),'B':str(B),'precision_bits':ctx.prec,
             'log_primitive_value':str(logvalue),
             'log_primitive_value_per_n2':str(logvalue/(n*n)),
             'log_integer_content_per_n2':str(arb(g).log()/(n*n)),
             'absolute_value_above_one':bool(abs(value)>1)}
        out.append(row)
        if n in [3,6,9,12]:
            print(n,name,float(logvalue/(n*n)),float(arb(g).log()/(n*n)),row['absolute_value_above_one'])
(here/'primitive-check.json').write_text(json.dumps({'records':out,
 'scope':'Exact full gcd normalization and Arb enclosures of the 24 already constructed forms; finite evidence only.'},indent=2)+'\n')
