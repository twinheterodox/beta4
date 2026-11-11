"""Rebuild the arithmetic partition and evaluate its maxima directly.

The C++ evaluator tests all discontinuities and intervening open cells
using the defining closed interval/odd-integer counts. It does not use
the supplied jump-accumulation algorithm. All rounding and logarithm
bounds here use exact fractions. No manuscript checker is imported.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import hashlib,json,subprocess,time
HERE=Path(__file__).resolve().parent
G=F(179,10);R=F(7,6);V=F(-744,25)
BLOCKS=[(F(0),8)]+[(F(k,160),1) for k in [*range(1,9),*range(10,47,2)]]

def partition():
    lines={(0,240),(-560,-240),(1040,-240)}
    for d,m in BLOCKS:
        s=int(480*d);lines.update(((s,0),(480-s,0)))
    result={1+F(j,20) for j in range(40)}|{3+F(j,2) for j in range(795)}
    for (s,t),(ss,tt) in combinations(sorted(lines),2):
        if ss==s:continue
        ds=ss-s;dt=t-tt
        lower=(ds-dt)//480+1
        # ds < dt+480*k < 400*ds; upper end is strict.
        upper=(400*ds-dt-1)//480
        for k in range(lower,upper+1):
            result.add(F(dt+480*k,ds))
    assert min(result)==1 and max(result)==400
    return sorted(result)

def log_interval(x,terms=70):
    x=F(x);assert x>0
    if x<1:
        lo,hi=log_interval(1/x,terms);return -hi,-lo
    k=0
    while x>2:x/=2;k+=1
    def unit_log(t):
        y=(t-1)/(t+1);power=y;total=F(0)
        for j in range(terms):
            total+=2*power/(2*j+1);power*=y*y
        return total,total+2*power/((2*terms+1)*(1-y*y))
    lo,hi=unit_log(x)
    if k:
        l2,h2=unit_log(F(2));lo+=k*l2;hi+=k*h2
    return lo,hi

def sum_logs(terms):
    lo=hi=F(0)
    for weight,x in terms:
        a,b=log_interval(x)
        lo+=weight*(a if weight>=0 else b)
        hi+=weight*(b if weight>=0 else a)
    return lo,hi

def real_bound(A):
    curvature=[]
    for a,b in [(R,F(6,5)),(F(6,5),F(5,4)),(F(5,4),F(4,3))]:
        upper=3/(a+1+R)-3/(b-R)+sum(m/(a+d)-m/(b+1-d) for d,m in BLOCKS)
        assert upper<0;curvature.append(str(upper))
    x=F(58911,50000)
    phase=[(3*(x+1+R),x+1+R),(-3*(x-R),x-R)]
    slope=[(F(3),x+1+R),(F(-3),x-R)]
    for d,m in BLOCKS:
        phase += [(-m*(x+1-d),x+1-d),(m*(x+d),x+d)]
        slope += [(-m,x+1-d),(m,x+d)]
    pl,ph=sum_logs(phase);sl,sh=sum_logs(slope)
    assert ph<F(-297613,10000) and F(-1,500)<sl and sh<0
    tangent=F(-297613,10000)+F(1,500)*(x-R)
    assert tangent<V
    sb2=sum(m*(1-2*d)/(3-d) for d,m in BLOCKS)
    sb3=sum(m*(1-2*d)/(4-d) for d,m in BLOCKS)
    assert sb2>F(191,20) and 3*log_interval(24)[1]<F(477,50)
    assert sb3>7 and 3*log_interval(6)[1]<F(27,5)
    assert F(180,11)-F(837,40)<-4
    assert 7*log_interval(3)[1]+35*log_interval(2)[1]<34
    l2,h2=log_interval(2);ql=V+A+20*l2;qh=V+A+20*h2
    assert qh<F(-357,5000)
    scale=10**70
    def rounded_pair(lo,hi):
        return [str(F((lo*scale).__floor__(),scale)),str(F((hi*scale).__ceil__(),scale))]
    return dict(curvature_upper=curvature,phi_at_x0=rounded_pair(pl,ph),
        derivative_at_x0=rounded_pair(sl,sh),tangent_upper=str(tangent),
        Q_enclosure=rounded_pair(ql,qh),Q_approx=float((ql+qh)/2),
        Q_strict_upper=str(F(-357,5000)),all_global_phase_and_tail_checks=True)

if __name__=='__main__':
    start=time.monotonic();events=partition();assert len(events)==62132
    mids=[(a+b)/2 for a,b in zip(events,events[1:])]
    data=''.join(f'{x.numerator} {x.denominator}\n' for x in mids)
    assert max(abs(x.numerator)*4000+x.denominator*1000 for x in mids)<2**62
    run=subprocess.run([str(HERE/'direct_estars')],input=data,text=True,capture_output=True,check=True)
    values=list(map(int,run.stdout.split()));assert len(values)==62131
    total=0;fingerprint=hashlib.sha256()
    for a,b,v10 in zip(events,events[1:],values):
        value=F(v10,10);term=(1/a-1/b)*value*2**96
        rounded=term.__ceil__();total+=rounded
        fingerprint.update(f'{a}|{b}|{value}|{rounded}\n'.encode())
    assert total==1243345293752876677540040193981
    A=F(total,2**96)+F(44,400)+G/798
    assert A==F(12507029707866288048422972493354331,790300921079786767495600904601600)
    phase=real_bound(A)
    output=dict(events=len(events),intervals=len(values),rounded_sum=str(total),
        A=str(A),event_value_sha256=fingerprint.hexdigest(),phase=phase,
        seconds=time.monotonic()-start,method='Direct closed-interval counts at all y-cells and boundaries')
    (HERE/'fresh-budget.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:v for k,v in output.items() if k!='phase'},indent=2),flush=True)
    print('Q =',phase['Q_approx'],'strict upper',phase['Q_strict_upper'],flush=True)
