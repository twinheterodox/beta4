from fractions import Fraction as F
from ctypes import CDLL, c_int
from itertools import combinations
from pathlib import Path
from collections import Counter

f = CDLL(str(Path(__file__).with_name('density_max.so').resolve())).envelope
f.argtypes = (c_int,c_int)
f.restype = c_int
offsets = [(j,0) for j in range(31)]+[(96-j,0) for j in range(31)]+[(-72,-48),(168,-48)]
classes=set()
for (s,t),(u,v) in combinations(offsets,2):
    delta=s-u
    residue=v-t
    if delta<0:delta,residue=-delta,-residue
    classes.add((delta,residue % 96))
events=set(map(F,range(1,401)))
for delta,r in classes:
    for a in range(delta+(r-delta)%96,400*delta+1,96):
        events.add(F(a,delta))
events=sorted(events)
tot=0
hist=Counter()
for a,b in zip(events,events[1:]):
    num=a.numerator+b.numerator
    den=a.denominator+b.denominator
    e16=f(num,den)
    hist[e16]+=1
    # Sum a rational upper bound for each cell, including negative ones.
    value=(1/a-1/b)*F(e16,16)*(1<<64)
    tot+=-(-value.numerator//value.denominator)
body=F(tot,1<<64)
tail=F(36,400)+F(221,32*399)
print('events',len(events),'cells',len(events)-1,'max/min e(x)',max(hist)/16,min(hist)/16)
print('body',body, float(body),'tail',tail,'with tail',body+tail, float(body+tail))
assert body == F(73794151590486879261,9223372036854775808)
assert body+tail == F(745969417189431202995707,92003136067626388684800)
