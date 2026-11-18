from math import prod

def strip(z,p):
    assert z
    v=0
    while z%p==0:
        v+=1
        z//=p
    return v,z%p

def factor_value(u,n,p,trim):
    # Evaluate the actual rational function, storing valuation and leading unit.
    x=u if trim else p*u
    halfwidth=239 if trim else 5*n//2+1
    v,top=strip(x,p)
    for r in range(-halfwidth,halfwidth+1,2):
        e,t=strip(x-r,p)
        top=top*pow(t,3,p)%p
        v+=3*e
    bot=1
    block_size=n//96
    for a in range(n+1):
        order=1+min(a//block_size,(n-a)//block_size,30)
        e,t=strip(x-n+2*a,p)
        v-=order*e
        bot=bot*pow(t,order,p)%p
    return v,top*pow(bot,-1,p)%p

def dilation(p):
    n=96*p
    zero=[0]*p
    pole=[0]*p
    for r in range(-240*p-1,240*p+2,2):
        zero[r%p]+=3
    for a in range(n+1):
        pole[(n-2*a)%p]+=1+min(a//p,(n-a)//p,30)
    A=pole[0]-zero[0]-1
    other_deficits={d:e for d,e in enumerate(pole) if d}
    othermax=max(e-zero[d] for d,e in other_deficits.items())
    assert (A,othermax)==(1356,1326)
    within=set(v%p for v in range(-239,240,2))|set((96-2*a)%p for a in range(97))
    samples=[a+p for a in range(p) if a not in within]
    assert len(samples)==p-337
    # Units outside all selected roots and poles; u>p avoids exact rational zeros.
    for u in samples:
        ve,ue=factor_value(u,n,p,False)
        vf,uf=factor_value(u,96,p,True)
        assert (ve+A,ue)==(vf,(-uf)%p)
    # Inspect all denominators of the full scaled harmonic sums at the selected poles.
    # Compute the surviving numerators using the actual h, without substituting h=pv+(p-1)/2.
    sums=[0]*100
    base=[0]*100
    chi=(-1)**((p-1)//2)
    for a in range(97):
        assert all((x-chi*y)%p==0 for x,y in zip(sums,base))
        if a==96: break
        for h in range(a*p,(a+1)*p):
            if (2*h+1)%p: continue
            v=(2*h+1)//p
            assert v%p
            iv=pow(v,-1,p)
            for s in range(1,100):
                sums[s]=(sums[s]+(-1)**h*pow(iv,s,p))%p
        iv=pow(2*a+1,-1,p)
        for s in range(1,100):
            base[s]=(base[s]+(-1)**a*pow(iv,s,p))%p
    print('prime',p,'full kernel samples',len(samples),'selected excess',A,'max unselected',othermax,'harmonic sign',chi,flush=True)

for p in (347,349):dilation(p)
