"""Selected exact Loewner construction functions; numerical searches omitted."""
from fractions import Fraction as Q

def family(name,lam):
    lam=Q(lam)
    if name=='three':
        return [(Q(1,4),Q(1)),(Q(1,2),24*lam-Q(5,2)),(Q(1),Q(3,2)-24*lam)]
    if name=='two':
        return [(Q(1,4),Q(1)),(Q(1),48*lam-6)]
    if name=='four_antiperiodic':
        return [(Q(1,4),Q(1,2)),(Q(1,2),24*lam),(Q(3,4),Q(-1,2)),(Q(1),-24*lam)]
    if name=='four_cancel':
        return [(Q(1,4),8*lam+Q(5,6)),(Q(1,2),-16*lam-Q(5,3)),
                (Q(3,4),8*lam-Q(1,6)),(Q(1),Q(1))]
    raise ValueError(name)

def target_coefficients(terms):
    # Coefficients of G and pi^2, respectively, in F'(0).
    values={Q(1,4):(8,Q(1)),Q(1,2):(0,Q(1,2)),
            Q(3,4):(-8,Q(1)),Q(1):(0,Q(1,6))}
    return tuple(sum((w*values[a][i] for a,w in terms),Q(0)) for i in (0,1))

def nodes(r,m,prefix=0,shift=0):
    """prefix is the number of border PAIRS preceding the common block."""
    assert 0<=prefix<=m
    common=list(range(2*prefix,2*prefix+r))
    left=list(range(2*prefix))
    right=list(range(2*prefix+r,r+2*m))
    border=left+right
    R=sorted(common+border[::2])
    C=sorted(common+border[1::2])
    return [t+shift for t in R],[t+shift for t in C]

def rational_data(terms,N):
    f=[Q(0)]*N
    d=[Q(0)]*N
    for t in range(1,N):
        f[t]=f[t-1]+sum((w/(t-1+a) for a,w in terms),Q(0))
        d[t]=d[t-1]-sum((w/(t-1+a)**2 for a,w in terms),Q(0))
    return f,d

def rational_matrix(terms,R,C,x=Q(0)):
    f,d=rational_data(terms,max(R+C)+1)
    return [[x+d[t] if t==u else (f[u]-f[t])/(u-t) for u in C] for t in R]
