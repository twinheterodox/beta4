from math import isqrt, prod
import numpy as np

roots = [96-2*a for a in range(97)]
orders = [sum(j <= a <= 96-j for j in range(31)) for a in range(97)]
odd = list(range(-239, 240, 2))

def polynomial_from_roots(multiset, p):
    result = np.array([1], dtype=np.int64)
    for r in multiset:
        result = np.convolve(result, [-r, 1]) % p
    return result

def taylor_at_roots(poly, powers, K, p):
    # Columns are the coefficients of poly(r+w) for every r, using global monomials.
    result = []
    degree = len(poly)-1
    binom = np.ones(degree+1, dtype=np.int64)
    for k in range(K):
        terms = poly[k:] * binom % p
        result.append((powers[:, :len(terms)] @ terms) % p)
        binom = (binom[:-1] * np.arange(k+1, degree+1) % p) * pow(k+1, -1, p) % p
    return np.stack(result, axis=1)

def horner(coefs, x, p):
    y = np.zeros_like(x)
    for c in reversed(coefs):
        y = (y*x + c)%p
    return y

def one(p):
    assert all(p % j for j in range(2,isqrt(p)+1))
    P = polynomial_from_roots([0]+odd*3, p)
    Q = polynomial_from_roots([r for r,b in zip(roots,orders) for _ in range(b)],p)
    assert (len(P)-1,len(Q)-1)==(721,2077)
    residues = np.array(roots) % p
    powers = np.ones((97,len(Q)),dtype=np.int64)
    for j in range(1,len(Q)):
        powers[:,j] = powers[:,j-1]*residues%p
    tp = taylor_at_roots(P,powers,3*max(orders),p)
    tq = taylor_at_roots(Q,powers,4*max(orders),p)
    principal = {1:[],2:[],3:[]}
    for a,b in enumerate(orders):
        assert np.all(tq[a,:b] == 0) and tq[a,b]
        div = tq[a,b:4*b]
        inv0 = pow(int(div[0]), -1, p)
        # Formal division of the two global Taylor polynomials, after removing w^b.
        value = np.zeros(3*b,dtype=np.int64)
        for j in range(3*b):
            value[j] = (int(tp[a,j])-int(div[1:j+1] @ value[:j][::-1]))*inv0%p
        sq = np.convolve(value,value)%(p)
        cu = np.convolve(sq[:3*b],value)%p
        for q,poly in [(1,value),(2,sq),(3,cu)]:
            principal[q].append(list(map(int,poly[:q*b][::-1])))
    inv = np.array([0]+[pow(j,-1,p) for j in range(1,p)],dtype=np.int64)
    output = []
    for q in range(1,4):
        coeffs = principal[q]
        for a in range(97):
            assert all(c == ((-1)**i*d)%p for i,(c,d) in enumerate(zip(coeffs[a],coeffs[96-a]),1))
        # The recombined fraction must agree with the original at q*deg Q points.
        # The numerator of the difference after multiplying by Q^q has smaller degree.
        x = np.arange(97,97+q*2077,dtype=np.int64)
        truth = (horner(P,x,p) * inv[horner(Q,x,p)])%p
        rec = np.zeros_like(x)
        for r,cs in zip(roots,coeffs):
            recip = inv[(x-r)%p]
            rec = (rec + recip*horner(cs,recip,p))%p
        assert np.all(rec == truth**q%p)
        sig = [0]*(33*q+1)
        k = 0
        beta = [0]*(31*q+1)
        for a,cs in enumerate(coeffs):
            sign = (-1)**a
            for i,c in enumerate(cs,1):
                t = sign * prod(range(i,i+2*q)) * c
                k = (k-t*sig[i+2*q])%p
                beta[i] = (beta[i]+t)%p
            inva = pow(2*a+1,-1,p)
            for s in range(1,len(sig)):
                sig[s] = (sig[s]+sign*pow(inva,s,p))%p
        assert not any(beta[1::2])
        output.append((q, sum(map(len,coeffs)), len(x), k))
    print(p,output,flush=True)
    return [entry[-1] for entry in output]

if __name__ == "__main__":
    printed_rows = {10007: [3738, 8279, 9978], 10009: [5333, 5481, 6289]}
    for modulus, expected in printed_rows.items():
        assert one(modulus) == expected
    # Additional modulus independent of the manuscript's table.
    assert one(10037) == [938, 5660, 4500]
