"""Exact modular rank gate for inhomogeneous polynomial multipliers.

Family: B=m=r=1,a=4. This tests cancellation of the beta(2)
coefficient. Full modular rank rules out a rational polynomial relation
on the explicitly tested indices; it does not by itself exclude a relation
restricted to an arbitrary untested eventual subsequence.
"""
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent
P=1000003


def inv(x):return pow(x%P,-1,P)


def multiply(left,right):
    out={}
    for (u,j),x in left.items():
        for (v,k),y in right.items():
            key=(u+v,j+k)
            out[key]=out.get(key,0)+x*y
    return {key:x for key,x in out.items() if x}


def power(poly,exponent):
    out={(0,0):1}
    for _ in range(exponent):out=multiply(out,poly)
    return out


def telescope(u,j):
    """Polynomial W_F for F(n,y)=n^u*y^(2j+1); exact integer coefficients."""
    first=multiply({(0,1):1,(1,0):-3,(0,0):-1},power({(0,1):1,(1,0):1},4))
    first=multiply(first,power({(0,1):1,(0,0):-1},2*j+1))
    second=multiply({(0,1):1,(1,0):3,(0,0):1},power({(0,1):1,(1,0):-1},4))
    second=multiply(second,power({(0,1):1,(0,0):1},2*j+1))
    out={}
    for key in first.keys()|second.keys():
        coefficient=first.get(key,0)+second.get(key,0)
        assert coefficient%2==0
        if coefficient:
            assert key[1]%2==0
            out[(key[0]+u,key[1])]=coefficient//2
    assert min(sum(key) for key in out)>=4
    return out


def beta2_rows(n,J):
    H=[0];H2=[0];O=[0];O2=[0]
    for j in range(1,2*n+1):
        q=inv(j);v=inv(2*j-1)
        H.append((H[-1]+q)%P);H2.append((H2[-1]+q*q)%P)
        O.append((O[-1]+v)%P);O2.append((O2[-1]+v*v)%P)
    weights=[1]
    for h in range(n):
        weights.append(weights[-1]*(2*(n+h)+1)*inv(2*(2*n-h)-1)
                       *pow((n-h)*inv(h+1)%P,4,P)%P)
    out=[0]*(J+1)
    target=[0]*(J+1)
    for h,g in enumerate(weights):
        x=n-2*h
        ell1=(O[2*n-h]-O[n+h]-2*(H[n-h]-H[h]))%P
        ell2=(-O2[2*n-h]-O2[n+h]+H2[n-h]+H2[h])%P
        A=(ell1*ell1+ell2)*inv(2)%P
        for j in range(J+1):
            term=pow(x,2*j,P)*A
            if j:term+=2*j*pow(x,2*j-1,P)*ell1+j*(2*j-1)*pow(x,2*j-2,P)
            out[j]=(out[j]+g*term)%P
            target[j]=(target[j]+g*pow(x,2*j,P))%P
    return out,target


def rank_and_kernel(matrix):
    mat=[row[:] for row in matrix];pivots=[];i=0
    for j in range(len(mat[0])):
        k=next((k for k in range(i,len(mat)) if mat[k][j]%P),None)
        if k is None:continue
        mat[i],mat[k]=mat[k],mat[i]
        z=inv(mat[i][j]);mat[i]=[x*z%P for x in mat[i]]
        for k in range(len(mat)):
            if k!=i and mat[k][j]:
                z=mat[k][j]
                mat[k]=[(x-z*y)%P for x,y in zip(mat[k],mat[i])]
        pivots.append(j);i+=1
        if i==len(mat):break
    free=[j for j in range(len(mat[0])) if j not in pivots]
    basis=[]
    for j in free:
        v=[0]*len(mat[0]);v[j]=1
        for k,col in enumerate(pivots):v[col]=-mat[k][j]%P
        basis.append(v)
    return len(pivots),basis


def main():
    assert all(P%q for q in range(2,1001))
    records=[]
    for E in (2,4,6,8,10,12,16):
        columns=[(u,j) for j in range(E//2+1) for u in range(E-2*j+1)]
        first=E+1
        indices=list(range(first,first+2*(len(columns)+8),2))
        mat=[];targets=[]
        for n in indices:
            rho,target=beta2_rows(n,E//2)
            mat.append([pow(n,u,P)*rho[j]%P for u,j in columns])
            targets.append([pow(n,u,P)*target[j]%P for u,j in columns])
        rank,kernel=rank_and_kernel(mat)
        target_survivors=sum(any(sum(x*y for x,y in zip(w,t))%P for t in targets) for w in kernel)
        generators=[]
        if E>=6:
            for j in range((E-6)//2+1):
                for u in range(E-6-2*j+1):
                    poly=telescope(u,j)
                    vector=[poly.get((v,2*k),0)%P for v,k in columns]
                    assert all(sum(x*y for x,y in zip(vector,row))%P==0 for row in mat)
                    assert all(sum(x*y for x,y in zip(vector,row))%P==0 for row in targets)
                    generators.append(vector)
            assert rank_and_kernel(generators)[0]==len(generators)
        assert len(generators)==len(kernel)
        assert target_survivors==0
        record={'total_degree':E,'dimension':len(columns),'rank_mod_p':rank,
                'kernel_dimension_mod_p':len(kernel),'kernel_vectors_with_nonzero_target':target_survivors,
                'independent_exact_telescoping_generators':len(generators),
                'indices':indices,'columns':columns,'kernel_mod_p':kernel}
        records.append(record)
        print(E,'columns',len(columns),'rank',rank,'kernel',len(kernel),'target_survivors',target_survivors,flush=True)
    first=telescope(0,0)
    result={'status':'Every modular kernel has exactly the dimension of the explicit rational telescoping subspace.',
            'scope':'Complete rational-kernel classification on each stated finite index set; the explicit telescopers are valid in general.',
            'modulus':P,'first_telescoper':{f'n^{u}*x^{j}':c for (u,j),c in sorted(first.items())},'records':records}
    (HERE/'polynomial-gate.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
