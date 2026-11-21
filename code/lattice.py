import sys, itertools; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import numpy as np
from fractions import Fraction as Fr
import sympy as sp
def theta_rows(G, wstar, nmax):
    """counts #{v: v^T G v/2 = n, v.wstar = l} for n<=nmax"""
    G=np.array(G); L=sp.Matrix(G)
    # enumeration box via eigenvalue bound
    lam=min(np.linalg.eigvalsh(G))/2
    R=int((nmax/lam)**0.5)+1
    rng=range(-R,R+1)
    out={}
    for v in itertools.product(rng,repeat=4):
        v=np.array(v); n=int(v@G@v)//2
        if n<=nmax:
            l=int(v@np.array(wstar)); out[(n,l)]=out.get((n,l),0)+1
    return out
def search(N,nmax=6):
    found=[]
    diag=[2,4,6]; off=[-2,-1,0,1,2]
    seen=set()
    for d in itertools.product(diag,repeat=4):
        if list(d)!=sorted(d): continue
        for o in itertools.product(off,repeat=6):
            G=[[d[0],o[0],o[1],o[2]],[o[0],d[1],o[3],o[4]],[o[1],o[3],d[2],o[5]],[o[2],o[4],o[5],d[3]]]
            M=sp.Matrix(G)
            if M.det()!=N*N: continue
            if not all(M[:k,:k].det()>0 for k in range(1,5)): continue
            # count norm-1 vectors quickly
            Gi=M.inv()
            # dual vectors wstar small with w^T G^-1 w /2 = 1/N
            for ws in itertools.product(range(-2,3),repeat=4):
                w=sp.Matrix(ws)
                if (w.T*Gi*w)[0]/2==sp.Rational(1,N):
                    found.append((G,ws)); break
            if len(found)>=40: return found
    return found
if __name__=='__main__':
    N=int(sys.argv[1])
    print(len(search(N)))
