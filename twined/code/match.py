import sys
from itertools import combinations
from collections import defaultdict
from g2theta import theta, mul, EVEN
from tt import phi
def product_series(t, C0num, C0den, Qmax8, Smax8, Rmin4):
    """Phi = q r s^{C0} prod (1-q^{tn} r^l s^m)^{c(nm,l)}, keyed like g2theta (q^/8, s^/8, r^/4)."""
    nmax=Qmax8//8; mmax=Smax8//8
    c={k:int(v) for k,v in phi(t,(nmax//t+1)*(mmax+1)+1).items()}
    S={(8, 8*C0num//C0den, 4):1}
    def fac(S,A,B,L,e):
        for _ in range(abs(e)):
            T=defaultdict(int); T.update(S)
            if e>0:
                for (a,b,r),v in S.items():
                    if a+8*A<=Qmax8 and b+8*B<=Smax8: T[(a+8*A,b+8*B,r+4*L)]-=v
            else:
                for (a,b,r),v in S.items():
                    k=1
                    while a+8*A*k<=Qmax8 and b+8*B*k<=Smax8 and r+4*L*k>=Rmin4 and k<300:
                        T[(a+8*A*k,b+8*B*k,r+4*L*k)]+=v; k+=1
            S={k:v for k,v in T.items() if v}
        return S
    S=fac(S,0,0,-1,c[(0,-1)])
    for n in range(0,nmax//t+1):
        for m in range(0,mmax+1):
            if n==0 and m==0: continue
            for (N,l),v in c.items():
                if N==n*m and v: S=fac(S,t*n,m,l,v)
    return S
if __name__=='__main__':
    # Phi_2^g = 2^{-8} (theta[01;00] theta[01;10] theta[11;00] theta[11;11])^2, through q^5 s^5
    Qm,Sm=40,40
    P=product_series(2,1,2,Qm,Sm,-200)
    P={k:v for k,v in P.items() if k[2]>=-40}
    th={ch:theta(ch[0],ch[1],Qm,Sm) for ch in EVEN}
    found=[]
    for four in combinations(EVEN,4):
        X={(0,0,0):1}
        for ch in four: X=mul(X,th[ch],Qm,Sm)
        X=mul(X,X,Qm,Sm)
        ks=[k for k in set(P)|set(X) if k[2]>=-40]
        if all(256*P.get(k,0)==X.get(k,0) for k in ks): found.append(four)
    print('nonzero coefficients compared',sum(1 for k in P if P[k]))
    print('matching sets of four even characteristics (a;b):',found)
