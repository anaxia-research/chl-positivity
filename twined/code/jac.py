# Weak Jacobi forms phi = 4(f2^2 f3^2 + f3^2 f4^2 + f4^2 f2^2)  (index 2)  and  8 f2^2 f3^2 f4^2 (index 3)
# computed exactly as (1/eta^6) * theta products.  Coefficients c[(n,l)] , n integer.
from fractions import Fraction as Fr
from collections import defaultdict
def mul(A,B,qmax):
    C=defaultdict(int)
    for (a,b),x in A.items():
        for (c,d),y in B.items():
            if a+c<=qmax: C[(a+c,b+d)]+=x*y
    return {k:v for k,v in C.items() if v}
# exponents in units: q^(e/8), r^(f/2)
def theta(i,z,E):
    T={}
    R=int((2*E)**0.5)+3
    for n in range(-R,R+1):
        if i==2: e=(2*n+1)**2; f=(2*n+1); sg=1
        elif i==3: e=4*n*n; f=2*n; sg=1
        elif i==4: e=4*n*n; f=2*n; sg=1 if n%2==0 else -1
        if e<=E:
            k=(e, f if z else 0); T[k]=T.get(k,0)+sg
    return T
def eta_inv6(E):  # 1/eta^6 = q^{-1/4} prod (1-q^n)^-6 ; in units q^{1/8}: shift -2
    p=[0]*(E//8+3); p[0]=1
    for _ in range(6):
        for n in range(1,len(p)):
            for i in range(n,len(p)): p[i]+=p[i-n]
    return {(8*i-2,0):v for i,v in enumerate(p) if 8*i-2<=E}
def phi(t,nmax):
    E=8*nmax+8
    th={(i,z):theta(i,z,E) for i in (2,3,4) for z in (0,1)}
    sq=lambda T: mul(T,T,E)
    S=defaultdict(int)
    if t==2:
        for (i,j,k) in ((2,3,4),(3,4,2),(4,2,3)):
            P=mul(mul(sq(th[(i,1)]),sq(th[(j,1)]),E),sq(th[(k,0)]),E)
            for kk,v in P.items(): S[kk]+=v
        S=mul(dict(S),eta_inv6(E),E)
        out={}
        for (e,f),v in S.items():
            assert e%8==0 and f%2==0
            if e//8<=nmax: out[(e//8,f//2)]=v
        return out
    if t==3:
        # 8 f2^2f3^2f4^2 = 8 th2^2th3^2th4^2(z)/(2 eta^3)^2 = 2 th(z)^6 / eta^6
        P=mul(mul(sq(th[(2,1)]),sq(th[(3,1)]),E),sq(th[(4,1)]),E)
        S=mul(P,eta_inv6(E),E)
        out={}
        for (e,f),v in S.items():
            assert e%8==0 and f%2==0
            if e//8<=nmax: out[(e//8,f//2)]=2*v
        return out
if __name__=='__main__':
    for t in (2,3):
        c=phi(t,3)
        for n in range(4): print(t,n,sorted((l,v) for (a,l),v in c.items() if a==n))
