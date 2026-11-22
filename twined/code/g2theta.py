# genus-two theta constants  theta[a;b](Z), Z=[[tau,z],[z,sigma]] with tau=sigma~, z=v~, sigma=rho~
# monomials keyed (Q,S,R): q^(Q/8) s^(S/8) r^(R/4)   (q=e(rho~), s=e(sigma~), r=e(v~))
from collections import defaultdict
from itertools import product
def theta(a,b,Qmax,Smax):
    T=defaultdict(int); L=12
    for l1 in range(-L,L+1):
        for l2 in range(-L,L+1):
            x1=2*l1+a[0]; x2=2*l2+a[1]          # 2(l+a/2)
            S=x1*x1; Q=x2*x2; R=x1*x2           # s^{x1^2/8} q^{x2^2/8} r^{x1x2/4}
            if S<=Smax and Q<=Qmax:
                ph=(l1*b[0]+l2*b[1])%2
                T[(Q,S,R)]+= -1 if ph else 1
    return {k:v for k,v in T.items() if v}
def mul(A,B,Qmax,Smax):
    C=defaultdict(int)
    for (a1,b1,c1),x in A.items():
        for (a2,b2,c2),y in B.items():
            if a1+a2<=Qmax and b1+b2<=Smax: C[(a1+a2,b1+b2,c1+c2)]+=x*y
    return {k:v for k,v in C.items() if v}
EVEN=[(a,b) for a in product((0,1),repeat=2) for b in product((0,1),repeat=2) if (a[0]*b[0]+a[1]*b[1])%2==0]
