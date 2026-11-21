from chl import CHL
from fractions import Fraction as F
for N in (3,5,7):
    C=CHL(N,40,2)
    seen={}
    for kap in range(0,3*N+1):
        for l in range(0,3*N+1):
            for j in range(-8,9):
                D=F(4*kap*l,N)-j*j
                if D<=0 and not (kap==0 and l==0 and j>=0 and D==0 and j==0):
                    v=C.E(kap,l,j)
                    key=(D,j%2,kap%N==0,l%N==0)
                    if v: seen.setdefault(key,set()).add(v)
    print(N, sorted(seen.items()))
    print(' E(0,0,0)=',C.E(0,0,0),' E(1,0,0)=',C.E(1,0,0), 'E(0,1,0)=',C.E(0,1,0))
