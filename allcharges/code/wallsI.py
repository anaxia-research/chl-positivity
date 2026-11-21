import sys
from fractions import Fraction as Fr
from exact import ginv, Exact
from walls import walls
N=int(sys.argv[1]); M=int(sys.argv[2]); kp2=24//(N+1)
E=Exact(N,N*(M+3),M); d=ginv(N,N*(4*M*M+20),kp2)
def D(k):
    k=Fr(k)
    return d[int(k)+1] if k.denominator==1 and k>=-1 else 0
bad=tot=0
for m in range(1,M+1):
    T=E.table(m)
    for (n,lt),v in T.items():
        nu=Fr(n,N)
        if lt>=2*m or 4*m*nu-lt*lt>=0: continue
        s=0
        for p,q,r,ss,mG,nG,lG in walls(m,nu,lt,R=4*m*N+4):
            s+=lG*((D(mG)*D(N*nG) if r%N==0 else 0)+(D(N*mG)*D(nG) if ss%N==0 else 0))
        tot+=1
        if s!=v: bad+=1; print('mismatch',m,n,lt,v,s) if bad<6 else None
print('N',N,'polar terms',tot,'mismatches',bad)
