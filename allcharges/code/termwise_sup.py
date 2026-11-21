# sup over c >= c1 of sum_{pairs} |d(A)d(B)| |rho_AB - e^{-cS}|, using
#   a_lo e^{-c D'} <= R <= a_hi e^{-c D'},  R = rho e^{cS},  D' = sqrt t0 - sqrt t - S,
#   a_hi = (k/(m+2)) (t/t0)^{nu/2} (t0/t)^{1/4},  a_lo = (k/(m+2)) (t/t0)^{nu}
# (from: sqrt z e^{-z} I_nu increasing, and I'/I <= 1+nu/z).  Usage: termwise_sup.py N m0 m1 x1 [D0]
import sys, math, numpy as np
N=int(sys.argv[1]); m0=int(sys.argv[2]); m1=int(sys.argv[3]); X1=float(sys.argv[4]); D0=float(sys.argv[5]) if len(sys.argv)>5 else 400
sys.argv=[sys.argv[0],str(N)]
from exact import ginv
kp2=24//(N+1); k=kp2-2; nu=k+1.5
d=np.array([float(v) for v in ginv(N,m1+400,kp2)])
def supf(al,Sv,Dp,c1):
    # f(c)=e^{-cS}(al e^{-c Dp}-1); return sup_{c>=c1}|f|
    f=lambda c,Sv,al,Dp: np.exp(-c*Sv)*(al*np.exp(-c*Dp)-1)
    v=np.abs(f(c1,Sv,al,Dp))
    with np.errstate(divide='ignore',invalid='ignore'):
        cs=np.log(al*(Sv+Dp)/Sv)/Dp
    ok=np.isfinite(cs)&(cs>c1)
    v[ok]=np.maximum(v[ok],np.abs(f(cs[ok],Sv[ok],al[ok],Dp[ok])))
    return v
for m in range(m0,m1+1):
    x1=max(X1,math.sqrt(D0)/m); c1=math.pi*x1
    A,B=np.meshgrid(np.arange(m+3),np.arange(m+3),indexing='ij'); A=A.ravel();B=B.ravel()
    a=A-1;b=B-1; kk=m-a-b; t=kk*kk-4*a*b
    ok=(kk>=1)&(t>0)&(b-a>=-m)&(b-a<m)&((A+B)>0)
    A,B,kk,t=A[ok],B[ok],kk[ok].astype(float),t[ok].astype(float)
    t0=m*m+4*m; S=(A+B).astype(float)
    Dp=math.sqrt(t0)-np.sqrt(t)-S
    ahi=kk/(m+2)*(t/t0)**(nu/2)*(t0/t)**0.25; alo=kk/(m+2)*(t/t0)**nu
    v=np.maximum(supf(ahi,S,Dp,c1),supf(alo,S,Dp,c1))
    dev=float(np.sum(d[A]*d[B]*v))
    # terms of |P|^2 that are absent from L: pairs outside P_m (all A,B < m+400)
    AA,BB=np.meshgrid(np.arange(m+400),np.arange(m+400),indexing='ij'); AA=AA.ravel();BB=BB.ravel()
    aa=AA-1;bb=BB-1;k2=m-aa-bb;t2=k2*k2-4*aa*bb
    inP=(k2>=1)&(t2>0)&(bb-aa>=-m)&(bb-aa<m)
    miss=float(np.sum((d[AA]*d[BB]*np.exp(-c1*(AA+BB)))[~inP]))
    dev+=miss
    print(N,m,'x1',round(x1,4),'dev_sup',dev,'miss',miss,flush=True)
