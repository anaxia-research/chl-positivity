# exact termwise tail  sum_{S>S1} |d~(A) d~(B)| rho_AB(m,x_min) for each m (rho decreasing in x)
import sys, math, numpy as np
N=int(sys.argv[1]); S1=int(sys.argv[2]); m0=int(sys.argv[3]); m1=int(sys.argv[4])
from cornerL import pairs, d, nu
from errterms import logI
x0=math.sqrt(3)/N
for m in range(m0,m1+1):
    x=max(x0,math.sqrt(float(__import__("os").environ.get("D0","400")))/m); c=math.pi*x
    A,B,kk,t=pairs(m,m+3); sel=(A+B>S1); A,B,kk,t=A[sel],B[sel],kk[sel],t[sel]
    t0=m*m+4*m
    lr=np.log(kk/(m+2))+(nu/2)*np.log(t/t0)+logI(nu,c*np.sqrt(t.astype(float)))-logI(nu,[c*math.sqrt(t0)])[0]
    print(m,'x_min',round(x,5),'tail',float(np.sum(np.abs(d[A]*d[B])*np.exp(lr))),flush=True)
