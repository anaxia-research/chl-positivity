# Upper bound for the leading-family sum L, uniform over reduced charges with D > D0:
#  (i)  m <= 400: |L| <= sum_{(A,B) in P_m} d~(A)d~(B) rho_AB(x_min(m)), since each rho_AB decreases in x;
#  (ii) m >= 300: rho_AB <= exp(-c S phi), phi=(2m-S)/(2m+4) >= phi_*(S), giving a bound independent of m.
# Usage: upper.py N D0
import sys, math, numpy as np
N=int(sys.argv[1]); D0=float(sys.argv[2])
sys.argv=[sys.argv[0],str(N)]
from cornerL import pairs, d, nu, ive_halfint
x0=math.sqrt(3)/N
worst=(0,None)
for m in range(1,401):
    x=max(x0,math.sqrt(D0)/m); c=math.pi*x
    A,B,kk,t=pairs(m,m+3); t0=m*m+4*m
    st=np.sqrt(t.astype(float)); s0=math.sqrt(t0)
    ut,inv=np.unique(t,return_inverse=True)
    iv=ive_halfint(nu,c*np.sqrt(ut.astype(float)))[inv]; i0=float(ive_halfint(nu,[c*s0])[0])
    rho=kk/(m+2)*(t/t0)**(nu/2)*iv/i0*np.exp(c*(st-s0))
    U=float(np.sum(np.abs(d[A]*d[B])*rho))
    if U>worst[0]: worst=(U,m)
print(N,'m<=400 upper',worst,flush=True)
from exact import ginv
kp2=24//(N+1); dd=[float(v) for v in ginv(N,1200,kp2)]
c0=math.pi*x0; tot=0.0
for S in range(0,1200):
    W=sum(dd[A]*dd[S-A] for A in range(S+1))
    phi=(600-S)/604 if S<=300 else 0.49
    tot+=W*math.exp(-c0*S*phi)
print(N,'m>=300 upper',tot)
