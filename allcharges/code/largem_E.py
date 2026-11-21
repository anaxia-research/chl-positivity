# E1 via Lagrange remainder with ball-arithmetic jets on real h-intervals.
# |F0(h)-T_J(h)| <= sup_{xi in [0,h*]} |F0^{(J+1)}(xi)/(J+1)!| h^{J+1}
import sys, math, numpy as np
from flint import arb, arb_series, ctx
from math import factorial
N=int(sys.argv[1]); Ms=int(sys.argv[2]); S1=int(sys.argv[3]); J=int(sys.argv[4]); XA=float(sys.argv[5]); XB=float(sys.argv[6]); W=float(sys.argv[7]); nsub=int(sys.argv[8])
sys.argv=[sys.argv[0],str(N)]
from exact import ginv
kp2=24//(N+1); k=kp2-2; nu2=2*k+3; n=(nu2-1)//2; al=arb(nu2)/4-arb(1)/4
dd=[int(v) for v in ginv(N,3000,kp2)]
pc=[arb(factorial(n+i))/(factorial(i)*factorial(n-i))*arb(-0.5)**i for i in range(n+1)]
ctx.cap=J+2
def p(w):
    s=0*w+1; wk=0*w+1
    for i in range(1,n+1): wk=wk*w; s=s+pc[i]*wk
    return s
def coefJ1(A,B,xi):
    S=A+B; y2=(A-B)**2
    h=arb_series([xi,1])
    T=1+(4-2*S)*h+y2*h*h; T0=1+4*h
    sT=T.sqrt(); sT0=T0.sqrt()
    E=(y2*h-2*S)/(sT+sT0)+S
    F=(1-S*h/(1+2*h))*(al*(T/T0).log()).exp()*(cb*E).exp()*p(h/(cb*sT))/p(h/(cb*sT0))
    return F.coeffs()[J+1]
import json
ad=np.abs(np.array([float(v) for v in dd])); conv=[float(np.sum(ad[:S+1]*ad[S::-1])) for S in range(len(dd))]
out=[]
xa=XA
while xa<XB-1e-12:
    xb=min(xa+W,XB); c0=math.pi*xa; c1=math.pi*xb; cb=arb((c0+c1)/2,(c1-c0)/2)
    hs=1.0/Ms
    E1=0.0; worst=[]
    for A in range(0,S1+1):
        for B in range(A,S1+1-A if False else S1+1):
            if A+B>S1: break
            w=abs(dd[A]*dd[B])*math.exp(-c0*(A+B))
            if w==0: continue
            sup=0.0
            for i in range(nsub):
                xi=arb((i+0.5)*hs/nsub, 0.5*hs/nsub*1.0001)
                sup=max(sup,float(abs(coefJ1(A,B,xi)).upper()))
            b=w*sup*hs**(J+1)*(1 if A==B else 2)
            E1+=b; worst.append((b,A,B))
    worst.sort(reverse=True)
    
    E2=sum(conv[S]*math.exp(-c0*S*((2*Ms-S)/(2*Ms+4) if S<=Ms else 0.49)) for S in range(S1+1,len(dd)))
    out.append(dict(xa=xa,xb=xb,E1=E1,E2=E2)); print(xa,xb,'E1',E1,'E2',E2,flush=True)
    json.dump(out,open(f'Ecells_{N}_{Ms}_{S1}_{J}_{XA:.3f}.json','w'))
    xa=xb
