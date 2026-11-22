# explicit upper bound U(T) for log|d_g| along the Y-family with uniform R_max, versus Part I lower bound
import math, sys
import numpy as np
from scipy.special import ive
KK={2:6,3:4}; CM={2:0.410,3:0.202}; C0={2:8,3:4}
def S(x):  # sum_{k>=1} -log(1-e^{-kx}) exact via eta, valid all x>0
    k=np.arange(1,200)
    if x<2*math.pi:
        return math.pi**2/(6*x)-x/24+0.5*math.log(x/(2*math.pi))+float(-np.log1p(-np.exp(-4*math.pi**2*k/x)).sum())
    return float(-np.log1p(-np.exp(-k*x)).sum())
def lam(y): return -math.log1p(-y)
def U(t,m,n,l,kappa,eps,Rmax,r):
    N=t; c0=C0[t]
    P=kappa/(4*t*t)/(1-eps*eps*t)
    Yr=math.sqrt(P*r); Ys=math.sqrt(P/r); v=eps*min(t*Yr,Ys)
    lin=2*math.pi*((m+1)*Yr+(n+1)/N*Ys-(l+1)*v)
    root=2*lam(math.exp(-2*math.pi*v))
    xr=2*math.pi*t*Yr; xs=2*math.pi*Ys; w=2*math.pi*v
    cusp=(c0+4)*(S(xr)+S(xs))+2*lam(math.exp(-(xr-w)))+2*lam(math.exp(-(xs-w)))
    return lin+root+cusp+Rmax
def Ubest(t,m,n,l,kappa,eps,Rmax):
    # golden search on log r
    a,b=-12.0,12.0
    f=lambda u: U(t,m,n,l,kappa,eps,Rmax,math.exp(u))
    for _ in range(80):
        c=a+(b-a)*0.382; d=a+(b-a)*0.618
        if f(c)<f(d): b=d
        else: a=c
    return f((a+b)/2)
def logLower(N,m,D):
    k=KK[N]; nu=k+1.5; t0=m*m+4*m; X=math.pi*math.sqrt(D*t0)/m
    return math.log(CM[N]*0.44)+X-0.5*math.log(2*math.pi*X)+math.log(2*math.pi*(m+2)/(N**((k+4)/2)*math.sqrt(2*m)))+(nu/2)*math.log(t0/D)
if __name__=='__main__':
    t=int(sys.argv[1]); kappa=float(sys.argv[2]); eps=float(sys.argv[3]); Rmax=float(sys.argv[4]); N=t
    extra=math.log(2) if t==3 else 0.0
    res={}
    for m in sorted(set([1,2,3,4,5,6,8,10,13,16,20,25,32,40,50,64,80,100,130,160,200,250,320,400,500,640,800,1000,1300,1700])):
        for l in sorted(set([0,m//2,m])):
            for D in [400,1000,2500,5000,10**4,3*10**4,10**5,3*10**5,10**6,10**7]:
                n=math.ceil(N*(D+l*l)/(4*m))
                if n<m: continue
                Dv=4*m*n/N-l*l
                if Dv<3*m*m/N**2: continue
                g=Ubest(t,m,n,l,kappa,eps,Rmax)+extra-logLower(N,m,Dv)
                key=D; res.setdefault(key,[]).append((g,(m,n,l)))
    for D in sorted(res):
        w=max(res[D]); print(D,'worst G=%.2f at'%w[0],w[1],'  (G<0 needed)')
