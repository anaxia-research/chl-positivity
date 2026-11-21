# Explicit upper bound U(m) for all non-leading pieces (walls rs>=2, frame-T walls, gamma>=2, K2/K3),
# relative to the extreme term, uniform in x >= x0 = sqrt3/N.  Usage: uniform.py N m0 m1 [step]
import sys, math, numpy as np
N=int(sys.argv[1]); m0=int(sys.argv[2]); m1=int(sys.argv[3]); step=int(sys.argv[4]) if len(sys.argv)>4 else 1
sys.argv=[sys.argv[0],str(N)]
from exact import ginv
import mpmath
kp2=24//(N+1); k=kp2-2; nu=k+1.5
X=N*m1+10
d=np.array([float(v) for v in ginv(N,X,kp2)]); cum=np.concatenate([[1.0],1.0+np.cumsum(d)])  # cum[i]=Db(i-1)
def Db(Xv):   # 1 + sum_{a=0}^{floor X} d(a)
    i=np.floor(np.asarray(Xv,dtype=float)).astype(int); return cum[np.clip(i,-1,X)+1]
kap0=kp2*(1+1/N)*math.pi**2/6
def Fb(lam):  # 1 + P(e^-lam) upper bound
    lam=np.asarray(lam,dtype=float); q=np.exp(-lam); return 1+np.exp(np.minimum(kap0*q/(1-q),700))
c=math.pi*math.sqrt(3)/N
def U(m,K=4):
    s3=3*math.sqrt(m); base=m+2-2/m; R=int(math.ceil(2*m*math.sqrt(N)))+2
    rs_=np.arange(1,2*m+1); RS=np.arange(1,R+1)
    SigS=float(np.sum(Db(m/rs_**2))); SigS2=float(np.sum(Db(m/rs_[1:]**2)))
    SigR=float(np.sum(Db(m/RS**2))); SigRN=float(np.sum(Db(N*m/RS[RS%N==0]**2)))
    small=[(r,s) for r in range(1,K*N+1) for s in range(1,K*N+1) if math.gcd(r,s)==1 and r*s<K*N]
    def Sfr(gam,first):   # frame S walls
        tot=0.0
        for r,s in small:
            if r*s>=K: continue
            if first and r*s<2: continue
            g=1/(r*s)
            tot+=math.exp(-c*(base-(m*g+s3)/gam))*min(Fb(c*s/(r*gam)),Db(m/s**2))*min(Fb(c*r/(s*gam)),Db(m/r**2))
        tot+=math.exp(-c*(base-(m/K+s3)/gam))*2*SigS*SigS2
        return tot
    def Tfr(gam):         # frame T walls (without the N^{kp2/2} weight)
        tot=0.0; KT=K*N
        for r,s in small:
            if r*s>=KT: continue
            g=1/(r*s)
            if r%N==0: tot+=math.exp(-c*(base-(m*g+s3)/gam))*min(Fb(c*s/(r*gam)),Db(m/s**2))*min(Fb(c*r/(N*s*gam)),Db(N*m/r**2))
            if s%N==0: tot+=math.exp(-c*(base-(m*g+s3)/gam))*min(Fb(c*s/(N*r*gam)),Db(N*m/s**2))*min(Fb(c*r/(s*gam)),Db(m/r**2))
        tot+=math.exp(-c*(base-(m/KT+s3)/gam))*2*SigR*SigRN
        return tot
    wN=N**(kp2/2)
    e1=Sfr(1,True)+(N-1)*wN*Tfr(1)
    tau=m+s3; Gam=int(math.ceil(c*tau))
    e2=sum(N*math.sqrt(2)*gam*(Sfr(gam,False)+wN*Tfr(gam)) for gam in range(2,Gam+1))
    y=c*math.sqrt(m*m+4*m)
    Wtot=SigS**2+wN*2*SigR*SigRN
    e3=N*math.sqrt(2)*math.exp(c*tau/Gam)*(c*tau/2)**nu*Gam**(1.5-nu)/(nu-1.5)*Wtot*math.sqrt(2*math.pi*y)*math.exp(-y)/(0.9*math.gamma(nu+1))
    return e1,e2,e3
def UK(m,G=400):
    # K2/K3 pieces as in errterms.eps_K, x=x0 (worst case), gamma<=G, plus power-law tail
    x0=math.sqrt(3)/N; D=(x0*m)**2
    t0=m*m+4*m
    E=(2/N)*math.pi/math.sqrt(2*m)*(m+2)*N**(-kp2/2)*(t0/D)**(nu/2)*float(mpmath.besseli(nu,math.pi*x0*math.sqrt(t0)))
    fh=float(d[m+1]); tot=0.0; last=0.0
    for g in range(1,G+1):
        Wm=(2*m)**-0.5 if g==1 else math.sqrt(math.gcd(m,g)/m)
        term=0.0
        for (cF,mu,cnt) in ((fh/N**(kp2/2),1.0,g),(fh,1.0/N,(N-1)*g)):
            z=2*math.pi/g*math.sqrt(mu*D/m)
            term+=cnt*math.sqrt(2*m)*g**-0.5*cF*Wm*(4*m*mu/D)**((k+2)/2)*float(mpmath.besseli(k+2,z))
            term+=cnt*(1/(2*math.pi))*g**-2*cF*Wm*(4*m*mu/D)**((2*k+5)/4)*(4*math.pi**2*m*m*g*g/3)*2*math.sqrt(mu/m)*float(mpmath.besseli(k+2.5,z))
        tot+=term; last=term
    # for g>G the g-dependence is at most g^{-(k+2)+1/2+... } ; bound tail by last*G*... (decay exponent >=1.5)
    tail=last*G/0.5
    return (tot+tail)/(N*E)
for m in range(m0,m1+1,step):
    e1,e2,e3=U(m); eK=UK(m); print(N,m,'gamma1',e1,'gamma>=2',e2,'tail',e3,'K',eK,'U',e1+e2+e3+eK,flush=True)
