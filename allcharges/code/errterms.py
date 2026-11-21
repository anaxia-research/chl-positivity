# Upper bounds (relative to the extreme term E(D,m)) for the non-leading pieces of the level-N expansion,
# uniform in x >= x0 and D >= D0.  Float evaluation; all quantities are sums of positive terms.
import math, numpy as np
from exact import ginv
from cornerL import ive_halfint
def setup(N,dmax=4000):
    kp2=24//(N+1); k=kp2-2; nu=k+1.5
    d=np.array([float(v) for v in ginv(N,dmax,kp2)])   # d[j]=coeff of q^{j-1}
    return kp2,k,nu,d
def dval(d,a):
    """d(a) for a in Z>=-1, array"""
    a=np.asarray(a); out=np.zeros(a.shape)
    ok=(a>=-1)&(np.abs(a-np.round(a))<1e-9)
    out[ok]=d[np.round(a[ok]).astype(int)+1]
    return out
from math import lgamma
def _logI_series(nu,z):
    z=np.asarray(z,dtype=float); K=int(np.max(z, initial=0)+60)
    k=np.arange(K)[:,None]; lz=np.log(z/2)[None,:]
    lt=(2*k)*lz-np.array([lgamma(i+1)+lgamma(nu+i+1) for i in range(K)])[:,None]
    mx=lt.max(0); return nu*np.log(z/2)+mx+np.log(np.exp(lt-mx).sum(0))
def logI(nu,z):   # log I_nu(z) for z>0, stable
    z=np.atleast_1d(np.asarray(z,dtype=float)); out=np.empty_like(z)
    sm=z<2*nu+20
    if sm.any(): out[sm]=_logI_series(nu,z[sm])
    if (~sm).any(): out[~sm]=np.log(ive_halfint(nu,z[~sm]))+z[~sm]
    return out
def Bgam(nu,t,t0,x,g):
    """(t/t0)^{nu/2} I(pi x sqrt t/g)/I(pi x sqrt t0)"""
    c=math.pi*x
    return np.exp((nu/2)*np.log(t/t0)+logI(nu,c*np.sqrt(t)/g)-logI(nu,np.array([c*math.sqrt(t0)]))[0])
def walls_S(m):
    """all frame-S walls (r,s,a,b,t): a,b in Z>=-1, a<=m/s^2, b<=m/r^2, r,s<=2m coprime"""
    out=[]
    for r in range(1,2*m+1):
        for s in range(1,2*m+1):
            if math.gcd(r,s)!=1: continue
            amax=m//(s*s); bmax=m//(r*r)
            a,b=np.meshgrid(np.arange(-1,amax+1),np.arange(-1,bmax+1),indexing='ij')
            a=a.ravel().astype(float); b=b.ravel().astype(float)
            out.append(np.stack([np.full(a.shape,r),np.full(a.shape,s),a,b]))
    W=np.concatenate(out,axis=1)
    return finish(W,m)
def walls_T(m,N):
    out=[]; R=int(math.ceil(2*m*math.sqrt(N)))+2
    for r in range(1,R+1):
        for s in range(1,R+1):
            if math.gcd(r,s)!=1 or (r%N and s%N): continue
            amax=m/(s*s); bmax=m/(r*r)
            for typ in ((0,) if r%N==0 else ())+((1,) if s%N==0 else ()):
                if typ==0:   # a integer, b in Z/N, weight d(a) d(N b)
                    A=np.arange(-1,math.floor(amax)+1); B=np.arange(-1,math.floor(N*bmax)+1)/N
                else:
                    A=np.arange(-1,math.floor(N*amax)+1)/N; B=np.arange(-1,math.floor(bmax)+1)
                a,b=np.meshgrid(A,B,indexing='ij'); a=a.ravel(); b=b.ravel()
                out.append(np.stack([np.full(a.shape,r),np.full(a.shape,s),a,b,np.full(a.shape,typ)]))
    W=np.concatenate(out,axis=1)
    return finish(W,m)
def finish(W,m):
    r,s,a,b=W[0],W[1],W[2],W[3]
    al=a/(m*r*r); be=b/(m*s*s); g=1/(r*s)
    p=(g+(al-be)/g)/2; q=g-p
    d2=p*p-al
    ok=(p>=0)&(q>=-1e-15)&(d2>0)&(a+b<m)
    t=4*m*m*d2
    return W[:,ok],t[ok]
def eps_gamma1(N,m,x0,WS,WT,kp2,nu,d,lead_excluded=True):
    t0=m*m+4*m
    W,t=WS; r,s,a,b=W[0],W[1],W[2],W[3]
    sel=(r*s>=2) if lead_excluded else np.ones(r.shape,bool)
    w=(m-a-b)*dval(d,a)*dval(d,b)
    eS=float(np.sum((w*Bgam(nu,t,t0,x0,1))[sel]))/(m+2)
    W,t=WT; r,s,a,b,typ=W
    w=(m-a-b)*np.where(typ==0,dval(d,a)*dval(d,N*b),dval(d,N*a)*dval(d,b))
    eT=(N-1)*N**(kp2/2)*float(np.sum(w*Bgam(nu,t,t0,x0,1)))/(m+2)
    return eS,eT
def eps_gamma2(N,m,x0,WS,WT,kp2,nu,d,G=40):
    t0=m*m+4*m; tot=0.0
    W,t=WS; a,b=W[2],W[3]; wS=(m-a-b)*dval(d,a)*dval(d,b)
    W2,t2=WT; a2,b2,typ=W2[2],W2[3],W2[4]
    wT=(m-a2-b2)*np.where(typ==0,dval(d,a2)*dval(d,N*b2),dval(d,N*a2)*dval(d,b2))
    for g in range(2,G+1):
        SS=float(np.sum(wS*Bgam(nu,t,t0,x0,g)))/(m+2)
        ST=N**(kp2/2)*float(np.sum(wT*Bgam(nu,t2,t0,x0,g)))/(m+2)
        # N_S + N_T <= N g ; |sqrt(2m) W| <= sqrt(2 gcd(m,g))
        tot+=math.sqrt(2*math.gcd(m,g))/g*N*g*max(SS,ST)
    # tail g>G: B_g <= g^{-nu} B_1 ... use B_g<=B_G*(G/g)^nu
    SS=float(np.sum(wS*Bgam(nu,t,t0,x0,G)))/(m+2); ST=N**(kp2/2)*float(np.sum(wT*Bgam(nu,t2,t0,x0,G)))/(m+2)
    tail=sum(math.sqrt(2*g)*N*max(SS,ST)*(G/g)**nu for g in range(G+1,20000))
    return tot+tail
def Eext(N,m,D,x,kp2,nu):
    t0=m*m+4*m
    return (2/N)*math.pi/math.sqrt(2*m)*(m+2)*N**(-kp2/2)*(t0/D)**(nu/2)*math.exp(logI(nu,np.array([math.pi*x*math.sqrt(t0)]))[0])
def eps_K(N,m,x0,kp2,k,nu,d,G=40):
    """bound for |K2|+|K3| /(N E) at x=x0 (worst case), D=(x0 m)^2"""
    D=(x0*m)**2; E=Eext(N,m,D,x0,kp2,nu); fh=d[m+1]
    tot=0.0
    for g in range(1,G+1):
        Wm=(2*m)**-0.5 if g==1 else math.sqrt(math.gcd(m,g)/m)
        for (cF,mu,cnt) in ((fh/N**(kp2/2),1.0,g),(fh,1.0/N,(N-1)*g+ (0))):
            # K2
            tot+=cnt*math.sqrt(2*m)*g**-0.5*cF*Wm*(4*m*mu/D)**((k+2)/2)*besI(k+2,2*math.pi/g*math.sqrt(mu*D/m))
            # K3
            tot+=cnt*(1/(2*math.pi))*g**-2*cF*Wm*(4*m*mu/D)**((2*k+5)/4)*(4*math.pi**2*m*m*g*g/3)*2*math.sqrt(mu/m)*besI(k+2.5,2*math.pi/g*math.sqrt(mu*D/m))
    return tot/(N*E)
import mpmath
def besI(n,z): return float(mpmath.besseli(n,z))
