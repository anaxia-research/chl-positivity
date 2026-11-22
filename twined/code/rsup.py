# Uniform bound for R(Y(r)) over an r-interval [r1,r2] (r2 may be inf) of the family Yr=sqrt(P r), Ys=sqrt(P/r),
# v = eps*sqrt(P)*min(t sqrt r, 1/sqrt r).  Per (n,m): tnYr+mYs >= sqrt(P)*min_{r in I}(tn sqrt r + m/sqrt r).
import math, numpy as np
from lam import table, ACONST, lam
def vrange(P,eps,t,r1,r2):
    f=lambda r: eps*math.sqrt(P)*min(t*math.sqrt(r),1/math.sqrt(r)) if r<math.inf else 0.0
    rc=1.0/t  # peak of min(t sqrt r, 1/sqrt r)
    vmax=f(rc) if r1<=rc<=r2 else max(f(r1),f(r2))
    vmin=min(f(r1),f(r2) if r2<math.inf else 0.0)
    return vmin,vmax
def Rsup(t,P,eps,r1,r2,G=70.0):
    arr,sg,Dt=table(t); Ac=ACONST[t]+4
    vmin,vmax=vrange(P,eps,t,r1,r2)
    beta=math.sqrt(1/t**2+4*vmax*vmax); sP=math.sqrt(P)
    delta=2*sP-beta
    assert delta>0
    R=0.0; n=1; prev=None
    while True:
        mstar=beta*beta*t*n/(4*P/ (math.sqrt(P/r2) if r2<math.inf else 1e-9)**2 if False else 1)  # unused
        # minimal tr per m: h(m) = sqrt(P)*min_{r in [r1,r2]} (t n sqrt r + m/sqrt r)
        Mcap=int(1e5)
        mm=np.arange(1,Mcap+1,dtype=float)
        ropt=mm/(t*n)
        rc=np.clip(ropt,r1,r2)
        h=sP*(t*n*np.sqrt(rc)+mm/np.sqrt(rc))
        g=np.log(2*Ac*(2*np.sqrt(4*t*n*mm+1)+1))+2*np.pi*beta*np.sqrt(t*n*mm)-2*np.pi*h
        assert g[-1]<-G-20, (n,g[-1])
        sel=np.nonzero(g>-G)[0]; rowmax=g.max()
        if sel.size:
            ms=(sel+1).astype(np.int64); hs=h[sel]
            Lms=np.array([math.isqrt(int(4*t*n*mv+1)) for mv in ms],dtype=np.int64)
            cnt=2*Lms+1
            m=np.repeat(ms,cnt); H=np.repeat(hs,cnt)
            starts=np.cumsum(cnt)-cnt
            l=np.arange(cnt.sum())-np.repeat(starts,cnt)-np.repeat(Lms,cnt)
            D=4*t*n*m-l*l; ok=D>=-1
            m=m[ok]; l=l[ok]; D=D[ok]; H=H[ok]
            e=2*np.pi*(H-np.where(l>0,l*vmax,l*vmin))
            assert np.all(e>0)
            x=np.exp(-e); inT=D<=Dt; Dc=np.minimum(D,Dt)
            lc=np.where(inT,np.log(np.maximum(arr[Dc+1],1e-300)),0.0)
            lc=np.where(inT&(arr[Dc+1]==0),-np.inf,lc)
            big=~inT
            if big.any():
                Db=np.maximum(D,1).astype(float)
                lb=np.log(ACONST[t]*np.exp(-np.pi*np.sqrt(Db)/(2*t))+4)+np.pi*np.sqrt(Db)/t
                lc=np.where(big,lb,lc)
            sgn=np.where(inT,sg[Dc+1],1.0)
            val=np.where(sgn>0,lam(x),np.log1p(x))
            lv=np.where(e<25,np.log(np.maximum(val,1e-300)),-e+1e-10)
            R+=float(np.exp(lc+lv).sum())
        if rowmax<-G and n>3 and prev is not None and rowmax<prev: break
        prev=rowmax; n+=1
    K=np.arange(1,4000000,dtype=float)
    S=float((2*np.sqrt(K)*np.sqrt(2*Ac*(2*np.sqrt(4*t*K+1)+1))*np.exp(-np.pi*delta*np.sqrt(t*K))).sum())
    return R+math.exp(-G/2)*S+1e-12, delta
if __name__=='__main__':
    import sys
    t=int(sys.argv[1]); kappa=float(sys.argv[2]); eps=float(sys.argv[3])
    P=kappa/(4*t*t)/(1-eps*eps*t)
    edges=[0.0]+[math.exp(x/8) for x in range(-40,41)]+[math.inf]
    worst=0; out=[]
    for r1,r2 in zip(edges,edges[1:]):
        Rv,dl=Rsup(t,P,eps,max(r1,1e-12),r2)
        out.append((r1,r2,Rv)); worst=max(worst,Rv)
    print('t',t,'kappa',kappa,'eps',eps,'P',P,'Rmax over all r',worst)
    print([ (round(a,3),round(c,3)) for a,b,c in out[::8]])
