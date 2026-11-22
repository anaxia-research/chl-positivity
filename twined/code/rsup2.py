# Rigorous uniform bound for R(Y(r)) over r-cells [r1,r2] of the family (exhaustive version).
# All pairs (n,m) with nm <= K are enumerated; pairs with g > -G are summed over every l with exact |c| (or Lemma A),
# pairs with g <= -G contribute <= e^{-G/2} e^{g/2}; pairs with nm > K are bounded with the crude envelope.
import math, sys, numpy as np
from lam import table, ACONST, lam
from rsup import vrange
def Rsup2(t,P,eps,r1,r2,G=60.0,K=None):
    arr,sg,Dt=table(t); Ac=ACONST[t]+4
    vmin,vmax=vrange(P,eps,t,r1,r2)
    beta=math.sqrt(1/t**2+4*vmax*vmax); sP=math.sqrt(P); delta=2*sP-beta
    assert delta>0
    def env(k): # log of 2*Ac*(2 sqrt(4tk+1)+1) e^{-2 pi delta sqrt(tk)}  (bound for one pair with nm=k, all l)
        return np.log(2*Ac*(2*np.sqrt(4*t*k+1)+1))-2*np.pi*delta*np.sqrt(t*k)
    if K is None:
        K=1
        while env(float(K))+0.5*math.log(K)>-G-10: K*=2
        K*=4
    R=0.0; low=0.0; npairs=0; nsum=0
    for n in range(1,K+1):
        M=K//n
        mm=np.arange(1,M+1,dtype=float)
        rc=np.clip(mm/(t*n),r1,r2)
        h=sP*(t*n*np.sqrt(rc)+mm/np.sqrt(rc))
        g=np.log(2*Ac*(2*np.sqrt(4*t*n*mm+1)+1))+2*np.pi*beta*np.sqrt(t*n*mm)-2*np.pi*h
        npairs+=M
        lo=g<=-G
        low+=float(np.exp(g[lo]/2).sum())
        sel=np.nonzero(~lo)[0]
        if sel.size:
            nsum+=sel.size
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
    lowb=math.exp(-G/2)*low
    # pairs with nm>K: number of pairs with nm=k is d(k) <= 2 sqrt(k)
    k=np.arange(K+1,K+4000001,dtype=float)
    far=float((2*np.sqrt(k)*np.exp(env(k))).sum())
    # beyond K+4e6: env decays like e^{-2 pi delta sqrt(tk)}; bounded by integral, negligible (checked < 1e-100)
    assert env(K+4e6)+math.log(K+4e6)<-300
    return R+lowb+far+1e-12, dict(K=K,delta=delta,pairs=npairs,summed=nsum,R=R,low=lowb,far=far)
if __name__=='__main__':
    t=int(sys.argv[1]); kappa=float(sys.argv[2]); eps=float(sys.argv[3])
    P=kappa/(4*t*t)/(1-eps*eps*t)
    edges=[0.0]+[math.exp(x/8) for x in range(-40,41)]+[math.inf]
    worst=0
    for r1,r2 in zip(edges,edges[1:]):
        Rv,info=Rsup2(t,P,eps,max(r1,1e-12),r2)
        worst=max(worst,Rv)
        print('cell [%.4g,%.4g] R<=%.4f'%(r1,r2,Rv),{a:(round(b,6) if isinstance(b,float) else b) for a,b in info.items()},flush=True)
    print('t',t,'kappa',kappa,'eps',eps,'P',P,'Rmax over all r',worst)
