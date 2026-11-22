# Rigorous upper bound  Lambda(Y) >= log sup_X |1/Phi_g(X+iY)|  for Y=(Yr,Ys,-v) in W0 with det Y > 1/(4t^2).
# Float64; every truncated sum gets an explicit tail bound; result rounded up by 1e-9 relative + 1e-9.
import math, pickle, os
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
ACONST={2:11.5,3:6.7}; C0={2:8,3:4}
_T={}
def table(t):
    if t in _T: return _T[t]
    N={2:400,3:300}[t]
    c=pickle.load(open(os.path.join(HERE,'cD_%d_%d.pkl'%(t,N)),'rb'))
    Dt=max(c); arr=np.zeros(Dt+2); sg=np.ones(Dt+2)
    for d,v in c.items():
        if d>=-1: arr[d+1]=abs(v); sg[d+1]=1.0 if v>0 else -1.0
    _T[t]=(arr,sg,Dt); return _T[t]
def cbound(t,D):  # |c_t(D)| for D>=1 (array ok)
    return ACONST[t]*np.exp(np.pi*np.sqrt(D)/(2*t))+4*np.exp(np.pi*np.sqrt(D)/t)
def lam(x): return -np.log1p(-x)
def cusp(x,w,c0):
    """sum_{k>=1} [c0 lam(e^{-kx}) + 2 lam(e^{-kx+w}) + 2 lam(e^{-kx-w})], 0<=w<x"""
    K=int(60/x)+5
    k=np.arange(1,K+1)
    s=(c0*lam(np.exp(-k*x))+2*lam(np.exp(-k*x+w))+2*lam(np.exp(-k*x-w))).sum()
    # tail k>K: lam(y)<=y/(1-y); each <= 2 e^{-kx+w}; geometric
    y=math.exp(-(K+1)*x+w); s+=(c0+4)*2*y/(1-math.exp(-x))
    return s
def Lambda(Y,t,G=70.0,detail=False):
    Yr,Ys,v=Y
    assert v>0 and t*Yr>v and Ys>v and Yr*Ys-v*v>1/(4*t*t)
    arr,sg,Dt=table(t); c0=C0[t]; Ac=ACONST[t]+4
    W=2*math.pi*(Yr+Ys/t-v)
    root=2*lam(math.exp(-2*math.pi*v))
    cr=cusp(2*math.pi*t*Yr,2*math.pi*v,c0); cs=cusp(2*math.pi*Ys,2*math.pi*v,c0)
    beta=math.sqrt(1/t**2+4*v*v); delta=2*math.sqrt(Yr*Ys)-beta
    assert delta>0
    R=0.0; n=1; rowmax_prev=None
    while True:
        # envelope over m: g(m) = log(2 Ac (2 sqrt(4tnm+1)+1)) + 2 pi beta sqrt(t n m) - 2 pi (t n Yr + m Ys)
        mstar=beta*beta*t*n/(4*Ys*Ys); Mcap=int(6*mstar+300/Ys+200)
        mm=np.arange(1,Mcap+1,dtype=float)
        g=np.log(2*Ac*(2*np.sqrt(4*t*n*mm+1)+1))+2*np.pi*beta*np.sqrt(t*n*mm)-2*np.pi*(t*n*Yr+mm*Ys)
        assert g[-1]<-G-50
        sel=np.nonzero(g>-G)[0]
        rowmax=g.max()
        if sel.size:
            ms=(sel+1).astype(np.int64)
            Lms=np.array([math.isqrt(int(4*t*n*mv+1)) for mv in ms],dtype=np.int64)
            cnt=2*Lms+1
            m=np.repeat(ms,cnt)
            starts=np.cumsum(cnt)-cnt
            l=np.arange(cnt.sum())-np.repeat(starts,cnt)-np.repeat(Lms,cnt)
            D=4*t*n*m-l*l
            ok=D>=-1
            m=m[ok]; l=l[ok]; D=D[ok]
            e=2*np.pi*(t*n*Yr+m*Ys-l*v)
            assert np.all(e>0)
            x=np.exp(-e)
            inT=D<=Dt
            Dc=np.minimum(D,Dt)
            lc=np.where(inT,np.log(np.maximum(arr[Dc+1],1e-300)),0.0)
            lc=np.where(inT&(arr[Dc+1]==0),-np.inf,lc)
            big=~inT
            if big.any():
                Db=np.maximum(D,1).astype(float)
                lb=np.log(ACONST[t]*np.exp(np.pi*np.sqrt(Db)/(2*t)-np.pi*np.sqrt(Db)/t)+4)+np.pi*np.sqrt(Db)/t
                lc=np.where(big,lb,lc)
            sgn=np.where(inT,sg[Dc+1],1.0)
            val=np.where(sgn>0,lam(x),np.log1p(x))
            lv=np.where(e<25,np.log(np.maximum(val,1e-300)),-e+1e-10)
            R+=float(np.exp(lc+lv).sum())
        if rowmax<-G and n>3 and (rowmax_prev is not None and rowmax<rowmax_prev): break
        rowmax_prev=rowmax; n+=1
        assert n<100000
    # tail: all (n,m) with g<=-G (not summed) contribute <= e^{-G/2} * S, S = sum_{n,m} sqrt(2Ac(2sqrt(4tk+1)+1)) e^{-pi delta sqrt(t k)} (k=nm), using lam(x)<=2x for x<=1/2
    K=np.arange(1,2000000,dtype=float)
    # number of (n,m) with nm=k <= 2 sqrt(k)
    S=float((2*np.sqrt(K)*np.sqrt(2*Ac*(2*np.sqrt(4*t*K+1)+1))*np.exp(-np.pi*delta*np.sqrt(t*K))).sum())
    tail=math.exp(-G/2)*S+1e-30
    tot=W+root+cr+cs+R+tail
    tot=tot*(1+1e-9)+1e-9 if tot>0 else tot*(1-1e-9)+1e-9
    if detail: return tot,dict(W=W,root=root,cr=cr,cs=cs,R=R,tail=tail,delta=delta)
    return tot
if __name__=='__main__':
    import time
    for Y in [(0.515,0.1609,0.073),(0.2211,0.4424,0.1475),(0.2304,0.3608,0.1046)]:
        t0=time.time(); v,d=Lambda(Y,2,detail=True); print(Y,v,{k:round(x,4) for k,x in d.items()},'%.2fs'%(time.time()-t0))
