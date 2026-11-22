# Exact twined single-centred indices d_{*,g}(m,n,l) for the Z2xZ2 (t=2) and Z3xZ3 (t=3) models:
#   1/Phi_g = q^-1 r^-1 s^{-1/t} (1-r^-1)^-2 prod_{m>=1}(1-s^m)^-c0 (1-s^m r)^-2 (1-s^m/r)^-2 * exp(sum_L q^{tL} H_L)
# coefficient of q^m s^{n/t} r^l (m=tL-1, n=tB-1), expanded in the chamber W0 (powers of r^-1).
import sys, os, pickle, json
from math import isqrt
from flint import fmpz_poly
from cD import build
HERE=os.path.dirname(os.path.abspath(__file__))
class Tw:
    def __init__(s, t, Lmax, Bmax):
        s.t=t; s.Lmax=Lmax; s.Bmax=Bmax
        J=2*isqrt(4*t*Lmax*Bmax+4)+Lmax+Bmax+8
        s.J=J; s.W=4*J+1; s.off=2*J
        need=4*t*Lmax*Bmax+4
        Nm=need//(4*t)+2
        fn=os.path.join(HERE,'cD_%d_%d.pkl'%(t,Nm))
        s.c=pickle.load(open(fn,'rb')) if os.path.exists(fn) else build(t,Nm)
        if not os.path.exists(fn): pickle.dump(s.c,open(fn,'wb'))
    def cc(s,D): return s.c.get(D,0) if D>=-1 else 0
    def poly(s,d):
        # d: dict (B,j)->int  -> Kronecker poly  (j in [-off,off])
        arr=[0]*((s.Bmax+1)*s.W)
        for (B,j),v in d.items():
            if B<=s.Bmax and -s.off<=j<=s.off: arr[B*s.W+j+s.off]+=v
        return fmpz_poly(arr)
    def unpoly(s,p,shift=1):
        # p was a product of `shift` polys: offsets add
        out={}
        for i,v in enumerate(p.coeffs()):
            if v==0: continue
            B,rem=divmod(i-(shift-1)*s.off, s.W)
            if B>s.Bmax: continue
            j=rem-s.off
            out[(B,j)]=int(v)
        return out
    def mul(s,a,b):
        return s.unpoly(s.poly(a)*s.poly(b),2)
    def run(s):
        t=s.t; L=s.Lmax; Bm=s.Bmax
        # L*H_L as dict (B,j)
        LW={}
        for Lk in range(1,L+1):
            d={}
            for k in range(1,Lk+1):
                if Lk%k: continue
                Lq=Lk//k
                for m in range(0,Bm//k+1):
                    lm=isqrt(4*t*Lq*m+1)
                    for l in range(-lm,lm+1):
                        v=s.cc(4*t*Lq*m-l*l)
                        if v and m*k<=Bm and abs(l*k)<=s.off:
                            d[(m*k,l*k)]=d.get((m*k,l*k),0)+Lq*v
            LW[Lk]=d
        E={0:{(0,0):1}}
        for Lk in range(1,L+1):
            acc=fmpz_poly([])
            for k in range(1,Lk+1):
                acc+=s.poly(LW[k])*s.poly(E[Lk-k])
            d=s.unpoly(acc,2)
            for key in d: 
                assert d[key]%Lk==0; d[key]//=Lk
            E[Lk]={k:v for k,v in d.items() if v}
        # P = prod_{m>=1} (1-s^m)^-c0 (1-s^m r)^-2 (1-s^m/r)^-2  (positive-power series in s)
        c0=s.cc(0); c1=s.cc(-1)
        P={(0,0):1}
        def facinv(P,B,j,e):
            for _ in range(e):
                Q=dict(P)
                for (b,jj) in sorted(P, key=lambda x:x[0]):
                    pass
                # multiply by 1/(1-x): Q[k] = P[k] + Q[k-x] processed in increasing B
                keys=sorted(set(P)|set(), key=lambda x:x[0])
                R={}
                for b in range(0,Bm+1):
                    for jj in range(-s.off,s.off+1):
                        v=P.get((b,jj),0)+R.get((b-B,jj-j),0) if b-B>=0 else P.get((b,jj),0)
                        if v: R[(b,jj)]=v
                P=R
            return P
        # faster: build P via flint as 1d series per needed j using log/exp is overkill; do direct with dict-of-rows
        P=s.Pseries(c0)
        self_res={}
        for Lk in range(1,L+1):
            G=s.mul(E[Lk],P)
            self_res[Lk]=G
        s.G=self_res
        return self_res
    def Pseries(s,c0):
        # prod_{m>=1} (1-s^m)^-c0 (1-s^m r)^-2 (1-s^m/r)^-2 : rows B -> numpy object array over j
        import numpy as np
        Bm=s.Bmax; off=s.off; Wd=2*off+1
        rows=np.zeros((Bm+1,Wd),dtype=object); rows[0,off]=1
        for m in range(1,Bm+1):
            for (j,e) in ((0,c0),(1,2),(-1,2)):
                for _ in range(e):
                    for B in range(m,Bm+1):
                        src=rows[B-m]
                        if j==0: rows[B]+=src
                        elif j==1: rows[B,1:]+=src[:-1]
                        else: rows[B,:-1]+=src[1:]
        out={}
        nz=np.nonzero(rows)
        for B,i in zip(*nz): out[(int(B),int(i)-off)]=int(rows[B,i])
        return out
    def d(s,m,n,l):
        t=s.t
        if (m+1)%t or (n+1)%t: return 0
        Lk=(m+1)//t; B=(n+1)//t
        G=s.G[Lk]
        tot=0
        for (b,j),v in G.items():
            if b==B and j>=l+1: tot+=(j-l)*v
        return (-1)**(l+1)*tot
if __name__=='__main__':
    t=int(sys.argv[1])
    T=Tw(t,2,3); T.run()
    for (m,n) in [(1,1),(1,3),(3,3),(1,2)] if t==2 else [(2,2),(2,5),(5,5)]:
        print((m,n),[T.d(m,n,l) for l in range(0,m+2)])
