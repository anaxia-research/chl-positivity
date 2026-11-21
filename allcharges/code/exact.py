# exact single-centred indices c_R(m,n,l), 0<=l<=2m, via the strip expansion of phi_m
import sys, math, time
from flint import fmpz_poly, fmpz
from fractions import Fraction
from fast import Fast
def series_1d_PA(N,X):
    """1/P_A in Kronecker layout handled by caller: returns dict (a,b)->int for prod (1-p^n y)^-2 (1-p^n/y)^-2 (1-p^n)^4, p=x^N"""
    from collections import defaultdict
    Tp=X//N; P={(0,0):1}
    def mulfac(P,n,ye,power):
        R=defaultdict(int); terms=[];c=1;k=0
        while k*n<=Tp:
            terms.append((k*n,k*ye,c)); k+=1; c=c*(power-k+1)//k*(-1)
        for (a,b),v in P.items():
            for da,db,cc in terms:
                if a+da<=Tp: R[(a+da,b+db)]+=v*cc
        return dict(R)
    for n in range(1,Tp+1):
        P=mulfac(P,n,1,-2);P=mulfac(P,n,-1,-2);P=mulfac(P,n,0,4)
    return {(N*a,b):v for (a,b),v in P.items() if v}
def ginv(N,X,kp2):
    """coefficients of prod (1-x^n)^-kp2 (1-x^{Nn})^-kp2 up to x^X (1/g = x^-1 * this)"""
    a=fmpz_poly([1])
    from fastinit import eta_pow_inv
    e=eta_pow_inv(kp2,X)
    p1=fmpz_poly(e); eN=[0]*(X+1)
    for i,v in enumerate(e):
        if N*i<=X: eN[N*i]=v
    p=(p1*fmpz_poly(eN)).coeffs()[:X+1]
    return [int(v) for v in p]+[0]*(X+1-len(p))
class Exact:
    def __init__(s,N,X,M):
        s.F=Fast(N,X+1,M); s.N=N; s.X=X; s.M=M
        F=s.F
        F.YB+=2*(X+1)//N+4; F.W=2*F.YB+1; F.K=F.YB   # room for 1/P_A, no wrap-around
        t=time.time(); s.Ev=F.numerators(); s.t_num=time.time()-t
        PA=series_1d_PA(N,X+1)
        s.PA=F.enc({k:Fraction(v) for k,v in PA.items()})
        gi=ginv(N,X+1,F.C.kp2)
        s.G=F.enc({(a,0):Fraction(v) for a,v in enumerate(gi) if v})
    def table(s,m):
        """returns dict (n,l)->Fraction c_R(m,n,l) for 0<=l<=2m, n<=X"""
        F=s.F
        R=F.mul(F.mul(s.Ev[m+1],s.PA),s.G)
        p,d=R; c=p.coeffs(); c=c+[0]*((F.X+1)*F.W-len(c))
        out={}
        for a in range(F.X+1):
            row=c[a*F.W:(a+1)*F.W]   # index j <-> y^{j-K}
            n=a-1                     # 1/g carries x^{-1}
            if n<-1 or n>s.X: continue
            for l in range(0,2*m+1):
                tot=0
                for j in range(l+1+F.K,F.W):
                    v=row[j]
                    if v: tot+=(j-F.K-l)*int(v)
                out[(n,l)]=Fraction(tot,int(d))
        return out
if __name__=='__main__':
    N,X,M=map(int,sys.argv[1:4])
    t=time.time(); E=Exact(N,X,M); print('build',time.time()-t)
    from fastinit import FastCHL
    from chl import finite_parts, dstar
    C=FastCHL(N,X,M); Fp=finite_parts(C)
    for m in range(1,M+1):
        T=E.table(m)
        ok=all(T[(n,l)]==Fraction(int(dstar(Fp,N,m,n,l).p),int(dstar(Fp,N,m,n,l).q))*(1 if l%2 else -1) for n in range(1,X-3) for l in range(0,m+1) if 4*m*n-N*l*l>0)
        print(m,ok)
