# frame S (cusp 0): second product Phi^S, its Fourier-Jacobi coefficients N^{(k+2)/2} phi_m|S
import sys, math, time
from fractions import Fraction
from flint import fmpz_poly
from fast import Fast
from exact import series_1d_PA, ginv, Exact
from fastinit import theta_tables, c01_table
from chl import EN_coeffs
class FastS(Fast):
    def __init__(s,N,X,M):
        s.N,s.X,s.M=N,X,M
        s.YB=int(math.isqrt(4*(M+2)*X+1))+M+4 + 2*(X+1)+4
        s.W=2*s.YB+1; s.K=s.YB
        Tq=X*(M+1)+3
        al=theta_tables(Tq); c01=c01_table(al,Tq)
        e=[int(v) for v in EN_coeffs(N,Tq+1)]
        # H(D): coefficient of phi_{-2,1} E_N at D=4n-j^2, read at j in {0,1}
        H={}
        for b in (0,1):
            a=[al.get(4*n-b,0) for n in range(Tq+1)]
            h=fmpz_poly(a)*fmpz_poly(e[:Tq+1])
            for n in range(Tq+1):
                H[4*n-b]=int(h[n]) if n<=h.degree() else 0
        s.C0={};s.C1={}
        for D in set(H)|set(c01):
            c=Fraction(c01.get(D,0)); h=Fraction(H.get(D,0))
            s.C0[D]=Fraction(4,N+1)*c+Fraction(2*(N-1),N+1)*h
            s.C1[D]=Fraction(2,N+1)*(c-h)
        s.kp2=24//(N+1)
    def L(s,LL):
        N=s.N; X=s.X; terms={}
        for t in range(1,LL+1):
            if LL%t: continue
            l=LL//t
            for k in range(0,X//t+1):
                jm=int(math.isqrt(4*k*l+1))+1
                for j in range(-jm,jm+1):
                    D=4*k*l-j*j
                    if D< -1: continue
                    v=s.C0.get(D,0)-s.C1.get(D,0)+(N*s.C1.get(D,0) if t%N==0 else 0)
                    if v: terms[(k*t,j*t)]=terms.get((k*t,j*t),0)+Fraction(v)/t
        return s.enc(terms)
class ExactS(Exact):
    def __init__(s,N,X,M):
        s.F=FastS(N,X+1,M); s.N=N; s.X=X; s.M=M
        F=s.F
        t=time.time(); s.Ev=F.numerators(); s.t_num=time.time()-t
        PA=series_1d_PA(1,X+1)
        s.PA=F.enc({k:Fraction(v) for k,v in PA.items()})
        gi=ginv(N,X+1,F.kp2)
        s.G=F.enc({(a,0):Fraction(v) for a,v in enumerate(gi) if v})
if __name__=='__main__':
    N,X,M=map(int,sys.argv[1:4])
    E=ExactS(N,X,M)
    for m in range(1,M+1):
        T=E.table(m)
        print(m,{k:str(v) for k,v in T.items() if k[0]<=2 and 4*m*k[0]-k[1]**2<0})
