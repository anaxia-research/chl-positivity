# Kronecker-substitution engine for numerators E_{m+1} = g A phi_m
import sys, math, time; sys.path.insert(0,'/tmp/claude-0/-home-user/e761e78e-04a9-5d22-97e1-0a0fe4aaf27f/scratchpad/code')
from chl import CHL, S, fmpq_poly, fmpq
from fastinit import FastCHL
from flint import fmpz_poly
from fractions import Fraction
class Fast:
    def __init__(s,N,X,M):
        s.C=FastCHL(N,X,M); s.N,s.X,s.M=N,X,M
        s.YB=int(math.isqrt(4*(M+2)*X//N+1))+M+4
        s.W=2*s.YB+1; s.K=s.YB
    def enc(s,terms):  # dict (a,b)->Fraction
        den=1
        for v in terms.values(): den=den*v.denominator//math.gcd(den,v.denominator)
        L=[0]*((s.X+1)*s.W)
        for (a,b),v in terms.items():
            if a<=s.X:
                assert abs(b)<=s.YB
                L[a*s.W+b+s.K]+=int(v*den)
        return (fmpz_poly(L),den)
    def mul(s,A,B):
        p=A[0].mul_low(B[0],(s.X+1)*s.W+s.K) if hasattr(A[0],'mul_low') else A[0]*B[0]
        # shift down by K
        c=p.coeffs()[s.K:(s.X+1)*s.W+s.K]
        return (fmpz_poly(c),A[1]*B[1])
    def L(s,LL):
        C=s.C; X=s.X; N=s.N; terms={}
        for t in range(1,LL+1):
            if LL%t: continue
            l=LL//t
            for kap in range(0,X//t+1):
                jm=int(math.isqrt(4*kap*l//N+1))+1
                for j in range(-jm,jm+1):
                    v=C.E(kap,l,j)
                    if v: terms[(kap*t,j*t)]=terms.get((kap*t,j*t),0)+Fraction(v)/t
        return s.enc(terms)
    def numerators(s):
        Ls=[None]+[s.L(l) for l in range(1,s.M+2)]
        one=(fmpz_poly([0]*s.K+[1]),1)
        Ev=[one]
        for n in range(1,s.M+2):
            acc=None; den=None
            parts=[]
            for k in range(1,n+1):
                p,d=s.mul(Ls[k],Ev[n-k]); parts.append((p*k,d))
            D=1
            for p,d in parts: D=D*d//math.gcd(D,d)
            tot=fmpz_poly(0)
            for p,d in parts: tot+=p*(D//d)
            D*=n
            g=tot.content() if hasattr(tot,'content') else 1
            g=math.gcd(int(g),D) if g else D
            Ev.append((fmpz_poly([c//g for c in tot.coeffs()]),D//g))
        return Ev
    def to_S(s,obj):
        p,d=obj; c=p.coeffs(); rows=[]
        for a in range(s.X+1):
            seg=c[a*s.W:(a+1)*s.W]
            rows.append(fmpq_poly([fmpq(int(v),d) for v in seg]) if any(seg) else fmpq_poly(0))
        return S(rows,s.K)
if __name__=='__main__':
    N,X,M=map(int,sys.argv[1:4])
    t=time.time(); F=Fast(N,X,M); print('init',time.time()-t)
    t=time.time(); Ev=F.numerators(); print('bell',time.time()-t, 'W',F.W)
    # compare with slow engine for small case
    if X<=40:
        from chl import coeff
        slow=F.C.numerators()
        for m in range(1,M+1):
            sm=slow[m]; fm=F.to_S(Ev[m+1])
            ok=all(coeff(sm,a,b)==coeff(fm,a,b) for a in range(X+1) for b in range(-F.YB,F.YB+1))
            print('m',m,'match',ok)
