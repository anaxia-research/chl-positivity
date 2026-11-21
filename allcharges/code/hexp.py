# Taylor expansion in h=1/m of R(A,B;h)=rho_AB e^{c(A+B)} for half-integer nu (e^{-2z} parts dropped).
# Truncated power series in h with coefficients in QQ[A,B,c,ic], ic=1/c.
import sympy as sp, pickle, sys
from math import factorial
A,B,c,ic=sp.symbols('A B c ic'); G=(A,B,c,ic)
def P(e): return sp.Poly(e,*G,domain='QQ')
Z=P(0)
def mul(f,g,J): return [sum((f[i]*g[j-i] for i in range(j+1)),Z) for j in range(J+1)]
def slog(f,J):  # f0=1
    g=[Z]*(J+1)
    for j in range(1,J+1):
        g[j]=(f[j]*j - sum((g[i]*f[j-i]*i for i in range(1,j)),Z))*sp.Rational(1,j)
    return g
def sexp(g,J):  # g0=0
    e=[P(1)]+[Z]*J
    for j in range(1,J+1):
        e[j]=sum((g[i]*e[j-i]*i for i in range(1,j+1)),Z)*sp.Rational(1,j)
    return e
def scal(f,a): return [x*a for x in f]
def add(f,g): return [x+y for x,y in zip(f,g)]
def pad(lst,J): return [P(v) for v in lst]+[Z]*(J+1-len(lst))
def expansion(nu2,J):
    n=(nu2-1)//2; nu=sp.Rational(nu2,2); S=A+B; y2=(A-B)**2
    K=J+2
    T=pad([1,4-2*S,y2],K); T0=pad([1,4],K)
    lT=slog(T,K); lT0=slog(T0,K)
    l1=add(slog(pad([1,2-S],K),K),scal(slog(pad([1,2],K),K),-1))
    l2=scal(add(lT,scal(lT0,-1)),nu/2-sp.Rational(1,4))
    sT=sexp(scal(lT,sp.Rational(1,2)),K); sT0=sexp(scal(lT0,sp.Rational(1,2)),K)
    dif=add(sT,scal(sT0,-1))      # coefficient 0 is 0
    l3=[ (dif[j+1]*P(c) if j+1<=K else Z) for j in range(K)]+[Z]
    l3[0]=l3[0]+P(c*S)
    def logp(isq):   # isq = series of T^{-1/2}; w = h*ic*isq
        w=[Z]+[x*P(ic) for x in isq[:K]]
        pw=pad([1],K); wk=pad([1],K)
        for kk in range(1,n+1):
            wk=mul(wk,w,K)
            pw=add(pw,scal(wk,sp.Rational(factorial(n+kk),factorial(kk)*factorial(n-kk))*sp.Rational(-1,2)**kk))
        return slog(pw,K)
    l4=add(logp(sexp(scal(lT,sp.Rational(-1,2)),K)),scal(logp(sexp(scal(lT0,sp.Rational(-1,2)),K)),-1))
    L=add(add(l1,l2),add(l3,l4))
    assert L[0].is_zero, L[0]
    R=sexp(L[:J+1],J)
    return [r.as_expr() for r in R]
if __name__=='__main__':
    nu2=int(sys.argv[1]); J=int(sys.argv[2])
    co=expansion(nu2,J)
    print(0,co[0]); print(1,sp.expand(co[1].subs(ic,1/c)))
    pickle.dump((co,(A,B,c,ic)),open(f'hexp_{nu2}_{J}.pkl','wb'))
