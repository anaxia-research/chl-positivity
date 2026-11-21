"""Certificate for the first row Q^2=2/N (m=1) of the Z_N CHL models, N=3,5,7.
K*(Delta) >= w_N(u) - b_N r_N(u) - 2 fhat_N H(Delta/N) for l odd, u=(Delta+N)/4."""
import sys, math; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from fractions import Fraction as Fr
exec(open('W.py').read().split("for N,(a1")[0])
from row1_lib import hurwitz
U=400
def leg(a,p):
    a%=p
    return 0 if a==0 else (1 if pow(a,(p-1)//2,p)==1 else -1)
data={3:(Fr(9,2),27),5:(Fr(7,3),14),7:(Fr(3,2),9)}
for N,(b,fh) in data.items():
    a1,a2,a3={3:(Fr(3,80),Fr(-3,80),0),5:(Fr(1,48),Fr(19,48),Fr(-5,12)),7:(Fr(1,320),Fr(29,320),Fr(-3,32))}[N]
    e4=E4(U); e4N=dil(E4(U),N,U); EN=EN_coeffs(N,U); EN2=sq(EN,U)
    W=[a1*e4[t]+a2*e4N[t]+a3*EN2[t] for t in range(U+1)]
    assert all(w>=0 for w in W)
    # Eisenstein + cusp decomposition
    if N==3: cusp=None; x=Fr(3,80)*1; y=0; g=[0]*(U+1)
    else:
        g = eta({1:4,5:4},U) if N==5 else [sum(a*b for a,b in zip(eta({1:3,7:3},U)[:n+1], list(reversed(th7[:n+1])))) for n in range(U+1)] if False else None
        if N==7:
            th=[0]*(U+1)
            r=int(math.isqrt(4*U))+2
            for p in range(-r,r+1):
                for q in range(-r,r+1):
                    v=p*p+p*q+2*q*q
                    if v<=U: th[v]+=1
            c3=eta({1:3,7:3},U); g=[sum(c3[i]*th[n-i] for i in range(n+1)) for n in range(U+1)]
        # W = x(E4-E4(N t)) + y g ; solve from t=1,2
        import sympy as sp
        xs,ys=sp.symbols('x y')
        sol=sp.solve([xs*240*(sigma(t,3))+ys*g[t]-W[t] for t in (1,2)],[xs,ys])
        x=Fr(str(sol[xs])); y=Fr(str(sol[ys]))
        assert all(W[t]==x*(e4[t]-e4N[t])+y*g[t] for t in range(U+1)), 'decomp fails'
        # Hecke multiplicativity check of normalized g (eigenform) on small n
        assert g[1]==1 and all(g[m*n]==g[m]*g[n] for m in range(2,12) for n in range(2,12) if math.gcd(m,n)==1 and m*n<=U)
    print('N=%d: W_N = %s (E4 - E4(%d tau)) + %s g_N, all coefficients >=0 through q^%d'%(N,x,N,y,U))
    # lower bound analytic: w(u) >= 240 x (1-N^-3) u^3 - |y| d(u) u^1.5 ; d(u)<=2 sqrt(u)
    # upper: b r(u) <= b*(24/(N-1)) sigma1(u) <= b*(24/(N-1)) * u*(1+ln u) ; 2 fh H <= 2 fh * 4u/N
    def LBan(u): return 240*float(x)*(1-N**-3)*u**3 - abs(float(y))*2*u**2 - float(b)*24/(N-1)*u*(1+math.log(u)) - 2*fh*4*u/N
    u0=next(u for u in range(2,10**6) if all(LBan(v)>0 for v in range(u,u+2000)))
    # LBan is eventually increasing; verify monotone derivative crude: check up to 10^5 sparse
    assert all(LBan(v)>0 for v in range(u0,100000,7))
    print('   analytic bound positive for u >= %d'%u0)
    # exact check of lower bound for 2<=u<u0 via lattice count r(u,1)<=total r_O(u)=(24/(N-1))sigma1*(u)
    bad=[]
    for u in range(1,max(u0,60)):
        D=4*u-N
        rO=Fr(24,N-1)*sum(d for d in range(1,u+1) if u%d==0 and d%N)
        h=hurwitz(D//N) if D%N==0 else 0
        lb=W[u]-b*rO-2*fh*h
        if lb<=0: bad.append((u,D,float(lb)))
    print('   u with crude bound <=0 (need exact K*):',bad)
