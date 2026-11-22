# Check: psi_L from the product exponents equals the twisted Hecke sum
#  (1/L) sum_{ad=L} sum_{b mod d} sum_{r,s} e(-(rb+sd)/N) F^{(r,s)}((a s+b)/d, a v)
import sys, cmath
from fractions import Fraction as Fr
N=int(sys.argv[1]); PREC=int(sys.argv[2]) if len(sys.argv)>2 else 12   # q-precision (integer powers)
def e(x): return cmath.exp(2j*cmath.pi*x)
# series: dict (n:Fraction, j:int) -> complex
def mul(A,B,cap):
    C={}
    for (n1,j1),a in A.items():
        for (n2,j2),b in B.items():
            n=n1+n2
            if n<=cap: C[(n,j1+j2)]=C.get((n,j1+j2),0)+a*b
    return C
def add(*args):
    C={}
    for c,A in args:
        for k,v in A.items(): C[k]=C.get(k,0)+c*v
    return C
cap=Fr(PREC)
# theta-type products for phi_{-2,1} = (y-2+1/y) prod (1-q^n y)^2(1-q^n/y)^2/(1-q^n)^4, phi_{0,1}
def prod_series():
    P={(Fr(0),1):1,(Fr(0),0):-2,(Fr(0),-1):1}
    for n in range(1,PREC+1):
        for f in [{(Fr(0),0):1,(Fr(n),1):-1},{(Fr(0),0):1,(Fr(n),-1):-1}]:
            P=mul(P,f,cap);P=mul(P,f,cap)
        inv={(Fr(n*k),0):(k+1)*(k+2)*(k+3)//6 for k in range(0,PREC//n+1)}  # (1-q^n)^-4
        P=mul(P,inv,cap)
    return P
phim2=prod_series()
# phi_{0,1} = 12 * phi_{-2,1} * wp/(2pi i)^2 ... use phi01 = (y+10+1/y)+q(10y^-2-64y^-1+108-64y+10y^2)+... compute via phi_{0,1}=4 sum_i theta_i^2/theta_i(0)^2 ; easier: phi01 = phim2 * (12 wp) with wp series
# Use known relation: phi_{0,1}(tau,z) = (3/pi^2) wp(tau,z) phi_{-2,1}; (3/pi^2)wp = 12*[ y/(1-y)^2 + 1/12 + sum_{n} sum_{d|n} d (y^d+y^-d-2) q^n ] / (-1)?? 
# (1/(2pi i)^2) wp = 1/12 + y/(1-y)^2 + sum_{n>=1} sum_{d|n} d (y^d - 2 + y^{-d}) q^n ; phi01 = 12 * phim2 * that.
# y/(1-y)^2 * (y-2+1/y) = -1, so handle: phim2*y/(1-y)^2 = -prod part
Q={(Fr(0),0):1}
for n in range(1,PREC+1):
    for f in [{(Fr(0),0):1,(Fr(n),1):-1},{(Fr(0),0):1,(Fr(n),-1):-1}]:
        Q=mul(Q,f,cap);Q=mul(Q,f,cap)
    inv={(Fr(n*k),0):(k+1)*(k+2)*(k+3)//6 for k in range(0,PREC//n+1)}
    Q=mul(Q,inv,cap)
W={(Fr(0),0):Fr(1,12)}
for n in range(1,PREC+1):
    for d in range(1,n+1):
        if n%d==0:
            for jj,cc in ((d,1),(-d,1),(0,-2)): W[(Fr(n),jj)]=W.get((Fr(n),jj),0)+d*cc
phi01=add((12,mul(phim2,W,cap)),(-12,Q))
def sigN(n): return sum(d for d in range(1,n+1) if n%d==0 and d%N)
def EN(scale_t=None):
    pass
# E_N(tau)= 1+24/(N-1) sum sigN(n) q^n ; E_N((tau+t)/N): q^{n/N} e(nt/N)
def EN_series(t=None):
    if t is None: return {(Fr(n),0):(1 if n==0 else Fr(24,N-1)*sigN(n)) for n in range(0,PREC+1)}
    return {(Fr(n,N),0):(1 if n==0 else Fr(24,N-1)*sigN(n))*e(Fr(n*t,N)) for n in range(0,N*PREC+1)}
F={}
aN=Fr(2,N*(N+1))
F[(0,0)]=add((Fr(2,N),phi01))
for s in range(1,N): F[(0,s)]=add((aN,phi01),(Fr(2,N+1),mul(phim2,EN_series(),cap)))
for r in range(1,N):
    for t in range(N):
        F[(r,(r*t)%N)]=add((aN,phi01),(-aN,mul(phim2,EN_series(t),cap)))
def coef(r,s,n,j): return F[(r%N,s%N)].get((n,j),0)
def E(k,l,j):
    return sum(e(Fr(-s*l,N))*coef(k,s,Fr(k*l,N),j) for s in range(N))
for L in range(1,4):
    lhs={}
    for t in range(1,L+1):
        if L%t: continue
        l=L//t
        for k in range(0,N*PREC//l+1 if l else 0):
            if Fr(k*l,N)>cap: break
            for j in range(-12,13):
                v=E(k,l,j)
                if abs(v)>1e-9:
                    key=(k*t,j*t); lhs[key]=lhs.get(key,0)+v/t
    rhs={}
    for a in range(1,L+1):
        if L%a: continue
        d=L//a
        for b in range(d):
            for r in range(N):
                for s in range(N):
                    w=e(Fr(-(r*b+s*d),N))/L
                    for (n,j),c in F[(r,s)].items():
                        # q^n at tau=(a sig+b)/d -> e(n b/d) x^{N n a/d}
                        ex=N*n*a/d
                        key=(ex,j*a)
                        rhs[key]=rhs.get(key,0)+w*c*e(n*b/d)
    # compare on integral x-exponents below cutoff
    bad=0;cnt=0
    keys=set(k for k in lhs)|set(k for k,v in rhs.items() if abs(v)>1e-8)
    for key in keys:
        ex,j=key
        if Fr(ex)>N*(PREC-1)*1 or abs(j)>8: continue
        if Fr(ex)>N*PREC//L-N: continue
        lv=lhs.get((int(ex),j),0) if Fr(ex).denominator==1 else 0
        rv=rhs.get(key,0)
        cnt+=1
        if abs(lv-rv)>1e-6: bad+=1; print('mismatch',L,key,lv,rv) if bad<5 else None
    print('N',N,'L',L,'compared',cnt,'mismatches',bad)
