# Z_N CHL dyon partition function: Fourier-Jacobi coefficients and finite parts.
# Series in x = e(sigma/N) truncated at order X; coefficients Laurent polys in y (offset k).
from flint import fmpq_poly, fmpq, fmpz_poly
from fractions import Fraction
import sys, math

class S:
    """sum_a rows[a](y) * y^{-k} x^{a+e}  (rows fmpq_poly), a=0..X"""
    __slots__=('rows','k','e')
    def __init__(s, rows, k=0, e=0): s.rows=rows; s.k=k; s.e=e
    @property
    def X(s): return len(s.rows)-1

def zero(X,k=0,e=0): return S([fmpq_poly(0) for _ in range(X+1)],k,e)
def shift_y(p, d):  # multiply poly by y^d, d>=0
    return p * fmpq_poly([0]*d+[1]) if d else p
def align(a,b):
    k=max(a.k,b.k)
    ra=[shift_y(p,k-a.k) for p in a.rows]; rb=[shift_y(p,k-b.k) for p in b.rows]
    return ra,rb,k
def add(a,b,cb=1):
    assert a.e==b.e
    ra,rb,k=align(a,b); X=min(len(ra),len(rb))-1
    return S([ra[i]+rb[i]*cb for i in range(X+1)],k,a.e)
def mul(a,b,X=None):
    if X is None: X=min(a.X,b.X)
    r=[fmpq_poly(0) for _ in range(X+1)]
    nzb=[j for j in range(min(b.X,X)+1) if b.rows[j]!=0]
    for i in range(min(a.X,X)+1):
        pa=a.rows[i]
        if pa==0: continue
        for j in nzb:
            if i+j>X: break
            r[i+j]+=pa*b.rows[j]
    return S(r,a.k+b.k,a.e+b.e)
def scal(a,c): return S([p*c for p in a.rows],a.k,a.e)
def coeff(a,n,l):
    i=n-a.e
    if i<0 or i>a.X: return fmpq(0)
    j=l+a.k
    p=a.rows[i]
    if j<0 or j>p.degree(): return fmpq(0)
    return p[j]

# ---------- q-series helpers (one variable, list of Fractions) ----------

def sigma(n,k=1): return sum(d**k for d in range(1,n+1) if n%d==0)

def jac_tables(T):
    """phi_{-2,1}, phi_{0,1} coefficient dicts {(n,j):int}, integer q-powers n<=T"""
    # P = prod (1-q^n y)^2(1-q^n/y)^2/(1-q^n)^4 as dict
    from collections import defaultdict
    P={(0,0):1}
    def mulfac(P, n, yexp, power):  # multiply by (1 - q^n y^yexp)^power (power may be negative)
        # binomial series
        terms=[]
        c=1; k=0
        while k*n<=T:
            terms.append((k*n,k*yexp,c))
            k+=1; c=c*(power-k+1)//k*(-1) if True else 0
        R=defaultdict(int)
        for (a,b),v in P.items():
            for (da,db,c) in terms:
                if a+da<=T: R[(a+da,b+db)]+=v*c
        return dict(R)
    for n in range(1,T+1):
        P=mulfac(P,n,1,2); P=mulfac(P,n,-1,2); P=mulfac(P,n,0,-4)
    m21={}
    for (a,b),v in P.items():
        for db,c in ((1,1),(0,-2),(-1,1)):
            m21[(a,b+db)]=m21.get((a,b+db),0)+v*c
    # wp/(2pi i)^2 regular part: 1/12 + sum_{n} sum_{d|n} d (y^d -2 + y^-d) q^n
    W={(0,0):Fraction(1,12)}
    for n in range(1,T+1):
        for d in range(1,n+1):
            if n%d==0:
                for db,c in ((d,1),(0,-2),(-d,1)): W[(n,db)]=W.get((n,db),0)+d*c
    m01={}
    for (a,b),v in m21.items():
        for (c,e),w in W.items():
            if a+c<=T: m01[(a+c,b+e)]=m01.get((a+c,b+e),0)+12*v*w
    for (a,b),v in P.items(): m01[(a,b)]=m01.get((a,b),0)+12*v
    m01={k:int(v) for k,v in m01.items() if v!=0}
    m21={k:v for k,v in m21.items() if v!=0}
    return m21,m01

def EN_coeffs(N,T):
    return [Fraction(1)]+[Fraction(24,N-1)*(sigma(n)-(N*sigma(n//N) if n%N==0 else 0)) for n in range(1,T+1)]

class CHL:
    def __init__(s,N,X,M):
        """N prime, X: x-order truncation, M: max m"""
        s.N=N; s.X=X; s.M=M
        kp2=24//(N+1); s.kp2=kp2; s.a=Fraction(2,N*(N+1))
        Tq=(X*(M+1))//N+3   # integer q-power range needed for genera
        s.m21,s.m01=jac_tables(Tq)
        e=EN_coeffs(N,N*Tq+N)
        # tables in t = N*n units
        Phi01={}; H={}; G={}
        for (n,j),v in s.m01.items(): Phi01[(N*n,j)]=v
        for (n,j),v in s.m21.items():
            for u in range(0,Tq+1-n):
                H[(N*(n+u),j)]=H.get((N*(n+u),j),0)+v*e[u]
            for u in range(0,N*(Tq-n)+1):
                G[(N*n+u,j)]=G.get((N*n+u,j),0)+v*e[u]
        s.Phi01,s.H,s.G=Phi01,H,G
    def E(s,kap,l,j):
        N=s.N; b=abs(j)%2; t=kap*l-N*(j*j-b)//4
        if t<0: return Fraction(0)
        r=kap%N; a=s.a
        p01=s.Phi01.get((t,b),0)
        if r==0:
            T00=Fraction(2,N)*p01
            T01=a*p01+Fraction(2,N+1)*s.H.get((t,b),0)
            return T00+T01*((N if l%N==0 else 0)-1)
        return a*N*((p01 if l%N==0 else 0)-s.G.get((t,b),0))
    def L(s,LL):
        X=s.X; N=s.N
        terms={}
        for t in range(1,LL+1):
            if LL%t: continue
            l=LL//t
            for kap in range(0,X//t+1):
                jm=int(math.isqrt(4*kap*l//N+1))+1
                for j in range(-jm,jm+1):
                    if kap==0 and False: pass
                    v=s.E(kap,l,j)
                    if v: terms[(kap*t,j*t)]=terms.get((kap*t,j*t),0)+v/t
        K=max([-b for (a,b) in terms]+[0])
        rows=[[] for _ in range(X+1)]
        D=[dict() for _ in range(X+1)]
        for (a,b),v in terms.items():
            if a<=X: D[a][b+K]=D[a].get(b+K,0)+v
        R=[]
        for a in range(X+1):
            if D[a]:
                deg=max(D[a]); lst=[fmpq(0)]*(deg+1)
                for b,v in D[a].items(): lst[b]=fmpq(v.numerator,v.denominator)
                R.append(fmpq_poly(lst))
            else: R.append(fmpq_poly(0))
        return S(R,K)
    def numerators(s):
        """Nm[m] = g(sigma/N) A phi_m for m=-1..M  (Bell recursion)"""
        Ls=[None]+[s.L(l) for l in range(1,s.M+2)]
        one=S([fmpq_poly(1)]+[fmpq_poly(0)]*s.X,0)
        Ev=[one]
        for n in range(1,s.M+2):
            acc=None
            for k in range(1,n+1):
                t=scal(mul(Ls[k],Ev[n-k]),k)
                acc=t if acc is None else add(acc,t)
            Ev.append(scal(acc,fmpq(1,n)))
        s.Ls=Ls
        return {m:Ev[m+1] for m in range(-1,s.M+1)}

def series_from_dict(D,X):
    K=max([-b for (a,b) in D]+[0])
    rows=[]
    for a in range(X+1):
        items={b+K:v for (aa,b),v in D.items() if aa==a}
        if items:
            lst=[fmpq(0)]*(max(items)+1)
            for b,v in items.items():
                v=Fraction(v); lst[b]=fmpq(v.numerator,v.denominator)
            rows.append(fmpq_poly(lst))
        else: rows.append(fmpq_poly(0))
    return S(rows,K)

def finite_parts(C):
    N,X,M=C.N,C.X,C.M
    Nm=C.numerators()
    # 1/P_A in p=x^N
    from collections import defaultdict
    Tp=X//N
    P={(0,0):Fraction(1)}
    def mulfac(P,n,ye,power):
        R=defaultdict(Fraction)
        terms=[];c=1;k=0
        while k*n<=Tp:
            terms.append((k*n,k*ye,c)); k+=1; c=c*(power-k+1)//k*(-1)
        for (a,b),v in P.items():
            for da,db,cc in terms:
                if a+da<=Tp: R[(a+da,b+db)]+=v*cc
        return dict(R)
    for n in range(1,Tp+1):
        P=mulfac(P,n,1,-2);P=mulfac(P,n,-1,-2);P=mulfac(P,n,0,4)
    invPA=series_from_dict({(N*a,b):v for (a,b),v in P.items()},X)
    # 1/g(sigma/N) = x^{-1} prod (1-x^n)^{-kp2}(1-x^{Nn})^{-kp2}
    gi=[Fraction(0)]*(X+1); gi[0]=Fraction(1)
    def mul1(f,n,power):
        out=[Fraction(0)]*(X+1); c=1;k=0;terms=[]
        while k*n<=X: terms.append((k*n,c)); k+=1; c=c*(power-k+1)//k*(-1)
        for i,v in enumerate(f):
            if v:
                for d,cc in terms:
                    if i+d<=X: out[i+d]+=v*cc
        return out
    for n in range(1,X+1):
        gi=mul1(gi,n,-C.kp2)
        if N*n<=X: gi=mul1(gi,N*n,-C.kp2)
    C.ginv=gi   # coefficients of x^{j} * x^{-1}
    # fhat(m): coeff q^m of 1/(eta^{k+2}(tau) eta^{k+2}(N tau)) = q^{-1} * same product in q
    fh={m:gi[m+1] for m in range(-1,min(M,X-1)+1)}
    C.fhat=fh
    out={}
    for m in range(1,M+1):
        R=mul(Nm[m],invPA)
        R=add(R,S([fmpq_poly(fmpq(fh[m].numerator,fh[m].denominator))]+[fmpq_poly(0)]*X,0),-1)
        # divide each row by y^{-1}(y-1)^2 : multiply by y, divide by (y-1)^2
        d2=fmpq_poly([1,-2,1])
        rows=[]
        for p in R.rows:
            q,r=divmod(p,d2)
            assert r==0, "non-divisible row"
            rows.append(q)
        Q=S(rows,R.k-1)
        # minus fhat * sum_{s!=0}
        D={}
        for u in range(1,X+1):
            for k in range(0,X+1):
                a=N*(m*u*u+u+u*k)
                if a>X: break
                D[(a,2*m*u+1+k)]=D.get((a,2*m*u+1+k),0)+(k+1)
                D[(a,-2*m*u-1-k)]=D.get((a,-2*m*u-1-k),0)+(k+1)
        if D:
            As=series_from_dict(D,X)
            Q=add(Q,scal(As,fmpq(fh[m].numerator,fh[m].denominator)),-1)
        G=S([fmpq_poly(fmpq(v.numerator,v.denominator)) for v in gi],0,-1)
        out[m]=mul(Q,G)
    return out

def dstar(F,N,m,n,l):
    return (1 if l%2 else -1)*coeff(F[m],n,l)
