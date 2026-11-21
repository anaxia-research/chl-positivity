# Weyl vectors of the CHL product and the Levi group Gamma_0(N).
# A root is T(l,kappa,j), matrix [[l, j/2],[j/2, kappa/N]], with E(T)!=0 and det T<0:
#   Z-type: N|l, N|kappa, j odd, 4 det = -1, multiplicity 2
#   P-type (N=5,7): kappa*l = 1 mod N, j odd, 4 det = 4/N-1, multiplicity -48/(N^2-1)
# For A in Gamma_0(N) (lower left divisible by N), Phi(A Om A^T) is the product of the
# chamber W''=A^{-1} C0 A^{-T} with Weyl vector A^T rho0 A.  We check
#   A^T rho0 A = rho0 + sum_{mu sep. C0,W''} c(mu) mu   and   sum c(mu) even.
import sys
from fractions import Fraction as Fr
from math import isqrt, floor, ceil
N=int(sys.argv[1])
def mat(l,k,j): return ((Fr(l),Fr(j,2)),(Fr(j,2),Fr(k,N)))
def tr(T,Y): return T[0][0]*Y[0][0]+2*T[0][1]*Y[0][1]+T[1][1]*Y[1][1]
def mult(T):
    l,k,j=T
    D4=Fr(4*k*l,N)-j*j
    if j%2==0: return 0
    if D4==-1 and l%N==0 and k%N==0: return 2
    if N>4 and D4==Fr(4,N)-1 and (k*l)%N==1: return Fr(-48,N*N-1)
    return 0
def roots_near(xlo,xhi,ymin):
    # all roots whose wall meets {xlo<=Re tau<=xhi, Im tau>=ymin}
    out=[]
    # l=0: walls kappa sigma/N + j v = 0 -> Re tau = kappa/(N j)  (j=+-1, N|kappa)
    for kk in range(floor(xlo)-1,ceil(xhi)+2):
        for j in (1,-1):
            T=(0,N*kk*j,j) if False else (0,N*kk,j)
            if mult(T): out.append(T)
    lmax=int(1/(2*ymin))+1
    for l in range(1,lmax+1):
        jlo=floor(2*l*(xlo-1))-1; jhi=ceil(2*l*(xhi+1))+1
        for j in range(jlo,jhi+1):
            if j%2==0: continue
            for D4 in ([-1] + ([Fr(4,N)-1] if N>4 else [])):
                # 4 k l /N - j^2 = D4  -> k = N (D4 + j^2)/(4l)
                kf=Fr(N,1)*(D4+j*j)/(4*l)
                if kf.denominator!=1: continue
                T=(l,int(kf),j)
                if mult(T): out.append(T)
    return out
def tau(Y):
    det=Y[0][0]*Y[1][1]-Y[0][1]**2
    return (-Y[0][1]/Y[1][1], (float(det)**0.5)/float(Y[1][1]))
def conj(A,Y):  # A Y A^T
    a,b,c,d=A
    M=((a,b),(c,d))
    R=[[sum(M[i][p]*Y[p][q]*M[j][q] for p in range(2) for q in range(2)) for j in range(2)] for i in range(2)]
    return ((R[0][0],R[0][1]),(R[1][0],R[1][1]))
def inv(A):
    a,b,c,d=A; D=a*d-b*c; assert D in (1,-1); return (d*D,-b*D,-c*D,a*D)
def actT(A,T):  # A^T T A as triple
    a,b,c,d=A
    M=mat(*T)
    At=(a,c,b,d)
    R=conj(At,M)
    l=R[0][0]; j=2*R[0][1]; k=N*R[1][1]
    assert l.denominator==1 and j.denominator==1 and k.denominator==1
    return (int(l),int(k),int(j))
rho0=(1,1,1)
Ya=((Fr(5,4),Fr(-1,2)),(Fr(-1,2),Fr(1)))   # tau = 1/2 + i
def check(A,verbose=False):
    Yb=conj(inv(A),Ya)
    ta=tau(Ya); tb=tau(Yb)
    xlo=min(ta[0],tb[0]); xhi=max(ta[0],tb[0]); ymin=min(ta[1],tb[1])
    S=[]
    for T in roots_near(xlo,xhi,ymin):
        sa=tr(mat(*T),Ya); sb=tr(mat(*T),Yb)
        assert sa!=0 and sb!=0
        if sa>0 and sb<0: S.append(T)
        if sa<0 and sb>0: S.append(tuple(-t for t in T))
    S=sorted(set(S))
    tot=[Fr(rho0[i]) for i in range(3)]
    csum=Fr(0)
    for T in S:
        c=mult(T); csum+=c
        for i in range(3): tot[i]+=c*T[i]
    lhs=actT(A,rho0)
    ok=(tuple(tot)==tuple(Fr(x) for x in lhs)) and csum.denominator==1 and int(csum)%2==0
    if verbose: print(A,'walls',len(S),'sum c',csum,'A^T rho A',lhs,'rho_W',tuple(int(x) for x in tot),'OK' if ok else 'FAIL')
    return ok,S
def mul(X,Y):
    a,b,c,d=X; e,f,g,h=Y
    return (a*e+b*g,a*f+b*h,c*e+d*g,c*f+d*h)
# coset reps of Gamma_0(N)\SL2(Z): Gamma_0(N) g ; line of (0,1) g = (c,d)
S_=(0,-1,1,0); T_=(1,1,0,1)
reps={}
def line(g):
    c,d=g[2]%N,g[3]%N
    if c==0 or N==1: return ('inf',)
    return ('pt',(d*pow(c,-1,N))%N)
reps[line((1,0,0,1))]=(1,0,0,1)
for j in range(N):
    g=mul(S_,(1,j,0,1)) ; reps.setdefault(line(g),g)
assert len(reps)==(N+1 if N>1 else 1)
gens=[]
for L,g in reps.items():
    for x in (S_,T_):
        gx=mul(g,x); h=mul(gx,inv(reps[line(gx)]))
        assert h[2]%N==0
        if h not in [(1,0,0,1),(-1,0,0,-1)]: gens.append(h)
print('N',N,'Schreier generators',len(gens))
allok=True
for h in gens:
    ok,S=check(h,verbose=True); allok&=ok
# also elliptic elements used in the paper and U_N
extra=[(1,0,0,-1)]+([(1,0,N,1),(2,1,-5,-2) if N==5 else (3,1,-7,-2) if N==7 else (1,-1,3,-2)] if N>1 else [(0,1,1,0)])
for h in extra:
    ok,S=check(h,verbose=True); allok&=ok
print('ALL OK' if allok else 'SOME FAIL')
