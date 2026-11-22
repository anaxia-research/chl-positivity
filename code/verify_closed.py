# Uniform first-row closed form, N=2,3,5,7:
#   g(sigma/N) phi^F_1 = W_N(sigma/N) phi_{-2,1} - b_N Theta_N - 2 d(1) H,
# Theta_N = theta series of a maximal order O_N (ramified at N, infinity) with characteristic w, nrd(w)=1/N.
import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from chl import *
from fractions import Fraction as Fr
from row1_lib import hurwitz
from lattice import theta_rows
import numpy as np, sympy as sp
exec(open('quat.py').read().split("for N in (2,3,5,7)")[0])
N=int(sys.argv[1]); X=int(sys.argv[2])
def E4c(T): return [Fr(1)]+[Fr(240*sigma(n,3)) for n in range(1,T+1)]
def etaq(exps,T):
    lead=sum(d*e for d,e in exps.items())//24
    f=[Fr(0)]*(T+1); f[0]=Fr(1)
    for d,e in exps.items():
        for n in range(1,T//d+1):
            g=[Fr(0)]*(T+1)
            c=1;k=0;terms=[]
            while k*d*n<=T: terms.append((k*d*n,c)); k+=1; c=c*(e-k+1)//k*(-1)
            for i,v in enumerate(f):
                if v:
                    for s_,cc in terms:
                        if i+s_<=T: g[i+s_]+=v*cc
            f=g
    return [Fr(0)]*lead+f[:T+1-lead]
e4=E4c(X); e4N=[e4[n//N] if n%N==0 else Fr(0) for n in range(X+1)]
if N==2: W=[(a-b)/5 for a,b in zip(e4,e4N)]
if N==3: W=[Fr(3,80)*(a-b) for a,b in zip(e4,e4N)]
if N==5: W=[15*c for c in etaq({5:10,1:-2},X)]
if N==7:
    th=[0]*(X+1); r=int((4*X)**0.5)+2
    for a in range(-r,r+1):
        for b in range(-r,r+1):
            v=a*a+a*b+2*b*b
            if v<=X: th[v]+=1
    e=etaq({7:7,1:-1},X); W=[3*sum(th[i]*e[n-i] for i in range(n+1)) for n in range(X+1)]
b={2:Fr(20,3),3:Fr(9,2),5:Fr(7,3),7:Fr(3,2)}[N]
d1={2:52,3:27,5:14,7:9}[N]
ws={2:(1,1,0,1),3:(-1,-1,0,0),5:(0,-1,-1,-1),7:(0,-1,0,0)}[N]
C=CHL(N,X,1); F=finite_parts(C)
assert C.fhat[1]==d1
A=series_from_dict({(N*a,bb):v for (a,bb),v in C.m21.items() if N*a<=X},X)
WA=mul(series_from_dict({(i,0):c for i,c in enumerate(W) if c},X),A)
th=theta_rows(np.array(order(N).tolist(),dtype=int),ws,X)
Hd={}
for n in range(0,X//N+1):
    for l in range(-2*int((n+1)**0.5)-2,2*int((n+1)**0.5)+3):
        h=hurwitz(4*n-l*l)
        if h: Hd[(N*n,l)]=h
K=add(WA,series_from_dict({k:v for k,v in th.items() if k[0]<=X},X),fmpq(-b.numerator,b.denominator))
K=add(K,series_from_dict(Hd,X),fmpq(-2*d1))
gi=C.ginv; inv=[Fr(0)]*(X+1); inv[0]=1/gi[0]
for i in range(1,X+1): inv[i]=-sum(gi[j]*inv[i-j] for j in range(1,i+1))/gi[0]
Kp=mul(F[1],S([fmpq_poly(fmpq(v.numerator,v.denominator)) for v in inv],0,1))
d=add(K,Kp,-1)
bad=[a for a in range(X-1) if d.rows[a]!=0]
if N==2:  # D4 description: r(n,l) = r_3(2n-l^2)
    r3=lambda D: sum(1 for a in range(-6,7) for bb in range(-6,7) for c in range(-6,7) if a*a+bb*bb+c*c==D)
    assert all(th.get((n,l),0)==r3(2*n-l*l) for n in range(0,13) for l in range(-5,6))
print('N=%d: closed form for g*phi^F_1 checked through x^%d; mismatches:'%(N,X-2),bad[:5])
