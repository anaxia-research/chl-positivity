import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from chl import *
from fractions import Fraction as Fr
import sympy as sp
N=int(sys.argv[1]); X=int(sys.argv[2])
C=CHL(N,X,1); Nm=C.numerators(); E2=Nm[1]
def jac_series(tab,scale): return series_from_dict({(scale*a,b):v for (a,b),v in tab.items() if scale*a<=X},X)
A=jac_series(C.m21,N); B=jac_series(C.m01,N)
def mod1(coeffs, scale=1): return series_from_dict({(scale*i,0):c for i,c in enumerate(coeffs) if scale*i<=X and c},X)
def E4c(T): return [Fr(1)]+[Fr(240*sigma(n,3)) for n in range(1,T+1)]
def etaq(exps,T):  # prod_d eta(d tau)^{e_d}, returns (leading power*24, coeff list) ; require integer leading
    lead=sum(d*e for d,e in exps.items()); assert lead%24==0; lead//=24
    f=[Fr(0)]*(T+1); f[0]=Fr(1)
    for d,e in exps.items():
        for n in range(1,T//d+1):
            # multiply by (1-q^{dn})^e
            c=1;k=0;terms=[]
            while k*d*n<=T: terms.append((k*d*n,c)); k+=1; c=c*(e-k+1)//k*(-1)
            g=[Fr(0)]*(T+1)
            for i,v in enumerate(f):
                if v:
                    for s,cc in terms:
                        if i+s<=T: g[i+s]+=v*cc
            f=g
    return [Fr(0)]*lead+f[:T+1-lead]
EN=EN_coeffs(N,X)
M2=[mod1(EN)]
M4=[mod1(E4c(X)), mod1(E4c(X//N+1),N), mod1([e*e for e in EN][:1]+[sum(EN[i]*EN[n-i] for i in range(n+1)) for n in range(1,X+1)])]
if N==5: M4.append(mod1(etaq({1:4,5:4},X)))
if N==7: M4.append(mod1(etaq({1:3,7:3},X))[:1] if False else mod1(etaq({1:3,7:3},X)))
basis=[('B2',mul(B,B))]+[('M2AB%d'%i,mul(mul(f,A),B)) for i,f in enumerate(M2)]+[('M4AA%d'%i,mul(mul(f,A),A)) for i,f in enumerate(M4)]
keys=[]
for a in range(min(X,60)+1):
  for b in range(-8,9): keys.append((a,b))
Mx=sp.Matrix([[sp.Rational(str(coeff(s,a,b))) for nm,s in basis] for (a,b) in keys])
v=sp.Matrix([sp.Rational(str(coeff(E2,a,b))) for (a,b) in keys])
try:
    sol,params=Mx.gauss_jordan_solve(v); print([nm for nm,s in basis]); print(sol.T, params)
except Exception as e: print('no solution',e)
