# Numerator identity and theta identity on the first row, N=2,3,5,7:
#  psi_2 + psi_1^2/2 = alpha pz^2 - beta E_N(s/N) pz pmi + gamma pmi^2,
#  E_N(s/N) pz + delta pmi = 12 Theta_N.
import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from chl import *
from fractions import Fraction as Fr
from lattice import theta_rows
import numpy as np, sympy as sp
exec(open('quat.py').read().split("for N in (2,3,5,7)")[0])
N=int(sys.argv[1]); X=int(sys.argv[2])
C=CHL(N,X,1); E2=C.numerators()[1]
A=series_from_dict({(N*a,b):v for (a,b),v in C.m21.items() if N*a<=X},X)
B=series_from_dict({(N*a,b):v for (a,b),v in C.m01.items() if N*a<=X},X)
def mod(c): return series_from_dict({(i,0):v for i,v in enumerate(c[:X+1]) if v},X)
E4=[Fr(1)]+[Fr(240*sigma(n,3)) for n in range(1,X+1)]
E4s=mod(E4)                                            # E_4(sigma/N)
E4N=mod([E4[n//N] if n%N==0 else 0 for n in range(X+1)])  # E_4(sigma)
EN=EN_coeffs(N,X); ENs=mod(EN)
ENsq=mod([sum(EN[i]*EN[n-i] for i in range(n+1)) for n in range(X+1)])
q=lambda f: fmpq(f.numerator,f.denominator)
def comb(terms):
    out=None
    for c,s in terms:
        t=scal(s,q(Fr(c)))
        out=t if out is None else add(out,t)
    return out
al={2:Fr(13,36),3:Fr(3,16),5:Fr(7,72),7:Fr(1,16)}[N]
be={2:Fr(5,9)}.get(N,2*al)
ga={2:[(Fr(4,45),E4s),(Fr(19,180),E4N)],
    3:[(Fr(3,160),E4s),(Fr(27,160),E4N)],
    5:[(Fr(-1,288),E4s),(Fr(-19,288),E4N),(Fr(1,6),ENsq)],
    7:[(Fr(-1,960),E4s),(Fr(-29,960),E4N),(Fr(3,32),ENsq)]}[N]
de={2:[(Fr(1,5),E4s),(Fr(-6,5),E4N)],
    3:[(Fr(1,20),E4s),(Fr(-21,20),E4N)],
    5:[(Fr(1,8),E4s),(Fr(15,8),E4N),(Fr(-3),ENsq)],
    7:[(Fr(1,30),E4s),(Fr(7,15),E4N),(Fr(-3,2),ENsq)]}[N]
gam=comb(ga); dl=comb(de)
rhs=comb([(al,mul(B,B)),(-be,mul(mul(ENs,B),A)),(1,mul(mul(gam,A),A))])
d=add(E2,rhs,-1); bad1=[a for a in range(X+1) if d.rows[a]!=0]
ws={2:(1,1,0,1),3:(-1,-1,0,0),5:(0,-1,-1,-1),7:(0,-1,0,0)}[N]
th=theta_rows(np.array(order(N).tolist(),dtype=int),ws,X)
T=series_from_dict({k:v for k,v in th.items() if k[0]<=X},X)
d2=add(add(mul(ENs,B),mul(dl,A)),scal(T,fmpq(12)),-1); bad2=[a for a in range(X+1) if d2.rows[a]!=0]
print('N=%d through x^%d: numerator identity mismatches %s, theta identity mismatches %s'%(N,X,bad1[:5],bad2[:5]))
