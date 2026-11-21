import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
N=int(sys.argv[1]); X=int(sys.argv[2])
sys.argv=['x',str(N),str(X)]
exec(open('hol3.py').read().split("mx=0")[0])   # gives T = 12 theta_O (verified below), A, B, M4...
from row1_lib import hurwitz
from quat import order
from lattice import theta_rows
import numpy as np
# lattice check
Gm=order(N); ws={3:(-1,-1,0,0),5:(0,-1,-1,-1),7:(0,-1,0,0)}[N]
nl=min(X,24)
th=theta_rows(np.array(Gm.tolist(),dtype=int),ws,nl)
assert all(12*th.get((n,l),0)==int(coeff(T,n,l)) for n in range(nl+1) for l in range(-10,11))
print('N=%d: E_N phi01 + c.A = 12 theta(O_N, w) through x^%d'%(N,nl))
b,fh={3:(Fr(-3,8),27),5:(Fr(-7,36),14),7:(Fr(-1,8),9)}[N]
Wc={3:(Fr(3,80),Fr(-3,80),0),5:(Fr(1,48),Fr(19,48),Fr(-5,12)),7:(Fr(1,320),Fr(29,320),Fr(-3,32))}[N]
q=lambda f: fmpq(f.numerator,f.denominator)
K=scal(T,q(b))
for c,f in zip(Wc,M4[:3]):
    if c: K=add(K,scal(mul(f,A),q(c)))
Hd={}
for n in range(0,X//N+1):
    for l in range(-2*int((n+1)**0.5)-2,2*int((n+1)**0.5)+3):
        h=hurwitz(4*n-l*l)
        if h: Hd[(N*n,l)]=h
K=add(K,scal(series_from_dict(Hd,X),fmpq(-2*fh)))
C2=CHL(N,X,1); F=finite_parts(C2)
gi=C2.ginv
inv=[Fr(0)]*(X+1); inv[0]=1/gi[0]
for i in range(1,X+1): inv[i]=-sum(gi[j]*inv[i-j] for j in range(1,i+1))/gi[0]
Gs=S([fmpq_poly(q(v)) for v in inv],0,1)
Kp=mul(F[1],Gs)
d=add(K,Kp,-1)
print('   g*phi^F_1 (from product) == W A + b*12theta - 2fhat H : mismatches through x^%d:'%(X-2),[a for a in range(X-1) if d.rows[a]!=0][:5])
