# Orbits of charges (m,n,l) under Gamma_1(N); membership in BSS set A
import math, sys
from fractions import Fraction as Fr
def q(m,n,l,N,a,b): return m*a*a+N*l*a*b+N*n*b*b
def inA(m,n,l,N):
    """BSS (1.3): (m1 l+2ms+m)^2 + m1^2 D >= m^2 for m1 in NZ, m1 | s(s+1); D=(4mn-N l^2)/N.
       Multiply by N: N(m1 l+2ms+m)^2 + m1^2 (4mn-N l^2) >= N m^2."""
    DN=4*m*n-N*l*l
    if DN<=0: return False
    k=1
    while (N*k)**2*DN < N*m*m:      # else m1^2 D >= m^2 already
        for m1 in (N*k,-N*k):
            # need |m1 l + m(2s+1)| < m  -> s in a small window
            c=-(m1*l)/(2*m)-0.5
            for s in range(math.floor(c)-2, math.ceil(c)+3):
                if (s*(s+1))%m1==0 and N*(m1*l+2*m*s+m)**2+m1*m1*DN < N*m*m:
                    return False
        k+=1
    return True
def minimal(m,n,l,N):
    """True if q(a,b) >= m for all a = +-1 mod N, gcd(a,b)=1, (a,b) != (+-1,0)."""
    # q(a,b) = m |a + b w|^2, Im w = N sqrt(D)/(2m); q<m needs |b| Im w < 1
    DN=4*m*n-N*l*l
    # |b| * N*sqrt(D)/(2m) < 1  <=>  b^2 N^2 D < 4 m^2  <=> b^2 N DN < 4 m^2
    b=1
    while b*b*N*DN < 4*m*m:
        for bb in (b,-b):
            # a near -bb*Re w = -bb*N*l/(2m)
            c=-bb*N*l/(2*m)
            for a in range(math.floor(c)-2, math.ceil(c)+3):
                if a%N in (1,N-1) and math.gcd(a,bb)==1 and q(m,n,l,N,a,bb)<m: return False
        b+=1
    return True
def reduced_reps(N,DNmax):
    """all charges with |l|<=m, |l|<=n, D>0, minimal m in orbit; key = N*D"""
    out=[]
    mmax=int(math.isqrt(DNmax))+2
    for m in range(1,10**9):
        # minimality with n>=|l|: DN=4mn-Nl^2 >= 4m|l|-N l^2 ... loop until no rep
        found=False
        for l in range(-m,m+1):
            n0=max(abs(l),1)
            for n in range(n0,10**9):
                DN=4*m*n-N*l*l
                if DN>DNmax: break
                if DN<=0: continue
                found=True
                if minimal(m,n,l,N): out.append((m,n,l))
        if not found and m>mmax: break
        if m> DNmax//4+2: break
    return out
if __name__=='__main__':
    N=int(sys.argv[1]); DNmax=int(sys.argv[2])
    R=reduced_reps(N,DNmax)
    bad=[r for r in R if not inA(*r,N)]
    print('N',N,'reps',len(R),'not in A',len(bad),bad[:20])
