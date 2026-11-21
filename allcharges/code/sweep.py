import math, sys, time, json
from fractions import Fraction
from reduce import q, inA
def minimal0(m,n,l,N):
    """q(a,b) >= m for all coprime (a,b), N does not divide a, b != 0"""
    DN=4*m*n-N*l*l; b=1
    while b*b*N*DN < 4*m*m:
        for bb in (b,-b):
            c=-bb*N*l/(2*m)
            for a in range(math.floor(c)-2, math.ceil(c)+3):
                if a%N and math.gcd(a,bb)==1 and q(m,n,l,N,a,bb)<m: return False
        b+=1
    return True
def reps(N,D0):
    """Gamma_0(N)-minimal charges with 0<=l<=m, 0<D<=D0, oriented by Fricke: (min,max,l)"""
    out=set(); DNmax=N*D0
    for m in range(1,DNmax//4+3):
        for l in range(0,m+1):
            for n in range(max(l,1),10**9):
                DN=4*m*n-N*l*l
                if DN>DNmax: break
                if DN<=0: continue
                if minimal0(m,n,l,N):
                    assert inA(m,n,l,N)
                    out.add((min(m,n),max(m,n),l))
    return sorted(out)
if __name__=='__main__':
    N=int(sys.argv[1]); D0=int(sys.argv[2])
    R=reps(N,D0)
    print(N,D0,len(R),'max m',max(r[0] for r in R),'max n',max(r[1] for r in R))
    xs=min(math.sqrt(Fraction(4*m*n-N*l*l,N))/m for m,n,l in R); print('min x',xs)
    json.dump(R,open(f'reps_{N}_{D0}.json','w'))
