# 1/Phi for the Borcherds product  Phi = q^1 r^1 s^{C0} prod_{(n,l,m)>0} (1 - q^{t n} r^l s^m)^{c(nm,l)}
# series in q (int), s (int; prefactor handled separately), r (int, bounded below by Cmin)
import sys
from collections import defaultdict
from jac import phi
def inv_coeffs(t, Amax, Bmax, Cmin):
    """coefficients of prod (1-X)^{-c} truncated: q^A (A<=Amax), s^B (B<=Bmax), r^C (C>=Cmin)."""
    cmax_N = (Amax//t)*Bmax+1
    c = phi(t, cmax_N)
    c = {k:int(v) for k,v in c.items()}
    S = {(0,0,0):1}
    def mulfac(S, A, B, L, e):  # multiply by (1 - q^A s^B r^L)^(-e)
        # use repeated multiplication by 1/(1-X) or (1-X)
        for _ in range(abs(e)):
            if e>0:  # divide by (1-X): S += X*S iteratively (geometric)
                T=defaultdict(int)
                keys=sorted(S.keys(), key=lambda k:(k[0],k[1],-k[2]))
                T.update(S)
                # process in order of increasing A,B so X*T uses updated values
                order=sorted(set(S.keys()))
                # simple approach: T = S * sum_k X^k
                T=defaultdict(int)
                for (a,b,cc),v in S.items():
                    k=0
                    while True:
                        aa,bb,ll=a+k*A,b+k*B,cc+k*L
                        if aa>Amax or bb>Bmax or ll<Cmin or (A==0 and B==0 and k>400): break
                        T[(aa,bb,ll)]+=v; k+=1
                        if A==0 and B==0 and L==0: break
                S={k:v for k,v in T.items() if v}
            else:
                T=defaultdict(int); T.update(S)
                for (a,b,cc),v in S.items():
                    aa,bb,ll=a+A,b+B,cc+L
                    if aa<=Amax and bb<=Bmax and ll>=Cmin: T[(aa,bb,ll)]-=v
                S={k:v for k,v in T.items() if v}
        return S
    # (0,-1,0)
    S = mulfac(S,0,0,-1,c[(0,-1)])
    for n in range(0, Amax//t+1):
        for m in range(0, Bmax+1):
            if n==0 and m==0: continue
            for (N,l),v in c.items():
                if N==n*m and v: S = mulfac(S, t*n, m, l, v)
    return S
if __name__=='__main__':
    t=int(sys.argv[1])
    S=inv_coeffs(t, 4, 2, -12)
    # 1/Phi = q^-1 r^-1 s^-C0 * S
    for A in range(0,5):
        for B in range(0,3):
            row=[(C-1,S.get((A,B,C),0)) for C in range(-3,7)]
            print('q^%d s^(%d-C0):'%(A-1,B), row)
