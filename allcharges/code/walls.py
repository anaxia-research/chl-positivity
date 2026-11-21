import sys, math
from fractions import Fraction as Fr
from exact import ginv
from exactS import ExactS
def walls(m,nu,lt,R=None):
    """G=(p q; r s) in SL2(Z), r,s>=1, q/s <= w0 < p/r, with mG,nG >= -1"""
    w0=Fr(lt,2*m); out=[]
    R=R or 2*m+2
    for r in range(1,R+1):
        for s in range(1,R+1):
            if math.gcd(r,s)!=1: continue
            # p s - q r = 1; q/s <= w0 < p/r ; p/r - q/s = 1/(rs)
            # q = floor(w0*s) candidates
            for q in range(math.floor(w0*s)-2, math.floor(w0*s)+3):
                if (1+q*r)%s: continue
                p=(1+q*r)//s
                if not (Fr(q,s)<=w0<Fr(p,r)): continue
                mG=r*r*nu+p*p*m-p*r*lt; nG=s*s*nu+q*q*m-q*s*lt
                lG=-2*r*s*nu-2*p*q*m+(p*s+q*r)*lt
                if mG>=-1 and nG>=-1: out.append((p,q,r,s,mG,nG,lG))
    return out
if __name__=='__main__':
    N=int(sys.argv[1]); M=int(sys.argv[2]); kp2=24//(N+1)
    E=ExactS(N,M+3,M); d=ginv(N,4*M*M+20,kp2)   # d[k+1] = coeff of q^k
    D=lambda k: d[k+1] if k>=-1 else 0
    bad=0;tot=0
    for m in range(1,M+1):
        T=E.table(m)
        for (nu,lt),v in T.items():
            if lt>=2*m or 4*m*nu-lt*lt>=0: continue
            s=sum(lG*D(mG)*D(nG) for (_,_,_,_,mG,nG,lG) in walls(m,nu,lt))
            tot+=1
            if s!=v: bad+=1; print('mismatch',m,nu,lt,v,s)
    print('N',N,'polar terms',tot,'mismatches',bad)
