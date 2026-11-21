# numerical level-N Rademacher sums (polar part only) for testing normalisations
import sys, math
from fractions import Fraction as Fr
from mpmath import mp, mpf, mpc, besseli, sqrt, pi, exp, expj
mp.dps=30
def e(x): return expj(2*pi*x)
def polar_data(N,X,M):
    from exact import Exact
    from exactS import ExactS
    E=Exact(N,X,M); ES=ExactS(N,X//N+2,M)
    kp2=24//(N+1)
    out={}
    for m in range(1,M+1):
        T=E.table(m); TS=ES.table(m)
        inf=[(Fr(n,N),l,v) for (n,l),v in T.items() if l<2*m and Fr(4*m*n,N)-l*l<0 and v]
        S=[(Fr(n),l,v) for (n,l),v in TS.items() if l<2*m and 4*m*n-l*l<0 and v]
        out[m]=(inf,S,T)
    return out
def W(m,g,a,d,lt,l):
    s=mpc(0)
    for T in range(g):
        x=lt-2*m*T
        s+=e(Fr(a*x*x-2*l*x+d*l*l,4*m*g))
    return s/sqrt(2*m*g*mpc(0,1))
def Msum(N,m,n,l,inf,S,G=30):
    k=24//(N+1)-2; nu=k+mpf(3)/2
    D=Fr(4*m*n,N)-l*l; Dm=mpf(D.numerator)/D.denominator
    Minf=mpc(0); MS=mpc(0)
    for g in range(1,G+1):
        for md in range(0,N*g):
            d=-md
            if math.gcd(g,d)!=1: continue
            a=0 if g==1 else pow(d,-1,g)
            b=(a*d-1)//g
            if d%N==0: frame='S'; data=S; jf=None
            else: frame='T'; jf=(-b*pow(d,-1,N))%N; data=inf
            for (nt,lt,v) in data:
                Dt=4*m*nt-lt*lt
                ph=e(Fr(a*Dt.numerator,Dt.denominator*4*m*g)+Fr(d*D.numerator,D.denominator*4*m*g))
                if frame=='T': ph*=e(jf*nt)
                Dtm=mpf(-Dt.numerator)/Dt.denominator
                term=mpf(v.numerator)/v.denominator*ph*W(m,g,a,d,lt,l)*(2*pi/g)*(Dtm/Dm)**(nu/2)*besseli(nu,pi*sqrt(Dtm*Dm)/(m*g))
                if frame=='S': MS+=term
                else: Minf+=term
    return Minf,MS
if __name__=='__main__':
    N=int(sys.argv[1]); kp2=24//(N+1)
    PD=polar_data(N,60,3)
    rows=[]
    for m in (1,2,3):
        inf,S,T=PD[m]
        for n in (6,10,20,40):
            for l in range(0,m+1):
                if 4*m*n-N*l*l<=0: continue
                tgt=(1 if l%2 else -1)*T[(n,l)]
                Mi,Ms=Msum(N,m,n,l,inf,S,G=12)
                Ms=Ms/mpf(N)**(mpf(kp2)/2)
                pre=(-1)**l/mpf(N)
                rows.append((m,n,l,tgt,pre*Mi,pre*Ms))
    for r in rows: 
        t=mpf(r[3].numerator)/r[3].denominator
        print(r[0],r[1],r[2],float(t), complex(r[4]/t), complex(r[5]/t))
