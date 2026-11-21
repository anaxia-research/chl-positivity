# independent check of the circle-method formula derived from the meromorphic Jacobi form
import math, sys
from fractions import Fraction as Fr
from mpmath import mp, mpf, mpc, besseli, besselj, sqrt, pi, sinh, quad, expj, inf
from rad import polar_data
from exact import ginv
mp.dps=20
def e(x): return expj(2*pi*x)
def chi(N,d):
    if N!=7: return 1
    d%=7
    return 0 if d==0 else (1 if pow(d,3,7)==1 else -1)
def W(m,g,a,d,lt,l):
    s=mpc(0)
    for T in range(g):
        x=lt-2*m*T
        s+=e(Fr(a*x*x-2*l*x+d*l*l,4*m*g))
    return s/sqrt(2*m*g*mpc(0,1))
def pairs(N,g):
    for md in range(0,N*g):
        d=-md
        if math.gcd(g,d)!=1: continue
        a=0 if g==1 else pow(d,-1,g)
        b=(a*d-1)//g
        yield a,b,d
def Bfun(A,B,nu):
    # (A/B)^{nu/2} I_nu(2 sqrt(AB)) as entire function of B
    AB=A*B
    if AB>0:
        y=sqrt(AB); return A**nu*y**(-nu)*besseli(nu,2*y)
    if AB<0:
        y=sqrt(-AB); return A**nu*y**(-nu)*besselj(nu,2*y)
    return A**nu/mp.gamma(nu+1)
def cF(N,m,n,l,inf_data,S_data,fh,G=8,GK=5):
    k=24//(N+1)-2; nu=k+mpf(3)/2
    D=Fr(4*m*n,N)-l*l; Dm=mpf(D.numerator)/D.denominator
    tot=mpc(0); parts={'MT':mpc(0),'MS':mpc(0),'KT':mpc(0),'KS':mpc(0)}
    for g in range(1,G+1):
        for a,b,d in pairs(N,g):
            if d%N:
                jp=(b*pow(d,-1,N))%N
                eps=chi(N,d)
                data=[(nt,lt,mpf(v.numerator)/v.denominator*e(-jp*nt)) for (nt,lt,v) in inf_data]
                cst=eps*mpf(fh)*e(Fr(jp,N)); mu=Fr(1,N); key='T'
            else:
                eps=chi(N,-b)*(-1)**k*mpc(0,1)**(k+2)
                data=[(nt,lt,mpf(v.numerator)/v.denominator/mpf(N)**(mpf(k+2)/2)) for (nt,lt,v) in S_data]
                cst=eps*mpf(fh)/mpf(N)**(mpf(k+2)/2); mu=Fr(1); key='S'
            # M-terms
            for (nt,lt,v) in data:
                Dt=4*m*nt-lt*lt; Dtm=mpf(-Dt.numerator)/Dt.denominator
                ph=e(Fr(a*Dt.numerator,Dt.denominator*4*m*g)+Fr(d*D.numerator,D.denominator*4*m*g))
                t=(2*pi/(N*g))*mpc(0,1)**(k+mpf(1)/2)*eps*v*ph*W(m,g,a,d,lt,l)*(Dtm/Dm)**(nu/2)*besseli(nu,pi*sqrt(Dtm*Dm)/(m*g))
                parts['M'+key]+=t
            # K-terms
            if g>GK: continue
            A=2*pi*mpf(mu.numerator)/mu.denominator/g
            pre=mpc(0,1)**k/(2*pi*N*g)*cst*e(-mu*a/g+Fr(d*D.numerator,D.denominator*4*m*g))
            sK=mpc(0)
            for h0 in range(2*m*g):
                if (h0+l)%(2*m): continue
                ph=e(Fr(d*h0*h0,4*m*g))
                Bx=lambda x: Bfun(A,pi*(Dm/(2*m*g)-2*m*g*x*x),k+2)
                if h0%(2*m*g):
                    c=mpf(h0)/(2*m*g)
                    f=lambda x: pi**2/sinh(pi*(x+1j*c))**2*Bx(x)
                    I=quad(f,[-inf,-2,0,2,inf])
                else:
                    def reg(x):
                        if abs(x)<mpf('1e-3'): return -pi**2/3+pi**4*x*x/15
                        return pi**2/sinh(pi*x)**2-1/x**2
                    f=lambda x: reg(x)*Bx(x)
                    B1=lambda x: Bfun(A,pi*(Dm/(2*m*g)-2*m*g*x*x),k+3)
                    I=quad(f,[-inf,-2,0,2,inf])-8*pi*m*g*quad(B1,[0,2,20,inf])
                sK+=ph*I
            parts['K'+key]+=pre*sK
    return parts
if __name__=='__main__':
    N=int(sys.argv[1]); kp2=24//(N+1)
    PD=polar_data(N,60,3); gi=ginv(N,10,kp2)
    for (m,n) in [tuple(map(int,s.split(','))) for s in sys.argv[2:]]:
        inf_data,S_data,T=PD[m]; fh=gi[m+1]
        for l in range(0,m+1):
            if 4*m*n-N*l*l<=0: continue
            tgt=T[(n,l)]; t=mpf(tgt.numerator)/tgt.denominator
            P=cF(N,m,n,l,inf_data,S_data,fh)
            tot=sum(P.values())
            print(N,m,n,l,float(t),'ratio',complex(tot/t),{k:complex(v/t) for k,v in P.items()},flush=True)
