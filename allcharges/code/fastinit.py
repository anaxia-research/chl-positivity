# fast construction of the twisted-genus tables used by CHL.E
import math
from flint import fmpz_poly
from fractions import Fraction
from chl import CHL, sigma, EN_coeffs
def eta_pow_inv(e,T):  # prod (1-q^n)^{-e} up to q^T
    p=fmpz_poly([1])
    # use log-derivative recursion: a(n) = (e/n) sum_{k=1}^n sigma(k) a(n-k)
    sig=[0]+[sigma(k) for k in range(1,T+1)]
    a=[1]+[0]*T
    for n in range(1,T+1):
        a[n]=e*sum(sig[k]*a[n-k] for k in range(1,n+1))//n
    return a
def theta_tables(T):
    """alpha(D), c01(D) for D <= 4T (D=4n-r^2)"""
    # in variable Q=q^{1/4}: theta_{1,mu}(tau,0)=sum_{r=mu mod 2} Q^{r^2}
    TT=4*T+4
    th0=[0]*(TT+1); th1=[0]*(TT+1)
    r=0
    while r*r<=TT:
        (th0 if r%2==0 else th1)[r*r]+= (1 if r==0 else 2); r+=1
    eta6=eta_pow_inv(6,T+1)          # q-series of prod(1-q^n)^-6
    E6=[0]*(TT+1)
    for i,v in enumerate(eta6):
        if 4*i<=TT: E6[4*i]=v
    # 1/eta^6 = q^{-1/4} * prod^-6 ; in Q: Q^{-1} * E6
    P0=fmpz_poly(th0)*fmpz_poly(E6); P1=fmpz_poly(th1)*fmpz_poly(E6)
    # phi_{-2,1}: h0 = -theta11(0)/eta^6, h1 = theta10(0)/eta^6 ; coefficient of q^n zeta^r is h_{r%2} at Q^{4n-r^2}
    # h has Q-expansion Q^{-1}*(...) so alpha(D) = coeff of Q^{D+1} in P
    def coef(P,k):
        return int(P[k]) if 0<=k<=P.degree() else 0
    alpha={}
    for D in range(-1,4*T+1):
        if D%4==0: alpha[D]=-coef(P1,D+1)
        elif D%4==3: alpha[D]=coef(P0,D+1)
    return alpha

def c01_table(alpha,T):
    """phi_{0,1} coefficients from the heat-operator identity
       c01(D) = -6 D alpha(D) - 5 sum_k e2(k) alpha(D-4k),  E2 = 1 - 24 sum sigma(k) q^k"""
    e2=[1]+[-24*sigma(k) for k in range(1,T+1)]
    A=[0]*(4*T+2); E=[0]*(4*T+2)
    for D,v in alpha.items():
        if D+1<len(A): A[D+1]=v
    for k,v in enumerate(e2):
        if 4*k<len(E): E[4*k]=v
    P=fmpz_poly(A)*fmpz_poly(E)
    c={}
    for D,v in alpha.items():
        conv=int(P[D+1]) if D+1<=P.degree() else 0
        c[D]=-6*D*v-5*conv
    return c

class FastCHL(CHL):
    def __init__(s,N,X,M):
        s.N=N; s.X=X; s.M=M
        kp2=24//(N+1); s.kp2=kp2; s.a=Fraction(2,N*(N+1))
        Tq=(X*(M+1))//N+3
        al=theta_tables(Tq); c01=c01_table(al,Tq)
        e=[int(v) for v in EN_coeffs(N,N*Tq+N)]  # EN integral
        s.Phi01={};s.H={};s.G={}
        for b in (0,1):
            a=[al.get(4*n-b,0) for n in range(Tq+1)]
            for n in range(Tq+1):
                v=c01.get(4*n-b,0)
                if v: s.Phi01[(N*n,b)]=v
            # H(N t, b) = sum_{n+u=t} a(n) e(u)
            h=fmpz_poly(a)*fmpz_poly(e[:Tq+1])
            for t in range(Tq+1):
                v=int(h[t]) if t<=h.degree() else 0
                if v: s.H[(N*t,b)]=v
            # G(t,b) = sum_{N n+u=t} a(n) e(u)
            aN=[0]*(N*Tq+1)
            for n in range(Tq+1): aN[N*n]=a[n]
            g=fmpz_poly(aN)*fmpz_poly(e[:N*Tq+1])
            for t in range(N*Tq+1):
                v=int(g[t]) if t<=g.degree() else 0
                if v: s.G[(t,b)]=v
