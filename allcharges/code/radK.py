import sys, math
from fractions import Fraction as Fr
from mpmath import mp, mpf, mpc, besseli, sqrt, pi, sinh, quad, expj
from rad import e, W, Msum, polar_data
from exact import ginv
mp.dps=25
def frames(N,g):
    for md in range(0,N*g):
        d=-md
        if math.gcd(g,d)!=1: continue
        a=0 if g==1 else pow(d,-1,g)
        b=(a*d-1)//g
        if d%N==0: yield a,b,d,'S',None
        else: yield a,b,d,'T',(-b*pow(d,-1,N))%N
def cmu(N,k,fh,fr,j):
    if fr=='S': return mpf(fh)/mpf(N)**(mpf(k+2)/2), mpf(1)
    return fh*e(Fr(-j,N)), mpf(1)/N
def fgg(m,g,gg,u):
    if gg%(2*m*g): return pi**2/sinh(pi*(u/g-1j*gg/(2*m*g)))**2
    if u==0: return -pi**2/3
    return pi**2/sinh(pi*u/g)**2-g*g/u**2
def K23(N,m,n,l,fh,G=6,do3=True):
    k=24//(N+1)-2
    D=Fr(4*m*n,N)-l*l; Dm=mpf(D.numerator)/D.denominator
    K2=mpc(0);K3=mpc(0)
    for g in range(1,G+1):
        for a,b,d,fr,j in frames(N,g):
            c,mu=cmu(N,k,fh,fr,j)
            ph=e(Fr(d*D.numerator,D.denominator*4*m*g))
            K2+= -sqrt(2*m)*g**-0.5*c*e(-a*mu*4*m/(4*m*g))*ph*W(m,g,a,d,0,l)*(4*m*mu/Dm)**(mpf(k+2)/2)*besseli(k+2,2*pi/g*sqrt(mu*Dm/m))
            if do3:
                for gg in range(-m*g,m*g):
                    pre=c*e(Fr(d*D.numerator,D.denominator*4*m*g))*e(-a*(4*m*mu+gg*gg)/(4*m*g))*W(m,g,a,d,gg,l)*(4*m*mu/Dm)**(mpf(2*k+5)/4)
                    lim=sqrt(mu/m)
                    I=quad(lambda u: fgg(m,g,gg,u)*(1-m*u*u/mu)**(mpf(2*k+5)/4)*besseli(k+mpf(5)/2,2*pi/g*sqrt((mu-m*u*u)*Dm/m)),[-lim,0,lim])
                    K3+=pre*I/(2*pi)*g**-2
    return K2,K3
if __name__=='__main__':
    N=int(sys.argv[1]); kp2=24//(N+1)
    PD=polar_data(N,60,3); gi=ginv(N,10,kp2)
    for m in (2,3):
        inf,S,T=PD[m]; fh=gi[m+1]
        for n in (12,40):
            for l in range(0,m+1):
                if 4*m*n-N*l*l<=0: continue
                tgt=(1 if l%2 else -1)*T[(n,l)]; t=mpf(tgt.numerator)/tgt.denominator
                Mi,Ms=Msum(N,m,n,l,inf,S,G=10); Ms=Ms/mpf(N)**(mpf(kp2)/2)
                K2,K3=K23(N,m,n,l,fh,G=5)
                pre=(-1)**l/mpf(N)*e(Fr(1,8))
                tot=pre*(Mi+Ms+K2+K3)
                print(m,n,l,float(t),'M only',float((pre*(Mi+Ms)).real/t),'with K',complex(tot/t))
