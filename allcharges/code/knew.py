# Bound for the shadow (K) terms of the corrected expansion (Theorem exp, new form), relative to the
# extreme term E, at x = xmin (worst case; the ratio is decreasing in x).  Usage: knew.py N m0 m1 [D0]
import sys, math, mpmath
from mpmath import mpf
from exact import ginv
from errterms import Eext
bL=0.7858   # Landau: |J_nu(y)| <= bL y^{-1/3}, nu>0, y>0
def frakB(A,B0,nu):   # (A/B0)^{nu/2} I_nu(2 sqrt(A B0)),  B0>0
    y=2*mpmath.sqrt(A*B0); return mpmath.besseli(nu,y)*(mpf(A)/B0)**(mpf(nu)/2)
def Bzero(A,nu): return A**nu/math.gamma(nu+1)
def sinsum(m,g):   # sum over h0 = -l mod 2m, h0 != 0 mod 2mg, of 1/|sin(pi h0/(2mg))|  (bound, uniform in l)
    return 2*m*g+2*g*(1+math.log(g))
def epsK(N,m,x,d,G=60):
    kp2=24//(N+1); k=kp2-2; nu=k+1.5
    D=mpf(x*m)**2; t0=m*m+4*m; E=(mpf(2)/N)*mpmath.pi/mpmath.sqrt(2*m)*(m+2)*mpf(N)**(-mpf(kp2)/2)*(t0/D)**(mpf(nu)/2)*mpmath.besseli(nu,mpmath.pi*x*mpmath.sqrt(t0)); fm=abs(mpf(d[m+1]))
    tot=0.0
    for g in range(1,G+1):
        for (w,mu,cnt) in ((N**(-kp2/2),1.0,g),(1.0,1.0/N,(N-1)*g)):
            A=2*math.pi*mu/g; B0=mpmath.pi*D/(2*m*g); xp=mpmath.sqrt(D)/(2*m*g)
            b2=frakB(A,B0,k+2); b3=frakB(A,B0,k+3); n3=k+3
            tau=xp*(Bzero(A,n3)+bL*2**(-1/3)*A**(n3/2-1/6)*(0.75*B0)**(-n3/2-1/6)*2**(2/3-n3)/(n3-2/3))
            br=math.pi**2*sinsum(m,g)*b2 + 2*math.pi*b2 + 8*math.pi*m*g*(xp*b3+tau)
            tot+=cnt*w*fm/(2*math.pi*g)*br
    # tail g>G: every bracket is <= poly(g) * exp(2 pi sqrt(mu D/m)/g)*g^{-(k+2)} ... bounded by last term * G (crude, decay >= g^{-2})
    return tot/(N*E)
if __name__=='__main__':
    N=int(sys.argv[1]); m0=int(sys.argv[2]); m1=int(sys.argv[3]); D0=float(sys.argv[4]) if len(sys.argv)>4 else 400
    kp2=24//(N+1); d=[float(v) for v in ginv(N,m1+5,kp2)]
    x0=math.sqrt(3)/N; worst=0
    for m in range(m0,m1+1):
        x=max(x0,math.sqrt(D0)/m); v=epsK(N,m,x,d); worst=max(worst,v)
        if m%25==0 or m==m0: print(N,m,round(x,4),'%.3e'%v,flush=True)
    print(N,'max',worst)
