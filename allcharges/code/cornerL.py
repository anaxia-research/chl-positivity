# per-charge evaluation of the leading-family sum L_m(x,theta) for corner charges (numpy, float)
import sys, math, numpy as np, mpmath
from exact import ginv
from sweep import minimal0
N=int(sys.argv[1]); kp2=24//(N+1); k=kp2-2; nu=k+1.5
d=np.array([float(v) for v in ginv(N,400,kp2)])  # d[A] = d~(A) = d(A-1)
def ive_halfint(nu,z):
    """I_nu(z) e^{-z}, nu half-integer, by upward recurrence from nu=-1/2, 1/2"""
    z=np.asarray(z,dtype=float); e=np.exp(-2*z); f=np.sqrt(2/(np.pi*z))
    im=f*(1+e)/2; i0=f*(1-e)/2; v=0.5
    while v<nu-1e-9:
        im,i0=i0,im-(2*v/z)*i0; v+=1
    return i0
def pairs(m,Smax):
    A,B=np.meshgrid(np.arange(Smax),np.arange(Smax),indexing='ij')
    A=A.ravel();B=B.ravel(); a=A-1;b=B-1; kk=m-a-b; t=kk*kk-4*a*b
    ok=(A+B<Smax)&(kk>=1)&(t>0)&(b-a>=-m)&(b-a<m)
    return A[ok],B[ok],kk[ok],t[ok]
def Lval(m,x,th,P):
    A,B,kk,t=P; c=math.pi*x; t0=m*m+4*m
    st=np.sqrt(t.astype(float)); s0=math.sqrt(t0)
    # rho = kk/(m+2) (t/t0)^{nu/2} I(c st)/I(c s0)
    zz=c*st
    # half-integer Bessel ratio via exact elementary form is faster; use mpmath on unique t
    ut,inv=np.unique(t,return_inverse=True)
    iv=ive_halfint(nu,c*np.sqrt(ut.astype(float)))[inv]
    i0=float(ive_halfint(nu,[c*s0])[0])
    rho=kk/(m+2)*(t/t0)**(nu/2)*iv/i0*np.exp(c*(st-s0))
    return float(np.sum(d[A]*d[B]*rho*np.cos(th*(A-B))))
if __name__=='__main__':
    from lead_min import P as Pf
    import os; Smax=int(sys.argv[2]); D0=float(os.environ.get("D0","400")); x1=float(sys.argv[3])
    for m in map(int,sys.argv[4:]):
        PP=pairs(m,Smax); worst=(9,None)
        cnt=0
        for l in range(0,m+1):
            for n in range(max(l,m),10**6):
                DN=4*m*n-N*l*l
                if DN<=N*D0: continue
                D=DN/N; x=math.sqrt(D)/m
                if x>x1: break
                if not minimal0(m,n,l,N): continue
                cnt+=1
                th=math.pi*l/m; L=Lval(m,x,th,PP)
                p2=abs(Pf(complex(math.exp(-math.pi*x)*math.cos(th),math.exp(-math.pi*x)*math.sin(th)),N,kp2))**2
                if L<worst[0]: worst=(L,(m,n,l,round(x,3),round(p2,4)))
        print(m,'charges',cnt,'min L',worst,flush=True)
