# Certify |d_g| < d_*/(M-1) for every reduced(-superset) charge with D_lo < D <= D_big, via boxes in (n,l) for each m.
# Twined bound: U(T;r) (Cauchy estimate along the family, uniform R_max); untwined: Part I, d_* >= c^- E(m,D), D > D0.
import math, sys, time, json
from tailG import U, logLower
def Sopt(t,m,n,l,kappa,eps,Rmax):
    a,b=-14.0,14.0
    f=lambda u: U(t,m,n,l,kappa,eps,Rmax,math.exp(u))
    for _ in range(50):
        c=a+(b-a)*0.382; d=a+(b-a)*0.618
        if f(c)<f(d): b=d
        else: a=c
    return math.exp((a+b)/2)
def run(N,kappa,eps,Rmax,Dlo,Dbig,extra,safety=1.0,log=None):
    t=N; boxes=0; worst=1e9; cert=1e9; perl=0; t0=time.time()
    mmax=int(N*math.sqrt(Dbig/3))+1
    P=kappa/(4*t*t)/(1-eps*eps*t)
    for m in range(1,mmax+1):
        # n range: reduced needs n>=m and D>=3m^2/N^2 ; we need D in (Dlo, Dbig] for some l in [0,m]
        nlo=m; nhi=math.floor(N*(Dbig+m*m)/(4*m))
        n1=nlo
        while n1<=nhi:
            # skip n where even l=0 gives D<=Dlo
            if 4*m*n1/N<=Dlo: n1=math.floor(N*Dlo/(4*m))+1; continue
            l1,l2=0,m
            # effective minimal D for lower bound in box: D(n1,l2) but at least max(Dlo, 3m^2/N^2)
            def lowerD(n,l): return max(4*m*n/N-l*l, Dlo+1e-9, 3*m*m/N**2)
            r=Sopt(t,m,n1,l1,kappa,eps,Rmax)
            Ys=math.sqrt(P/r)
            slope=2*math.pi*Ys/N
            Lw=logLower(N,m,lowerD(n1,l2))
            u1=U(t,m,n1,l1,kappa,eps,Rmax,r)+extra
            gap=Lw-u1-safety
            if gap<0:
                # split l-range: handle each l separately at this n1 with step 1
                ok=True
                for l in range(0,m+1):
                    Dv=4*m*n1/N-l*l
                    if Dv<=Dlo or Dv<3*m*m/N**2 or Dv>Dbig: continue
                    rr=Sopt(t,m,n1,l,kappa,eps,Rmax)
                    g=logLower(N,m,Dv)-U(t,m,n1,l,kappa,eps,Rmax,rr)-extra
                    worst=min(worst,g); cert=min(cert,g); boxes+=1; perl+=1
                    if g<=0: print('FAIL',m,n1,l,g,flush=True); ok=False
                n1+=1; continue
            step=max(1,int(gap/slope))
            n2=n1+step-1
            worst=min(worst,gap+safety); cert=min(cert,gap+safety-(n2-n1)*slope); boxes+=1
            n1=n2+1
        if log and m%200==0: print('m',m,'boxes',boxes,'corner margin %.2f'%worst,'certified %.3f'%cert,'per-l',perl,'%.0fs'%(time.time()-t0),flush=True)
    return boxes,worst,cert,perl
if __name__=='__main__':
    N=int(sys.argv[1]); kappa=float(sys.argv[2]); eps=float(sys.argv[3]); Rmax=float(sys.argv[4])
    Dlo=float(sys.argv[5]); Dbig=float(sys.argv[6])
    extra=math.log(2) if N==3 else 0.0
    safety=float(sys.argv[7]) if len(sys.argv)>7 else 1.0
    b,w,c,p=run(N,kappa,eps,Rmax,Dlo,Dbig,extra,safety=safety,log=True)
    print('N',N,'Dlo',Dlo,'Dbig',Dbig,'safety',safety,'boxes',b,'corner margin',w,'certified margin',c,'per-l checks',p)
