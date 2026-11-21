# Lower bound for G(x,theta,h)=sum_{S<=S1} d~(A)d~(B) e^{-cS} cos(theta(A-B)) T_J(A,B,h;c) over
# x in [XA,XB], theta in [0,pi], h in [0,1/Ms]; cells of width W/nx in x (aligned with E-cells of width W).
import sys, math, pickle, json, numpy as np, sympy as sp
N=int(sys.argv[1]); Ms=int(sys.argv[2]); S1=int(sys.argv[3]); J=int(sys.argv[4])
XA=float(sys.argv[5]); XB=float(sys.argv[6]); W=float(sys.argv[7]); nx=int(sys.argv[8]); dth=float(sys.argv[9])
sys.argv=[sys.argv[0],str(N)]
from exact import ginv
kp2=24//(N+1); nu2=2*(kp2-2)+3
dd=np.array([float(v) for v in ginv(N,S1+5,kp2)])
co,(As,Bs,cs,ics)=pickle.load(open(f'hexp_{nu2}_{J}.pkl','rb'))
co=[sp.Poly(sp.sympify(e),As,Bs,cs,ics) for e in co]
def arr(p):
    t=p.terms(); E=np.array([m for m,_ in t],dtype=float).reshape(-1,4); K=np.array([float(v) for _,v in t]); return E,K
def dcarr(EK):   # d/dc with ic=1/c :  d(c^k ic^l)/dc = k c^{k-1} ic^l - l c^k ic^{l+1}
    E,K=EK; E1=E.copy(); E1[:,2]-=1; E2=E.copy(); E2[:,3]+=1
    En=np.vstack([E1,E2]); Kn=np.concatenate([K*E[:,2],-K*E[:,3]]); ok=Kn!=0
    return En[ok],Kn[ok]
def ev(EK,A,B,c,ic,absval=False):
    E,K=EK
    if len(K)==0: return np.zeros_like(A)
    K=np.abs(K) if absval else K
    pa=[A**k for k in range(int(E[:,0].max())+1)]; pb=[B**k for k in range(int(E[:,1].max())+1)]
    out=np.zeros_like(A)
    for (i,j,k,l),kv in zip(E.astype(int),K):
        out+=kv*(c**k)*(ic**l)*pa[i]*pb[j]
    return out
P0=[arr(p) for p in co]; P1=[dcarr(e) for e in P0]; P2=[dcarr(e) for e in P1]
f =[lambda A,B,c,ic,e=e: ev(e,A,B,c,ic) for e in P0]
f1=[lambda A,B,c,ic,e=e: ev(e,A,B,c,ic) for e in P1]
g0=[lambda A,B,c,ic,e=e: ev(e,A,B,c,ic,True) for e in P0]
g1=[lambda A,B,c,ic,e=e: ev(e,A,B,c,ic,True) for e in P1]
g2=[lambda A,B,c,ic,e=e: ev(e,A,B,c,ic,True) for e in P2]
A,B=np.meshgrid(np.arange(S1+1),np.arange(S1+1),indexing='ij'); ok=(A+B<=S1); A=A[ok].astype(float); B=B[ok].astype(float)
S=A+B; Y=A-B; one=np.ones_like(A); D2=dd[A.astype(int)]*dd[B.astype(int)]
hs=1.0/Ms; hp=np.array([hs**j for j in range(J+1)])
nth=int(math.ceil(math.pi/(2*dth))); th=(np.arange(nth)+0.5)*2*dth   # centres, half-width dth
C=np.cos(np.outer(Y,th)); Sn=np.sin(np.outer(Y,th))
res=[]
xa=XA
while xa<XB-1e-12:
    xb=min(xa+W,XB); hw=(xb-xa)/(2*nx); cellmin=9e9
    for ix in range(nx):
        xc=xa+(2*ix+1)*hw; c=math.pi*xc; clo=math.pi*(xc-hw); chi=math.pi*(xc+hw)
        w=D2*np.exp(-c*S)
        R=np.array([f[j](A,B,c,1/c)*one for j in range(J+1)])          # (J+1,P)
        R1=np.array([f1[j](A,B,c,1/c)*one for j in range(J+1)])
        Q=(w*R)@C                                   # (J+1,nth): coefficients of h^j of G
        Qx=math.pi*((w*(-S*R+R1))@C)
        Qt=-((w*Y*R)@Sn)
        def lo(V): return V[0]-np.abs(V[1:]).T@hp[1:]
        def hi(V): return np.abs(V[0])+np.abs(V[1:]).T@hp[1:]
        # termwise second-derivative sups over the cell, all h
        wl=np.abs(D2)*np.exp(-clo*S)
        a0=sum(hp[j]*g0[j](A,B,chi,1/clo)*one for j in range(J+1))
        a1=sum(hp[j]*g1[j](A,B,chi,1/clo)*one for j in range(J+1))
        a2=sum(hp[j]*g2[j](A,B,chi,1/clo)*one for j in range(J+1))
        Htt=float(np.sum(wl*Y*Y*a0)); Hxt=math.pi*float(np.sum(wl*np.abs(Y)*(S*a0+a1))); Hxx=math.pi**2*float(np.sum(wl*(S*S*a0+2*S*a1+a2)))
        quad=0.5*(Hxx*hw*hw+2*Hxt*hw*dth+Htt*dth*dth)
        LB=lo(Q)-hi(Qx)*hw-hi(Qt)*dth-quad
        i=int(np.argmin(LB)); cellmin=min(cellmin,float(LB[i]))
    res.append(dict(xa=xa,xb=xb,G=cellmin,quad=quad))
    print(round(xa,5),round(xb,5),'G_lb',cellmin,'quad',quad,flush=True)
    json.dump(res,open(f'Gcells_{N}_{Ms}_{S1}_{J}_{XA:.3f}.json','w'))
    xa=xb
