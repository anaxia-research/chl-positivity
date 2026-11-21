# Fixed-curve check for elliptic elements h of Gamma0(N): the curve {z Y0}, z in H, of points fixed by
# Omega -> M^T Omega M (M=(a,Nb;c/N,d)) lies in no divisor of Phi~.  A divisor containing it is T^perp with
# tr(T Y0)=0, T=(m,l/2;l/2,n/N), D=4mn/N-l^2 in {-1, 4/N-1} (polar classes).  Y0^perp is negative definite.
import sympy as sp, itertools, sys
def check(N,h):
    a,b,c,d=h; M=sp.Matrix([[a,N*b],[sp.Rational(c,N),d]])
    y1,y2,y3=sp.symbols('y1 y2 y3'); Y=sp.Matrix([[y1,y2],[y2,y3]])
    sol=sp.solve(list(M.T*Y*M-Y),[y1,y2,y3],dict=True)[0]
    Y0=Y.subs(sol); free=list(Y0.free_symbols); Y0=Y0.subs(free[0],1) if free else Y0
    assert Y0.det()>0
    Y0=Y0*sp.sign(Y0[0,0])
    # tr(T Y0) = m Y0[0,0] + l Y0[0,1] + n/N Y0[1,1] = 0 ; enumerate with |D| small: Y0^perp neg def, so bounded
    found=[]
    B=60
    for m in range(-B,B+1):
        for l in range(-B,B+1):
            # solve for n
            if Y0[1,1]==0: continue
            n=-(m*Y0[0,0]+l*Y0[0,1])*N/Y0[1,1]
            if n!=int(n): continue
            n=int(n); D=sp.Rational(4*m*n,N)-l*l
            if D in (-1,sp.Rational(4,N)-1) and l%2==1:
                cls=((n%N),(m%N))
                ok=(D==-1 and cls==(0,0)) or (D==sp.Rational(4,N)-1 and (m*n)%N==1)
                if ok: found.append((m,n,l,D))
    return Y0,found
for N,h in [(5,(2,-1,5,-2)),(7,(3,-1,7,-2)),(3,(-2,1,-3,1))]:
    Y0,f=check(N,h); print(N,h,'Y0=',list(Y0),'divisor vectors on the curve:',f)
for h in [(-2,-1,5,2),(3,-2,5,-3),(-3,-2,5,3),(2,1,-5,-2)]:
    Y0,f=check(5,h); print(5,h,list(Y0),f)
