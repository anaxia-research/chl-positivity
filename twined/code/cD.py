# Coefficients c_t(Delta) of phi_t = 2*phi_{0,t} (t=2,3), Delta = 4tN - l^2, exact, via theta functions.
#   phi_2 = eta^-6 * sum_cyc th_i(z)^2 th_j(z)^2 th_k(0)^2 ,  phi_3 = 2 eta^-6 th_2(z)^2 th_3(z)^2 th_4(z)^2
import sys, os, pickle
from math import isqrt
from flint import fmpz_poly
HERE=os.path.dirname(os.path.abspath(__file__))
def build(t, Nmax):
    E=8*Nmax+8                      # q exponent units 1/8
    Rr=2*isqrt(E)+8                  # half-units of r: |f| <= Rr per theta... product bound
    F=6*Rr+4; W=2*F+1                # Kronecker width for r (half-units), offset F
    def th(i,z):
        c=[0]*((E+1)*W)
        R=isqrt(2*E)+3
        for n in range(-R,R+1):
            if i==2: e=(2*n+1)**2; f=2*n+1; sg=1
            else: e=4*n*n; f=2*n; sg=1 if (i==3 or n%2==0) else -1
            if e<=E:
                ff=f if z else 0
                c[e*W+ff+F]+=sg
        return fmpz_poly(c)
    def trunc(p):
        c=p.coeffs(); 
        # remove offset accumulation: product of k polys has offset k*F ; handled by caller
        return c
    # build product with offset tracking: each factor shifts by F
    def prod(fs):
        p=fmpz_poly([1]); k=0
        for f in fs:
            p=p*f; k+=1
            c=p.coeffs()
            # truncate q-degree > E
            c=c[:(E+1)*W + k*F + W]
            p=fmpz_poly(c)
        return p,k
    T={(i,z):th(i,z) for i in (2,3,4) for z in (0,1)}
    if t==2:
        terms=[]
        for (i,j,k) in ((2,3,4),(3,4,2),(4,2,3)):
            terms.append(prod([T[(i,1)],T[(i,1)],T[(j,1)],T[(j,1)],T[(k,0)],T[(k,0)]]))
        mult=1
    else:
        terms=[prod([T[(2,1)],T[(2,1)],T[(3,1)],T[(3,1)],T[(4,1)],T[(4,1)]])]
        mult=2
    # extract 2D coefficients: index = e*W + f + kF  (k=6 factors) -> careful: offsets add as F per factor
    num={}
    for p,k in terms:
        c=p.coeffs()
        for idx,v in enumerate(c):
            if v==0: continue
            g=idx-k*F
            e,f=divmod(g+F,W); f-=F
            # g = e*W + f  with f in [-F,F]
            if e<=E: num[(e,f)]=num.get((e,f),0)+int(v)
    # eta^-6 = q^{-1/4} prod (1-q^n)^-6 ; units 1/8 -> shift -2
    P=[0]*(Nmax+3); P[0]=1
    for _ in range(6):
        for n in range(1,len(P)):
            for i in range(n,len(P)): P[i]+=P[i-n]
    out={}
    for l in range(0,t+1):
        # coefficient of r^l (half-units 2l) as series in q^(1/8): num[(e,2l)], then times eta^-6
        ser={}
        for (e,f),v in num.items():
            if f==2*l: ser[e]=v
        for N in range(0,Nmax+1):
            # total exponent 8N = e + 8j - 2  -> e = 8N+2-8j
            s=0
            for j in range(0,N+1):
                s+=ser.get(8*N+2-8*j,0)*P[j]
            out[(N,l)]=mult*s
    cD={}
    for (N,l),v in out.items():
        D=4*t*N-l*l
        if D in cD: assert cD[D]==v,(D,cD[D],v)
        cD[D]=v
    return cD
if __name__=='__main__':
    t=int(sys.argv[1]); Nmax=int(sys.argv[2])
    cD=build(t,Nmax)
    pickle.dump(cD,open(os.path.join(HERE,'cD_%d_%d.pkl'%(t,Nmax)),'wb'))
    print(sorted(cD.items())[:14])
    import math
    for D in [100,400,1000,2000,4*t*Nmax-1]:
        if D in cD and cD[D]: print(D, cD[D], math.log(abs(cD[D]))/(math.pi*math.sqrt(D)/t))
