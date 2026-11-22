# Duality check for the twined index: for charges whose attractor point lies in the closure of C0,
# compare the C0-chamber coefficient with that of a Gamma_0(N)-equivalent reduced charge.
import sys, math
sys.path.insert(0,'.')
from exact_tw import Tw
def act(N,m,n,l,g):
    (a,b),(c,d)=g
    m2=a*a*m+N*a*b*l+N*b*b*n; n2=c*c*m//N+c*d*l+d*d*n; l2=(2*a*c*m)//N+(a*d+b*c)*l+2*b*d*n
    assert (c*c*m)%N==0 and (2*a*c*m)%N==0
    return m2,n2,l2
def inC0(N,m,n,l):
    D=4*m*n/N-l*l
    if D<=0 or not (0<=l<=2*m): return False
    x=l/(2*m); y=math.sqrt(D)/(2*m)
    for L in range(N,4*N*m+N,N):          # walls |tau - j/(2L)| = 1/(2L), L in N Z
        for j in range(-1,4*L+2,2):
            if (x-j/(2*L))**2+y*y < 1/(4*L*L)-1e-12: return False
    return True
def reduce(N,m,n,l):
    # greedy: lower m using (a,b) with q(a,b)<m, a=+-1 mod N; then shift l, flip sign, exchange
    for _ in range(200):
        changed=False
        if n<m: m,n=n,m; changed=True
        k=round(-l/(2*m)); m,n,l=act(N,m,n,l,((1,0),(N*k,1))) if k else (m,n,l)
        if l<0: l=-l
        best=None; B=1
        while B<=40:
            for b in (B,-B):
                c=-b*N*l/(2*m)
                for a in range(math.floor(c)-2,math.ceil(c)+3):
                    if a%N and math.gcd(a,b)==1:
                        v=m*a*a+N*l*a*b+N*n*b*b
                        if v<m and (best is None or v<best[0]): best=(v,a,b)
            B+=1
        if best:
            v,a,b=best
            # complete (a b; c d) with a d - b c = 1, N | c
            for cp in range(-200,201):
                c=N*cp
                if (1+b*c)%a==0:
                    d=(1+b*c)//a; break
            m,n,l=act(N,m,n,l,((a,b),(c,d))); changed=True
        if not changed: return m,n,l
    raise RuntimeError
for N,M,NM in ((2,9,25),(3,8,26)):
    T=Tw(N,(M+1)//N+1,(NM+1)//N+12); T.run()
    pairs=0; bad=0; nz=0
    for m in range(N-1,M+1,N):
        for n in range(N-1,NM+1,N):
            for l in range(0,2*m+1):
                if not inC0(N,m,n,l): continue
                r=reduce(N,m,n,l)
                if r==(m,n,l): continue
                if (r[0]+1)//N>(M+1)//N+1 or (r[1]+1)//N>(NM+1)//N+12: continue
                a=T.d(m,n,l); b=T.d(*r); pairs+=1; nz+=(a!=0)
                if abs(a)!=abs(b) or a!=b: bad+=1; print('MISMATCH',N,(m,n,l),r,a,b)
    print('N',N,'pairs',pairs,'nonzero',nz,'mismatches',bad,flush=True)
