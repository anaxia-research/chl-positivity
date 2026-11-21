import sys, math
from fastinit import FastCHL
from chl import finite_parts, dstar
from reduce import q, minimal, inA
def act(N,m,n,l,g):
    (a,b),(c,d)=g
    m2=a*a*m+N*a*b*l+N*b*b*n
    assert (c*c*m)%N==0
    n2=c*c*m//N+c*d*l+d*d*n
    l2=(2*a*c*m)//N+(a*d+b*c)*l+2*b*d*n
    assert (2*a*c*m)%N==0
    return m2,n2,l2
def canon1(N,m,n,l):
    DN=4*m*n-N*l*l
    best=(m,1,0); B=1
    while B*B*N*DN<4*m*m*  (1 if True else 1) *1+ 4*m*m:  # generous
        for bb in range(-B,B+1):
            if bb==0: continue
            c=-bb*N*l/(2*m)
            for a in range(math.floor(c)-3,math.ceil(c)+4):
                if a%N in (1,N-1) and math.gcd(a,bb)==1:
                    v=q(m,n,l,N,a,bb)
                    if v<best[0]: best=(v,a,bb)
        B+=1
        if B>200: break
    v,a,b=best
    s=1
    if a%N==N-1: a,b,s=-a,-b,-1
    # complete to Gamma1(N): a d - b c = 1, c = N c'
    # solve a d - b N c' = 1
    g,x,y=egcd(a,-b*N) if b!=0 else (a,1,0)
    if b==0: d,c=1,0
    else:
        d,cp=x,y; c=N*cp
        # adjust d = 1 mod N automatically
    assert a*d-b*c==1
    m2,n2,l2=act(N,m,n,l,((a,b),(c,d)))
    # shift l into [-m,m] by (1 0; N k 1)? l -> l + 2 m2 k
    k=round(-l2/(2*m2)); 
    m3,n3,l3=act(N,m2,n2,l2,((1,0),(N*k,1)))
    return m3,n3,abs(l3)
def egcd(a,b):
    if b==0: return (a,1,0) if a>0 else (-a,-1,0)
    g,x,y=egcd(b,a%b); return g,y,x-(a//b)*y
N=int(sys.argv[1]); X=int(sys.argv[2]); M=int(sys.argv[3])
C=FastCHL(N,X,M); F=finite_parts(C)
def ds(m,n,l):
    if m>M or n>X-2: return None
    return dstar(F,N,m,n,l)
gens={5:[((2,-1),(5,-2)),((-2,-1),(5,2))],7:[((3,-1),(7,-2)),((-2,-1),(7,3))],3:[((1,-1),(3,-2))]}[N]
res=[]
for m in range(1,M+1):
    for n in range(1,X-2):
        for l in range(0,m+1):
            if 4*m*n-N*l*l<=0 or not minimal(m,n,l,N) or not inA(m,n,l,N): continue
            v=ds(m,n,l)
            for g in gens:
                im=canon1(N,*act(N,m,n,l,g))
                if inA(*im,N) and im[0]<=M:
                    w=ds(*im)
                    if w is not None: res.append(((m,n,l),im,int(v.p) if hasattr(v,'p') else v,int(w.p) if hasattr(w,'p') else w))
same=sum(1 for r in res if r[2]==r[3]); print(N,'pairs',len(res),'equal',same)
for r in res[:12]: print(r)
