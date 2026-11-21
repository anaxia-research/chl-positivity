import sys,time; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from chl import *
from fractions import Fraction as Fr
def inA(N,m,n,l):
    D=Fr(4*m*n,N)-l*l
    if not(m>0 and n>0 and D>0): return False
    # m1 in N Z \ {0}; need (m1 l+2ms+m)^2 + m1^2 D >= m^2
    m1=N
    while m1*m1*D < m*m:
        for sg in (1,-1):
            M1=sg*m1
            # s with M1 | s(s+1): bracket min over s; search s near minimizer
            s0=-(M1*l+m)/(2*m)
            for s in range(int(s0)-abs(M1)-2,int(s0)+abs(M1)+3):
                if (s*(s+1))%M1==0 and (M1*l+2*m*s+m)**2+M1*M1*D < m*m: return False
        m1+=N
    return True
N=int(sys.argv[1]); X=int(sys.argv[2]); M=int(sys.argv[3])
t=time.time(); C=CHL(N,X,M); F=finite_parts(C); print('built',time.time()-t,flush=True)
cnt=0;bad=[];notA=0;minrat=None
for m in range(1,M+1):
    for n in range(m,X):  # Fricke: m<=n
        for l in range(0, m+1):
            if Fr(4*m*n,N)-l*l<=0: continue
            if not inA(N,m,n,l): notA+=1; continue
            d=dstar(F,N,m,n,l); cnt+=1
            if d<=0: bad.append((m,n,l,d))
print(N,'checked',cnt,'charges in A; not in A:',notA,'; nonpositive:',bad[:20])
