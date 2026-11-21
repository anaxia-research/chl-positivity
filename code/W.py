import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from chl import sigma, EN_coeffs
from fractions import Fraction as Fr
T=60
def E4(T): return [Fr(1)]+[Fr(240*sigma(n,3)) for n in range(1,T+1)]
def dil(f,N,T): g=[Fr(0)]*(T+1); 
def dil(f,N,T):
    g=[Fr(0)]*(T+1)
    for i,v in enumerate(f):
        if N*i<=T: g[N*i]=v
    return g
def sq(f,T): return [sum(f[i]*f[n-i] for i in range(n+1)) for n in range(T+1)]
def eta(exps,T):
    f=[Fr(0)]*(T+1); f[0]=Fr(1)
    for d,e in exps.items():
        for n in range(1,T//d+1):
            c=1;k=0;terms=[]
            while k*d*n<=T: terms.append((k*d*n,c)); k+=1; c=c*(e-k+1)//k*(-1)
            g=[Fr(0)]*(T+1)
            for i,v in enumerate(f):
                if v:
                    for s,cc in terms:
                        if i+s<=T: g[i+s]+=v*cc
            f=g
    lead=sum(d*e for d,e in exps.items())//24
    return [Fr(0)]*lead+f[:T+1-lead]
for N,(a1,a2,a3) in {3:(Fr(3,80),Fr(-3,80),0),5:(Fr(1,48),Fr(19,48),Fr(-5,12)),7:(Fr(1,320),Fr(29,320),Fr(-3,32))}.items():
    e4=E4(T); e4N=dil(E4(T),N,T); EN=EN_coeffs(N,T); EN2=sq(EN,T)
    W=[a1*e4[t]+a2*e4N[t]+a3*EN2[t] for t in range(T+1)]
    # decompose: W = a*(E4 - E4(N tau)) + cusp part
    a=W[1]/240
    cusp=[W[t]-a*(e4[t]-e4N[t]) for t in range(T+1)]
    print(N,'W[:12]=',[str(w) for w in W[:12]], 'min over 1..T:',min(W[1:]))
    print('   a=',a,' cusp part[:8]=',[str(c) for c in cusp[:8]])
    if N==5:
        f=eta({1:4,5:4},T); print('   eta^4eta5^4 ratio',[str(cusp[t]/f[t]) if f[t] else '-' for t in range(2,8)])
    if N==7:
        # weight 4 cusp forms on Gamma0(7): dim 1: eta(t)^4 eta(7t)^4 has weight 4, level 7? check
        f=eta({1:4,7:4},T); print('   ratio to eta^4eta7^4',[str(cusp[t]/f[t]) if f[t] else '-' for t in range(1,8)])
print('----')
def leg(a,p):
    a%=p
    if a==0: return 0
    return 1 if pow(a,(p-1)//2,p)==1 else -1
f5=eta({5:5,1:-1},T)
print('eta5^5/eta',[str(x) for x in f5[:10]])
print('guess5    ',[sum(leg(n//d,5)*d for d in range(1,n+1) if n%d==0) if n else 0 for n in range(10)])
sq5=sq(f5,T); W5=[Fr(0)]*0
f7=eta({7:7,1:-1},T)
print('eta7^7/eta',[str(x) for x in f7[:12]])
for shift in (0,1,2):
  print('guess7',shift,[ (sum(leg(d,7)*(n//d)**2 for d in range(1,n+1) if n%d==0) if n>0 else 0) for n in range(0,12)])
th=[0]*(T+1)
for x in range(-20,21):
    for y in range(-20,21):
        v=x*x+x*y+2*y*y
        if v<=T: th[v]+=1
prod=[sum(th[i]*f7[n-i] for i in range(n+1)) for n in range(T+1)]
print('theta_Q7*eta7^7/eta',prod[:12])
print('====')
import sympy as sp
chi=lambda a: leg(a,7)
EA=[0]+[sum(chi(d)*(n//d)**2 for d in range(1,n+1) if n%d==0) for n in range(1,T+1)]
EB=[0]+[sum(chi(n//d)*d**2 for d in range(1,n+1) if n%d==0) for n in range(1,T+1)]
E0=[Fr(-8,7) if False else 0]  # placeholder
cf=eta({1:3,7:3},T)
M=sp.Matrix([[EA[n],EB[n],cf[n],1 if n==0 else 0] for n in range(0,T+1)])
v=sp.Matrix([f7[n] for n in range(T+1)])
try:
    s,p=M.gauss_jordan_solve(v); print('fit',s.T)
except Exception as e: print('nofit',e)
