import sys, json, time, os
from fractions import Fraction
from exact import Exact
from sweep import reps
N=int(sys.argv[1]); D0=int(sys.argv[2])
t0=time.time()
fn=f'reps_{N}_{D0}.json'
R=[tuple(r) for r in json.load(open(fn))] if os.path.exists(fn) else reps(N,D0)
json.dump(R,open(fn,'w'))
print('reps',len(R),'enum time',time.time()-t0,flush=True)
Mmax=max(r[0] for r in R)
bounds=[0]; b=1
while bounds[-1]<Mmax:
    bounds.append(min(b,Mmax)); b*=2
res=[]; bad=[]
MLO=int(sys.argv[3]) if len(sys.argv)>3 else 0
for lo,hi in zip(bounds,bounds[1:]):
    if hi<=MLO: continue
    sub=[r for r in R if lo<r[0]<=hi]
    if not sub: continue
    X=max(r[1] for r in sub)+1
    t=time.time(); E=Exact(N,X,hi)
    for m in range(lo+1,hi+1):
        S=[r for r in sub if r[0]==m]
        if not S: continue
        T=E.table(m)
        for (mm,n,l) in S:
            v=T[(n,l)]*(1 if l%2 else -1); assert v.denominator==1
            res.append((mm,n,l,str(Fraction(4*mm*n-N*l*l,N)),int(v)))
            if v<=0: bad.append(res[-1])
    del E
    print(f'block m in ({lo},{hi}] X={X} charges={len(sub)} time={time.time()-t:.1f}s bad so far={len(bad)}',flush=True)
print('N',N,'D0',D0,'charges',len(res),'nonpositive',len(bad),bad[:10],'total',time.time()-t0,flush=True)
json.dump(res,open(f'sweep_{N}_{D0}_{MLO}.json','w'))
