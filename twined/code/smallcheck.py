import json, sys, time
from exact_tw import Tw
N=int(sys.argv[1]); D0=int(sys.argv[2])
S=json.load(open('../../allcharges/code/sweep_%d_%d_0.json'%(N,D0)))
t0=time.time()
Lmax=max((m+1)//N for m,n,l,D,v in S); Bmax=max((n+1)//N for m,n,l,D,v in S)
print('Lmax',Lmax,'Bmax',Bmax,flush=True)
T=Tw(N,Lmax,Bmax); T.run(); print('built',time.time()-t0,flush=True)
worst=None; out=[]; nz=0
for m,n,l,D,v in S:
    g=T.d(m,n,l)
    out.append((m,n,l,D,v,g))
    if g: nz+=1
    r=abs(g)/v
    if worst is None or r>worst[0]: worst=(r,(m,n,l,D,v,g))
print('charges',len(S),'twined nonzero',nz,'max |d_g|/d_*',worst)
json.dump(out,open('../data/tw_%d_%d.json'%(N,D0),'w'))
