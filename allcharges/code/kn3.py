import sys; N=int(sys.argv[1]); sys.argv=[sys.argv[0],str(N)]
import knew,math
from exact import ginv
D0={3:2500,5:400,7:400}[N]; mU={3:123,5:192,7:260}[N]
d=[int(v) for v in ginv(N,1500,24//(N+1))]
rows=[l.split() for l in open('../errtab_%d.txt'%N)]
best=(0,0)
for m,x,Ko,S in rows:
    m=int(m); 
    if m>=mU: continue
    x=max(math.sqrt(3)/N,math.sqrt(D0)/m)
    new=float(S)-float(Ko)+knew.epsK(N,m,x,d)
    best=max(best,(new,m))
print(N,'per-charge max new',best)
w=(0,0)
for m in range(mU,1461):
    w=max(w,(knew.epsK(N,m,math.sqrt(3)/N,d,G=40),m))
print(N,'uniform K new max',w)
