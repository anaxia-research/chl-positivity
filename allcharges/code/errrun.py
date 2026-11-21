import sys,math,time
N=int(sys.argv[1]); sys.argv=[sys.argv[0],str(N)]+sys.argv[2:]
import numpy as np
from errterms import *
ms=[int(v) for v in sys.argv[2].split(',')]
X0=math.sqrt(3)/N
kp2,k,nu,d=setup(N)
for m in ms:
    t1=time.time(); x0=max(X0,math.sqrt(float(__import__("os").environ.get("D0","400")))/m)
    WS=walls_S(m); WT=walls_T(m,N)
    eS,eT=eps_gamma1(N,m,x0,WS,WT,kp2,nu,d)
    e2=eps_gamma2(N,m,x0,WS,WT,kp2,nu,d)
    eK=eps_K(N,m,x0,kp2,k,nu,d)
    print(N,m,'x0',round(x0,4),'rs',eS,'T',eT,'g2',e2,'K',eK,'sum',eS+eT+e2+eK,'%.0fs'%(time.time()-t1),flush=True)
