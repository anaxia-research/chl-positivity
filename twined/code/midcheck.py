# N=3, 400 < D <= 2500: exact d_* (Part I sweep) against the Cauchy bound for the twined index: need d_* > 2|d_g|
import json, math, sys
from fractions import Fraction
from boxcheck import Sopt
from tailG import U
N=3; kappa,eps,Rmax=1.3,0.15,2.08
S=json.load(open(sys.argv[1]))
worst=(1e9,None); cnt=0; bad=0
for m,n,l,D,d in S:
    D=Fraction(D)
    if not (400<D<=2500): continue
    if d<=0: bad+=1; continue
    r=Sopt(N,m,n,l,kappa,eps,Rmax)
    g=math.log(d)-math.log(2)-U(N,m,n,l,kappa,eps,Rmax,r)
    cnt+=1
    if g<worst[0]: worst=(g,(m,n,l,float(D)))
    if g<=0: bad+=1; print('FAIL',m,n,l,g)
print('charges',cnt,'bad',bad,'worst log-margin',worst)
