import sys, itertools; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import sympy as sp, numpy as np
from lattice import theta_rows
R=sp.Rational
def order(N):
    if N==3: a,b=-1,-3; basis=[(1,0,0,0),(R(1,2),0,R(1,2),0),(0,1,0,0),(0,R(1,2),0,R(1,2))]
    if N==7: a,b=-1,-7; basis=[(1,0,0,0),(R(1,2),0,R(1,2),0),(0,1,0,0),(0,R(1,2),0,R(1,2))]
    if N==5: a,b=-2,-5; basis=[(1,0,0,0),(R(1,2),0,R(1,2),R(1,2)),(0,R(1,4),R(1,2),R(1,4)),(0,0,0,1)]
    if N==2: a,b=-1,-1; basis=[(1,0,0,0),(0,1,0,0),(0,0,1,0),(R(1,2),R(1,2),R(1,2),R(1,2))]
    wts=[1,-a,-b,a*b]
    G=sp.Matrix(4,4,lambda r,c: 2*sum(wts[t]*basis[r][t]*basis[c][t] for t in range(4)))
    return G
for N in (2,3,5,7):
    G=order(N); print(N,'Gram',G.tolist(),'det',G.det(), 'integral even:', all(x.q==1 for x in G) and all(G[i,i]%2==0 for i in range(4)))
