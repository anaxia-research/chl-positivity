import sys, itertools; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import numpy as np, sympy as sp
from quat import order
def reduced(N):
    G=np.array(order(N).tolist(),dtype=np.int64)
    lam=min(np.linalg.eigvalsh(G))/2; R=int((4/lam)**0.5)+1
    vs=[]
    for v in itertools.product(range(-R,R+1),repeat=4):
        v=np.array(v); n=v@G@v//2
        if 0<n<=4: vs.append((n,tuple(v)))
    vs.sort()
    basis=[]
    for n,v in vs:
        M=sp.Matrix(basis+[list(v)])
        if M.rank()==len(basis)+1: basis.append(list(v))
        if len(basis)==4: break
    U=sp.Matrix(basis)
    assert abs(U.det())==1, U.det()
    return U, sp.Matrix(G)
def theta_jac(N, wstar, nmax):
    U,G=reduced(N)
    Gr=np.array((U*G*U.T).tolist(),dtype=np.int64)
    wr=np.array((U*sp.Matrix(wstar)).T.tolist()[0],dtype=np.int64)  # pairing v_red . (U w*)
    lam=min(np.linalg.eigvalsh(Gr))/2; R=int((nmax/lam)**0.5)+1
    ax=np.arange(-R,R+1)
    V=np.stack(np.meshgrid(ax,ax,ax,ax,indexing='ij'),-1).reshape(-1,4)
    n=np.einsum('ij,jk,ik->i',V,Gr,V)//2; l=V@wr
    m=n<=nmax
    out={}
    for a,b in zip(n[m],l[m]): out[(int(a),int(b))]=out.get((int(a),int(b)),0)+1
    return out, Gr
if __name__=='__main__':
    for N in (3,5,7):
        U,G=reduced(N); print(N, (U*G*U.T).tolist())
