import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
N=int(sys.argv[1]); X=int(sys.argv[2])
sys.argv=['x',str(N),str(X)]
exec(open('hol3.py').read().split("mx=0")[0])
from lat2 import theta_jac
ws={3:(-1,-1,0,0),5:(0,-1,-1,-1),7:(0,-1,0,0)}[N]
th,Gr=theta_jac(N,ws,X)
ok=all(12*th.get((n,l),0)==int(coeff(T,n,l)) for n in range(X+1) for l in range(-12,13))
print('N=%d reduced Gram %s : E_N phi01 + c.phi-2 == 12 theta(O,w) through x^%d : %s'%(N,Gr.tolist(),X,ok))
