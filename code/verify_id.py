import sys,time; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
N=int(sys.argv[1]); X=int(sys.argv[2])
sys.argv=['x',str(N),str(X)]
t0=time.time()
exec(open('fitN.py').read().split("basis=")[0])
coef={3:(Fr(3,16),Fr(-3,8),[Fr(3,160),Fr(27,160),0]),5:(Fr(7,72),Fr(-7,36),[Fr(-1,288),Fr(-19,288),Fr(1,6)]),7:(Fr(1,16),Fr(-1,8),[Fr(-1,960),Fr(-29,960),Fr(3,32)])}[N]
al,be,ga=coef
q=lambda f: fmpq(f.numerator,f.denominator)
G=scal(mul(B,B),q(al))
G=add(G,scal(mul(mul(M2[0],A),B),q(be)))
AA=mul(A,A)
for c,f in zip(ga,M4[:3]):
    if c: G=add(G,scal(mul(f,AA),q(c)))
diff=add(G,E2,-1)
bad=[(a) for a in range(X+1) if diff.rows[a]!=0]
print('N=%d numerator identity E_2 = alpha phi01^2 + beta E_N phi01 phi-2 + gamma phi-2^2: mismatching x-orders through %d:'%(N,X),bad[:5],' time %.0fs'%(time.time()-t0))
