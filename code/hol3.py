import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
exec(open('fitN.py').read().split("basis=")[0])
P=mul(M2[0],B); Qs=[mul(f,A) for f in M4[:(2 if N==3 else 3)]]
R=lambda s,n,l: sp.Rational(str(coeff(s,n,l)))
eqs={3:[(0,1,0),(1,1,24)],5:[(0,1,0),(1,1,0),(2,1,72)],7:[(0,1,0),(1,1,0),(2,1,24)]}[N]
Mm=sp.Matrix([[R(q,n,l) for q in Qs] for n,l,_ in eqs]); r=sp.Matrix([v-R(P,n,l) for n,l,v in eqs])
cs=list(Mm.LUsolve(r)); print('c =',cs)
T=P
for c,q in zip(cs,Qs): T=add(T,scal(q,fmpq(int(c.p),int(c.q))))
mx=0
for n in range(0,X+1):
    vals=[coeff(T,n,l) for l in range(-12,13)]
    assert all(v>=0 and v.q==1 for v in vals),(n,vals)
    assert Fr(int(sum(vals).p),int(sum(vals).q))==12*EN[n],(n,sum(vals),12*EN[n])
print('Theta_N = E_N*phi01 + c.A : nonneg integer coeffs and row sums 12 E_N(n) through x^%d'%X)
for n in range(0,8): print(n,[str(coeff(T,n,l)//12) for l in range(0,4)])
