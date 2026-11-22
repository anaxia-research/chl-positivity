# W_N = alpha E4(N t) + gamma_N + beta delta_N equals the closed forms, and the Eisenstein + cusp decompositions.
import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from fractions import Fraction as Fr
exec(open('verify_closed.py').read().split("e4=E4c(X)")[0].replace("N=int(sys.argv[1]); X=int(sys.argv[2])","X=int(sys.argv[1])"))
from chl import EN_coeffs, sigma
for N in (2,3,5,7):
    e4=E4c(X); e4N=[e4[n//N] if n%N==0 else Fr(0) for n in range(X+1)]
    EN=EN_coeffs(N,X); EN2=[sum(EN[i]*EN[n-i] for i in range(n+1)) for n in range(X+1)]
    al={2:Fr(13,36),3:Fr(3,16),5:Fr(7,72),7:Fr(1,16)}[N]; be={2:Fr(5,9)}.get(N,2*al)
    ga={2:(Fr(4,45),Fr(19,180),0),3:(Fr(3,160),Fr(27,160),0),5:(Fr(-1,288),Fr(-19,288),Fr(1,6)),7:(Fr(-1,960),Fr(-29,960),Fr(3,32))}[N]
    de={2:(Fr(1,5),Fr(-6,5),0),3:(Fr(1,20),Fr(-21,20),0),5:(Fr(1,8),Fr(15,8),Fr(-3)),7:(Fr(1,30),Fr(7,15),Fr(-3,2))}[N]
    W=[al*e4N[n]+ga[0]*e4[n]+ga[1]*e4N[n]+ga[2]*EN2[n]+be*(de[0]*e4[n]+de[1]*e4N[n]+de[2]*EN2[n]) for n in range(X+1)]
    if N==2: Wc=[(a-b)/5 for a,b in zip(e4,e4N)]; dec=Wc
    if N==3: Wc=[Fr(3,80)*(a-b) for a,b in zip(e4,e4N)]; dec=Wc
    if N==5:
        Wc=[15*c for c in etaq({5:10,1:-2},X)]; g5=etaq({1:4,5:4},X)
        dec=[Fr(1,208)*(a-b)-Fr(15,13)*c for a,b,c in zip(e4,e4N,g5)]
    if N==7:
        th=[0]*(X+1); r=int((4*X)**0.5)+2
        for a in range(-r,r+1):
            for b in range(-r,r+1):
                v=a*a+a*b+2*b*b
                if v<=X: th[v]+=1
        e=etaq({7:7,1:-1},X); Wc=[3*sum(th[i]*e[n-i] for i in range(n+1)) for n in range(X+1)]
        e3=etaq({1:3,7:3},X); g7=[sum(th[i]*e3[n-i] for i in range(n+1)) for n in range(X+1)]
        dec=[Fr(1,800)*(a-b)-Fr(3,10)*c for a,b,c in zip(e4,e4N,g7)]
    print('N=%d through q^%d: W formula %s, decomposition %s, W>=0 %s, first coefficients %s'%(N,X,W==Wc,Wc==dec,all(w>=0 for w in W),[str(w) for w in W[:6]]))
