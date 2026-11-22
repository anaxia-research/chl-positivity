# Constants of the closed-form estimate (Lemmas closed and theta, Proposition tail) for D > 1e8.
import math
CM={2:0.410,3:0.202}; KW={2:6,3:4}; RMAX={2:1.23,3:2.08}; C0={2:8,3:4}
for N in (2,3):
    eps=0.15; P=1.3/(4*N*N)/(1-eps*eps*N); s=math.sqrt(P)
    cp=(C0[N]+4)*math.pi/12; al=cp/(2*math.pi*P)
    K=math.sqrt(2+al)/(4*math.pi*eps*s); b=6*math.log(1+K)+RMAX[N]
    th2a=16*P/(N*(4/N-1))*(101+al)/104*(101+al)/100
    th2b=4*P/(1-1e-4)*(2+al)/5*1.00002
    th=math.sqrt(max(th2a,th2b))
    k=KW[N]; nu=k+1.5
    g=math.log(0.44*CM[N])+math.log(4*math.pi/N**((k+4)/2))-0.5*math.log(2*math.pi**2*math.sqrt(5))+nu/2*math.log(5)
    D=1e8
    thN={2:0.865,3:0.847}[N]; assert th<=thN
    first=(1-thN)*math.pi*math.sqrt(D)-(nu/2+13/4)*math.log(D+1)-b+g-math.log(N-1)
    second=(0.9-thN)*math.pi*math.sqrt(D)-(3+0.9*(nu/2+0.25))*math.log(D+1)-b+0.9*g
    print('N',N,'P %.7f s %.7f alpha %.4f K %.4f b %.4f'%(P,s,al,K,b))
    print('  Theta^2: m>=100 %.4f, m<=99 %.4f -> Theta %.4f <= theta_N %.3f'%(th2a,th2b,th,thN))
    print('  gamma %.4f; at D=1e8: positivity margin %.1f, 0.9-power margin %.1f (both increase with D)'%(g,first,second))
    print('  Theorem C constants: 4 pi s/sqrt N = %.5f, 1+alpha = %.4f, b = %.4f'%(4*math.pi*s/math.sqrt(N),1+al,b))
