# Rigorous constant A_N in Lemma A:
#   A_N >= A* = sum_{(n,l): Delta>=0} |c(n,l)| exp(-2 pi y0 [Delta/(4N) + N g_l^2]),  y0 = sqrt3/2,  g_l = (|l|/(2N)-1/2)_+ .
# Exact |c| for n <= n0; for n > n0 the a priori bound |c(Delta)| <= (A* + 4 e^{pi/N}) e^{pi Delta/(4N)} (Cauchy at height 1/2)
# gives A* <= (A_{<=n0} + 4 e^{pi/N} eps)/(1-eps), eps = sum_{n>n0,l} exp(-(2 y0 - 1) pi Delta/(4N) - 2 pi y0 N g_l^2).
import pickle, math
y0=math.sqrt(3)/2
for N,n0 in ((2,400),(3,300)):
    c=pickle.load(open('cD_%d_%d.pkl'%(N,n0),'rb'))
    A0=0.0
    for n in range(0,n0+1):
        lm=math.isqrt(4*N*n)
        for l in range(-lm,lm+1):
            D=4*N*n-l*l; g=max(abs(l)/(2*N)-0.5,0.0)
            A0+=abs(c[D])*math.exp(-2*math.pi*y0*(D/(4*N)+N*g*g))
    eps=0.0
    for n in range(n0+1,n0+20001):
        lm=math.isqrt(4*N*n)
        for l in range(-lm,lm+1):
            D=4*N*n-l*l; g=max(abs(l)/(2*N)-0.5,0.0)
            e=(2*y0-1)*math.pi*D/(4*N)+2*math.pi*y0*N*g*g
            if e<745: eps+=math.exp(-e)
    # n > n0+20000: every term has exponent >= min over l of ... >= 2 pi y0 N (sqrt(n/N)-1/2)^2 *0 + ... ; bounded crudely
    Astar=(A0+4*math.exp(math.pi/N)*eps)/(1-eps)
    print('N',N,'A_{<=n0}=%.6f'%A0,'eps=%.3e'%eps,'A* <= %.6f'%Astar)
