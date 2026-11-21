# closed-form bound for all remaining terms when m >= M (decreasing in m), uniform in x >= x0
import math
for N in (3,5,7):
    kp2=24//(N+1); nu=kp2-2+1.5
    c=math.pi*math.sqrt(3)/N; kap=kp2*(1+1/N)*math.pi**2/6
    def logDb(X):   # log(1+ sqrt(X/kap) e^kap e^{2 sqrt(kap X)})  (X>=4kap)
        return math.log(math.sqrt(X/kap))+kap+2*math.sqrt(kap*X)+math.log(2)
    def logU(m):
        R=2*m*math.sqrt(N)+2; Gam=c*(m+3*math.sqrt(m))+1
        # every wall term of gamma>=2 or rs>=2: exp(-c(m+2-2/m-(m/2+3 sqrt m))) at worst (g/gamma<=1/2)
        lw=-c*(m+2-2/m-(m/2+3*math.sqrt(m)))
        lcount=math.log(2*R*R)+math.log(N*math.sqrt(2)*Gam**2.5+1)+math.log(math.sqrt(2)*N**(kp2/2)*N)
        lwall=lw+lcount+2*logDb(m)
        # shadow terms: <= (number of (gamma,delta) up to Gam)*poly * d(m) * exp(-c(m+2-2/m-2 sqrt m)), poly=(2m)^3 (4 pi^2 m^2 Gam^2/3)
        # shadow terms of the corrected expansion: per pair (gamma,delta) the bracket is at most
        # (2pi)^(k+3) [pi^2(2m g+2g(1+log g))+2pi+8 pi m g (x0/(2g)+1)] e^{y}, y<=2c sqrt m; summed over <=N g pairs, g<=Gam,
        # and over g>Gam with decay g^{-(k+2)}: total <= 2 x (sum up to Gam)
        k=kp2-2; Gm=Gam
        poly=(2*math.pi)**(k+3)*(math.pi**2*(2*m*Gm+2*Gm*(1+math.log(Gm)))+2*math.pi+8*math.pi*m*Gm*(math.sqrt(3)/N/2+1))
        lK=math.log(2*N*Gm*Gm)+math.log(poly)+logDb(m)-c*(m+2-2/m-2*math.sqrt(m))
        return max(lwall,lK)+math.log(2)
    M=1500
    print(N,'logU(1500)=%.1f'%logU(M), 'slope', (logU(M+1)-logU(M)), 'slope at 1e4', logU(10001)-logU(10000))
