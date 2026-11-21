import cmath, math
def P(z,N,e,T=80):
    v=1
    for n in range(1,T):
        v*= (1-z**n)**(-e)*(1-z**(N*n))**(-e)
    return v
def inF(w,N):
    if abs(w.real)>N/2+1e-12 or abs(w)**2<N-1e-12: return False
    for b in range(1,8):
        for a in range(-40,41):
            if a%N and math.gcd(a,b)==1 and abs(a+b*w)<1-1e-12: return False
    return True
if __name__=="__main__":
    for N,e in ((3,6),(5,4),(7,3)):
      best=(9,None)
      for i in range(0,401):
          re=-N/2+N*i/400
          # lowest admissible Im w along this vertical line
          lo,hi=0.0,3.0
          for _ in range(50):
              mid=(lo+hi)/2
              if inF(complex(re,mid),N): hi=mid
              else: lo=mid
          w=complex(re,hi)
          v=abs(P(cmath.exp(2j*math.pi*w/N),N,e))**2
          if v<best[0]: best=(v,w)
      print(N,best)
