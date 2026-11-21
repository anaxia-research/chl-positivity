from fractions import Fraction as Fr
def hurwitz(D):
    if D==0: return Fr(-1,12)
    if D<0 or D%4 not in (0,3): return Fr(0)
    h=Fr(0); b=D%2
    while 3*b*b<=D:
        ac=(b*b+D)//4; a=max(b,1)
        while a*a<=ac:
            if ac%a==0:
                c=ac//a
                if a==b==c: h+=Fr(1,3)
                elif b==0 and a==c: h+=Fr(1,2)
                elif a==b or a==c or b==0: h+=1
                else: h+=2
            a+=1
        b+=2
    return h
