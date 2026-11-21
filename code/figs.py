import sys, math, json; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from chl import *
from fractions import Fraction as Fr
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'serif','font.size':10})
out=__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)),'..','paper','figures')+'/'
rows={}; tables={}
for N in (3,5,7):
    X=4*N*12//4+8
    X=int(N*40/4)+6
    C=CHL(N,X,1); F=finite_parts(C)
    pts=[]
    for n in range(1,X-2):
        for l in range(0,2):
            D=Fr(4*n,N)-l*l
            if D>0 and D<=40:
                v=dstar(F,N,1,n,l); pts.append((float(D),int(v.p)))
    pts=sorted(set(pts)); rows[N]=pts
    tables[N]=[(str(Fr(4*n,N)-l*l),n,l,int(dstar(F,N,1,n,l).p)) for n in range(1,7) for l in range(0,3) if Fr(4*n,N)-l*l>0]
json.dump({str(k):v for k,v in tables.items()},open(out+'../table_row1.json','w'))
fig,ax=plt.subplots(figsize=(5.2,3.4))
mk={3:('o','none'),5:('s','none'),7:('^','k')}
for N,pts in rows.items():
    ax.plot([p[0] for p in pts],[math.log10(p[1]) for p in pts],linestyle='none',marker=mk[N][0],mfc=mk[N][1],mec='k',ms=4,label=r'$N=%d$'%N)
ax.set_xlabel(r'$D=Q^2P^2-(Q\cdot P)^2$'); ax.set_ylabel(r'$\log_{10} d_*(Q,P)$')
ax.legend(frameon=False); ax.set_xlim(0,40.5); ax.grid(False)
for s in ax.spines.values(): s.set_color('k')
fig.tight_layout(); fig.savefig(out+'fig_row1.pdf'); fig.savefig(out+'fig_row1.png',dpi=200)
# figure 3: positive vs negative pieces for l odd, N=7 and N=5 and N=3
exec(open('W.py').read().split("for N,(a1")[0])
from row1_lib import hurwitz
fig,axs=plt.subplots(1,3,figsize=(7.0,2.6),sharey=False)
for ax,(N,b,fh,ab) in zip(axs,[(3,Fr(9,2),27,(Fr(3,80),Fr(-3,80),0)),(5,Fr(7,3),14,(Fr(1,48),Fr(19,48),Fr(-5,12))),(7,Fr(3,2),9,(Fr(1,320),Fr(29,320),Fr(-3,32)))]):
    U=30; e4=E4(U); e4N=dil(E4(U),N,U); EN=EN_coeffs(N,U); EN2=sq(EN,U)
    W=[ab[0]*e4[t]+ab[1]*e4N[t]+ab[2]*EN2[t] for t in range(U+1)]
    us=list(range(2,U+1))
    pos=[float(W[u]) for u in us]
    neg=[float(b*Fr(24,N-1)*sum(d for d in range(1,u+1) if u%d==0 and d%N)+2*fh*(hurwitz((4*u-N)//N) if (4*u-N)%N==0 else 0)) for u in us]
    ax.semilogy(us,pos,'k-',lw=1.2,label=r'$w_N(u)$')
    ax.semilogy(us,neg,'k--',lw=1.0,label=r'$b_N\, r_O(u)+2\hat{f}\, H$')
    ax.set_title(r'$N=%d$'%N,fontsize=10); ax.set_xlabel(r'$u$')
    for s in ax.spines.values(): s.set_color('k')
axs[0].legend(frameon=False,fontsize=8,loc='lower right')
fig.tight_layout(); fig.savefig(out+'fig_bounds.pdf'); fig.savefig(out+'fig_bounds.png',dpi=200)
print(tables)
