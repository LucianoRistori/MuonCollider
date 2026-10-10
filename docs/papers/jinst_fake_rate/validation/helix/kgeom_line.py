import math, itertools, json, numpy as np
def K_line(Y, W, s, nq=200001):
    N=len(Y); d=2*N-4
    qmax=1.5*(W/Y).max(); q=np.linspace(-qmax,qmax,nq)
    hi=np.min(W/2-q[:,None]*Y,1); lo=np.max(-W/2-q[:,None]*Y,1)
    area=np.trapezoid(np.maximum(hi-lo,0),q)
    detD=math.sqrt(N*(Y**2).sum()-Y.sum()**2)/s**2
    Vd=math.pi**(d/2)/math.gamma(d/2+1)
    return Vd*(area*detD)**2/np.prod((W/s)**2)
s=0.1
cal=json.load(open('../calib_sets.json')); ref=json.load(open('../ref_corrections.json'))
for N in (4,6):
    Y=np.linspace(100,100*N,N); W=2*s*N*Y/(100*N)
    print(N,"geom",K_line(Y,W,s),"ref",ref[str(N)]['K'])
Y=np.linspace(100,600,6); W=2*s*6*Y/600
for S in [tuple(range(6))]+list(itertools.combinations(range(6),5)):
    S=list(S); c=cal[",".join(map(str,S))]
    kg=K_line(Y[S],W[S],s); print(S, f"geom {kg:.5g} calib {c['K']:.5g}+-{c['K_err']:.2g} ratio {c['K']/kg:.4f}")
