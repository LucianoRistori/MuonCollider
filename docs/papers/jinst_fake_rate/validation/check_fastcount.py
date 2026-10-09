import itertools, time, numpy as np
from validate_line import projective, cut_for
from fastcount import Plane, search
rng=np.random.default_rng(7)
geo=projective(6,100,600,7*0.1*6,0.1)
planes=[Plane(y,w,s) for y,w,s in zip(geo.Y,geo.W,geo.sigma)]
n=5; tot_b=tot_f=0; mism=0
for e in range(300):
    hu=[rng.uniform(-W/2,W/2,n) for W in geo.W]; hv=[rng.uniform(-W/2,W/2,n) for W in geo.W]
    for s in [tuple(range(6))]+list(itertools.combinations(range(6),5)):
        g=geo.subset(s); c=cut_for(g.ndof)
        ix=np.stack([a.ravel() for a in np.meshgrid(*[np.arange(n)]*len(s),indexing="ij")],1)
        u=np.stack([hu[i][ix[:,j]]/geo.sigma[i] for j,i in enumerate(s)],1); v=np.stack([hv[i][ix[:,j]]/geo.sigma[i] for j,i in enumerate(s)],1)
        bf={tuple(r) for r in ix[g.chi2(u,v)<c]}
        res,ch=search([planes[i] for i in s],[(hu[i],hv[i]) for i in s],c)
        fs={tuple(r) for r in res}
        tot_b+=len(bf); tot_f+=len(fs); mism+= (bf!=fs)
print("brute",tot_b,"fast",tot_f,"events*sets with any difference:",mism)
