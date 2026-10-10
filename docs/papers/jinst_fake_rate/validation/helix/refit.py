import json,sys,numpy as np
from scipy.optimize import minimize
d=json.load(open(sys.argv[1])); nd=d['ndof']; n=d['trials']
g=np.array([0]+d['grid']); C=np.array([0]+d['counts'])
obs=np.diff(C)
for cmax in (0.5,1.0,2.0,4.0):
  m=g<=cmax*1.0001; e=g[m]; o=obs[:m.sum()-1]
  for order in (1,2,3):
    def nll(p):
      K=np.exp(p[0]); poly=1+sum(p[i]*e**i for i in range(1,order+1))
      mu=np.diff(n*K*e**(nd/2)*poly)
      return 1e30 if np.any(mu<=0) else float(np.sum(mu-o*np.log(mu)))
    x=[np.log(d['K'])]+[0]*order
    for _ in range(4): x=minimize(nll,x,method='Nelder-Mead',options=dict(xatol=1e-10,fatol=1e-10,maxiter=50000,maxfev=50000)).x
    print(f"cmax {cmax} order {order}: K/K0={np.exp(x[0])/d['K']:.4f} coeffs {np.round(x[1:],4)} nll {nll(x):.1f}")
