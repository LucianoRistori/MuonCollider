import numpy as np, math, time, gsurf as g
SH=1000/1486; s=0.1
def tower(N, wmin):
    Y1 = wmin*s/SH; Y = Y1*np.arange(1,N+1); W = SH*Y
    planes=[g.Plane([y,0,0],[1,0,0],[0,1,0],-w/2,w/2,-w/2,w/2,s) for y,w in zip(Y,W)]
    return planes, 0.2/Y[-1], Y, W
if __name__=="__main__":
  rng=np.random.default_rng(1)
  for N,Kref in ((4,5.8329e-3),):
      P,Y,W = None,None,None
      planes,kmax,Y,W = tower(N,2.0)
      box=[(-0.8,0.8),(-kmax,kmax),(-3*W[-1],3*W[-1]),(-1.5,1.5)]
      t0=time.time()
      Pc,err,edge = g.tube_P(planes, box, kmax, 1.0, 4_000_000, rng)
      print(N, "K_general", Pc, "+-", err, "edge frac", edge, "K_paper1", Kref, "ratio", Pc/Kref, f"{time.time()-t0:.0f}s")
