"""
  SpiceyPy-based class implementing Vallado and Alfano,
  "Curvilinear coordinate transformations for relative motion,"
  Celest Mech Dyn Astr (2014) 118:253-271, DOI 10.1007/s10569-014-9531-1

  Test usage:

    python va2014.py [--mu=398600.4418] [ChiefStateECI DeputyStateECI]

    - State in ECI is six-element vector:  X; Y; Z; Xdot; Ydot; Zdot.

"""
import os
import sys
import math
import numpy
import SpiceyPy as sp

class VA2014:
  def __init__(self,stateChf,StateDep,mu):
    """
    mu = GM, km**3 s**-2
"""

    ### Extract ECI Positions and ECI velocities from 6-element ECI states
    self.rChfEci,self.velChfEci = stateChf[:3], stateChf[3:]
    self.rDepEci,self.velDepEci = StateDep[:3], StateDep[3:]

    ### Convert chief position and velocity to RSW1 unit vectors
    self.rHatEci = sp.vhat(self.rChfEci)
    self.wHatEci = sp.ucrss(self.rChfEci,self.velChfEci)
    ### W cross R, not R cross W which is in paper
    self.sHatEci = sp.ucrss(self.wHatEci,self.rHatEci)

    ### Combine RSW1 unit vector into [Rhat|Shat|What]1 matrix
    self.mtxEci2Rsw1 = [self.rHatEci,self.sHatEci,self.wHatEci]

    ### Transform chief position and velocity to RSW1 frame
    self.rChfRsw1,self.velChfRsw1,self.rDepRsw1,self.velDepRsw1 = [
       sp.mxv(self.mtxEci2Rsw1,v)
       for v in 
       [self.rChfEci,self.velChfEci,self.rDepEci,self.velDepEci]
       ]

    ### Get delta-lambda to deputy as RA from (Radius,RA,DEC)
    ###   returned by recrad.
    self.deltaLambdaDep = sp.recrad(self.rDepRsw1)[1]

    ### Semi-major axis and eccentricity
    vHChf = sp.vcrss(self.rChfRsw1,self.velChfRsw1)
    pChf = sp.vdot(vHChf,vHChf) / mu
    vEccChf = sp.vsub(sp.vscl(1./mu,sp.vcrss(self.velChfRsw1,vHChf)),sp.vhat(self.rChfRsw1))
    eccChfSquared = sp.vdot(vEccChf,vEccChf)
    self.eccChf = math.sqrt(eccChfSquared)
    self.aChf = pChf / (1. - eccChfSquared)
    if eccChfSquared>0: self.vEccHatChf = sp.vhat(vEccChf)
    else              : self.vEccHatChf = sp.vhat(self.rChfRsw1)

    ### Use SPICE toolkit routine to calculate perfocal distance (rp), eccentricity, and semi-major axis
    self.rpChf,self.eccChfSpice,inc,lnode,argp,m0,t0,mu = sp.oscelt(stateChf,0.,mu)
    self.aChfSpice = self.rpChf / (1 - self.eccChfSpice)

########################################################################
if "__main__"==__name__:

  arg12 = [s for s in sys.argv[1:] if s[:5]!='--mu=']

  muArg = float(([398600.4418] + [s[5:] for s in sys.argv[1:] if s[:5]=='--mu='])[-1])

  if len(sys.argv[1:])!=12:
    arg12 = '1e4 1e4 1e4  .34e1 -.34e1 0   1.1e4 1.2e4 1.3e4  .29e1 -.29e1 0'.split()
    while arg12: sys.argv.insert(1,arg12.pop())

  stateChief = numpy.array(map(float,sys.argv[1:7]))
  stateDeputy = numpy.array(map(float,sys.argv[7:13]))

  print(dict(stateChief=stateChief,stateDeputy=stateDeputy))

  va = VA2014(stateChief,stateDeputy,muArg)
  import pprint
  pprint.pprint(vars(va))
