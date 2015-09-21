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
from scipy.special import ellipeinc as E
import SpiceyPy as sp

halfpi = sp.halfpi()

### Indices:  ECI (xyz); RSW (rsw); SEZ; R,RA,DEC (from sp.recrad(...))

iX, iY, iZ = iR, iSrsw, iW = iSsez, iE, iZsez = iRadius, iRA, iDEC = xrange(3)


########################################################################
class VA2014:
  def __init__(self,stateChf,mu,stateDep=None):
    """Make chief-only calculations.  Save deputy-based calculations for
Deputy() method below; call it now if stateDep is provided

mu = GM, km**3 s**-2

"""
    ####################################################################
    ### Extract ECI Positions and ECI velocities from 6-element ECI states
    ### Save mu in self object
    self.rChfEci,self.velChfEci = stateChf[:3], stateChf[3:]
    self.mu = mu

    ####################################################################
    ### Equations (1) and (2)
    ### Convert chief position and velocity to [Rhat|Shat|What]1 matrix
    self.mtxEci2Rsw1 = RVtoRSW(self.rChfEci,self.velChfEci)

    ####################################################################
    ### Transform chief position and velocity to RSW1 frame
    (self.rChfRsw1
    ,self.velChfRsw1
    ) = [ sp.mxv(self.mtxEci2Rsw1,v)
          for v in 
          [self.rChfEci,self.velChfEci]
        ]

    ####################################################################
    ### Equations (3) - postponed to deputy calculation

    ####################################################################
    ### Equations (4)
    ### Eccentricity and semi-major axis
    ### - vHChf = momentum vector
    ### - pChf = semiparameter or semilatus rectum = distance from the
    ###            primary focus to the orbit, measured perpendicular to
    ###            the semi-major axis
    ### - self.eccChf = Eccentricity (also Elliptic Modulus, k)
    ### - self.aChf = length of semi-major axis
    ### - self.bChf = length of semi-minor axis
    vHChf = sp.vcrss(self.rChfRsw1,self.velChfRsw1)
    self.pChf = sp.vdot(vHChf,vHChf) / self.mu
    vEccChf = sp.vsub(sp.vscl(1./self.mu,sp.vcrss(self.velChfRsw1,vHChf)),sp.vhat(self.rChfRsw1))
    self.eccChfSquared = sp.vdot(vEccChf,vEccChf)
    self.eccChf = math.sqrt(self.eccChfSquared)
    self.bOverA = math.sqrt(1 - self.eccChfSquared)
    self.aChf = self.pChf / (1. - self.eccChfSquared)
    self.bChf = self.pChf / self.bOverA
    if self.eccChfSquared>0.: self.vEccHatChf = sp.vhat(vEccChf)
    else                    : self.vEccHatChf = sp.vhat(self.rChfRsw1)

    ### Use SPICE toolkit routine to calculate perfocal distance (rp), eccentricity, and semi-major axis
    self.rpChf,self.eccChfSpice,inc,lnode,argp,m0,t0,muTmp = sp.oscelt(stateChf,0.,self.mu)
    self.aChfSpice = self.rpChf / (1. - self.eccChfSpice)

    ####################################################################
    ### Equations (5)
    ### Chief True anomaly posiiton
    ### Deputy True anomaly postponed to deputy calculation
    self.lambdaPerigee = sp.recrad(self.vEccHatChf)[1]
    self.trueAnom1 = (sp.twopi() - self.lambdaPerigee) % sp.twopi()

    ####################################################################
    ### Equations (6 through 9) - postponed to deputy calculation

    ####################################################################
    ### Equations (9) - chief only

    ### 9.1) Convert True Anomaly 1 (chief) to Eccentric Anomaly
    ### - https://en.wikipedia.org/wiki/Eccentric_anomaly#From_the_true_anomaly
    ### - Possible exceptions:  divBy0; domain error.
    self.eccAnom1 = math.atan2(self.bOverA * math.sin(self.trueAnom1)
                              , self.eccChf + math.cos(self.trueAnom1)
                              )

    ####################################################################
    ### Complete conversion to EQCM for Deputy if Deputy state is present

    if not (stateDep is None): self.Deputy(stateDep)


########################################################################
  def Deputy(self,stateDep):
    """Make depyty-based calculations.  Chief-only calculations were
performed in __init__() method above

"""
    ####################################################################
    ### Extract Deputy ECI Position and velocity from 6-element deputy
    ###   ECI state
    self.rDepEci,self.velDepEci = stateDep[:3], stateDep[3:]

    ####################################################################
    ### Equations (1) and (2)
    ### Convert deputy position and velocity to RSW using
    ###   [Rhat|Shat|What]1 matrix
    (self.rDepRsw1
    ,self.velDepRsw1
    ) = [ sp.mxv(self.mtxEci2Rsw1,v)
          for v in 
          [self.rDepEci,self.velDepEci]
        ]

    ####################################################################
    ### Equations (3)
    ### Get delta-lambda to deputy as RA from (Radius,RA,DEC)
    ###   returned by recrad.
    RrdDepRsw1 = sp.recrad(self.rDepRsw1)
    self.deltaLambdaDep = sp.recrad(self.rDepRsw1)[iRA]

    ####################################################################
    ### Equations (5)
    ### Deputy True anomaly
    self.trueAnom2 = (sp.twopi() + self.deltaLambdaDep - self.lambdaPerigee) % sp.twopi()

    ####################################################################
    ### Equations (6)
    ### Equivalent chief position at point 2, using PQW2
    rChf2 = self.pChf / (1. + (self.eccChf * math.cos(self.trueAnom2)))
    pHat = self.vEccHatChf
    qHat = sp.ucrss([0.,0.,1.],pHat)
    self.rChfPqw2 = sp.vscl(rChf2,sp.vadd(sp.vscl(math.cos(self.trueAnom2),pHat),sp.vscl(math.sin(self.trueAnom2),qHat)))
    self.vChfPqw2 = sp.vscl(-math.sqrt(self.mu/self.pChf)
                           ,sp.vadd(sp.vscl(math.sin(self.trueAnom2),pHat)
                                   ,sp.vscl(self.eccChf+math.cos(self.trueAnom2),qHat)))

    ####################################################################
    ### Equations (7)
    ### Convert from PQW2 to RSW2
    self.mtxPqw2toRsw2 = RVtoRSW(self.rChfPqw2,self.rChfPqw2)
    self.rChfRsw2 = sp.mxv(self.mtxPqw2toRsw2,self.rChfPqw2)
    self.velChfRsw2 = sp.mxv(self.mtxPqw2toRsw2,self.vChfPqw2)

    ####################################################################
    ### Equations (8)
    ### Transform Deputy vectors to SEZ frame
    ### - deltaPhiDep is DEC from [Radius,RA,DEC] of rDepRsw1
    ### - Matrix is ROT2[90-deltaPhiDep] ROT3[deltaLambdaDep]
    self.deltaPhiDep = RrdDepRsw1[iDEC]
    self.mtxRswToSez = sp.eul2m(halfpi-self.deltaPhiDep,self.deltaLambdaDep,0.,2,3,1)
    self.rDepSez = sp.mxv(self.mtxRswToSez,self.rDepRsw1)
    self.velDepSez = sp.mxv(self.mtxRswToSez,self.velDepRsw1)

    ####################################################################
    ### Equations (9)
    ### Transform Deputy vectors from SEZ to EQCM frame

    ### 9.1) Convert True Anomalies 2 (Deputy) to Eccentric Anomaly 2
    ### - https://en.wikipedia.org/wiki/Eccentric_anomaly#From_the_true_anomaly
    ### - Possible exceptions:  divBy0; domain error.
    ### - Chief done in __init__() method above
    self.eccAnom2 = math.atan2(self.bOverA * math.sin(self.trueAnom2), self.eccChf + math.cos(self.trueAnom2))

    ### 9.2) Get arc length between those Eccentric Anomalies
    arcLength1to2 = self.bChf * orbit_IEISK(self.eccChfSquared, self.eccAnom1, self.eccAnom2)

    ### 9.3) Relate the deputy relative to chief at point 2
    self.rDepEqcm = sp.vpack( self.rDepSez[iZsez] - self.rChfRsw2[iR]
                            , arcLength1to2
                            , self.deltaPhiDep * rChf2
                            )
    self.velDepEqcm = sp.vpack( self.velDepSez[iZsez] - self.velChfRsw2[iR]
                              , (self.velDepSez[iE] * rChf2 / (RrdDepRsw1[iRadius] * self.deltaPhiDep)) - self.velChfRsw2[iE]
                              , -self.velDepSez[iSsez] * rChf2 / RrdDepRsw1[iRadius]
                              )


########################################################################
def RVtoRSW(R,V):
    """Convert R and V (Postion and Velocity) into RWS matrix"""

    ####################################################################
    ### Equations (1)
    ### Convert chief position and velocity to RSW1 unit vectors
    rHatEci = sp.vhat(R)
    wHatEci = sp.ucrss(R,V)
    ### W cross R, not R cross W which is in paper
    sHatEci = sp.ucrss(wHatEci,R)

    ####################################################################
    ### Equations (2)
    ### Combine RSW1 unit vector into [Rhat|Shat|What] matrix
    return [rHatEci,sHatEci,wHatEci]


########################################################################
def orbit_IEISK(eccSquared,ea0,ea1):
  """Calculate the arc length from one position on an elliptical orbit
to a second position on that same orbit, using E, the Incomplete
Elliptic Integral of the Second Kind (IEISK).

The ellipse semi-major axis is 1 (= a); the semi-minor axis length is b,
and b <= a.

Arguments:

  eccSquared:  the ellipse eccentricity squared = (1 - (b/a)^2)
               N.B. b <= a; 0 <= eccSquared <= 1
               N.B. eccSquared equals ecc2 in the notes below
               N.B. eccSquared is also known as the paramete", m = k^2
               N.B. sqrt(eccSquared) is also known as the elliptic
                    modulus or eccentricity, e = k = sin(alpha), where
                    alpha is the modular angle

  ea0 = the first position on the orbit, expressed as an eccentric
        anomaly (angle, radians) positive from zero at the semi-major
        axis toward the direction of the orbiting body

  ea1 = the second position on the orbit, expressed as an eccentric
        anomaly (angle, radians) positive from zero at the semi-major
        axis toward the direction of the orbiting body

Implementation:

1) The arc length from the semi-major axis is definite integral of
sqrt(1 - ecc2 cos(ea)^2) evaluated from eccentric anomaly ea0 to ea1.

2) However, the integral available, E, the Incomplete Elliptic Integral
of the Second Kind, is the integral of sqrt(1 - ecc2 sin(ea)^2), from
eccentric anomaly zero to ea, which evalutates angles relative to a
semi-minor axis a of length < 1, with semi-major axis b of length 1,
and with ecc2 = (1 - (a/b)^2).

2.1) So E is off by PI/2 relative to the desired integral in (1) above.

3) The solution is to use E to calculate the arc length from eccentric
anomaly (PI/2 - ea) to eccentric anomaly (PI/2 - 0).  This is a variable
substitution that effectively replaces cos(ea) in (1) with its
equivalent sin(PI/2 - ea), and converts the integral in (1) to the form
in (2), at which point E can be used.

3.1) So

  arcLength0 = integral(1-ecc2 cos(ea)^2) from ea=0 to ea0
             = integral(1-ecc2 sin(PI/2-ea)^2) from ea=PI/2-ea0 to PI/2
             = E(ecc2,PI/2) - E(ecc2,PI/2-ea0)

3.2) and similarly

  arcLength1 = E(ecc2,PI/2) - E(ecc2,PI/2-ea1)

3.3) and finally the desired arc length from ea0 to ea1 in (1) above

  = arcLength1 - arcLength0
  = (E(ecc2,PI/2)-E(ecc2,PI/2-ea1)) - (E(ecc2,PI/2)-E(ecc2,PI/2-ea0)) 
  = E(ecc2,PI/2-ea0) - E(ecc2,PI/2-ea1)

3.4) Note the ea0 and ea1 terms reverse because d(PI/2-ea)/d(ea) = -1

  """

  ### With all that, this function is reduced to a one-liner
  ### scipy.special.ellipeinc provides E
  ### - Handles wraparound of eccentric anomaly
  ### - Routine uses the parameter, m = eccentricity^2

  return E(eccSquared,halfpi-ea0) - E(eccSquared,halfpi-ea1)


########################################################################
if "__main__"==__name__:

  arg12 = [s for s in sys.argv[1:] if s[:5]!='--mu=']

  muArg = float(([398600.4418] + [s[5:] for s in sys.argv[1:] if s[:5]=='--mu='])[-1])

  if len(sys.argv[1:])!=12:
    arg12 = '1e4 1e4 1e4  .34e1 -.34e1 0   1.1e4 1.2e4 1.3e4  .29e1 -.29e1 0'.split()
    while arg12: sys.argv.insert(1,arg12.pop())

  stateChief = numpy.array(map(float,sys.argv[1:7]))
  stateDeputy = numpy.array(map(float,sys.argv[7:13]))

  print(dict(stateChf=stateChief,stateDep=stateDeputy))

  va = VA2014(stateChief,muArg,stateDep=stateDeputy)
  import pprint
  pprint.pprint(vars(va))
