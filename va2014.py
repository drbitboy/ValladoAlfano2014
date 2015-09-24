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
twopi = sp.twopi()
dpr = sp.dpr()

### Indices:  ECI (xyz); RSW (rsw); SEZ; R,RA,DEC (from sp.recrad(...))

iX, iY, iZ = \
iRrsw, iSrsw, iW = \
iSsez, iE, iZsez = \
iRadius, iRA, iDEC = \
iReq, iAeq, iZeq = \
xrange(3)


########################################################################
class VA2014:


  ######################################################################
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
    self.mtxEciToRsw1 = RVtoRSW(self.rChfEci,self.velChfEci)

    ####################################################################
    ### Transform chief position and velocity to RSW1 frame
    (self.rChfRsw1
    ,self.velChfRsw1
    ) = [ sp.mxv(self.mtxEciToRsw1,v)
          for v in 
          [self.rChfEci,self.velChfEci]
        ]

    ####################################################################
    ### Equations (3) - postponed to deputy calculation

    ####################################################################
    ### Equations (4)
    ### Eccentricity and semi-major axis
    ### - vHChf = momentum vector
    ### - self.pChf = semiparameter or semilatus rectum = distance from
    ###               the primary focus to the orbit, measured
    ###               perpendicular to the semi-major axis
    ### - self.eccChfSquared = Square of Eccentricity
    ###   - derivation, dot product of vEccChf with itself, will always
    ###       be non-negative
    ###   - square root of (1-self.eccChfSquared) will throw math domain
    ###       exception if eccentricity exceeds unity.
    ###   - reciprocal of (1-self.eccChfSquared) will throw division by
    ###       zero exception if eccentricity is unity
    ### - self.eccChf = Eccentricity (also Elliptic Modulus, k)
    ### - self.bOverA = ratio of semi-minor to semi-major axes' lengths
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
    self.cChf = self.aChf * self.eccChf
    if self.eccChfSquared>0.: self.vEccHatChf = sp.vhat(vEccChf)
    else                    : self.vEccHatChf = sp.vhat(self.rChfRsw1)

    ### Use SPICE toolkit routine to calculate perfocal distance (rp), eccentricity, and semi-major axis
    self.rpChf,self.eccChfSpice,inc,lnode,argp,m0,t0,muTmp = sp.oscelt(stateChf,0.,self.mu)
    self.aChfSpice = self.rpChf / (1. - self.eccChfSpice)

    ### Quadrant arc length; used later
    self.quadArc = self.aChf * E(halfpi,self.eccChfSquared)

    ####################################################################
    ### Equations (5)
    ### Chief True anomaly posiiton
    ### Deputy True anomaly postponed to deputy calculation
    self.lambdaPerigee = sp.recrad(self.vEccHatChf)[iRA]
    self.trueAnom1 = (twopi - self.lambdaPerigee) % twopi

    ####################################################################
    ### Equations (6 through 9) - postponed to deputy calculation

    ####################################################################
    ### Equations (9) - chief only

    ### 9.1) Convert True Anomaly 1 (chief) to Eccentric Anomaly
    self.eccAnom1 = self.TrueToEccAnom(self.trueAnom1)

    ### Get the arc length from periapse to the chief
    self.arcChf = self.E(self.eccAnom1)

    ####################################################################
    ### Complete conversion to EQCM for Deputy if Deputy state is present

    if not (stateDep is None): self.Deputy(stateDep)


  ######################################################################
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
    ) = [ sp.mxv(self.mtxEciToRsw1,v)
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
    self.trueAnom2 = (twopi + self.deltaLambdaDep - self.lambdaPerigee) % twopi

    ####################################################################
    ### Equations (6)
    ### Equivalent chief position at point 2, using PQW2
    rChf2 = self.pChf / (1. + (self.eccChf * math.cos(self.trueAnom2)))
    pHat = self.vEccHatChf
    qHat = sp.ucrss([0.,0.,1.],pHat)
    self.rChfPqw2 = sp.vscl(rChf2,sp.vadd(sp.vscl(math.cos(self.trueAnom2),pHat),sp.vscl(math.sin(self.trueAnom2),qHat)))
    self.velChfPqw2 = sp.vscl(-math.sqrt(self.mu/self.pChf)
                             ,sp.vadd(sp.vscl(math.sin(self.trueAnom2),pHat)
                                     ,sp.vscl(self.eccChf+math.cos(self.trueAnom2),qHat)
                                     )
                             )

    ####################################################################
    ### Equations (7)
    ### Convert from PQW2 to RSW2
    self.mtxPqw2toRsw2 = RVtoRSW(self.rChfPqw2,self.velChfPqw2)
    self.rChfRsw2 = sp.mxv(self.mtxPqw2toRsw2,self.rChfPqw2)
    self.velChfRsw2 = sp.mxv(self.mtxPqw2toRsw2,self.velChfPqw2)

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
    self.eccAnom2 = self.TrueToEccAnom(self.trueAnom2)

    ### Ensure .eccAnom2 is within PI/2 of .eccAnom1
    while self.eccAnom2 >  (self.eccAnom1+math.pi): self.eccAnom2 -= (2 * math.pi)
    while self.eccAnom2 <= (self.eccAnom1-math.pi): self.eccAnom2 += (2 * math.pi)

    ### 9.2) Get arc length from positive semi-major axis to deputy
    self.arcDep = self.E(self.eccAnom2)

    ### 9.3) Relate the deputy relative to chief at point 2
    self.rDepEqcm = sp.vpack( self.rDepSez[iZsez] - self.rChfRsw2[iRrsw]
                            , self.arcDep - self.arcChf
                            , self.deltaPhiDep * rChf2
                            )
    self.velDepEqcm = sp.vpack( self.velDepSez[iZsez] - self.velChfRsw2[iRrsw]
                             #, (self.velDepSez[iE] * rChf2 / (RrdDepRsw1[iRadius] *          self.deltaPhiDep )) - self.velChfRsw1[iE]
                              , (self.velDepSez[iE] * rChf2 / (RrdDepRsw1[iRadius] * math.cos(self.deltaPhiDep))) - self.velChfRsw1[iE]
                              , -self.velDepSez[iSsez] * rChf2 / RrdDepRsw1[iRadius]
                              )

    return self.rDepEqcm, self.velDepEqcm


  ######################################################################
  def inverseDeputy(self,stateDepEqcm=False):
    """Invert depyty-based calculations.  Chief-only calculations were
performed in __init__() method above
"""

    ####################################################################
    ### Process argmument ...
    if stateDepEqcm:
      ### Use deputy state argument (6-vector), if supplied
      rDepEqcm = stateDepEqcm[:3]
      velDepEqcm = stateDepEqcm[3:]
    else:
      ### Use result of .Deputy() method if no state is supplied
      rDepEqcm = self.rDepEqcm
      velDepEqcm = self.velDepEqcm

    ####################################################################
    ### Equations (11 and 12) are already done in Equations 1 through 5
    ### in the constructor __init__() above
    ### - .mtxEciToRsw1
    ### - .{r,vel}ChfRsw1
    ### - .lambdaPerigee
    ### - .trueAnom1

    ####################################################################
    ### Equations (13)
    ### in the constructor __init__() above
    ### - .mtxEciToRsw1

    self.zinvArcDep = self.arcChf + rDepEqcm[iAeq]
    self.zinvEccAnom2 = self.invE(arcLengthArg=self.zinvArcDep)

    self.zinvTrueAnom2 = self.XYtoTrueAnom(self.aChf * math.cos(self.zinvEccAnom2)
                                          ,self.bChf * math.sin(self.zinvEccAnom2)
                                          )

    ####################################################################
    ### Equation (14)

    self.zinvDeltaLambdaDep = self.zinvTrueAnom2 - self.trueAnom1

    ####################################################################
    ### Equations (6) repeat from Deputy method, as 2 may be different
    ### Equivalent chief position at point 2, using PQW2
    rChf2 = self.pChf / (1. + (self.eccChf * math.cos(self.zinvTrueAnom2)))
    pHat = self.vEccHatChf
    qHat = sp.ucrss([0.,0.,1.],pHat)
    self.zinvRChfPqw2 = sp.vscl(rChf2,sp.vadd(sp.vscl(math.cos(self.zinvTrueAnom2),pHat),sp.vscl(math.sin(self.zinvTrueAnom2),qHat)))
    self.zinvVelChfPqw2 = sp.vscl(-math.sqrt(self.mu/self.pChf)
                                 ,sp.vadd(sp.vscl(math.sin(self.zinvTrueAnom2),pHat)
                                         ,sp.vscl(self.eccChf+math.cos(self.zinvTrueAnom2),qHat)
                                         )
                                 )

    ####################################################################
    ### Equations (7) repeat from Deputy method, as 2 may be different
    ### Convert from PQW2 to RSW2
    self.zinvMtxPqw2toRsw2 = RVtoRSW(self.zinvRChfPqw2,self.zinvVelChfPqw2)
    self.zinvRChfRsw2 = sp.mxv(self.zinvMtxPqw2toRsw2,self.zinvRChfPqw2)
    self.zinvVelChfRsw2 = sp.mxv(self.zinvMtxPqw2toRsw2,self.zinvVelChfPqw2)

    ####################################################################
    ### Equations (15)
    self.zinvDeltaPhiDep = rDepEqcm[iZeq] / rChf2

    rDepHatRsw1 = sp.radrec(1., self.zinvDeltaLambdaDep, self.zinvDeltaPhiDep)

    ####################################################################
    ### - Equations (8) repeat from Deputy method, as deputy may differ
    ###   - Transformation matrix from RSW to SEZ frame
    ###   - Matrix is ROT2[90-deltaPhiDep] ROT3[deltaLambdaDep]
    self.zinvMtxRswToSez = sp.eul2m(halfpi-self.zinvDeltaPhiDep,self.zinvDeltaLambdaDep,0.,2,3,1)

    ####################################################################
    ### Equations (16)
    ### - Transform deputy unit vector from RSW1 to SEZ
    ### - Scale deputy unit vectors in RSW1 and in SEZ
    ###   - factor is X (R) component of deputy EQCM X (R), plus R
    ###       component of chief RSW2
    ### - Transform deputy vector from RSW1 to ECI (transpose matrix)
    rDepHatSez = sp.mxv(self.zinvMtxRswToSez,rDepHatRsw1)
    sezScale = (rDepEqcm[iReq] + self.zinvRChfRsw2[iRrsw]) / rDepHatSez[iZsez]
    self.zinvRDepSez = sp.vscl(sezScale, rDepHatSez)
    self.zinvRDepRsw1 = sp.vscl(sezScale, rDepHatRsw1)
    self.zinvRDepEci = sp.mtxv(self.mtxEciToRsw1, self.zinvRDepRsw1)

    ####################################################################
    ### Equations (17)
    ### - final velocity displacements
    ### - Not working yet:  what is rChf?
    self.zinvVelDepSez = sp.vpack(-velDepEqcm[iZeq] * sezScale / rChf2
                                #,(velDepEqcm[iAeq] + self.velChfRsw1[iSrsw]) * sezScale *          self.zinvDeltaPhiDep  / rChf2
                                 ,(velDepEqcm[iAeq] + self.velChfRsw1[iSrsw]) * sezScale * math.cos(self.zinvDeltaPhiDep) / rChf2
                                 ,velDepEqcm[iReq] + self.zinvVelChfRsw2[iRrsw]
                                 )
    self.zinvVelDepRsw1 = sp.mtxv(self.zinvMtxRswToSez, self.zinvVelDepSez)
    self.zinvVelDepEci = sp.mtxv(self.mtxEciToRsw1, self.zinvVelDepRsw1)

    return self.zinvRDepEci,self.zinvVelDepEci


  ######################################################################
  def invE(self,offsetArcLength=None,arcLengthArg=None,rtnIter=False):
    """Find eccentric anomaly corresponding to the input arc length"""

    ### Parse arguments
    assert not (offsetArcLength is not None and arcLengthArg is not None)
    if arcLengthArg is None: arcLengthTarget = self.arcChf + offsetArcLength
    else                   : arcLengthTarget = float(arcLengthArg)

    ### Initial estimate of Eccentric Anomaly
    ### - highest positive or negative multiple of halfpi that gives
    ###   an arc length less than or equal to targetLength
    nQuad = int(arcLengthTarget / self.quadArc)
    arcLengthGuess = nQuad * self.quadArc
    if arcLengthGuess > arcLengthTarget:
      nQuad -= 1
      arcLengthGuess -= self.quadArc
    eaGuess0 = eaGuess = nQuad * halfpi
    eaGuess0PlusTwoPi = eaGuess0 + twopi

    deltaArc =  arcLengthTarget - arcLengthGuess
    tol = self.aChf * 1e-9
    itter = 0
    while abs(deltaArc) > tol:

      itter += 1
      ### Vector from primary focus to Ecc. Anom. guess
      cosea = math.cos(eaGuess)
      sinea = math.sin(eaGuess)
      x = self.aChf * cosea
      y = self.bChf * sinea

      ### Delta-[x,y] of tangent at (x,y) for distance of deltaArc
      delx = -self.aChf * sinea
      dely =  self.bChf * cosea
      fac = deltaArc / sp.vnorm(sp.vpack(delx,dely,0.))
      delx *= fac
      dely *= fac
   
      ### Delta-vector tangent at (x,y) for distance of deltaArc
      ### - [X,Y,Z] = [x+delx, y+dely, 0]
      ### - scale X by 1/a and Y by 1/b to get Eccentric Anomaly
      eaGuess = sp.recrad(sp.vpack((x+delx)/self.aChf
                                   ,(y+dely)/self.bChf
                                   ,0.
                                   )
                          )[iRA]
      while eaGuess0 > eaGuess: eaGuess += twopi
      while (eaGuess0PlusTwoPi) <= eaGuess: eaGuess -= twopi
      arcLengthGuess = self.E(eaGuess)
      deltaArc = arcLengthTarget - arcLengthGuess

    return rtnIter and (eaGuess,itter,) or eaGuess


  ######################################################################
  def TrueToEccAnom(self,ta):
    """Convert True Anomaly to Eccentric Anomaly

https://en.wikipedia.org/wiki/Eccentric_anomaly#From_the_true_anomaly

"""
    return sp.recrad(sp.vpack(self.eccChf + math.cos(ta)
                             ,self.bOverA * math.sin(ta)
                             ,0.
                             )
                    )[iRA]

  ######################################################################
  def E(self,ea):
    """Arc length from peripse to Eccentric Anomaly ea

=====================
Note about arc length
=====================

Calculate the arc length from one position on an elliptical orbit
to a second position on that same orbit, using E, the Incomplete
Elliptic Integral of the Second Kind (IEISK).

The ellipse semi-major axis length is a, the semi-minor axis length is
b, and b <= a.

Arguments to IEISK:

  eccSquared:  the ellipse eccentricity squared = (1 - (b/a)^2)
               N.B. b <= a; 0 <= eccSquared <= 1
               N.B. eccSquared equals ecc2 in the notes below
               N.B. eccSquared is also known as the paramete", m = k^2
               N.B. sqrt(eccSquared) is also known as the elliptic
                    modulus or eccentricity, e = k = sin(alpha), where
                    alpha is the modular angle

  ea1 = the first position on the orbit, expressed as an eccentric
        anomaly (angle, radians) positive from zero at the semi-major
        axis toward the direction of the orbiting body

  ea2 = the second position on the orbit, expressed as an eccentric
        anomaly (angle, radians) positive from zero at the semi-major
        axis toward the direction of the orbiting body

Implementation:

1) The arc length from a, the semi-major axis, is the definite integral
of

  a * sqrt(1 - ecc2 cos(ea)^2)

evaluated from eccentric anomaly ea1 to ea2.

2) However, the integral available, E, the Incomplete Elliptic Integral
of the Second Kind, is the integral of

  a * sqrt(1 - ecc2 sin(ea)^2)

from eccentric anomaly zero to ea, which evalutates angles relative to a
semi-MINOR axis a, with semi-MAJOR axis b, so a < b, and
ecc2 = (1 - (a/b)^2).

2.1) So E is off by PI/2 relative to the desired integral in (1) above.

3) The solution is to use E to calculate the arc length from eccentric
anomaly (PI/2 - ea) to eccentric anomaly (PI/2 - 0).  This is a variable
substitution that effectively replaces cos(ea) in (1) with its
equivalent sin(PI/2 - ea), and converts the integral in (1) to a form
siimlar to that in (2), at which point E can be used.

3.1) So

  arcLength0 = a * integral(1-ecc2 cos(ea)^2,ea=0..ea1)
             = a * integral(1-ecc2 sin(PI/2-ea)^2,ea=PI/2-ea1..PI/2)
             = a * (E(PI/2,ecc2) - E(PI/2-ea1,ecc2))

3.2) and similarly

  arcLength1 = a (E(PI/2,ecc2) - E(PI/2-ea2,ecc2))

3.3) and finally the desired arc length from ea1 to ea2 in (1) above

  = arcLength1 - arcLength0
  = a * ((E(PI/2,ecc2)-E(PI/2-ea2,ecc2)) - (E(PI/2,ecc2)-E(PI/2-ea1,ecc2)))
  = a * (E(PI/2-ea1,ecc2) - E(PI/2-ea2,ecc2))

3.4) Note the difference between the ea1 and ea2 terms are reverses0

============================
End of note about arc length
============================

"""
    return self.quadArc - (self.aChf * E(halfpi-ea,self.eccChfSquared))


  ######################################################################
  def XYtoTrueAnom(self,x,y):
    return sp.recrad(sp.vpack(x-self.cChf,y,0.))[iRA]

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
    return numpy.vstack((rHatEci,sHatEci,wHatEci,))


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
  zinvStateDeputy = va.inverseDeputy()

  import pprint
  d = dict()
  d.update(vars(va))
  for k in d:
    if k[:4]!='zinv': continue
    knoz = k[4].lower() + k[5:]
    if not (knoz in d): continue
    try:
      d[k] = (d[k],d[k]-d[knoz],)
    except:
      print(dict(k=k,knoz=knoz,dk=d[k],dknoz=d[knoz]))

  pprint.pprint(d)

  testEa = -13 * halfpi / 3.
  testArc = va.E(testEa)
  testEaOut,itter = va.invE(testArc,rtnIter=True)
  testEaErr = testEaOut - testEa
  testArcOut = va.E(testEaOut)
  testArcErr = testArcOut - testArc
  pprint.pprint(dict(testEa=testEa,testDegEa=testEa*dpr
                    ,testEaOut=testEaOut,testDegEaOut=testEaOut*dpr
                    ,testArc=testArc,testArcOut=testArcOut
                    ,eaErr=testEaErr,eaFracErr=testEaErr/testEa
                    ,arcErr=testArcErr,arcFracErr=testArcErr/testArc
                    ,iter=itter
                    )
               )
########################################################################
