"""
  SpiceyPy-based class implementing Vallado and Alfano,
  "Curvilinear coordinate transformations for relative motion,"
  Celest Mech Dyn Astr (2014) 118:253-271, DOI 10.1007/s10569-014-9531-1
"""
import os
import sys
import numpy
import SpiceyPy as sp

class VA2014:
  def __init__(self,stateChf,StateDep):

    ### Extract ECI Positions and ECI velocities from 6-element ECI states
    self.rChfEci,self.velChfEci = stateChf[:3], stateChf[3:]
    self.rDepEci,self.velDepEci = StateDep[:3], StateDep[3:]

    ### Convert chief position and velocity to [Rhat|Shat|What]1 matrix
    self.rHatEci = sp.vhat(self.rChfEci)
    self.wHatEci = sp.vcrss(self.rChfEci,self.velChfEci)
    self.sHatEci = sp.vcrss(self.rHatEci,self.wHatEci)
    self.mtxEci2Rsw1 = [self.rHatEci,self.wHatEci,self.sHatEci]

    self.rChfRsw1,self.velChfRsw1,self.rDepRsw1,self.velDepRsw1 = [sp.mxv(self.mtxEci2Rsw1,v)
       for v in 
       [self.rChfEci,self.velChfEci,self.rDepEci,self.velDepEci]
       ]

if "__main__"==__name__:
 if len(sys.argv[1:])!=12:
    arg12 = '1 1 1  .1 .1 -.1   1 1.1 1.1  .1 -.1 +.1'.split()
    while arg12: sys.argv.insert(1,arg12.pop())

 stateChief = numpy.array(map(float,sys.argv[1:7]))
 stateDeputy = numpy.array(map(float,sys.argv[7:13]))

 print(dict(stateChief=stateChief,stateDeputy=stateDeputy))

 va = VA2014(stateChief,stateDeputy)
 import pprint
 pprint.pprint(vars(va))
