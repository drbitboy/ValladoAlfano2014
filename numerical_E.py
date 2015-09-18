"""Numerical approximation of integral of ellipse arc length

- from Eccentric Anomaly of 0, to
- 1%, 5%, 50%, 95%, 99% and 100% of PI/2

- Compare results to scipy.special.ellipeinc, the
  Incomplete Elliptic Integral of the Second Kind

"""
import os
import sys
import math
import numpy
from scipy.special import ellipeinc

### Numerical steps per 1% of PI/2
steppct = 100000

### Total number of steps from 0 to PI/2
steps = steppct * 100

### Upper limit of definite integral
subpcts = (1,5,50,95,99,100,)

### Eccentric Anomalies at each step
thetas = math.pi * (numpy.arange(steps+1))/(2.*steps)

### Cosines and Sines of those Eccentric Anomalies at each step
coss = numpy.cos(thetas)
sins = numpy.sin(thetas)

########################################################################
"""Obsolete function eli2(b2)

- Numerical integral for arc length along ellipse
  = integral(sqrt(a**2 sin(Theta)**2 + b**2 cos(Theta)**2) dTheta)

- Assumes a = 1 = a**2
- b2 = b**2


### Step size
dtheta = thetas[1] - thetas[0]

### Eccentric Anomalies at the midpoint of each step
thetas2 = thetas + (dtheta/2)

### Cosines and Sines squared at those midpoints
coss2 = numpy.cos(thetas2)**2
sins2 = numpy.sin(thetas2)**2


########################################################################
def eli2(b2): 
  return numpy.sum(numpy.sqrt(sins2+(coss2*b2))) * dtheta
"""

########################################################################
def eli(a,b,substep): 
  """Numerical integral for arc length along ellipse from axis a

- sum of sqrt(deltaXi**2 + deltaYi**2)
- deltaXi is (Xi+1 - Xi); deltaYi is (Yi+1 - Yi)

Input arguments

- a,b semi-axes' lengths
- substep - number of steps to use in coss and sins arrays

"""
  xs = a * coss[:substep]   ### X0 to Xn
  ys = b * sins[:substep]   ### Y0 to Yn
  dxs = xs[1:] - xs[:-1]    ### Count of n delta-Xs (Xi+1 - Xi)
  dys = ys[1:] - ys[:-1]    ### Count of n delta-Ys (Yi+1 - Yi)
  return numpy.sum(numpy.sqrt(dxs*dxs + dys*dys))  ### Sum of chord lengths


########################################################################
### Initialize empty list of results
elis = []

### Loop over precentages of PI/2
for subpct in subpcts:

  ### Number of steps:  (steps per %) * pct
  ### and upper limit angle of definite integral:  PI/2 * pct / 100
  substep = steppct * subpct
  theta = math.pi * subpct / 200

  ### Loop over eccentricities from 0 to .99999999 in steps of .11111111
  for ie in xrange(0,10):
    ### Only use 0, .11111111, .55555555, .99999999
    if ie and not ((ie%4)==1): continue
    e = .11111111 * ie
    e2 = e * e

    ### Case 1:  a is semi-major axis of length 1; b < a
    a1a  = a1a2 = 1.
    a1b2 = 1 - (e*e)
    a1b = math.sqrt(a1b2)

    ### Case 2:  b is semi-major axis of length 1; b = 1 > a
    b1b = b1b2 = a1a
    b1a2 = a1b2
    b1a = a1b

    ### Append Cases 1 and 2 results to elis list
    elis.append(
      (subpct                         ### Percent of PI/2 to use for upper limit
      ,theta                          ### Upper limit Eccentric Anomaly angle
      ,e                              ### Eccentricity a.k.a. Elliptic Modulus
      ,a1b                            ### semi-minor axis b for a = 1 > b
      ,eli(a1a,a1b,substep)           ### Case 1 arc length from a for b < a = 1
      ,ellipeinc(math.pi/2,e2)-ellipeinc(math.pi/2-theta,e2)  ### Case 1 actual
      ,eli(b1a,b1b,substep)           ### Case 2 arc length from a for b = 1 > a
      ,ellipeinc(theta,e2)            ### Case 2 actual
      ,)
    )

########################################################################
### Print out results
for subpct,theta,e,a1b,a1eli,Ecompl,b1eli,E in elis:
  print('%s%3d %16.16f %10.8f %12.10f %12.10f %12.10f %12.4e %12.10f %12.10f %12.4e'
  %(e==0. and '\n' or ''
   ,subpct                               ### Upper limit, percent of PI/2
   ,theta                                ### Upper limit Ecc. Anom. angle
   ,e                                    ### Eccentricity
   ,a1b                                  ### semi-minor axis
   ,a1eli                                ### Case 1 numerical result
   ,Ecompl                               ### Case 1 scipy...ellipeinc result
   ,abs(a1eli-Ecompl)*2/(a1eli+Ecompl)   ### Case 1 error
   ,b1eli                                ### Case 2 numerical result
   ,E                                    ### Case 2 scipy...ellipeinc result
   ,abs(b1eli-E)*2/(b1eli+E)             ### Case 2 error
   ,)
  )
