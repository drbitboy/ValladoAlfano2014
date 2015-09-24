import math
import numpy
import pprint
import SpiceyPy as sp

rp=6500.
ecc=.77
inc=0.
lnode=0.
argp=0.
m0=0.
t0=0.
mu=398600.4415

a=rp/(1-ecc)
b=a*math.sqrt(1-ecc**2)
p=b*b/a

t=2*sp.pi()*math.sqrt(a**3/mu)

elts=[rp,ecc,inc,lnode,argp,m0,t0,mu]

print(sp.conics(elts,0.))
print(sp.conics(elts,0.+t/2))
print(sp.conics(elts,0.+t))

print(dict(a=a,b=b,p=p,t=t,elts=elts,sqrtMuOverP=math.sqrt(mu/p),eccSqrtMuOverP=ecc*math.sqrt(mu/p)))

N = 100

def Asincos(i):
  epoch = t * i / N
  state = sp.conics(elts,epoch)
  trueAnom = sp.recrad(state[:3])[1]
  return state[3:5],numpy.array([math.cos(trueAnom),math.sin(trueAnom),1])

tup = map(Asincos,xrange(N))

Ys,sincoss = map(numpy.array,zip(*tup))

pprint.pprint(Ys[:4])
pprint.pprint(sincoss[:4])

print('='*72)

result = numpy.linalg.lstsq(sincoss,Ys)

pprint.pprint(result)
