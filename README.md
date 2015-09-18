# ValladoAlfano2014
SpiceyPy implementation of Vallado and Alfano, 2014, Curvilinear coordinate transformations for relative motion

cf. Celest Mech Dyn Astr (2014) 118:253–271, DOI 10.1007/s10569-014-9531-1

# va2014.py implements the paper
## Pre-requisites, python modules:
### SpiceyPy:  [https://github.com/AndrewAnnex/SpiceyPy]
### ScyPy

# numerical_E.py is a test of the usage of SciPy special function ellipeinc
## scipy.special.ellipeinc is the Incomplete Elliptic Integral of the Second Kind (IEISK)
## > integral(sqrt(1 - k**2 sin(t)**2)
## Used to calculate arc lengths along an ellipse

# elit.f is extracted from SciPy code
## implements IEISK
## testelit.f is command-line manual interface to elit.f
