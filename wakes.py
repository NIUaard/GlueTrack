#from scipy import *

def wakeTM(s): return -1e12*(344*exp(-sqrt(s/0.00174)))
def wakeLOLA(s): return -1e12*(257.6*exp(-sqrt(s/0.00396))+1.16*cos(1760*s**0.72)/(sqrt(s)+1600*s**1.23+1e-20))

def wake0(s): return wakeTM(s) 
wake0_delta=0

def wake1(s): return 2*wakeTM(s) 
wake1_delta=0

def wake2(s): return 2*wakeTM(s)+wakeLOLA(s) 
wake2_delta=-0.036e12*20/22.44
