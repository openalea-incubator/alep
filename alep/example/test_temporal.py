import numpy as np
from pylab import *

t1, t2 = 100,100
dt = 1
Dt1 = 10
Dt2 = 30 

ts  = (t1, t2)
Dts = (Dt1, Dt2)
n = 150
growth_rate = 0.23

def lesion(Dts=Dts,ts=ts, n=n, growth_rate=growth_rate): 

    t1, t2 = ts
    Dt1, Dt2 = Dts
    Dt1 = float(t1) / (float(t1)/Dt1)
    Dt2 = float(t2) / (float(t2)/Dt2)

    s1 = np.zeros(int(t1/Dt1))
    s2 = np.zeros(int(t2/Dt2))

    def growth(s, dt, Dt, growth_rate):
        """ 
        Parameters
        ==========

        - s : array
            state variable
        - dt : float
            step of the simulation
        - Dt : float
            step for grouping surfaces
        - growth_rate : float
            constant input = growth by time step

        Returns
        =======
        - output_growth_rate: float
            produced growth rate in the new state
        """
        ds = s * (float(dt) / float(Dt))
        s = s - ds
        s[1:] += ds[0:-1]
        s[0]+= growth_rate

        return s, ds[-1]

    in1 = growth_rate
    for i in range(n):

        s1, in2 = growth(s1, dt, Dt1, in1)
        s2, in3 = growth(s2, dt, Dt2, in2)


    n1 = len(s1)*int(Dt1)
    n2 = len(s2)*int(Dt2)
    s3 =  np.zeros(n1+n2)
    s3[0:n1] = np.repeat(s1/Dt1,int(Dt1))
    s3[n1:] = np.repeat(s2/Dt2,int(Dt2))

    plot(s3)
    
for i in range(1,21,4):
    lesion(Dts=(i,i))

