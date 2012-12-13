import numpy as np
from pylab import *

dt = 1 # Time step of the simulation
t1, t2 = 100, 200 # Times needed to complete state 1 and state 2
ts  = (t1, t2)
nbSteps = 10 # Number of fractions into a state
n = 300 # Number of time steps achieved
growth_rate = 0.06 # Growth rate by time step

def lesion(nbSteps=nbSteps,ts=ts, n=n, growth_rate=growth_rate): 

    t1, t2 = ts
    Dt1, Dt2 = t1/nbSteps, t2/nbSteps # Number of time steps into each fraction of state

    s1 = np.zeros(nbSteps) # Initialisation of each state
    s2 = np.zeros(nbSteps)

    def growth(iSim, iStep, t, s, dt, Dt, growth_rate):
        """ 
        Parameters
        ==========

        - iSim : int
            current step of simulation
        - iStep : int
            current step of simulation into the state
        - t : float
            time needed to achieve state
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
        
        # If no surface has completed this state
        if iSim < t:
            if iStep == 1:
                current_Dt = 0
                # current_Dt is the fraction of the state in which the simulation is currently running
            else:
                # if iStep%Dt == 0:
                    # current_Dt = (float(iStep)/Dt)-1
                # else:
                current_Dt = floor(float(iStep)/Dt)

            # Balance between input and output until current_Dt
            s[0:current_Dt+1] -= ds[0:current_Dt+1]
            s[1:current_Dt+1] += ds[0:current_Dt]    
        
        # If this state has been completed by at least 1 surface
        else:
            ds = s * (float(dt) / float(Dt))
            s = s - ds
            s[1:] += ds[0:-1]
        
        # In all cases the input is the same
        s[0]+= growth_rate
        
        return s, ds[-1]

    in1 = growth_rate
    for iSim in range(1, n):
    
        iStep=iSim
        s1, in2 = growth(iSim, iStep, t1, s1, dt, Dt1, in1)
        if iSim > t1:
            iStep = iSim-t1
            s2, in3 = growth(iSim, iStep, t2, s2, dt, Dt2, in2)

    n1 = len(s1)*int(Dt1)
    n2 = len(s2)*int(Dt2)
    s3 =  np.zeros(n1+n2)
    s3[0:n1] = np.repeat(s1/Dt1,int(Dt1))
    s3[n1:] = np.repeat(s2/Dt2,int(Dt2))

    plot(s3)
    show()
    
for i in range(1,101):
    if t1%i==0 & t2%i==0:
        lesion(nbSteps=i)

