import multiprocessing
import numpy as np
from poc_sensi import *

def sensi_step(num_step):
    print '\nstep %d' % num_step
    param_values = np.loadtxt('SGInput'+str(num_step)+'.txt', delimiter=' ')
    Y = evaluate(param_values)
    np.savetxt('SGOutput'+str(num_step)+'.txt', Y, delimiter=' ')
    return
    
if __name__ == '__main__':
    for i in range(8):
        p = multiprocessing.Process(target=sensi_step, args=(i+1,))
        p.start()
