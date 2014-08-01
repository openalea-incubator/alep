import multiprocessing
import numpy as np
from poc_sensi import *

param_values = np.loadtxt('SGInput.txt', delimiter=' ')

def sensi_step(num_step, nb_cpu):
    nb_params = int(len(param_values)/nb_cpu)
    Y = evaluate(param_values[nb_params*(num_step-1):nb_params*num_step])
    np.savetxt('SGOutput'+str(num_step)+'.txt', Y, delimiter=' ')
    return
    
if __name__ == '__main__':
    nb_cpu = multiprocessing.cpu_count()
    for i in range(nb_cpu):
        p = multiprocessing.Process(target=sensi_step, args=(i+1, nb_cpu))
        p.start()
