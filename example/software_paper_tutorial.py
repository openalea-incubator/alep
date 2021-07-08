#Tutorial by Marc LABADIE for software paper rebuild echap loop with adel protocol

# 1. Get echap reconstruction

from alinea.echap.architectural_reconstructions import EchapReconstructions

echap = EchapReconstructions()
adel = echap.get_reconstruction(name="Mercia", nplants=2,seed=0)

# 2. run simulation 
timestep = 30 #Day degree
steps = 3

canopy_age=300
g = adel.setup_canopy(age=canopy_age)

adel.plot(g)


for i in range(steps):
    canopy_age+=timestep
    g = adel.setup_canopy(age=canopy_age)
    adel.plot(g)
