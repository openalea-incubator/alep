#Tutorial by Marc LABADIE for software paper rebuild echap loop with adel protocol

# 1. Get echap reconstruction

from pathlib import Path

from alinea.echap.architectural_reconstructions import EchapReconstructions, echap_reconstructions

from alinea.alep.architecture import update_healthy_area

from alinea.adel.newmtg import move_properties


filename= 'echap_reconstruction.pckl'

#if not os.path.exists(filename):
if not Path(filename).exists():
    echap = EchapReconstructions()
    echap.save(filename=filename)
else:
    echap= EchapReconstructions.load(filename=filename)

adel = echap.get_reconstruction(name="Mercia", nplants=2,seed=0)

# 2. run simulation 

#init sumulation
timestep = 30 #Day degree
steps = 10



# 2.1 Grow wheat canopy and vizualized development

## init canopy
canopy_age=300
g = adel.setup_canopy(age=canopy_age)

## init alep
### Add the property 'healthy_area' on the leaves
update_healthy_area(g, label = 'LeafElement')

adel.plot(g)
for i in range(steps):
    canopy_age+=timestep
    
    # update canopy
    newg = adel.setup_canopy(age=canopy_age)
    adel.canopy_age=canopy_age
    move_properties(g, newg)
    g= newg
    update_healthy_area(g, label = 'LeafElement')
    adel.plot(g)
    
# 3. Disease development (septo3D)
    
## Initialize the models for septoria
septoria = plugin_septoria()

