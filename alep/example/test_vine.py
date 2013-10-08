import pandas
from alinea.astk.TimeControl import *
from datetime import datetime, timedelta
from alinea.alep.vine import Vine
from alinea.astk.plant_interface import new_canopy, grow_canopy
from alinea.alep.architecture import get_leaves
from alinea.adel.mtg_interpreter import plot3d
from openalea.plantgl.all import Viewer
vine = Vine()
g,_ = new_canopy(vine, age = 6)
start_date = datetime(2001, 03, 1, 1, 00, 00)
end_date = datetime(2001, 6, 1, 1, 00, 00)
date = start_date
nb_steps = len(pandas.date_range(start_date, end_date, freq='H'))
vine_timing = TimeControl(delay=24, steps=nb_steps)
plot_timing = TimeControl(delay=24, steps=nb_steps)
timer = TimeControler(vine=vine_timing, plotting=plot_timing)
for t in timer:
    print(start_date+timedelta(hours=timer.numiter))
    g,_ = grow_canopy(g, vine, t['vine'])
    if t['plotting'].dt > 0:
        g.node(188).color = (180,0,0)
        scene = plot3d(g)
        Viewer.display(scene)