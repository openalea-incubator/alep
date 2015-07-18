"""
Main steps to run a simulation of a complex epidemics of septoria & brown rust
"""
# General imports
import pandas as pd



# Imports for weather and scheduling of simulation
from septo_decomposed import get_weather
from alinea.alep.alep_time_control import CustomIterWithDelays
from alinea.astk.TimeControl import (time_filter, IterWithDelays,
                                     thermal_time_filter, DegreeDayModel,
                                     time_control)
                                     
# Imports for septoria
from alinea.alep.brown_rust import plugin_septoria
from alinea.septo3d.dispersion.alep_interfaces import SoilInoculum
from alinea.popdrops.alep_interface import (PopDropsSoilContamination, 
                                            PopDropsEmission,
                                            PopDropsTransport)

# Imports for rust
from alinea.alep.septoria_age_physio import BrownRustFungus
from alinea.alep.inoculation import AirborneContamination
from alinea.alep.dispersal_transport import BrownRustDispersal

# Imports for both diseases
from alinea.alep.protocol import infect, update, disperse, external_contamination
from alinea.alep.infection_control import BiotrophDUProbaModel
from alinea.alep.growth_control import SeptoRustCompetition
from alinea.alep.disease_outputs import SeptoRustRecorder


# Imports for wheat
from alinea.adel.newmtg import move_properties
from alinea.echap.architectural_reconstructions import echap_reconstructions
from alinea.alep.architecture import set_properties

# Imports for weather and scheduling of simulation
from septo_decomposed import get_weather
from alinea.alep.alep_time_control import CustomIterWithDelays
from alinea.astk.TimeControl import (time_filter, IterWithDelays,
                                     thermal_time_filter, DegreeDayModel,
                                     time_control)