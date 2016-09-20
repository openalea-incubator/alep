# General useful imports
import random as rd
import numpy as np

# Imports for weather
import pandas
import alinea.septo3d
import datetime
from openalea.deploy.shared_data import shared_data
from alinea.astk.Weather import Weather
from alinea.alep.alep_weather import wetness_rapilly, basic_degree_days
from alinea.astk.TimeControl import *
from alinea.alep.alep_time_control import *

# Imports for wheat
from alinea.alep.wheat import initialize_stand
from alinea.astk.plant_interface import grow_canopy
from alinea.alep.architecture import set_properties, update_healthy_area, get_leaves

# Imports for septoria
from alinea.alep.protocol import *
from alinea.alep.septo3d_v2 import plugin_septoria
from alinea.alep.disease_operation import DU_Generator, generate_stock_du
from alinea.alep.inoculation import InoculationLowerLeaves
from alinea.septo3d.dispersion.alep_interfaces import SoilInoculum, Septo3DSoilContamination, Septo3DTransport
from alinea.alep.dispersal_emission import SeptoriaRainEmission
from alinea.septo3d.alep_interfaces import Septo3DSplash
from alinea.popdrops.alep_interface import PopDropsEmission, PopDropsTransport, emission_csv, diameter_csv
from alinea.alep.washing import RapillyWashing
from alinea.alep.growth_control import PriorityGrowthControl
from alinea.alep.infection_control import BiotrophDUPositionModel
from alinea.alep.senescence import WheatSeptoriaPositionedSenescence
from alinea.alep.disease_outputs import LeafInspector, SeptoRecorder

