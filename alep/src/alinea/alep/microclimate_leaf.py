# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 10:58:05 2013

@author: lepse
"""
from alinea.caribu.CaribuScene import CaribuScene
import alinea.caribu.sky_tools.turtle as turtle


def runcaribu(sectors, scene, energy):
    """ 
    Calls Caribu for differents energy sources

    :Parameters:
    ----------

    - `sectors` (int)
    - `scene` - Scene containing the simulated system
    - `energy` (float) - Meteorological mean variable at the global scale. Could be:
        - 'PAR' : Quantum PAR (ppfd) in micromol.m-2.sec-1
        - 'Pluie' : Precipitation (mm)

    :Returns:
    ----------
        
    - 'id_out' (dict) - Meteorological variable at the leaf scale
    """
    energie, emission, direction, elevation, azimuth = turtle.turtle(sectors=sectors, energy=energy) 
    sources = zip(energie,direction)
    c_scene = CaribuScene()    
    idmap = c_scene.add_Shapes(scene)    
    c_scene.addSources(sources)
    output = c_scene.runCaribu(infinity=False)
    id_out = c_scene.output_by_id(output, idmap)
    return id_out


def microclimate_leaf(sectors, mean_globalclimate, scene):
    """ Calculate the radiation and rain interception with the Caribu interception model and the Tair, humidity and wind for each leaf and stem element id

        :Parameters:
        ----------

        - `sectors` (int) 
        - `mean_globalclimate` (Dict) - Dict of global climate averaged over an hourly defined time step
        expected variables are:
            - 'PAR' : Quantum PAR (ppfd) in micromol.m-2.sec-1
            - 'Pluie' : Precipitation (mm)
            - 'Tair' : Temperature of air (Celcius)
            - 'HR': Humidity of air (kPa)
            - 'Vent' : Wind speed (m.s-1)
        - `scene` - Scene containing the simulated system

        :Returns:
        -------

        - `microclimate` (Dict) - Dict of microclimate variables (PAR, radiation, rain, Tair, humidity, wind, Tleaf) for each id of leaf and stem element
            - `radiation`: The global radiation in kJ.m-2.h-1
            PAR=ppfd in micromol.m-2.sec-1 (1 ppfd = 0.2174 Watts.m-2.sec-1 PAR, 1 W.m-2.sec-1 global = 0.48 Watts.m-2.sec-1 PAR)
            - `Tleaf`: Temperature of air near the leaf (Celcius)
    """
    microclimate = {}  
    PAR_leaf = {}
# Temp
    id_out = runcaribu(sectors, scene, energy = mean_globalclimate['PAR'])
    EiInf = id_out['EiInf']
    EiSup = id_out['EiSup']
    for Infid, e in EiInf.iteritems():
        for Supid, a in EiSup.iteritems():
            if Infid == Supid:
                PAR_leaf[Infid] = {'PAR': e + a}
# Rain
    id_out = runcaribu(sectors, scene, energy = mean_globalclimate['Pluie'])
    rain_leaf = id_out['Einc']
# PAR
    id_out = runcaribu(sectors, scene, energy = (((mean_globalclimate['PAR']*0.2174)/0.48)/1000)*3600)
    EiInf = id_out['EiInf']
    EiSup = id_out['EiSup']
    for Infid, e in EiInf.iteritems():
        for Supid, a in EiSup.iteritems():
            if Infid == Supid:
                microclimate[Infid] = {'PAR': PAR_leaf[Infid]['PAR'], 'radiation': e + a, 'rain': rain_leaf[Infid], 'Tair': mean_globalclimate['Tair'], 'humidity':mean_globalclimate['HR'], 'wind':mean_globalclimate['Vent'], 'Tleaf':0} 
                if microclimate[Infid]['PAR'] == 0:
                    microclimate[Infid]['Tleaf'] = microclimate[Infid]['Tair']
                else:
                    microclimate[Infid]['Tleaf'] = microclimate[Infid]['Tair'] + (microclimate[Infid]['PAR'] / 300)
    return microclimate


class CaribuMicroclimModel(object):
    """ Adaptor for Caribu model compliying echap local_microclimate model protocol (climate_model)
    """
    def __init__(self, sectors='46'):
        self.sectors = sectors
    def microclim(self, mean_globalclimate, scene):
        """ Return the local meteo calculated with the microclimate_leaf function

        :Parameters:
        ----------

        - `mean_globalclimate` (Dict) - Dict of global climat averaged over an hourly defined time step 
        - `scene` - Scene containing the simulated system

        :Returns:
        -------

        - `local_meteo` (Dict) - Dict of microclimate variables (radiation, rain, Tair, humidity and wind) for each id of leaf and stem element
        """            
        local_meteo = microclimate_leaf(self.sectors, mean_globalclimate, scene)
        return local_meteo





