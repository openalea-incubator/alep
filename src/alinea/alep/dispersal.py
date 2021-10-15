""" Generic dispersal API interacting with Alep Fungus
"""
from functools import reduce
from typing import Tuple, Dict, Any
import random
from collections import Counter


def count_dus(dispersal_units: Dict[Any, list]) -> Dict[Any, int]:
    def _sum_dus(du_list):
        return reduce(lambda x, y: x + y.nb_dispersal_units, du_list, 0)

    return {vid: _sum_dus(du_list) for vid, du_list in dispersal_units.items()}


class Dispersal(object):
    """ Generic class for dispersal models interacting with alep Fungus objects"""

    def __init__(self, **parameters):
        """ Initialize the model with fixed parameters.

        Parameters
        ----------
        """
        self.parameters = parameters

    # generic APIs
    ##############

    def get_sporulating_lesions(self, lesions: Dict[Any, list], fungus_name: str = None) -> Dict[Any, list]:
        """Get sporulating lesions of a given fungus

        Parameters
        ----------
        lesions : dict
             a {vid:[lesion,...], ...} dict of list of lesion objects
        fungus_name : str, optional
            Name of fungus, by default None

        Returns
        -------
        Union[dict,int]
            a {vid:[lesion,...], ...} dict of list of sporulating lesion objects

        """
        if fungus_name is None:
            les = {k: [l for l in v if l.is_sporulating()]
                   for k, v in lesions.items()}
        else:
            les = {k: [l for l in v if l.fungus.name is fungus_name and l.is_sporulating()]
                   for k, v in lesions.items()}
        return les

    def get_dispersal_units(self, sporulating_lesions: Dict[Any, list], emission_demands: Dict[Any, list])->Dict[Any, list]:
        """Collect actual emissions of sporulating lesions

        Parameters
        ----------
        sporulating_lesions : dict
            a {vid:[lesion,...], ...} dict of list of lesion objects
        emission_demands:
            a {vid:[nb_du,...],...} dict of list of emission demands per lesion attached to vid

        Returns
        -------
        dict
            Dispersal units emitted by sources. {source_vid : [dispersal unit, ...], ...}
        """

        DU = {}
        for vid, lesions in sporulating_lesions.items():
            for il, lesion in enumerate(lesions):
                if vid not in DU:
                    DU[vid] = []
                DU[vid].append(lesion.emission(emission_demand=emission_demands[vid][il]))
        return DU

    def deposits(self, transport_map: Dict[Any, Dict[Any, int]], dispersal_units: Dict[Any, list])->Dict[Any, list]:
        """

        Parameters
        ----------
        transport_map : a {target_vid:{source_vid: nbDU, ...}, ...} dict of dict counting deposits on targets,
            indexed by source origin
        dispersal_units: Dispersal units emitted by sources. a {source_vid : [dispersal unit, ...], ...} dict

        Returns
        -------
            a {target_vid: [dispesal_unit, ...], ...} dict
        """
        deposits = {}
        emissions = {vid: reduce(lambda x, y: x + y,
                                 [[i] * du.nb_dispersal_units for i, du in enumerate(du_list)],
                                 [])
                     for vid, du_list in dispersal_units.items()}
        random.shuffle(emissions)
        for target, sources in transport_map.items():
            deposits[target] = []
            for source, ntot in sources.items():
                origins = Counter([emissions[source].pop() for i in range(ntot)])
                for idu, nb in origins.items():
                    mother_du = dispersal_units[source][idu]
                    deposits[target].append(mother_du.fungus.dispersal_unit(nb_dispersal_units=nb))
        return deposits

    def disperse(self, lesions: Dict[Any, list], fungus_name: str = None,
                 emission_args: dict = None, transport_args: dict = None) -> Tuple[Dict[Any, list], int]:
        if emission_args is None:
            emission_args = {}
        if transport_args is None:
            transport_args = {}
        sporulating_lesions = self.get_sporulating_lesions(lesions, fungus_name= fungus_name)
        emission_demands = self.emission_demands(sporulating_lesions, **emission_args)
        dispersal_units = self.get_dispersal_units(sporulating_lesions, emission_demands)
        sources = count_dus(dispersal_units)
        tmap = self.transport_map(sources, **transport_args)
        deposits = self.deposits(tmap, dispersal_units)
        loss = sum(sources.values()) - sum(count_dus(deposits).values())
        return deposits, loss

    # specific methods to be overwritten

    def emission_demands(self, sporulating_lesions: Dict[Any, list], **kwds) -> Dict[Any, list]:
        """ Emissions as driven by environmental and internal variables"""
        return {vid: [1 for les in lesions] for vid, lesions in sporulating_lesions.items()}

    def transport_map(self, sources: Dict[Any, int], targets: list=None, **kwds)-> Dict[Any, Dict[Any, int]]:
        """

        Parameters
        ----------
        sources: a {source_vid: nbDU, ...} dict of total number of DU associated to a vid
        targets:  a list of potential target_vids for dispersal units. if None (default), they just stay where they
            are
        kwds : other args to dispersal transport model

        Returns
        -------
        a {target_vid:{source_vid: nbDU, ...}, ...} dict of dict counting deposits on targets, indexed by source origin
        """
        deposits = {vid: {vid: emission} for vid, emission in sources.items()}
        return deposits

