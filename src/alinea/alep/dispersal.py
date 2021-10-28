""" Generic dispersal API interacting with Alep Fungus
"""
from functools import reduce
from typing import Tuple, Dict, Any
import random
from collections import Counter
from inspect import getfullargspec

class DotDict(dict):
    """
    a dictionary that supports dot notation
    as well as dictionary access notation
    usage: d = DotDict() or d = DotDict({'val1':'first'})
    set attributes: d.val2 = 'second' or d['val2'] = 'second'
    get attributes: d.val2 or d['val2']
    """
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __init__(self, dct = None):
        if dct is None:
            dct = {}
        for key, value in dct.items():
            if hasattr(value, 'keys'):
                value = DotDict(value)
            self[key] = value


class ParametricModel(object):

    _parameters = {}

    def __init__(self, **parameters):
        self._parameters.update(parameters)

    def parameters(self):
        return DotDict(self._parameters)


def count_dus(dispersal_units: Dict[Any, list]) -> Dict[Any, int]:
    def _sum_dus(du_list):
        return reduce(lambda x, y: x + y.nb_dispersal_units, du_list, 0)

    return {vid: _sum_dus(du_list) for vid, du_list in dispersal_units.items()}


# base dispersal model

class BasicEmission(ParametricModel):

    def __call__(self, sporulating_lesions: Dict[Any, list], demand = 1) -> Dict[Any, list]:
        """ Emissions demands as driven by environmental and/or internal variables"""
        return {vid: [demand for _ in lesions] for vid, lesions in sporulating_lesions.items()}


class BasicTransport(ParametricModel):

    def __call__(self, sources: Dict[Any, int], targets: list=None)-> Dict[Any, Dict[Any, int]]:
        """

        Parameters
        ----------
        sources: a {source_vid: nbDU, ...} dict of total number of DU associated to a vid
        targets:  a list of potential target_vids for dispersal units. if None (default), they just stay where they
            are
        kwds : other args to dispersal transport model

        Returns
        -------
        a {target_vid:{source_vid: nbDU, ...}, ...} dict of dict counting deposits o
        n targets, indexed by source origin
        """
        deposits = {vid: {vid: emission} for vid, emission in sources.items()}
        return deposits


class Dispersal(object):
    """ Generic class for dispersal models interacting with alep Fungus objects"""

    def __init__(self, emission=BasicEmission, transport=BasicTransport, name='basic', **parameters):
        """ Initialize the fungus.

        :Parameters:
         - 'Lesion' (Lesion class): Lesion class for a specific fungus.
         - 'DispersalUnit' (Dispersal class): Dispersal unit class for a specific fungus.
            operates with strictly individual lesions and DUs.
         - 'length_unit' (float): Unit of conversion for dimensions (uses meters as reference:
            thus 0.01 is centimeter)
        """
        #check emission
        try:
            assert isinstance(emission, ParametricModel)
            assert hasattr(emission, '__call__')
            spec = getfullargspec(emission.__call__)
            assert spec.args[1] == 'sporulating_lesions'
        except:
            raise ValueError('emission does not match generic signature of emission models')
        emission_parameters = {k:v for k, v in parameters.items() if k in emission._parameters}
        self.emission = emission(**emission_parameters)
        self.emission_args = spec.args
        # check transport
        try:
            assert isinstance(emission, ParametricModel)
            assert hasattr(emission, '__call__')
            spec = getfullargspec(emission.__call__)
            assert spec.args[1] == 'sources'
            assert spec.args[2] == 'targets'
        except:
            raise ValueError('transport does not match generic signature of transport models')
        transport_parameters = {k: v for k, v in parameters.items() if k in transport._parameters}
        self.transport = transport(**parameters)
        self.transport_args = spec.args
        self.name = name


    def parameters(self):
        return DotDict({'emission': self.emission.parameters(), 'transport': self.transport.parameters()})


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

        dus = {}
        for vid, lesions in sporulating_lesions.items():
            for il, lesion in enumerate(lesions):
                if vid not in dus:
                    dus[vid] = []
                du = lesion.emission(emission_demand=emission_demands[vid][il])
                if du.nb_dispersal_units > 0:
                    found = False
                    for d in dus[vid]:
                        if du.is_like(d):
                            d.merge(du)
                            found = True
                            break
                    if not found:
                        dus[vid].append(du)
        return dus

    def deposits(self, transport: Dict[Any, Dict[Any, int]], dispersal_units: Dict[Any, list])->Dict[Any, list]:
        """

        Parameters
        ----------
        transport : a {target_vid:{source_vid: nbDU, ...}, ...} dict of dict counting deposits on targets,
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
        for vid in emissions:
            random.shuffle(emissions[vid])
        for target, sources in transport.items():
            deposits[target] = []
            for source, ntot in sources.items():
                origins = Counter([emissions[source].pop() for _ in range(ntot)])
                for idu, nb in origins.items():
                    mother_du = dispersal_units[source][idu]
                    du = mother_du.fungus.dispersal_unit(nb_dispersal_units=nb)
                    found = False
                    for d in deposits[target]:
                        if du.is_like(d):
                            d.merge(du)
                            found = True
                            break
                    if not found:
                        deposits[target].append(du)
        return deposits

    def disperse(self, lesions: Dict[Any, list], fungus_name: str = None, targets=None, **kwds) -> Tuple[Dict[Any, list], int]:
        emission_args = {k: v for k, v in kwds.items() if k in self.emission_args}
        transport_args = {k: v for k, v in kwds.items() if k in self.transport_args}
        sporulating_lesions = self.get_sporulating_lesions(lesions, fungus_name= fungus_name)
        if len(emission_args) > 0:
            emissions = self.emission(sporulating_lesions, **emission_args)
        else:
            emissions = self.emission(sporulating_lesions)
        dispersal_units = self.get_dispersal_units(sporulating_lesions, emissions)
        sources = count_dus(dispersal_units)
        if len(transport_args) > 0:
            transport = self.transport(sources, targets, **transport_args)
        else:
            transport = self.transport(sources, targets)
        deposits = self.deposits(transport, dispersal_units)
        loss = sum(sources.values()) - sum(count_dus(deposits).values())
        return deposits, loss



