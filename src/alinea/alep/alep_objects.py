"""Base Alep objects"""
import copy

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


class FungalObject(object):
    """Mother class for Lesion and DU"""

    fungus = None

    def __init__(self, kind, n_objects=1, **mutations):
        assert kind in ('DispersalUnit', 'Lesion')
        self.kind = kind
        self.n_objects = n_objects
        self.mutations = mutations

    def parameters(self):
        """ Get parameters of parent fungus, updated with mutations
        """
        pars = {}
        if self.fungus:
            pars.update(self.fungus._parameters)
        pars.update(self.mutations)
        return DotDict(pars)

    def is_like(self, other):
        """ return True if foreign has the same origin and parameters as itself"""
        return type(self) is type(other) and self.parameters() == other.parameters()

    def merge(self, other):
        self.n_objects += other.n_objects

    def produce(self, n_objects=1, **mutations):
        if self.fungus is None:
            raise TypeError("fungus undefined, production unavailable")
        transmitted_mutations = copy.deepcopy(self.mutations)
        transmitted_mutations.update(mutations)
        if self.kind == 'Lesion':
            return self.fungus.dispersal_unit('DispersalUnit', n_objects, **transmitted_mutations)
        elif self.kind == 'DispersalUnit':
            return self.fungus.lesion('Lesion', n_objects, **transmitted_mutations)
        else:
            raise TypeError('Undefined production type : ' + self.kind)
