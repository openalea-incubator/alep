import random
from typing import Dict, Any
from alinea.alep.dispersal import ParametricModel, Dispersal


class RandomTransport(ParametricModel):

    _parameters = {'p_loss': 0.9}

    def __call__(self, sporulating_lesions: Dict[Any, list], **kwds) -> Dict[Any, list]:
        """ Emissions demands as driven by environmental and/or internal variables"""
        p = self.parameters()
        return {vid: [p.p_loss for _ in lesions] for vid, lesions in sporulating_lesions.items()}


RandomDispersal = Dispersal(transport = RandomTransport)


class RandomDispersal:
    """ Template class for a dispersal model that complies with the guidelines of Alep.

    A class for a model of dispersal must include a method 'disperse'. In this example,
    dispersal units are randomly distributed.

    """


    def __init__(self, p_loss=0.9):


    def disperse(self, g, dispersal_units, time_control = None):
        """ Example method for dispersal with random distribution.

        Parameters
        ----------
        g: MTG
            MTG representing the canopy (and the soil)
        dispersal_units : dict
            Dispersal units emitted by the lesions on leaves

        Returns
        -------
        deposits : dict
            Dispersal units deposited on new position on leaves
        """
        # vids = scene.todict().keys()
        try:
            dt = time_control.dt
        except:
            dt = 1

        vids = [id for id,v in g.property('geometry').items()]
        n = len(vids)
        deposits = {}
        if dt > 0:
            for vid, dlist in dispersal_units.items():
                for d in dlist:
                    if random.random() < 0.1:
                        if n>=1:
                            idx = random.randint(0,n-1)
                            v = vids[idx]
                            d.set_position([0, 0])
                            deposits.setdefault(v,[]).append(d)

        return deposits