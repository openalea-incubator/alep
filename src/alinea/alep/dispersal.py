""" Gather different strategies for modeling emission of fungus propagules, 
    first step of dispersal.
"""
from typing import Union

def get_sporulating_lesions(lesions:dict, fungus_name:str=None)-> Union[dict,int]:
    """Get and copy sporulating lesions

    Parameters
    ----------
    lesions : dict
         a {vid:[lesion,...], ...} dict of list of lesion objects
    fungus_name : str, optional
        Name of fungus, by default None

    Returns
    -------
    Union[dict,int]
        a {vid:[lesion,...], ...} dict of list of sporulating lesion objects or a number
    
    """    
    # get and copy sporulating lesions

    if fungus_name is None: 
        les = {k: [l for l in v if l.is_sporulating] 
                for k, v in lesions.items()}
    else:
        les = {k:[l for l in v if l.fungus.name is fungus_name and l.is_sporulating] 
                for k, v in g.property('lesions').items()}
    return les

def create_dus(sporulating_lesions:dict, du_per_lesions:Union[dict,int]=1, **kwds) -> dict:    
    """Call emission methods of sporulating_lesions with appropriate du_per_lesions

    Parameters
    ----------
    sporulating_lesions : dict
        dict with sporulating lesion
    du_per_lesions : Union(dict,int), 
        number of dispersal Unit by lesion either a number or dict, by default 1

    Returns
    -------
    dict
        Dispersal unit emitted by the lesion
    """  
 
    DU = {}
    for vid, l in sporulating_lesions.items():
        for il,lesion in enumerate(l):
            if vid not in DU:
                DU[vid] = []
            if isinstance(du_per_lesions, int):
                DU[vid].append(lesion.emission(nb_DU=du_per_lesions, **kwds))
            else:
                DU[vid].append(lesion.emission(nb_DU=du_per_lesions[vid][il], **kwds))
    return DU

# Template emission model ###########################################################
### TODO : check if/how 'group_dus' is interfering:

class Emission(object):
    
    """ Template class for a model of emission of dispersal units 
        that complies with the guidelines of Alep.
    
    Emission is the first step of dispersal, before transport. A class for a model 
    of emission must include a method 'get_dispersal_units'. In this example, each
    lesion that is found will emit.
    """
    
    def __init__(self):
        """ Initialize the model with fixed parameters.
        
        Parameters
        ----------
        """
        pass
        
    def get_dispersal_units(self, lesions:dict, fungus_name:str=None, du_per_lesion:int=1, **kwds)->dict:
        """Compute emission of dispersal units

        Parameters
        ----------
        lesions : dict
            a {vid:[lesion,...], ...} dict of list of lesion objects
        fungus_name : str, optional
            name of fungus, by default None
        du_per_lesion : int, 
            number of dispersal unit by lesion, by default 1

        Returns
        -------
        dict
            Dispersal units emitted by leaf. dict([leaf_id, list of dispersal units])
        """
        sporulating_lesions = get_sporulating_lesions(lesions, fungus_name)
        # eventualy compute here actual number of du_per_lesion
        DU = create_dus(sporulating_lesions, du_per_lesion, **kwds)
        return DU



