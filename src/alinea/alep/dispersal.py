""" Gather different strategies for modeling emission of fungus propagules, 
    first step of dispersal.
"""

# Template emission model ###########################################################
### TODO : check if/how 'group_dus' is interfering
class Emission(object):
    """ Template class for a model of emission of dispersal units 
        that complies with the guidelines of Alep.
    
    Emission is the first step of dispersal, before transport. A class for a model 
    of emission must include a method 'get_dispersal_units'. In this example, each
    lesion that is found will emit.
    """
    
    def __init__(self, nb_DU=1):
        """ Initialize the model with fixed parameters.
        
        Parameters
        ----------
        - nb_DU (int): number of DU emited (default 1)
        """
        self.nb_DU = nb_DU
        pass
        
    def get_dispersal_units(self, lesions, fungus_name=None, **kwds):
        """ Compute emission of dispersal units by rain splash on wheat.
        
        Parameters
        ----------
        lesions : a {vid:[lesion,...], ...} dict of list of lesion objects 
        fungus_name: str (optional)
            restrict emmision to lesions belonging to fungus_name
                    
        Returns
        -------
        dispersal_units : dict([leaf_id, list of dispersal units])
            Dispersal units emitted by leaf.
        """
        # get and copy sporulating lesions
        if fungus_name is None: 
            les = {k: [l for l in v if l.is_sporulating] 
                    for k, v in lesions.items()}
        else:
            les = {k:[l for l in v if l.fungus.name is fungus_name and l.is_sporulating] 
                    for k, v in g.property('lesions').items()} 

        DU = {}
        for vid, l in les.items():
            for lesion in l:
                # Compute number of dispersal units emitted by lesion
                if vid not in DU:
                    DU[vid] = []
                DU[vid] += lesion.emission(self.nb_DU, **kwds)
        return DU



