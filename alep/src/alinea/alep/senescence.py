""" Gather different strategies for response to senescence for fungal lesions on a MTG.

"""

# Imports #########################################################################
from alinea.alep.disease_outputs import compute_green_lesion_areas_by_leaf
# Random inoculation ##############################################################

class WheatSeptoriaPositionedSenescence:
    """ Template class for fungal response to senescence that complies with the 
        guidelines of Alep.
    
    A class for a model of response to senescence must include the following
    method that function at the lesion level:
        - 'find_senescent_lesions': find lesions on senescent tissues and 
        transmit them the information.
                                
    In this example, senescence of a wheat leaf is positionned on the blade axis.   
    """
    
    def __init__(self, g, label='LeafElement'):
        """ Initialize the model of senescence.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy
        label: str
            Label of the part of the MTG concerned by the calculation
        """
        self.position_senescence = {}
        labels = g.property('label')
        vids = (v for v,l in labels.iteritems() if l.startswith(label))
        for vid in vids:
            leaf = g.node(vid)
            self.position_senescence[vid] = leaf.position_senescence
    
    def find_senescent_lesions(self, g, label='LeafElement'):
        """ Find lesions affected by leaf senescence.
        
        Parameters
        ----------
        g: MTG
            MTG representing the canopy
        label: str
            Label of the part of the MTG concerned by the calculation
        
        lesion : Lesion instantiation
            A lesion of septoria with the properties: 'position',
            'dt_before_senescence', 'is_senescent'
        leaf: Leaf sector node of an MTG 
            A leaf sector with the property 'position_senescence'
            
        Returns
        -------
        True or False:
            The lesion is affected or not by senescence
        """
        # Select all the leaf elements
        labels = g.property('label')
        vids = (v for v,l in labels.iteritems() if l.startswith(label))
        
        # Find senescent lesions and transmit them useful data
        lesions = g.property('lesions')
        for v in vids:
            leaf = g.node(v)
            leaf_lesions = [l for l in lesions.get(v,[]) if l.is_active and not l.is_senescent]
            if leaf_lesions:
                for l in leaf_lesions:
                    if leaf.position_senescence!=None and l.position[0] >= leaf.position_senescence:
                        l.become_senescent(old_position_senescence = self.position_senescence[v])
            
            # Save senescence position for next time step
            self.position_senescence[v] = leaf.position_senescence