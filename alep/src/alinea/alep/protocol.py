""" Define the protocol between plant architecture and lesions """

import random

def initiate(g, 
             dispersal_units_stock, 
             initiation_model, 
             label="LeafElement"):
    """ 
    Allocates dispersal units (objects) on elements of the MTG according to initiation_model 

    :Parameters:
      - `g` : MTG representing the canopy (and the soil).
      - `dispersal_units_stock` : source of dispersal units to disperse in the scene.
      - `initiation_model` : model that allows positioning each DU in stock on g.

    :Example:
      >>> g = MTG()
      >>> stock = [SeptoriaDU(fungus = septoria(), nbSpores=random.randint(1,100), nature='emitted') for i in range(100)]
      >>> inoculator = RandomInoculation()
      >>> initiate(g, stock, inoculator)
    """
    # Allocation of stock of inoculum
    initiation_model.disperse(g, dispersal_units_stock)

    return g,

class RandomInoculation:
    """ Example of class created to allocate inoculum on MTG randomly.
    
    """
    
    def disperse(self, g, inoculum, label='LeafElement'):
        """ Select randomly elements of the MTG and allocate them a random part of the inoculum.

        :Parameters:
          - `g` : MTG representing the canopy (and the soil).
          - `inoculum` : source of dispersal units to disperse in the scene.
        """        
        import random
        vids = [n for n in g if g.label(n).startswith(label)]
        infected = random.sample(vids,len(inoculum))
        n = len(vids)
        for vid in vids:
            g.node(vid).dispersal_units = []
            
        for i in inoculum:
            idx = random.randint(0,n-1)
            v = vids[idx]
            # Deposit a DU from inoculum on node v of the MTG  
            i.deposited()
            g.node(v).dispersal_units.append(i)

def infect(g, dt):
    """ Compute infection success by dispersal units.
    
    Dispersal units are stored in the MTG as lesions property.
    """
    dispersal_units = g.property('dispersal_units')
    for vid, du in dispersal_units.iteritems():
        for dispersal_unit in du:
            leaf = g.node(vid)
            if dispersal_unit.active: # TODO : Condition here ?
                dispersal_unit.infect(dt, leaf)

    return g,
    
def update(g, dt):
    """ Update lesion status.
    
    Lesions are stored in the MTG as lesions property.
    """   
    lesions = g.property('lesions')
    for vid, l in lesions.iteritems():
        for lesion in l:
            # proposition 1 on fait ici une correspondance nom attendus dans g <-> noms caracterisant un environnement de lesion (classe a faire ??)
            leaf=g.node(vid)
            #proposition 2 : les conventions de nommage noms sont definies ailleurs (ou???) et on passe juste une leaf qui repond a cette convention
            lesion.update(dt, leaf)
          
    return g,
    
def disperse(g, dispersal_model, lesion_factory, label="LeafElement"):
    """ Disperse spores of the lesions of fungus identified by fungus_name.
    
    New infections occur only on nodes identified by label.
    """
    fungus_name = lesion_factory.fungus.name
    # arrachage
    dispersal_units = {} # dict of mtg id -> quantity of spores/udin emitted by lesions of a given type
    lesions = g.property('lesions')
    for vid, l in lesions.iteritems():
        for lesion in l:
            if lesion.fungus is fungus_name:
                leaf = g.node(vid)
                dispersal_units[vid] = lesion.emission(leaf) # other derterminant (microclimate...) are expected on leaf
    # dispersion
    deposits = dispersal_model.disperse(g, dispersal_units) # deposits is a list of aggregates of spores defined by a (mtg_id, relative_position, nbSpores_in_the_aggregate)
    # creation of new lesions
    for d in deposits:
        vid, pos, nbSpores = d
        if g.label(vid).startswith(label):
            leaf = g.node(vid)
            les = lesion_factory.instantiate(nbSpores)
            if not 'lesions' in vid.properties():
                vid.lesions=[]
            vid.lesions.append(les)

    return g,

def wash(g, washing_model): 
    """ compute spores loss by washing """
    dispersal_units = g.property('dispersal_units')
    for vid, du in dispersal_units.iteritems():
        for dispersal_unit in du:
            leaf = g.node(vid)
            if washing_model(g, dispersal_unit):
                dispersal_unit.inactive()
                        
    # A washing model could use the formalism of Rapilly :
    # washing_rate = rain_intensity / (healthy_surface + rain_intensity)
    # The implementation above must be wrong. Each du will not be treated individually.

def growth_control(g):
    pass
    
def nutrients_uptake(g):
    pass
    
