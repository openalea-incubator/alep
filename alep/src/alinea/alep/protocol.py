""" Define the protocol between plant architecture and lesions """

import random

def initiate(g, 
             dispersal_units_stock, 
             initiation_model, 
             label="LeafElement"):
    """ Allocates dispersal units (objects) on elements of the MTG according to initiation_model 

    :Parameters:
      - `g` : MTG representing the canopy (and the soil).
      - `dispersal_units_stock` : source of dispersal units to disperse in the scene.
      - `initiation_model` : model that allows positioning each DU in stock on g.

    :Example:
      >>> g = MTG()
      >>> stock = [SeptoriaDU(fungus = septoria(), nbSpores=random.randint(1,100), nature='emitted') for i in range(100)]
      >>> inoculator = RandomInoculation()
      >>> initiate(g, stock, inoculator)
      >>> return g
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
    
    :Parameters:
      - `g` : MTG representing the canopy (and the soil). 'dispersal_units' are stored in the MTG as a property.
      - `dt` : time step of the simulation.
    
    :Example:
      >>> g = MTG()
      >>> stock = [SeptoriaDU(fungus = septoria(), nbSpores=random.randint(1,100), nature='emitted') for i in range(100)]
      >>> inoculator = RandomInoculation()
      >>> initiate(g, stock, inoculator)
      >>> dt = 1
      >>> nb_steps = 50
      >>> for i in range(nb_steps):
      >>>     update_climate(g)
      >>>     infect(g, dt)
      >>> return g
      
    """
    dispersal_units = g.property('dispersal_units')
    for vid, du in dispersal_units.iteritems():
        for dispersal_unit in du:
            leaf = g.node(vid)
            if dispersal_unit.active: # TODO : Condition here ?
                dispersal_unit.infect(dt, leaf)

    return g,
    
def update(g, dt):
    """ Update the status of every lesion on the MTG.
    
    :Parameters:
      - `g` : MTG representing the canopy (and the soil). 'lesions' are stored in the MTG as a property.
      - `dt` : time step of the simulation.
    
    :Example:
      >>> g = MTG()
      >>> stock = [SeptoriaDU(fungus = septoria(), nbSpores=random.randint(1,100), nature='emitted') for i in range(100)]
      >>> inoculator = RandomInoculation()
      >>> initiate(g, stock, inoculator)
      >>> dt = 1
      >>> nb_steps = 1000
      >>> for i in range(nb_steps):
      >>>     update_climate(g)
      >>>     infect(g, dt)
      >>>     update(g,dt)
      >>> return g
    
    """   
    lesions = g.property('lesions')
    for vid, l in lesions.iteritems():
        for lesion in l:
            # proposition 1 on fait ici une correspondance nom attendus dans g <-> noms caracterisant un environnement de lesion (classe a faire ??)
            leaf=g.node(vid)
            #proposition 2 : les conventions de nommage noms sont definies ailleurs (ou???) et on passe juste une leaf qui repond a cette convention
            lesion.update(dt, leaf)
          
    return g,
    
def disperse(g, 
             scene,
             dispersal_model,
             fungus_name, 
             label="LeafElement"):
    """ Disperse spores of the lesions of fungus identified by fungus_name.
    
    New infections occur only on nodes identified by label.
    
    :Parameters:
      - `g` : MTG representing the canopy (and the soil).
      - `dispersal_model` : model that allows positioning each DU in stock on g.
      - `lesion_factory` : class that generates a lesion of a given fungus.

    :Example:
      >>> 
    
    """
  
    # arrachage
    lesions = g.property('lesions')
    DU = {}
    for vid, l in lesions.iteritems():
        for lesion in l:
            if lesion.fungus.name is fungus_name:
                leaf = g.node(vid)
                if vid not in DU:
                    DU[vid] = []
                DU[vid].append(lesion.emission(leaf)) # other derterminant (microclimate...) are expected on leaf
    # dispersion, position a sortir de du ??
    deposits = dispersal_model.disperse(scene, DU) # update DU in g , change position, nature
    # allocation of new dispersal units
    for vid,dlist in deposits.iteritems():
        if g.label(vid).startswith(label):
            leaf = g.node(vid)
            for d in dlist:
                d.deposited()
                if not 'dispersal_units' in vid.properties():
                    vid.dispersal_units=[]
                vid.dispersal_units.append(d)

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
    
