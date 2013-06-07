""" Define the protocol between plant architecture and lesions """

import random

def initiate(g, 
             fungal_objects_stock, 
             initiation_model, 
             label="LeafElement",
             activate=True):
    """ Allocates fungal objects (dispersal units OR lesions) on elements of the MTG
        according to initiation_model.

    Parameters
    ----------
    g: MTG
        MTG representing the canopy (and the soil)
    fungal_objects_stock: list of fungal objects (dispersal units OR lesions)
        Source of fungal objects to distribute in the scene
    initiation_model: model
        Model that sets the position of each DU/lesion in stock on g
        Requires a method named 'allocate' (see doc)
    label: str
        Label of the part of the MTG concerned by the calculation
    activate: bool
        True if computation is achieved, False otherwise
    
    Returns
    -------
    g: MTG
        Updated MTG representing the canopy (and the soil)
    
    Example
    ------- 
      >>> g = MTG()
      >>> stock = [SeptoriaDU(fungus = septoria(), nbSpores=random.randint(1,100), status='emitted') for i in range(100)]
      >>> inoculator = RandomInoculation()
      >>> initiate(g, stock, inoculator)
      >>> return g
    """
    if activate:
        vids = [n for n in g if g.label(n).startswith(label)]
        if vids:
            # Allocation of stock of inoculum
            initiation_model.allocate(g, fungal_objects_stock, label)

    return g

def infect(g, dt, 
          position_checker_model=None, 
          label="LeafElement",
          activate=True):
    """ Compute infection success by dispersal units.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy (and the soil)
        'dispersal_units' are stored in the MTG as a property
    position_checker_model: model
        Model that disable the DU if it falls on an existing lesion or senescent tissue
        Requires methods: 'check position' and 'disable' (see doc)
    dt: int
        Time step of the simulation
    label: str
        Label of the part of the MTG concerned by the calculation
    activate: bool
        True if computation is achieved, False otherwise
    
    Returns
    -------
    g: MTG
        Updated MTG representing the canopy (and the soil)
    
    Example
    -------
      >>> g = MTG()
      >>> stock = [SeptoriaDU(fungus = septoria(), nbSpores=random.randint(1,100), status='emitted') for i in range(100)]
      >>> inoculator = RandomInoculation()
      >>> initiate(g, stock, inoculator)
      >>> dt = 1
      >>> nb_steps = 50
      >>> for i in range(nb_steps):
      >>>     update_climate(g)
      >>>     infect(g, dt)
      >>> return g
      
    """
    if activate:
        # Check if its position prevent it from infecting (optional)
        if position_checker_model:
            position_checker_model.check_position(g)
        
        # Find dispersal units on MTG
        dispersal_units = g.property('dispersal_units')
        for vid, du in dispersal_units.iteritems():
            # By leaf element, keep only those which are deposited and active
            du = [d for d in du if d.is_active and d.status=="deposited"]
            leaf = g.node(vid)       
            for d in du:
                if not d.can_infect_at_position:
                    d.disable()     
                # If not compute infection success
                if d.is_active:
                    d.infect(dt, leaf)
            # Update the list of dispersal unit by leaf element
            dispersal_units[vid] = [d for d in du if d.is_active]

    return g
    
def update(g, dt,
           growth_control_model,
           senescence_model=None, 
           label="LeafElement",
           activate=True):
    """ Update the status of every lesion on the MTG.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy (and the soil)
        'lesions' are stored in the MTG as a property
    dt: int
        Time step of the simulation
    growth_control_model: model
        Model with rules of competition between the lesions
    senescence_model:
        Model that find lesions on senescent tissue and provoke their response
        Requires a method 'find_senescent_lesions' (see doc)
    label: str
        Label of the part of the MTG concerned by the calculation
    activate: bool
        True if computation is achieved, False otherwise
    
    Returns
    -------
    g: MTG
        Updated MTG representing the canopy (and the soil)
    
    Example
    -------
      >>> g = MTG()
      >>> stock = [SeptoriaDU(fungus = septoria(), nbSpores=random.randint(1,100), status='emitted') for i in range(100)]
      >>> inoculator = RandomInoculation()
      >>> initiate(g, stock, inoculator)
      >>> controler = GrowthControlModel()
      >>> senescence_model = SenescenceModel()
      >>> dt = 1
      >>> nb_steps = 1000
      >>> for i in range(nb_steps):
      >>>     update_climate(g)
      >>>     infect(g, dt)
      >>>     update(g, dt, controler, senescence_model)
      >>> return g
    
    """
    if activate:
        # 1. Determine which lesions will be affected by senescence (optional)
        if senescence_model:
            senescence_model.find_senescent_lesions(g, dt)
        
        # 2. Compute growth demand
        lesions = g.property('lesions')
        for vid, l in lesions.iteritems():
            # Remove inactive lesions
            l = [les for les in l if les.is_active]
            # lesions[vid] = l
            # Update active lesions
            for lesion in l:
                leaf=g.node(vid)
                lesion.update(dt, leaf)
        
        # 3. Allocate or not growth demand, and compute corresponding production of spores 
        growth_control_model.control(g, label=label)
        
        if senescence_model:
            # 4. Call a specific response if lesions are on senescent tissue
            lesions = g.property('lesions')
            l = [l for les in lesions.values() for l in les if l.is_senescent]
            for lesion in l:             
                lesion.senescence_response()
    
    return g
    
def disperse(g, 
             dispersal_model,
             fungus_name, 
             label="LeafElement",
             activate=True):
    """ Disperse spores of the lesions of fungus identified by fungus_name.
        
    Parameters
    ----------
    g: MTG
        MTG representing the canopy (and the soil)
    dispersal_model: model
        Model that is used to position each DU in stock on g
    fungus_name: str
        Name of the fungus
    label: str
        Label of the part of the MTG concerned by the calculation
    activate: bool
        True if computation is achieved, False otherwise
    
    Returns
    -------
    g: MTG
        Updated MTG representing the canopy (and the soil)
        
    Example
    -------
      >>> g = MTG()
      >>> stock = [SeptoriaDU(fungus = septoria(), nbSpores=random.randint(1,100), status='emitted') for i in range(100)]
      >>> inoculator = RandomInoculation()
      >>> initiate(g, stock, inoculator)
      >>> controler = GrowthControlModel()
      >>> dispersor = SplashModel()
      >>> dt = 1
      >>> nb_steps = 1000
      >>> for i in range(nb_steps):
      >>>     update_climate(g)
      >>>     infect(g, dt)
      >>>     update(g,dt, controler)
      >>>     if dispersal_event():
      >>>       scene = plot3d(g)
      >>>       disperse(g, dispersor, "septoria")
      >>> return g
    
    """
    if activate:
        # Release of dispersal units 
        lesions = g.property('lesions')   
        DU = {}
        for vid, l in lesions.iteritems():
            for lesion in l:
                if lesion.fungus.name is fungus_name and lesion.stock_spores > 0.:
                    leaf = g.node(vid)
                    if vid not in DU:
                        DU[vid] = []
                    DU[vid] += lesion.emission(leaf) # other derterminant (microclimate...) are expected on leaf
        # Transport of dispersal units
        if DU:
            deposits = dispersal_model.disperse(g, DU) # update DU in g , change position, status
            # Allocation of new dispersal units
            for vid,dlist in deposits.iteritems():
                if g.label(vid).startswith(label):
                    leaf = g.node(vid)
                    for d in dlist:
                        d.deposited()
                        if not 'dispersal_units' in leaf.properties():
                            leaf.dispersal_units=[]  
                        # /!\ Temp /!\ G.Garin 07/06/2013:
                        # Limit propagation
                        # if random.random()<0.05:                        
                            # leaf.dispersal_units.append(d)
                        leaf.dispersal_units.append(d)

    return g

def wash(g, 
         washing_model, 
         global_rain_intensity, 
         label="LeafElement", 
         activate=True): 
    """ Compute spores loss by washing.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy (and the soil)
    washing_model: model
        Model used to wash the DUs out of the leaf
    global_rain_intensity: float
        Rain intensity over the canopy to trigger washing
    label: str
        Label of the part of the MTG concerned by the calculation
    activate: bool
        True if computation is achieved, False otherwise
    
    Returns
    -------
    g: MTG
        Updated MTG representing the canopy (and the soil)
    
    Example
    -------
      >>> g = MTG()
      >>> stock = [SeptoriaDU(fungus = septoria(), nbSpores=random.randint(1,100), status='emitted') for i in range(100)]
      >>> inoculator = RandomInoculation()
      >>> initiate(g, stock, inoculator)
      >>> washor = WashingModel()
      >>> dt = 1
      >>> nb_steps = 1000
      >>> for i in range(nb_steps):
      >>>     update_climate(g)
      >>>     if global_rain_intensity > 0.:
      >>>       wash(g, washor, global_rain_intensity)
      >>> return g
    """
    if activate:
        # compute washing rate on each leaf
        washing_model.compute_washing_rate(g, global_rain_intensity)
        
        dispersal_units = g.property('dispersal_units')
        # TODO : sort DU with chosen status before the loop.
        for vid, du in dispersal_units.iteritems():
            if g.label(vid).startswith(label):
                leaf = g.node(vid)
                if du: # Sometimes, the list is created but is empty
                    nb_dus = len(du)
                    nb_washed = int(round(leaf.washing_rate*nb_dus))
                    for dispersal_unit in random.sample(du, nb_washed):
                        # disable the DUs according to washing_rate on the leaf
                        dispersal_unit.disable()
                           
            dispersal_units[vid] = [d for d in du if d.is_active]
    return g

def control_growth(g, control_model, label="LeafElement"):
    """ Regulate the growth of lesion according to their growth demand,
    the free space on leaves and a set of rules in the growth model.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy (and the soil)
    control_model: 
        Model with rules of competition between the lesions
    label: str
        Label of the part of the MTG concerned by the calculation
    
    Returns
    -------
    g: MTG
        Updated MTG representing the canopy (and the soil)
    
    """
    control_model.control(g, label)
    return g,

def nutrients_uptake(g):
    pass
