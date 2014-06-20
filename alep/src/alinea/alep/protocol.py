""" Define the protocol between plant architecture and lesions """

import random

def external_contamination(g, 
             contamination_source, 
             contamination_model,
             weather_data, 
             label = 'LeafElement', **kwds):
    """ Innoculates fungal objects (dispersal units) on elements of the MTG
        according to contamination_model.

    Parameters
    ----------
    g: MTG
        MTG representing the canopy (and the soil)
    contamination_source : model
        A model that return a stock of DUs objects of the contaminant
        Requires a method named 'emmission(g, weather_data)'
    contamination_model: model
        Model that set the position of each DU/lesion in stock on g
        Requires a method named 'contaminate(g, DU, weather_data)'
    weather_data: pandas DataFrame
        Weather data for the time step
    
    Returns
    -------
    g: MTG
        Updated MTG representing the canopy (and the soil)
    
    """
    
    stock = contamination_source.emission(g, weather_data, **kwds)
    if len(stock) > 0:
        # Allocation of stock of inoculum
        deposits = contamination_model.contaminate(g, stock, weather_data)
        stock = [] # stock has been used (avoid uncontrolled future re-use)
        # Allocation of new dispersal units
        for vid,dlist in deposits.iteritems():
            if g.label(vid).startswith(label):
                leaf = g.node(vid)
                for d in dlist:
                    if not 'dispersal_units' in leaf.properties():
                        leaf.dispersal_units=[]  
                    leaf.dispersal_units.append(d)
    return g


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
        fungal_objects_stock = [] # stock has been used (avoid uncontrolled future re-use)
    return g

def infect(g, dt, 
          infection_control_model=None, 
          label="LeafElement",
          activate=True):
    """ Compute infection success by dispersal units.
    
    Parameters
    ----------
    g: MTG
        MTG representing the canopy (and the soil)
        'dispersal_units' are stored in the MTG as a property
    infection_control_model: model
        Model that checks if the dispersal unit can infect at its position
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
        if infection_control_model:
            infection_control_model.control_position(g, label)
       
        # Find dispersal units on MTG
        dispersal_units = g.property('dispersal_units')

        for vid, du in dispersal_units.iteritems():
            # By leaf element, keep only those which are deposited and active
            dispersal_units[vid] = [d for d in du if d.is_active]
            leaf = g.node(vid)
            for du in dispersal_units[vid]:
                du.infect(dt, leaf)
    return g
    
def update(g, dt,
           growth_control_model,
           senescence_model=None, 
           label="LeafElement",
           activate=True,
           weather_data=None):
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
        # Get lesions with inactive lesions removed
        # les = {k:[l for l in v if l.is_active] 
                    # for k, v in g.property('lesions').iteritems()}
        lesions = g.property('lesions')
        
        # 1. Compute growth demand
        for vid, les in lesions.iteritems():
            les = [lesion for lesion in les if lesion.is_active]
            # Update active lesions
            for lesion in les:
                leaf=g.node(vid)
                if weather_data is None:
                    lesion.update(dt, leaf)
                else:
                    lesion.update(dt, leaf, weather_data)
        
        # 2. Allocate or not growth demand, and compute corresponding production of spores 
        growth_control_model.control(g, label=label)
            
            # if senescence_model:
                # # 4. Call a specific response if lesions are on senescent tissue
                # l = [l for v in les.values() for l in v if l.is_senescent]
                # for lesion in l:             
                    # lesion.senescence_response()
        
        
            # # 1. Determine which lesions will be affected by senescence (optional)
            # if senescence_model:
                # senescence_model.find_senescent_lesions(g)
            
            # # 2. Compute growth demand
            # for vid, l in les.iteritems():
                # # Update active lesions
                # for lesion in l:
                    # leaf=g.node(vid)
                    # if weather_data is None:
                        # lesion.update(dt, leaf)
                    # else:
                        # lesion.update(dt, leaf, weather_data)
            
            # # 3. Allocate or not growth demand, and compute corresponding production of spores 
            # growth_control_model.control(g, label=label)
            
            # if senescence_model:
                # # 4. Call a specific response if lesions are on senescent tissue
                # l = [l for v in les.values() for l in v if l.is_senescent]
                # for lesion in l:             
                    # lesion.senescence_response()     
    return g

def disperse(g,
             emission_model=None,
             transport_model=None,
             fungus_name='', 
             label="LeafElement",
             weather_data=None):
    """ Disperse spores of the lesions of fungus identified by fungus_name.
        
    Parameters
    ----------
    g: MTG
        MTG representing the canopy (and the soil)
    emission_model: model
        Model that is used to calculate how many DU are emitted and the spores they contain
    transport_model: model
        Model that is used to position each DU emitted on g
    fungus_name: str
        Name of the fungus
    label: str
        Label of the part of the MTG concerned by the calculation
    weather_data: panda dataframe
        Weather data for the timestep, or None
    
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

    # (C. Fournier, 22/11/2013 : keep compatibility with previous protocol and add a new one (used in septo3d/echap)
    if weather_data is None:
        DU = emission_model.get_dispersal_units(g, fungus_name, label)
    else: 
        DU = emission_model.get_dispersal_units(g, fungus_name, label, weather_data)
    
    # Transport of dispersal units
    if sum([len(v) for v in DU.itervalues()])>0:
        # (C. Fournier, 22/11/2013 : keep compatibility with previous protocol and add a new one (used in septo3d/echap)
        if weather_data is not None:
            deposits = transport_model.disperse(g, DU, weather_data)
        else:
            deposits = transport_model.disperse(g, DU) # update DU in g , change position, status 
        # Allocation of new dispersal units
        for vid,dlist in deposits.iteritems():
            if g.label(vid).startswith(label):
                leaf = g.node(vid)
                for d in dlist:
                    if not 'dispersal_units' in leaf.properties():
                        leaf.dispersal_units=[]  
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
        for vid, du in dispersal_units.iteritems():
            if g.label(vid).startswith(label):
                leaf = g.node(vid)
                # Sometimes, the list is created but is empty
                if du: 
                    nb_dus = len(du)
                    nb_washed = int(round(leaf.washing_rate*nb_dus))
                    for dispersal_unit in random.sample(du, nb_washed):
                        # disable the DUs according to washing_rate on the leaf
                        dispersal_unit.disable()
                           
            dispersal_units[vid] = [d for d in du if d.is_active]
    return g
