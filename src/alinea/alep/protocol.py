""" Define the protocol between plant architecture and lesions """

def external_contamination(g, 
             contamination_source, 
             contamination_model,
             weather_data = None, 
             label = 'LeafElement', **kwds):
    """ Inoculate fungal objects (DispersalUnit) on elements of the MTG by using a
    contamination_model.

    :Parameters:
     - 'g' (MTG): MTG representing the canopy (Cf Pradal and Godin 2016)
     - 'contamination_source' (Class): A model that returns a stock of DUs objects.
        Requires a method named 'emission(g, weather_data, **kwds)' (see doc)
     - 'contamination_model' (Class): A model that sets the position of each DU in stock on g
        Requires a method named 'contaminate(g, DU, weather_data, label, **kwds)' (see doc)
     - 'weather_data' (pandas DataFrame): Weather data for the time step
     - 'label' (str): Label of the part of the MTG concerned by the calculation
    
    :Returns:
     - 'g' (MTG): Updated MTG representing the canopy
    """
    stock = contamination_source.emission(g, weather_data, **kwds)
    labels = g.property('label')
    if stock > 0:
        # Allocation of stock of inoculum
        deposits = contamination_model.contaminate(g, stock, weather_data, label=label, **kwds)
        stock = 0 # stock has been used (avoid uncontrolled future re-use)

        # Allocation of new dispersal units
        for vid, dlist in deposits.iteritems():
            if len(dlist)>0 and labels[vid].startswith(label):
                leaf = g.node(vid)
                try:
                    leaf.dispersal_units += dlist
                except:
                    leaf.dispersal_units = dlist
    return g

def initiate(g, 
             fungal_objects_stock, 
             initiation_model, 
             label="LeafElement"):
    """ Allocate fungal objects (dispersal units OR lesions) on elements of the MTG
        according to initiation_model.

    :Parameters:
     - 'g' (MTG): MTG representing the canopy (Cf Pradal and Godin 2016)
     - 'fungal_objects_stock' (list): List of fungal objects (dispersal units OR lesions)
       Source of fungal objects to distribute in the scene
     - 'initiation_model' (Class): Model that sets the position of each DU/lesion in stock on g
          Requires a method named 'allocate' (see doc)
     - 'label' (str): Label of the part of the MTG concerned by the calculation

    :Returns:
     - 'g' (MTG): Updated MTG representing the canopy
    """
    vids = [n for n in g if g.label(n).startswith(label)]
    if vids:
        # Allocation of stock of inoculum
        initiation_model.allocate(g, fungal_objects_stock, label)
    fungal_objects_stock = [] # stock has been used (avoid uncontrolled future re-use)
    return g

def infect(g, dt, 
          infection_control_model=None, 
          label="LeafElement"):
    """ Coordinate infection process by dispersal units.
    
    :Parameters:
     - 'g' (MTG): MTG representing the canopy 
        'DispersalUnit' objects and climatic variables are stored in the MTG as properties
     - 'dt' (int): Time step of the simulation
     - 'infection_control_model' (Class): Model that controls if the dispersal unit can infect 
        Requires methods: 'control(g, label)' (see doc)
     - 'label' (str): Label of the part of the MTG concerned by the calculation

    :Returns:
     - 'g' (MTG): Updated MTG representing the canopy    
    """
    # Check if infection possible according to conditions on leaf
    if infection_control_model:
        infection_control_model.control(g, label)
   
    # Find dispersal units on MTG
    dispersal_units = g.property('dispersal_units')

    for vid, du in dispersal_units.iteritems():
        # By leaf element, keep only those which are deposited and active
        leaf = g.node(vid)
        dispersal_units[vid] = [d for d in du if d.is_active]
        for du in dispersal_units[vid]:
            du.infect(dt, leaf)
    return g
    
def update(g, dt,
           growth_control_model=None,
           weather_data=None,
           label="LeafElement"):
    """ Update the status of every lesion on the MTG.
    
    In the framework, the growth of lesions is calculated in two steps:
        - 1. Each lesion ages and calculates its potential growth ('growth_demand') according to 
        climatic conditions.
        - 2. A model at leaf scale gathers all the demands and limits the growth of each lesion
        ('growth_offer') depending on the space available on the leaf.
    
    :Parameters:
     - 'g' (MTG): MTG representing the canopy 
        'DispersalUnit' objects and climatic variables are stored in the MTG as properties
     - 'dt' (int): Time step of the simulation
     - 'growth_control_model' (Class): Model that coordinates thegrowth of lesions on the leaf. 
        Requires methods: 'control(g, label)' (see doc)
     - 'label' (str): Label of the part of the MTG concerned by the calculation

    :Returns:
     - 'g' (MTG): Updated MTG representing the canopy    
    """
    lesions = g.property('lesions')
    # 1. Compute growth demand
    for vid, les in lesions.iteritems():
        for lesion in lesions[vid]:
            if lesion.is_active:
                leaf=g.node(vid)
                if weather_data is None:
                    lesion.update(dt, leaf)
                else:
                    lesion.update(dt, leaf, weather_data)
    
    # 2. Allocate or not growth demand
    if growth_control_model:
        growth_control_model.control(g, label=label)  
    return g

def disperse(g,
             emission_model=None,
             transport_model=None,
             fungus_name='', 
             weather_data=None,
             label="LeafElement",
			 **kwds):
    """ Disperse spores of the lesions of fungus identified by fungus_name.
        
    In the framework, the dispersal of DUs is calculated in two steps:
        - 1. the 'emission_model' browses the MTG to find all emitted DUs
        - 2. the 'dispersal_model' browses the MTG to deposit intercepted DUs
        
    :Parameters:
     - 'g' (MTG): MTG representing the canopy 
        'DispersalUnit' objects and climatic variables are stored in the MTG as properties
    emission_model: model
        Model that is used to calculate how many DU are emitted and the spores they contain
    transport_model: model
        Model that is used to position each DU emitted on g
    fungus_name: str
        Name of the fungus
     - 'weather_data' (pandas DataFrame): Weather data for the time step
     - 'label' (str): Label of the part of the MTG concerned by the calculation
    
    :Returns:
     - 'g' (MTG): Updated MTG representing the canopy
    """
    # Emission of dispersal units
    if weather_data is None:
        DU = emission_model.get_dispersal_units(g, fungus_name, label, **kwds)
    else: 
        DU = emission_model.get_dispersal_units(g, fungus_name, label, weather_data, **kwds)
    # DU is in the following format: dict: {'leaf_id in the MTG': number of DU emitted}
 
    # Transport of dispersal units
    if sum(DU.values())>0:
        if weather_data is not None:
            deposits = transport_model.disperse(g, DU, weather_data, **kwds)
        else:
            deposits = transport_model.disperse(g, DU, **kwds)
        # deposits is in the following format: dict: {'leaf_id in the MTG': list of DU deposited}
            
        # Allocation of new dispersal units
        labels = g.property('label') 
        for vid, dlist in deposits.iteritems():
            if len(dlist)>0 and labels[vid].startswith(label):
                leaf = g.node(vid)
                try:
                    leaf.dispersal_units += dlist
                except:
                    leaf.dispersal_units = dlist
    return g