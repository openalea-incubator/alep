""" Define the protocol between plant architecture and lesions """


def initiate(g, lesion_factory, label="LeafElement"):
    """ Instantiate a Lesion on nodes of g identified by label """
    import random
    vids = [n for n in g if g.label(n).startswith(label)]
    nbinfected = int(round(random.random() * len(vids)))
    nbinfected = max(nbinfected,1)
    infected = random.sample(vids,nbinfected)
    for i in infected:
        n = g.node(i)
        les = lesion_factory.instantiate_at_stage(spores = 1)
        if not 'lesions' in n.properties():
            n.lesions=[]
        n.lesions.append(les)
    return g,

def update(g, dt):
    """ Update lesion status .
    
    Lesions are stored in the MTG as lesions property.
    """   
    lesions = g.property('lesions')
    for vid, l in lesions.iteritems():
        for lesion in l:
        # proposition 1 on fait ici une correspondance nom attendus dans g <-> noms caracterisant un environnement de lesion (classe a faire ??)
            leaf=g.node(vid)
            #proposition 2 : les conventions de nomage noms sont definies ailleurs (ou???) et on passe juste une leaf qui repond a cette convention
            lesion.update(dt, leaf)
          
    return g,
    
def disperse(g, dispersion_model, lesion_factory, label="LeafElement"):
    """ Disperse spores of the lesions of fungus identified by fungus_name.
    
    New infections occur only on nodes identified by label.
    """
    fungus_name = lesion_factory.fungus
    # arrachage
    dispersal_units = {} # dict of mtg id -> quantity of spores/udin emitted by lesions of a given type
    lesions = g.property('lesions')
    for vid, l in lesions.iteritems():
        for lesion in l:
            if l.fungus is fungus_name:
                leaf = g.node(vid)
                dispersal_units[vid] = l.emission(leaf) # other derterminant (microclimate...) are expected on leaf
    # dispersion
    deposits = dispersion_model.disperse(g, dispersal_units) # deposits is a list of aggregates of spores defined by a (mtg_id, relative_position, nbSpores_in_the_aggregate)
    # creation of new lesions
    for d in deposits:
        vid, pos, nbSpores = d
        if g.label(vid).startswith(label):
            leaf = g.node(vid)
            les = lesion_factory.instantiate(nbSpores)
            if not 'lesions' in n.properties():
                n.lesions=[]
            n.lesions.append(les)

    return g,

    
# def wash(g, washing_model): 
    # """ compute spores loss by washing """
    # lesions = g.property('lesions')
    # spores = {}
    # for vid, l in lesions.iteritems():
        # for lesion in l:
            # spores[vid] = l.getSpores()
    # losses = washing_model(g,spores)
    