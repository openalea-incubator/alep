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
        les = lesion_factory.instantiate()
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
    
def disperse(g, disperion_algo):
    """ Spores emit by the lesions """
    # TODO : envoi spores a l'algo + recuperation depots/creation nouvelle lesions
    return g,