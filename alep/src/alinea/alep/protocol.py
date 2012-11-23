""" Define the protocol between plant architecture and lesions """


def initiate(g, lesion_class, label="LeafElement"):
    """ Instantiate a Lesion on nodes of g identified by label """
    vids = [n for n in g if g.label(n).startsWith(label)
    infected = random.sample(vids)
    for i in infected:
        lesion_class.instantiate(g,i)
    return g,

def update(g, dt):
    """ Update lesion status .
    
    Lesions are stored in the MTG as lesions property.
    """   
    lesions = g.property('lesions')
    for vid, l in lesions.iteritems():
        for lesion in l:
            lesion.update(dt, g.node(vid))
          
    return g,
    
def disperse(g, disperion_algo):
    """ Spores emit by the lesions """
    # TODO : envoi spores a l'algo + recuperation depots/creation nouvelle lesions
    return g,