import matplotlib.pyplot as plt
import string

def int2str(integer):
    """ Convert an integer to a string.
    """
    return "%d" % integer

fig = plt.figure()
nb_graphs = 10
length = 5
width = nb_graphs/length
all_letters = string.ascii_uppercase
for i in range(nb_graphs):
    graph_handle = 'ax'+int2str(i+1)
    # if i < 4:
        # globals()[graph_handle] = fig.add_subplot(length,width,i+1,
                                                  # xticklabels=[], yscale='log')
    if i < 8:
        globals()[graph_handle] = fig.add_subplot(length,width,i+1, xticklabels=[])
    else:
        globals()[graph_handle] = fig.add_subplot(length,width,i+1)
    graph_letter = all_letters[i]
    globals()[graph_letter] = plt.text(0.05, 0.9, graph_letter, 
                                    fontsize=18, ha='center',
                                    va='center', transform=globals()[graph_handle].transAxes)

# Titles
droplets = plt.text(-0.15, 2.1, 'Number of Infectious Droplets', 
                    fontsize=18, rotation='vertical', transform=ax3.transAxes)
leaf_area = plt.text(-0.15, 1.4, 'Leaf Area', 
                    fontsize=18, rotation='vertical', transform=ax7.transAxes)
coverage = plt.text(-0.15, 0.75, '% Recovered', 
                    fontsize=18, rotation='vertical', transform=ax9.transAxes)
    
        
plt.show(False)
