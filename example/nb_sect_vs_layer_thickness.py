from multi_default_simus import *
from mpl_toolkits.mplot3d import Axes3D



def get_recorder_with_dh(yr, nplants, nsect, frac, dh, i_rep):
    stored_rec = '.\mercia\\stability\\recorder_'+str(yr)+'_'+str(nplants)+'pl_'+str(nsect)+'sect'+'_frac'+str(len(str(frac))-2)+'_rep'+str(i_rep)+'_dh'+str(dh)+'.pckl'
    f_rec = open(stored_rec)
    recorder = pickle.load(f_rec)
    f_rec.close()
    return recorder
    
fracs = [1e-3]
# years = [1998, 2003, 2004]
years = [2004]
nb_plants = [6]
nb_sects = [1, 3, 5, 7, 10, 15]
layer_thicknesses = [1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1]
scenarios = list(product(years, nb_plants, nb_sects, fracs, layer_thicknesses))

# get_norm_audpc_by_leaf(recorder_getter = get_recorder_with_dh, scenario = scenarios[0], num_leaves = [1], nb_rep = 5)

def find_missing_scenarios():
    missing_scenarios = []
    for scenario in scenarios:
        for i_rep in range(5):
            try:
                yr, nplants, nsect, frac, dh = scenario
                stored_rec = '.\mercia\\stability\\recorder_'+str(yr)+'_'+str(nplants)+'pl_'+str(nsect)+'sect'+'_frac'+str(len(str(frac))-2)+'_rep'+str(i_rep)+'_dh'+str(dh)+'.pckl'
                f_rec = open(stored_rec)
                f_rec.close()
                del f_rec
            except:
                missing_scenarios.append(scenario)
    return missing_scenarios
    
def plot_sect_vs_thick(num_leaves = [1, 3, 5, 7, 9]):
    fig = plt.figure()
    for num_yr, yr in enumerate(years):
        for num_lf, lf in enumerate(num_leaves):
        
            print(yr, lf)
        
            df = pandas.DataFrame(index = list(range(len(nb_sects))), columns = list(range(len(layer_thicknesses))))
            for i, nsect in enumerate(nb_sects):
                for j, dh  in enumerate(layer_thicknesses):
                    df.ix[i, j] = sum(get_norm_audpc_by_leaf(recorder_getter=get_recorder_with_dh, scenario = (yr, 6, nsect, 1e-3, dh), num_leaves = [lf], nb_rep = 5))
                    # df.ix[i, j] = sum(get_audpc_by_leaf(recorder_getter=get_recorder_with_dh, scenario = (yr, 6, nsect, 1e-3, dh), num_leaves = [lf], nb_rep = 5))
            ax = fig.add_subplot(len(num_leaves), len(years), len(years)*num_lf+num_yr+1, projection='3d')
            X, Y = numpy.meshgrid(df.columns, df.index)
            ax.plot_surface(X,Y, df.values, rstride=1, cstride=1, color='b', alpha=0.5)
            xs = list(range(len(layer_thicknesses)))
            ax.set_xticks(xs)
            ax.set_xticklabels(layer_thicknesses)
            
            ys = list(range(len(nb_sects)))
            ax.set_yticks(ys)
            ax.set_yticklabels(nb_sects)
            # ax.tick_params(axis='both', labelsize=14)
            ax.set_zlim([0,1])
            ax.annotate('F%d'%lf, xy=(0.05, 0.85), xycoords='axes fraction', fontsize=14)
            if num_lf == 0:
                ax.set_title(str(yr), fontsize=14)
            del df          
    