""" Generate figures of simulation outputs """

try:
    import cPickle as pickle
except:
    import pickle
import pandas
import numpy
from alinea.alep.septoria import is_iterable
from collections import OrderedDict
import matplotlib.pyplot as plt
from matplotlib import cm
plt.ion()

def load_out(param_short_name='', params=['1e-3', '1e-2', '1e-1'], nb_rep=15):
    out = {}
    for par in params:
        out[par] = {}
        for i_rep in range(nb_rep):
            print par, i_rep
            if len(param_short_name)>0 and param_short_name[-1]!='_':
                param_short_name += '_'
            # stored_rec = '.\mercia\\recorder_'+param_short_name+par+'_'+str(i_rep)+'.pckl'
            stored_rec = '.\mercia\\old_simu\\recorder_'+param_short_name+par+'_'+str(i_rep)+'.pckl'
            f_rec = open(stored_rec)
            out[par][i_rep] = pickle.load(f_rec)
            f_rec.close()
    return out

def load_out_stability(first_rep=0, nb_rep=200):
    out = {}
    for i_rep in range(nb_rep):
        print i_rep
        # stored_rec = '.\mercia\\recorder_'+str(first_rep+i_rep)+'.pckl'
        stored_rec = '.\mercia\\stability\\recorder_2004_6pl_5sect_frac3_rep'+str(i_rep)+'.pckl'
        f_rec = open(stored_rec)
        out[i_rep] = pickle.load(f_rec)
        f_rec.close()
        del f_rec
    return out
    
def plot_stability(leaf='F1'):
    fig, ax = plt.subplots(1, 1)
    # groups = [5, 10, 15, 20, 50, 100]
    groups = [1, 5, 10, 50, 100, 300, 500]
    # groups = [5]
    for ind in range(len(groups)):
        out = load_out_stability(first_rep=int(sum(groups[:ind])), nb_rep=groups[ind])
        plot_by_plant({'':out}, param_name='', variable = 'necrosis_percentage', group_by='plant',
                plants = ['P%d' % p for p in range(1, 4)], leaves = [leaf], ylabel='% Septoria necrotic', axs=ax)
        del out
        
    # Custom
    cmap=cm.jet
    indcolors = numpy.linspace(0, 300, len(groups))
    lines = ax.get_lines()
    i_line = -1
    for line in lines:
        i_line += 1
        line.set_color(cmap(int(indcolors[i_line])))

    ax.legend(['%d rep' % g for g in groups], bbox_to_anchor=(1.0, 0.5), loc=6, borderaxespad=0.5)
    ax.set_title('Mean of repetitions', fontsize = 18)
    
def plot_distribution(out, leaf='F1', variable='necrosis_percentage', ylabel='% Septoria necrotic'):
    fig, ax = plt.subplots(1,1)
    for i_rep, out_rep in out.iteritems():
        data = pandas.DataFrame({pl:v[leaf].data[variable].values for pl,v in out_rep.iteritems()},
                                 index=out_rep['P1'][leaf].data.degree_days)
        ax.plot(data.index, data.mean(axis=1), 'b')
    ax.set_ylabel(ylabel, fontsize=18)
    ax.set_xlabel('Degree days', fontsize=18)
    ax.set_title('Distribution of %d repetitions' % len(out), fontsize=18)
        
def plot_by_plant(out, param_name='Fraction spo', variable = 'leaf_disease_area', group_by=None,
                plants = ['P%d' % p for p in range(1, 4)], leaves = ['F%d' % lf for lf in range(1, 13)], 
                cmap=cm.jet, ylabel='Disease area (cm2)', axs=None, xlims=None, ylims=None, notations=pandas.DataFrame()):
    
    def customize_ax(ax, param, title):
        if od.keys().index(param)==len(od.keys())-1:
            ax.set_xlabel('Degree days', fontsize=18)
        ax.set_ylabel(ylabel, fontsize=18)
        ax.set_title(title, fontsize=18)
        ax.annotate(param_name+': '+param, xy=(0.05, 0.85), xycoords='axes fraction', fontsize=14)
    
    indcolors = numpy.linspace(0, 300, len(leaves))
    cmapdict = {leaves[i]:cmap(int(indcolors[i])) for i in range(len(indcolors))}
    od = OrderedDict(sorted(out.items(), key=lambda x: x[0]))
    h={}
    if group_by!='plant':
        if axs==None:
            fig, axs = plt.subplots(len(od.keys()),len(plants))
        if group_by==None:
            for param, out_par in od.iteritems():
                for i_rep, out_par_rep in out_par.iteritems():
                    for pl in plants:
                        for lf in leaves:
                            reco = out_par_rep[pl][lf]
                            if len(od.keys())>1:
                                ax = axs[od.keys().index(param)][plants.index(pl)]
                            else:
                                ax = axs[plants.index(pl)]
                            if i_rep==0 and ax==axs.flat[0]:
                                h[lf] = ax.plot(reco.data.degree_days, reco.data[variable], color=cmapdict[lf])
                            else:
                                ax.plot(reco.data.degree_days, reco.data[variable], color=cmapdict[lf])
                            customize_ax(ax, param, pl)
            xmax = max(reco.data.degree_days)

        if group_by=='repetition' or group_by=='leaf':
            for param, out_par in od.iteritems():
                for pl in plants:
                    ax = axs[od.keys().index(param)][plants.index(pl)]
                    for lf in leaves:
                        data = pandas.DataFrame({k:v[pl][lf].data[variable].values for k,v in out_par.iteritems()},
                                                index=out_par[out_par.keys()[0]]['P1']['F10'].data.degree_days)
                        if group_by=='repetition':
                            h[lf] = ax.plot(data.mean(axis=1).index, data.mean(axis=1).values, color=cmapdict[lf])
                        elif group_by=='leaf':
                            try:
                                df_mean[lf] = data.mean(axis=1).values
                            except:
                                df_mean = pandas.DataFrame({lf:data.mean(axis=1).values},
                                        index=out_par[out_par.keys()[0]]['P1']['F10'].data.degree_days)
                    if group_by=='leaf':
                        ax.plot(df_mean.mean(axis=1).index, df_mean.mean(axis=1).values, 'b')
                    
                    xmax = max(data.mean(axis=1).index)
                    customize_ax(ax, param, pl)
    else:
        if axs == None:
            fig, axs = plt.subplots(1, len(od.keys()))
        if not is_iterable(axs):
            axs = numpy.array([axs])
        for param, out_par in od.iteritems():
            ax = axs[od.keys().index(param)]
            for lf in leaves:
                for pl in plants:
                    data = pandas.DataFrame({k:v[pl][lf].data[variable].values for k,v in out_par.iteritems()},
                                            index=out_par[out_par.keys()[0]]['P1']['F10'].data.degree_days)
                    try:
                        df_mean[pl] = data.mean(axis=1).values
                    except:
                        df_mean = pandas.DataFrame({pl:data.mean(axis=1).values},
                                index=out_par[out_par.keys()[0]]['P1']['F10'].data.degree_days)
                h[lf] = ax.plot(df_mean.mean(axis=1).index, df_mean.mean(axis=1).values, color=cmapdict[lf])
                xmax = max(data.mean(axis=1).index)
                customize_ax(ax, param, '')
            
    # Adjust the limits of axes
    ymax = max([ax.get_ylim()[1] for ax in axs.flat])
    for ax in axs.flat:
        if xlims==None:
            xlims = [0, xmax]
        if ylims==None:
            ylims = [0, ymax]
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)
        ax.set_xticks(numpy.arange(xlims[0], xlims[-1], 100), minor=True)
        ax.xaxis.grid(True, which='minor')
        ax.yaxis.grid(True)
    
    # Set legend
    if group_by!='leaf':
        labels = sorted(h, key=lambda x: float(x[1:]))
        handles = sum([h[lb] for lb in labels],[])
        if len(od.keys())>1 and group_by!='plant':
            axs[len(axs)/2][-1].legend(handles, labels, title='Leaf number', bbox_to_anchor=(1.05, 0.5), loc=6, borderaxespad=0.)
        else:
            axs[-1].legend(handles, labels, title='Leaf number', bbox_to_anchor=(1.05, 0.5), loc=6, borderaxespad=0.)
    else:
        title = 'Mean on leaves\n'+leaves[0]+' to '+leaves[-1] if len(leaves)>1 else 'Leaf number\n'+leaves[0]
        if len(od.keys())>1 and group_by!='plant':
            axs[len(axs)/2][-1].legend([], [], title=title, bbox_to_anchor=(1.05, 0.5), loc=6, borderaxespad=0.)
        else:
            axs[-1].legend([], [], title=title, bbox_to_anchor=(1.05, 0.5), loc=6, borderaxespad=0.)
        
    # Add notations if given
    if len(notations)>0:
        for ax in axs.flat:
            ax.plot(notations.index, notations.values, 'r*', markersize=10)

def plot_by_state(out, param_name='Fraction spo', plant='P1',
                 leaves = ['F%d' % lf for lf in range(1, 13)], 
                 cmap=cm.jet, xlims=None):
    
    indcolors = numpy.linspace(0, 300, len(leaves))
    cmapdict = {leaves[i]:cmap(int(indcolors[i])) for i in range(len(indcolors))}
    od = OrderedDict(sorted(out.items(), key=lambda x: x[0]))
    fig, axs = plt.subplots(6, len(od.keys()))
    h = {}
    states = {0:'ratio_inc', 1:'ratio_chlo', 2:'ratio_nec', 3:'ratio_spo', 4:'ratio_empty', 5:'necrosis_percentage'}
    ylabels = {0:'% Incubating\nsurface', 1:'% Chlorotic\nsurface', 2:'% Necrotic\nsurface',
                3:'% Sporulating\nsurface', 4:'% Empty\nsurface', 5:'Total necrotic\nseptoria surface'}
    for param, out_par in od.iteritems():
        for state, variable in states.iteritems():
            ax = axs[state][od.keys().index(param)]
            for lf in leaves:
                data = pandas.DataFrame({k:v[plant][lf].data[variable].values for k,v in out_par.iteritems()},
                                        index=out_par[out_par.keys()[0]][plant]['F10'].data.degree_days)
                h[lf] = ax.plot(data.mean(axis=1).index, data.mean(axis=1).values, color=cmapdict[lf])
            xmax = max(data.mean(axis=1).index)
            
            if state==states.keys()[-1]:
                ax.set_xlabel('Degree days', fontsize=16)
            if param==od.keys()[0]:
                ax.set_ylabel(ylabels[state], fontsize=16)
            if state==states.keys()[0]:
                ax.set_title(param_name+': '+param, fontsize=16)
                                        
    # Adjust the limits of axes
    ymax = max([ax.get_ylim()[1] for ax in axs.flat])
    for ax in axs.flat:
        ax.set_ylim([0, ymax])
        if xlims==None:
            xlims = [0, xmax]
        ax.set_xlim(xlims)
        ax.set_xticks(numpy.arange(xlims[0], xlims[-1], 100), minor=True)
        ax.xaxis.grid(True, which='minor')
        ax.yaxis.grid(True)
    
    # Set legend
    labels = sorted(h, key=lambda x: float(x[1:]))
    handles = sum([h[lb] for lb in labels],[])
    axs[2][-1].legend(handles, labels, title='Leaf number', bbox_to_anchor=(1.05, 0.5), loc=6, borderaxespad=0.)
    
# fig.tight_layout(rect=(0,0,0.98,0.98))
# out = load_out(param_name=param_name, params=params, nb_rep=nb_rep)
# notations = pandas.DataFrame({'severity_on_green':0.01},index=[472])