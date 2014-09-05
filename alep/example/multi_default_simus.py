""" Run simulations with default parameters for:
    - varied years
    - varied number of plants in canopies
    - varied number of phyto-elements by leaf. """

from openalea.multiprocessing.parallel import pymap
from multiprocessing import cpu_count
from septo_decomposed import run_disease
from itertools import product
import numpy
import pandas
import random
import matplotlib.pyplot as plt
from matplotlib import colors, cm
plt.ion()
try:
    import cPickle as pickle
except:
    import pickle
    
years = [1998, 2003, 2004, 2010]
# nb_plants = [3, 6, 10, 20, 30]
nb_plants = [3, 6, 10, 20]
nb_sects = [1, 3, 5, 7]
fracs = [1e-2, 1e-3, 1e-4]

default_yr = 2004
default_nplants = 6
default_nsect = 5
default_frac = 1e-3

scen_default = [(default_yr, default_nplants, default_nsect, default_frac)]
scen_yr = [(yr, default_nplants, default_nsect, default_frac) for yr in years if yr!=default_yr]
scen_pl = [(default_yr, npl, default_nsect, default_frac) for npl in nb_plants if npl!=default_nplants]
scen_ns = [(default_yr, default_nplants, ns, default_frac) for ns in nb_sects if ns!=default_nsect]
scen_fr = [(default_yr, default_nplants, default_nsect, fr) for fr in fracs if fr!=default_frac]
scenarios = scen_default + scen_yr + scen_pl + scen_ns + scen_fr
# scenarios = scen_default

def annual_loop((yr, nplants, nsect, frac)):
    try:
        for i_rep in range(nrep):
            g, recorder = run_disease(start_date = str(yr)+"-10-15 12:00:00", 
                             end_date = str(yr+1)+"-08-01 00:00:00",
                             nplants = nplants,
                             nsect = nsect,
                             disc_level = 30, 
                             dir = './adel/adel_'+str(yr)+'_'+str(nplants)+'pl_'+str(nsect)+'sect',
                             sporulating_fraction = frac)
            stored_rec = '.\mercia\\default\\recorder_'+str(yr)+'_'+str(nplants)+'pl_'+str(nsect)+'sect'+'_frac'+str(len(str(frac))-2)+'_rep'+str(i_rep)+'.pckl'
            f_rec = open(stored_rec, 'w')
            pickle.dump(recorder, f_rec)
            f_rec.close()
            del recorder
        print 'evaluation successful'
    except:
        print 'evaluation failed'

def num_leaf_to_str(num_leaves=range(1,5)):
    return ['F%d' % lf for lf in num_leaves]
        
def get_recorder(yr, nplants, nsect, frac, i_rep):
    stored_rec = '.\mercia\\default\\recorder_'+str(yr)+'_'+str(nplants)+'pl_'+str(nsect)+'sect'+'_frac'+str(len(str(frac))-2)+'_rep'+str(i_rep)+'.pckl'
    f_rec = open(stored_rec)
    recorder = pickle.load(f_rec)
    f_rec.close()
    return recorder

def get_recorder_list(scenario = (default_yr, default_nplants, default_nsect, default_frac), nb_rep = 5):
    return [get_recorder(*scenario+(i_rep,)) for i_rep in range(nb_rep)]

def mean_by_leaf(recorder, variable='necrosis_percentage'):
    df_mean_by_leaf = pandas.DataFrame(index = recorder.values()[0].values()[0].degree_days)
    for leaf in range(1, max(len(v) for v in recorder.itervalues())+1):
        lf = 'F%d' % leaf
        df_mean_by_leaf[lf] = numpy.mean([v[lf].data[variable] for v in recorder.itervalues() if lf in v], axis=0)
    return df_mean_by_leaf

def get_mean_by_leaf(scenario = (default_yr, default_nplants, default_nsect, default_frac),
                     nb_rep = 5, 
                     variable='necrosis_percentage'):
    df_means = [mean_by_leaf(reco, variable=variable) for reco in get_recorder_list(scenario=scenario, nb_rep=nb_rep)]
    return glue_df_means(df_means, nb_rep)
    
def glue_df_means(df_means, nb_rep=5):
    glued = pandas.concat(df_means, axis=1, keys=range(nb_rep))
    glued.swaplevel(0, 1, axis=1).sortlevel(axis=1)
    return glued.groupby(level=1, axis=1).mean()

def plot_from_df(df,
                 cmap = cm.jet,
                 ax = None, 
                 xlabel = True,
                 ylabel = None,
                 legend = True,
                 title = None,
                 annotation = None,
                 xlim = None,
                 ylim = None):
    # Plot
    if ax==None:
        fig, ax = plt.subplots(1,1)
    indcolors = numpy.linspace(0, 300, len(df.columns))
    cmapdict = {col:cmap(int(indcolors[i])) for i,col in enumerate(df.columns)}
    h = []
    for lf in df.columns:
        h += ax.plot(df.index, df.ix[:, lf], color=cmapdict[lf])
    # Customize axes
    if xlabel == True:
        ax.set_xlabel('Degree-days', fontsize=18)
    if ylabel != None:
        ax.set_ylabel(ylabel, fontsize=18)
    ax.xaxis.grid(True)
    ax.yaxis.grid(True)
    if legend == True:
        ax.legend(h, df.columns.tolist(), title='Leaf number', bbox_to_anchor=(1.01, 0.5), loc=6, borderaxespad=0.) 
    if title != None:
        ax.set_title(title, fontsize = 20)
    if annotation != None:
        ax.annotate(annotation, xy=(0.05, 0.85), xycoords='axes fraction', fontsize=14)
    return ax

def plot_stability(scenario = (default_yr, default_nplants, default_nsect, default_frac), 
                    nb_reps = [1, 5, 10, 50, 100, 300, 500],
                    variable='necrosis_percentage',
                    ylabel = 'Necrosis percentage',
                    num_leaves = range(1,4)):
    yr, nplants, nsect, frac = scenario
    fig, axs = plt.subplots(len(num_leaves), 2)
    color_names = random.sample([c for c in colors.cnames if c.startswith('dark')], len(nb_reps))
    rep_artists = [plt.Line2D((0,1),(0,0), color = colors.colorConverter.to_rgb(col)) for col in color_names]
    for ind, nb_rep in enumerate(nb_reps):
        df_means = []
        for i_rep in range(nb_rep):
            print i_rep+sum(nb_reps[:ind])
            stored_rec = '.\mercia\\stability\\recorder_'+str(yr)+'_'+str(nplants)+'pl_'+str(nsect)+'sect'+'_frac'+str(len(str(frac))-2)+'_rep'+str(i_rep+sum(nb_reps[:ind]))+'.pckl'
            f_rec = open(stored_rec)
            recorder = pickle.load(f_rec)
            f_rec.close()
            del f_rec
            df_means.append(mean_by_leaf(recorder = recorder, variable = variable).ix[:, num_leaf_to_str(num_leaves)])
            del recorder
        df_mean = glue_df_means(df_means = df_means, nb_rep = nb_rep)
        
        for num_lf in num_leaves:
            lf = 'F%d' % num_lf
            df = df_mean.ix[:, lf].to_frame()
            ax = axs[num_lf-1]
            color = colors.ListedColormap(color_names[ind])
            ax[0] = plot_from_df(df, ax = ax[0], cmap = color,
                         xlabel = True if num_lf == num_leaves[-1] else False,
                         ylabel = ylabel, 
                         legend = False,
                         annotation = '\nyear: %d\nnb plants: %d\nnb sectors: %d\nfraction spo: %.4f\n' % scenario + lf if num_lf == num_leaves[0] else lf)
            ax[0].set_ylim([0,1])
            ax[1] = plot_from_df(df, ax = ax[1], cmap = color,
                         xlabel = True if num_lf == num_leaves[-1] else False,
                         ylabel = None, 
                         legend = False,
                         annotation = None)
            xlims = [0.99*df[lf].argmax(), 1.01*df[lf].argmax()]
            ylims = [0.99*df[lf].max(), 1.01*df[lf].max()]
            ax[1].set_xlim(xlims)
            ax[1].set_ylim(ylims)
    
    h, l = axs[0][1].get_legend_handles_labels()
    axs[0][1].legend(h+rep_artists, l+['%d rep' % nrep for nrep in nb_reps], bbox_to_anchor=(1.01, 0.5), loc=6, borderaxespad=0.)
    
def plot_mean_by_leaf(scenario = (default_yr, default_nplants, default_nsect, default_frac), 
                      nb_rep = 5, 
                      variable='necrosis_percentage',
                      ylabel = 'Necrosis percentage',
                      num_leaves = range(1, 13),
                      cmap=cm.jet):
    df_mean = get_mean_by_leaf(scenario = scenario, nb_rep = nb_rep, variable = variable)
    ax = plot_from_df(df_mean.ix[:, num_leaf_to_str(num_leaves)], cmap=cmap, ylabel=ylabel, 
                      title = 'Mean of %d repetitions' % nb_rep, 
                      annotation = 'year: %d\nnb plants: %d\nnb sectors: %d\nfraction spo: %.4f' % scenario)

def plot_states_by_leaf(scenario = (default_yr, default_nplants, default_nsect, default_frac), 
                        nb_rep = 5, 
                        num_leaves = range(1, 13),
                        xlim = None,
                        cmap=cm.jet):
    variables = ['ratio_inc', 'ratio_chlo', 'ratio_nec', 'ratio_spo', 'ratio_empty', 'necrosis_percentage']
    ylabels = ['% Incubating', '% Chlorotic', '% Necrotic', '% Sporulating', '% Empty', '% Necrotic + Sporulating + Empty']
    fig, axs = plt.subplots(2,3)
    for i_var, ax in enumerate(axs.flat):
        df_mean = get_mean_by_leaf(scenario = scenario, nb_rep = nb_rep, variable = variables[i_var])
        ax = plot_from_df(df_mean.ix[:, num_leaf_to_str(num_leaves)], ax=ax, cmap = cmap,
                          xlabel = False if ax in axs[0] else True,
                          ylabel = ylabels[i_var], 
                          legend = True if ax == axs[0][-1] else False,
                          annotation = '\nyear: %d\nnb plants: %d\nnb sectors: %d\nfraction spo: %.4f' % scenario if ax == axs[0][0] else None)
        ax.set_ylim([0,1])
        if xlim!=None:
            ax.set_xlim(xlim)
        
def plot_leaves_by_plant(scenario = (default_yr, default_nplants, default_nsect, default_frac), 
                         nb_rep = 5, 
                         variable='necrosis_percentage',
                         ylabel = 'Necrosis percentage',
                         num_leaves = range(1, 13),
                         cmap=cm.jet):   
    recorders = get_recorder_list(scenario = scenario, nb_rep = nb_rep)
    
    fig, axs = plt.subplots(int(numpy.ceil(scenario[1]/3.)),min(scenario[1],3))
    for i_plant, ax in enumerate(axs.flat):
        pl = 'P%d' % (i_plant+1)
        df_mean = glue_df_means([mean_by_leaf({pl:reco[pl]},variable=variable) for reco in recorders])
        ax = plot_from_df(df_mean, ax=ax,
                          xlabel=True if ax in axs[-1] else False,
                          ylabel=ylabel if i_plant%3==0 else None,
                          legend=True if ax == axs[0][-1] else False, 
                          annotation=pl+'\nyear: %d\nnb plants: %d\nnb sectors: %d\nfraction spo: %.4f' % scenario if ax == axs[0][0] else pl)
    
    ymin = min([ax.get_ylim()[0] for ax in axs.flat])
    ymax = max([ax.get_ylim()[1] for ax in axs.flat])
    for ax in axs.flat:
        ax.set_ylim([ymin, ymax])
    
def group_recorders(scenario = (default_yr, default_nplants, default_nsect, default_frac),
                    num_leaves = range(1,5), 
                    nb_rep = 5):
    recorders = []
    for i_rep in range(nb_rep):
        recorder = get_recorder(*scenario+(i_rep,))
        Fs = num_leaf_to_str(num_leaves)
        Ps = recorder.keys()
        recorders.append([recorder[p][f] for p,f in product(Ps,Fs)])
    return recorders

def get_audpc_by_leaf(scenario = (default_yr, default_nplants, default_nsect, default_frac), max_ddays=500., num_leaves = range(1,9,2), nb_rep = 5):
    df_audpc_by_leaf = pandas.DataFrame(index = num_leaf_to_str(num_leaves), columns=range(nb_rep))
    for i_rep in range(nb_rep):
        recorder = get_recorder(*scenario+(i_rep,))
        df_audpc_by_leaf[i_rep] = [numpy.mean([r[lf].get_audpc(max_ddays=max_ddays) for r in recorder.values()]) for lf in num_leaf_to_str(num_leaves)]
    return df_audpc_by_leaf.mean(axis=1)
    
def get_max_necrosis_by_leaf(scenario = (default_yr, default_nplants, default_nsect, default_frac), num_leaves = range(1,9,2), nb_rep = 5):
    df_max_nec_by_leaf = pandas.DataFrame(index = num_leaf_to_str(num_leaves), columns=range(nb_rep))
    for i_rep in range(nb_rep):
        recorder = get_recorder(*scenario+(i_rep,))
        df_max_nec_by_leaf[i_rep] = [numpy.mean([max(r[lf].data.necrosis_percentage) for r in recorder.values()]) for lf in num_leaf_to_str(num_leaves)]
    return df_max_nec_by_leaf.mean(axis=1)

def get_scenarios_and_labels(factor='year'):
    if factor == 'year':
        scen = [sc for sc in scenarios if sc[1]==default_nplants and sc[2]==default_nsect and sc[3]==default_frac]
        scen.sort(key=lambda x: x[0])
        labels = [str(x[0]) for x in scen]
    elif factor == 'nb_plants':
        scen = [sc for sc in scenarios if sc[0]==default_yr and sc[2]==default_nsect and sc[3]==default_frac]
        scen.sort(key=lambda x: x[1])
        labels = [str(x[1]) for x in scen]
    elif factor == 'nb_sects':
        scen = [sc for sc in scenarios if sc[0]==default_yr and sc[1]==default_nplants and sc[3]==default_frac]
        scen.sort(key=lambda x: x[2])
        labels = [str(x[2]) for x in scen]
    elif factor == 'fraction_spo':
        scen = [sc for sc in scenarios if sc[0]==default_yr and sc[1]==default_nplants and sc[2]==default_nsect]
        scen.sort(key=lambda x: x[3])
        labels = [str(x[3]) for x in scen]
    return scen, labels

def customize_ax(ax, factor, xs, labels, ylims=[0,350]):
    ax.set_xticks(xs)
    ax.set_xticklabels(labels)
    ax.tick_params(axis='both', labelsize=14)
    ax.set_ylim(ylims)
    xlims = ax.get_xlim()
    ax.set_xlim([min(xlims)-0.1*max(xlims), max(xlims)*1.1])    
    if factor == 'year':
        xlabel = 'Year'
    elif factor == 'nb_plants':
        xlabel = 'Number of plants in canopy'
    elif factor == 'nb_sects':
        xlabel = 'Number of sectors by leaf'
    elif factor == 'fraction_spo':
        xlabel = 'Initial infectious fraction of soil (m2)'
    ax.set_xlabel(xlabel, fontsize=18)
    ax.yaxis.grid(True)
       
def plot_audpc(factor='year', max_ddays=500., num_leaves=range(1,9,2), nb_rep=5, ylims=[0,150], cmap=cm.jet):
    scen, labels = get_scenarios_and_labels(factor=factor)
    xs = range(1, len(labels)+1)
    fig, ax = plt.subplots(1, 1)
    indcolors = numpy.linspace(0, 300, len(num_leaves))
    df_audpcs = pandas.concat([get_audpc_by_leaf(sc, max_ddays = max_ddays, num_leaves = num_leaves, nb_rep = nb_rep) for sc in scen], axis=1, keys=xs)
    h = []
    for i, lf in enumerate(df_audpcs.index.tolist()):
        h+=ax.plot(df_audpcs.columns, df_audpcs.ix[lf], color=cmap(int(indcolors[i])), linestyle='-', marker = 'o', markersize=10)
    customize_ax(ax, factor, xs, labels, ylims=ylims)
    ax.set_ylabel('AUDPC', fontsize=18)
    ax.legend(h, df_audpcs.index.tolist(), bbox_to_anchor=(1.01, 0.5), loc=6, borderaxespad=0.)
    return h
    
def plot_max_necrosis(factor='year', num_leaves=range(1,9,2), nb_rep=5, ylims=[0,1], cmap=cm.jet):
    scen, labels = get_scenarios_and_labels(factor=factor)
    xs = range(1, len(labels)+1)
    fig, ax = plt.subplots(1, 1)
    indcolors = numpy.linspace(0, 300, len(num_leaves))
    df_max_necs = pandas.concat([get_max_necrosis_by_leaf(sc, num_leaves = num_leaves, nb_rep = nb_rep) for sc in scen], axis=1, keys=xs)
    h = []
    for i, lf in enumerate(df_max_necs.index.tolist()):
        h+=ax.plot(df_max_necs.columns, df_max_necs.ix[lf], color=cmap(int(indcolors[i])), linestyle='-', marker = 'o', markersize=10)
    customize_ax(ax, factor, xs, labels, ylims=ylims)
    ax.set_ylabel('AUDPC', fontsize=18)
    ax.legend(h, df_max_necs.index.tolist(), bbox_to_anchor=(1.01, 0.5), loc=6, borderaxespad=0.)
    return h
    
def plot_comparisons(factor='year',
                     nb_rep=5, 
                     variable='necrosis_percentage',
                     ylabel = 'Necrosis percentage',
                     num_leaves = range(1, 13),
                     cmap=cm.jet):
    scen, labels = get_scenarios_and_labels(factor=factor)
    annotation = 'year: %d\nnb plants: %d\nnb sectors: %d\nfraction spo: %.4f' % scen[0]
    fig, axs = plt.subplots(1, len(scen))
    for i_scen, ax in enumerate(axs.flat):
        df_mean = get_mean_by_leaf(scenario = scen[i_scen], nb_rep = nb_rep, variable = variable)
        ax = plot_from_df(df_mean, ax = ax, ylabel = ylabel, title = labels[i_scen],
                          legend = True if ax==axs[-1] else False,
                          annotation = annotation if ax==axs[0] else None)
    
    ymin = min([ax.get_ylim()[0] for ax in axs.flat])
    ymax = max([ax.get_ylim()[1] for ax in axs.flat])
    for ax in axs.flat:
        ax.set_ylim([ymin, ymax])

def plot_all_audpcs(num_leaves=range(1,5)):
    factors = ['year', 'nb_plants', 'nb_sects', 'fraction_spo']
    for fc in factors:
        plot_audpc(fc, num_leaves)
        
def plot_all_max_necrosis(num_leaves=range(1,5)):
    factors = ['year', 'nb_plants', 'nb_sects', 'fraction_spo']
    for fc in factors:
        plot_max_necrosis(fc, num_leaves)
   
# if __name__ == '__main__':
    # nb_cpu = cpu_count()
    # pymap(annual_loop, scenarios, nb_cpu-1)