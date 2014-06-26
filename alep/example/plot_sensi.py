""" Plot outputs for first sensitivity test of the model. """
import pickle
import collections
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def is_iterable(obj):
    """ Test if object is iterable """
    return isinstance(obj, collections.Iterable)

def load_out(year, **kwds):
    if len(kwds)>0:
        stored_rec =  '.\sensitivity\recorder_'+str(year)+'_'+kwds.keys()[0][:3]+'-'+str(kwds.values()[0])+'.pckl'
    else:
        stored_rec = '.\sensitivity\recorder_'+str(year)+'.pckl'
    f_rec = open(stored_rec)
    out = pickle.load(f_rec)
    f_rec.close()
    return out

def read_outputs(start_years = [1998, 2001, 2002], **kwds):
    out = {}
    for year in start_years:
        if len(kwds)>0:
            for k in kwds.keys():
                if not is_iterable(kwds[k]):
                    kwds[k] = [kwds[k]]
                for v in kwds[k]:
                    out[str(year)+'_'+k[:3]+'-'+str(v)] = load_out(year, **{k:v})

        else:
            out[str(year)] = load_out(year)
    return out

def get_audpc(out, leaf='rosette', **kwds):
    years = np.unique([k[:4] for k in out.iterkeys()])
    params = np.unique([k[5:] for k in out.iterkeys()])
    if len(kwds.keys())>1:
        raise "can only extract one parameter at a time"
    rows = np.array(kwds.values()[0])
    df = pd.DataFrame()
    for yr in years:
        out_yr = {k:v for k,v in out.iteritems() if k.startswith(yr)}
        audpcs = {k[5:]:out_yr[k][leaf].audpc for k in out_yr.keys()}
        df[yr] = np.array([audpcs[par] for par in params])
    df.index = rows    
    return df
    
def get_mean(out, variable_getter=get_audpc, leaves=['F 1', 'F 2', 'F 3'], **kwds):
    dfs =[]
    for lf in leaves:
        dfs.append(variable_getter(out, leaf=lf, **kwds))
    return np.mean(dfs)
    
# TODO group in same function    
def get_necrosis(out, year=1998, leaf='rosette', **kwds):
    years = np.unique([k[:4] for k in out.iterkeys()])
    params = np.unique([k[5:] for k in out.iterkeys()])
    df = pd.DataFrame()
    out_yr = {k:v for k,v in out.iteritems() if k.startswith(str(year))}
    seq = pd.date_range(start = str(year)+"-10-01 01:00:00", end = str(year+1)+"-07-01 01:00:00", freq='H')
    seq = seq[1:] # TODO check if normal
    surf_necs = {k[5:]:out_yr[k][leaf].data.surface_total_nec for k in out_yr.keys()}
    for par in params:
        df[par] = surf_necs[par]
    df.index = seq
    return df
    
def get_necrosis_percentage(out, year=1998, leaf='rosette', **kwds):
    years = np.unique([k[:4] for k in out.iterkeys()])
    params = np.unique([k[5:] for k in out.iterkeys()])
    df = pd.DataFrame()
    out_yr = {k:v for k,v in out.iteritems() if k.startswith(str(year))}
    seq = pd.date_range(start = str(year)+"-10-01 01:00:00", end = str(year+1)+"-07-01 01:00:00", freq='H')
    seq = seq[1:] # TODO check if normal
    ratio_necs = {k[5:]:out_yr[k][leaf].data.ratio_total_nec for k in out_yr.keys()}
    for par in params:
        df[par] = ratio_necs[par]
    df.index = seq
    return df
    
def get_disease_area(out, year=1998, leaf='rosette', **kwds):
    years = np.unique([k[:4] for k in out.iterkeys()])
    params = np.unique([k[5:] for k in out.iterkeys()])
    df = pd.DataFrame()
    out_yr = {k:v for k,v in out.iteritems() if k.startswith(str(year))}
    seq = pd.date_range(start = str(year)+"-10-01 01:00:00", end = str(year+1)+"-07-01 01:00:00", freq='H')
    seq = seq[1:] # TODO check if normal
    areas = {k[5:]:out_yr[k][leaf].data.leaf_disease_area for k in out_yr.keys()}
    for par in params:
        df[par] = areas[par]
    df.index = seq
    return df
    
def get_green_area(out, year=1998, leaf='rosette', **kwds):
    years = np.unique([k[:4] for k in out.iterkeys()])
    params = np.unique([k[5:] for k in out.iterkeys()])
    df = pd.DataFrame()
    out_yr = {k:v for k,v in out.iteritems() if k.startswith(str(year))}
    seq = pd.date_range(start = str(year)+"-10-01 01:00:00", end = str(year+1)+"-07-01 01:00:00", freq='H')
    seq = seq[1:] # TODO check if normal
    areas = {k[5:]:out_yr[k][leaf].data.leaf_green_area for k in out_yr.keys()}
    for par in params:
        df[par] = areas[par]
    df.index = seq
    return df
    
def plot_necrosis(df_necrosis):
    """
    Dimension of df_necrosis:
            nec_sim1    nec_sim2    nec_sim3 ... nec_sim10
    date1
    date2
    ...
    """
    df_necrosis.plot()
    plt.show(False)
    
def plot_audpc(df_audpc):
    """
    Dimension of df_audpc:
                1998-1999       2001-2002       2002-2003
    audpc_sim1
    audpc_sim2
    audpc_sim3
    ...
    audpc_sim10
    """
    df_audpc.plot()
    plt.show(False)
    
def plot_disease_area(df_area):
    pass

def plot_disease_area_by_leaf(out, year=1998, **kwds):    
    fig, ax = plt.subplots(2,2)
    df_ros = get_disease_area(out, year=year, leaf='rosette', **kwds)
    df_3 = get_disease_area(out, year=year, leaf='F 3', **kwds)
    df_2 = get_disease_area(out, year=year, leaf='F 2', **kwds)
    df_1 = get_disease_area(out, year=year, leaf='F 1', **kwds)
    
    df_1.plot(ax=ax[0][0])
    ax[0][0].annotate('Leaf 1', xy=(0.05, 0.85), xycoords='axes fraction', fontsize=18)
    df_2.plot(ax=ax[0][1])
    ax[0][1].annotate('Leaf 2', xy=(0.05, 0.85), xycoords='axes fraction', fontsize=18)
    df_3.plot(ax=ax[1][0])
    ax[1][0].annotate('Leaf 3', xy=(0.05, 0.85), xycoords='axes fraction', fontsize=18)
    df_ros.plot(ax=ax[1][1])
    ax[1][1].annotate('Rosette', xy=(0.05, 0.85), xycoords='axes fraction', fontsize=18)
    
def plot_green_area_by_leaf(out, year=1998, **kwds):    
    fig, ax = plt.subplots(2,2)
    df_ros = get_green_area(out, year=year, leaf='rosette', **kwds)
    df_3 = get_green_area(out, year=year, leaf='F 3', **kwds)
    df_2 = get_green_area(out, year=year, leaf='F 2', **kwds)
    df_1 = get_green_area(out, year=year, leaf='F 1', **kwds)
    
    df_1.plot(ax=ax[0][0])
    ax[0][0].annotate('Leaf 1', xy=(0.05, 0.85), xycoords='axes fraction', fontsize=18)
    df_2.plot(ax=ax[0][1])
    ax[0][1].annotate('Leaf 2', xy=(0.05, 0.85), xycoords='axes fraction', fontsize=18)
    df_3.plot(ax=ax[1][0])
    ax[1][0].annotate('Leaf 3', xy=(0.05, 0.85), xycoords='axes fraction', fontsize=18)
    df_ros.plot(ax=ax[1][1])
    ax[1][1].annotate('Rosette', xy=(0.05, 0.85), xycoords='axes fraction', fontsize=18)