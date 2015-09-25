""" Functions used for sensitivity analysis of septoria model
"""
from disease_sensi_morris import *
from alinea.alep.disease_outputs import (AdelWheatRecorder,
                                        get_synthetic_outputs_by_leaf)
from alinea.alep.simulation_tools.septo_decomposed import annual_loop_septo

### Run and save simulation
def variety_code():
    return {'Mercia':1, 'Rht3':2, 'Tremie12':3, 'Tremie13':4, 'Custom':5}
    
def variety_decode():
    return {v:k for k,v in variety_code().iteritems()}

def get_stored_rec(variety, year, i_sample, i_boot):
    return './septoria/'+variety.lower()+'_'+ str(year)+ \
            '_'+str(int(i_sample))+'_boot'+str(int(i_boot))+'.csv'
    
def run_septoria(sample):
    i_sample = sample.pop('i_sample')
    print '------------------------------'
    print 'i_sample %d' %i_sample
    print '------------------------------'
    i_boot = sample.pop('i_boot')
    year = int(sample.pop('year'))
    variety = variety_decode()[sample.pop('variety')]
    sowing_date = '10-15'
    if year == 2012:
        sowing_date = '10-21'
    elif year == 2013:
        sowing_date = '10-29'
    output_file = get_stored_rec(variety, year, i_sample, i_boot)
    annual_loop_septo(year = year, variety = variety, sowing_date=sowing_date,
                        nplants = 15, output_file = output_file, **sample)
                        
def run_custom_septoria(sample):
    i_sample = sample.pop('i_sample')
    print '------------------------------'
    print 'i_sample %d' %i_sample
    print '------------------------------'
    i_boot = sample.pop('i_boot')
    year = int(sample.pop('year'))
    variety = variety_decode()[sample.pop('variety')]
    sowing_date = '10-15'
    if year == 2012:
        sowing_date = '10-21'
    elif year == 2013:
        sowing_date = '10-29'
    output_file = get_stored_rec(variety, year, i_sample, i_boot)
    annual_loop_septo(year = year, variety = variety, sowing_date=sowing_date,
                        nplants = 15, output_file = output_file,
                        proba_inf=1, Smin=0.02, growth_rate=0.0006,
                        degree_days_to_chlorosis=160., 
                        degree_days_to_necrosis=160.,
                        degree_days_to_sporulation=50.,
                        sporulating_fraction=0.01,
                        reduction_by_rain=0., 
                        rain_events_to_empty=5, 
                        leaf_duration=2.5, 
                        keep_leaves=True, 
                        rh_effect=True, apply_rh='chlorosis', 
                        rh_max=50., rh_min=50.,
                        **sample)

def get_septo_morris_path(year = 2012, variety = 'Tremie12'):
    return './septoria/septo_morris_output_'+variety.lower()+'_'+str(year)+'.csv'
        
def save_sensitivity_outputs(year = 2012, variety = 'Tremie12',
                             parameter_range_file = './septoria/septo_param_range.txt',
                             input_file = './septoria/septo_morris_input_full.txt',
                             nboots = 5):
    parameter_names = pd.read_csv(parameter_range_file, header=None, sep = ' ')[0].values.tolist()
    list_param_names = ['i_sample', 'i_boot', 'year', 'variety'] + parameter_names
    df_out = pd.DataFrame()
    for boot in range(nboots):
        input_boot = input_file[:-9]+'_boot'+str(boot)+input_file[-9:]
        df_in = pd.read_csv(input_boot, sep=' ',
                            index_col=0, names=list_param_names)
        vc = variety_code()
        df_in = df_in[(df_in['year']==year) & (df_in['variety']==vc[variety])]
        i_samples = df_in.index
        df_out_b = pd.DataFrame(columns = ['num_leaf_top'] + list(df_in.columns) + ['normalized_audpc', 'audpc', 'audpc_500'])
        for i_sample in i_samples:
            stored_rec = get_stored_rec(variety, year, i_sample, boot)
            df_reco = pd.read_csv(stored_rec)
            df_s = get_synthetic_outputs_by_leaf(df_reco)
            for lf in np.unique(df_reco['num_leaf_top']):
                df_reco_lf = df_reco[(df_reco['num_leaf_top']==lf)]
                output = {}
                output['num_leaf_top'] = lf
                for col in df_in.columns:
                    output[col] = df_in.loc[i_sample, col]
                output['normalized_audpc'] = df_reco_lf.normalized_audpc.mean()
                output['audpc'] = df_reco_lf.audpc.mean()
                output['audpc_400'] = df_reco_lf.audpc_400.mean()
                output['max_severity'] = df_s[df_s['num_leaf_top']==lf].max_severity.astype(float).values[0]
                df_out_b = df_out_b.append(output, ignore_index = True)
        df_out_b['i_boot'] = boot
        df_out = pd.concat([df_out, df_out_b])
    output_file = get_septo_morris_path(year=year, variety=variety)
    df_out.to_csv(output_file)
    
def plot_septo_morris_by_leaf(year = 2012, variety = 'Tremie12',
                             variable = 'normalized_audpc',
                             parameter_range_file = './septoria/septo_param_range.txt',
                             input_file = './septoria/septo_morris_input.txt',
                             nboots = 5, ylims=None):
    output_file = get_septo_morris_path(year=year, variety=variety)
    df_out = pd.read_csv(output_file)
    plot_morris_by_leaf(df_out, variable=variable,
                        parameter_range_file=parameter_range_file,
                        input_file=input_file, nboots=nboots, ylims=ylims)

def plot_septo_morris_3_leaves(year = 2012, variety = 'Tremie12', 
                               leaves = [10, 5, 1], variable = 'normalized_audpc',
                               parameter_range_file = './septoria/septo_param_range.txt',
                               input_file = './septoria/septo_morris_input.txt',
                               nboots=5, ylims=None):
    output_file = get_septo_morris_path(year=year, variety=variety)
    df_out = pd.read_csv(output_file)
    plot_morris_3_leaves(df_out, leaves=leaves, variable=variable,
                        parameter_range_file=parameter_range_file,
                        input_file=input_file, nboots=nboots, ylims=ylims)
                        
def septo_scatter_plot_by_leaf(year = 2012, variety = 'Tremie12',
                                variable = 'normalized_audpc',
                                parameter = 'Smax'):
    output_file = get_septo_morris_path(year=year, variety=variety)
    df_out = pd.read_csv(output_file)
    scatter_plot_by_leaf(df_out, variable=variable, parameter=parameter)
    
def septo_boxplot_by_leaf(year = 2012, variety = 'Tremie12',
                            variable = 'normalized_audpc',
                            parameter = 'Smax', ylims=None):
    output_file = get_septo_morris_path(year=year, variety=variety)
    df_out = pd.read_csv(output_file)
    boxplot_by_leaf(df_out, variable=variable, 
                    parameter=parameter, ylims=ylims)
    
def septo_boxplot_3_leaves(year = 2012, leaves = [10, 5, 1],
                            variety = 'Tremie12',
                            variable = 'normalized_audpc',
                            parameter = 'Smax', ylims=None):
    output_file = get_septo_morris_path(year=year, variety=variety)
    df_out = pd.read_csv(output_file)
    boxplot_3_leaves(df_out, leaves=leaves, variable=variable,
                    parameter=parameter, ylims=ylims)