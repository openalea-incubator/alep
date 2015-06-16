from collections import OrderedDict
from multiprocessing import cpu_count
from disease_sensi_morris import *
from openalea.multiprocessing.parallel import pymap
from septo_decomposed import run_disease, make_canopy
from alinea.alep.disease_outputs import get_recorder
try:
    import cPickle as pickle
except:
    import pickle
from alinea.alep.disease_outputs import *

### Run and save simulation
def wheat_path((year, variety, nplants, nsect)):
    if variety.lower().startswith('tremie'):
        variety = 'tremie'
    return '../adel/'+variety.lower()+'_'+str(int(year))+'_'+str(nplants)+'pl_'+str(nsect)+'sect'

def make_canopies((yr, day, variety, nplants, nsect, wheat_path)):
    make_canopy(start_date = str(int(yr-1))+"-10-"+str(int(day))+" 12:00:00", end_date = str(int(yr))+"-08-01 00:00:00",
            variety = variety, nplants = nplants, nsect = nsect, disc_level = 5, dir = wheat_path)
    
def reconstruct_wheat(qualitative_parameters, nb_plants = 6, nb_sects = 5):   
    combinations = list(itertools.product(*[qualitative_parameters['year']['values'], variety_code.values(), [nb_plants], [nb_sects]]))
    combinations = map(lambda x: x + (wheat_path(x),), combinations)
    make_canopies(combinations)

def annual_loop(sample):
    try:
        # Get indice of simulation
        i_sample = sample.pop('i_sample')

        # Get year of simulation
        if 'year' in sample:
            year = int(sample.pop('year'))
            if year == 2011:
                start_date = "2010-10-15"
                end_date = "2011-06-20"
            elif year == 2012:
                start_date = "2011-10-21"
                end_date = "2012-07-18"
                variety = 'Tremie12'
            elif year == 2013:
                start_date = "2012-10-29"
                end_date = "2013-08-01"
                variety = 'Tremie13'
            else:
                start_date = str(year-1)+"-10-15"
                end_date = str(year)+"-08-01"
        else:
            year = 2005
            start_date = str(year-1)+"-10-15"
            end_date = str(year)+"-08-01"
        start_date += " 12:00:00"
        end_date += " 00:00:00"
            
        # Get variety
        if 'variety' in sample:
            variety_code = {1:'Mercia', 2:'Rht3', 3:'Tremie12', 4:'Tremie13'}
            variety = variety_code[sample.pop('variety')]
        #else:
         #   variety = 'Mercia'

        # Get wheat path
        nplants = 10
        #nplants = 1
        nsect = 7
        w_path = wheat_path((year, variety, nplants, nsect))

        # Run and save simulation
        g, recorder = run_disease(start_date = start_date, 
                         end_date = end_date, 
                         variety = variety, nplants = nplants, nsect = nsect,
                         dir = w_path, reset_reconst = True, 
                         degree_days_to_necrosis = 70.,
                         degree_days_to_sporulation = 70., **sample)
        stored_rec = './'+variety.lower()+'/recorder_'+str(int(i_sample))+'.pckl'
        f_rec = open(stored_rec, 'w')
        pickle.dump(recorder, f_rec)
        f_rec.close()
        del recorder
    except:
        print 'evaluation failed'
        
def save_sensitivity_outputs(quantitative_parameters, qualitative_parameters, weather,
                             scenario = {'year':2012},
                             variety = 'Tremie12',
                             filename = 'output_morris_tremie_12.csv'):
    list_param_names = ['i_sample'] + qualitative_parameters.keys() + quantitative_parameters.keys()
    df_in = pd.read_csv('septo_morris_input_full.txt', sep=' ', index_col=0, names=list_param_names)
    df_in = df_in[df_in[scenario.keys()[0]]==scenario.values()[0]]
    for k,v in qualitative_parameters.iteritems():
        if k!=scenario.values()[0]:
            df_in = df_in[df_in[k]==v['default']]
    i_samples = df_in.index

    df_out = pd.DataFrame(columns = ['num_leaf_top'] + list(df_in.columns) + ['normalized_audpc_pycnidia_coverage', 
                                                                              'normalized_audpc_necrosis_percentage',
                                                                              'max_pycnidia_coverage', 
                                                                              'max_necrosis_percentage',
                                                                              'speed', 'date_thr'])
    for i_sample in i_samples:
        stored_rec = './'+variety.lower()+'/recorder_'+str(int(i_sample))+'.pckl'
        recorder = get_recorder(stored_rec)
        df_reco = recorder.data
        df_reco = df_reco.rename(columns = {'num_plant':'plant', 'date':'datetime'})
        df_dates = get_date_threshold(df_reco, weather, variable = 'necrosis_percentage', threshold = 0.2).mean()
        df_speed = get_speed(df_reco, weather, variable = 'necrosis_percentage').mean()
        for lf in np.unique(df_reco['num_leaf_top']):
            df_reco_lf = df_reco[(df_reco['num_leaf_top']==lf)]
            output = {}
            output['num_leaf_top'] = lf
            for col in df_in.columns:
                output[col] = df_in.loc[i_sample, col]
            output['normalized_audpc_pycnidia_coverage'] = df_reco_lf.normalized_audpc_pycnidia_coverage[df_reco_lf['normalized_audpc_pycnidia_coverage'].apply(np.isreal)].mean()
            output['normalized_audpc_necrosis_percentage'] = df_reco_lf.normalized_audpc_necrosis_percentage[df_reco_lf['normalized_audpc_necrosis_percentage'].apply(np.isreal)].mean()
            output['max_pycnidia_coverage'] = min(df_reco_lf.groupby('plant').max()['pycnidia_coverage'].mean(), 100.)
            output['max_necrosis_percentage'] = min(df_reco_lf.groupby('plant').max()['necrosis_percentage'].mean(), 100.)
            output['speed'] = df_speed[lf] if lf in df_speed.index else np.nan
            output['date_thr'] = df_dates[lf] if lf in df_dates.index else np.nan
            df_out = df_out.append(output, ignore_index = True)
        del recorder
    df_out.to_csv(filename)

variety_code = {1:'Mercia', 2:'Rht3', 3:'Tremie12', 4:'Tremie13'}

quantitative_parameters = OrderedDict([('sporulating_fraction', [0., 2.e-4]),
                                       ('degree_days_to_chlorosis', [130., 300.]),
                                       ('Smin', [0.01, 0.09]),
                                       ('Smax', [0.1, 1.0]),
                                       ('growth_rate', [0.0001, 0.001]),
                                       ('age_physio_switch_senescence', [0.01, 1.]),
                                       ('sporulating_capacity', [0., 1.]),
                                       ('density_dus_emitted', [1e3, 3e3]),
                                       ('reduction_by_rain', [0., 1.]), 
                                       ('proba_inf', [0., 1.]), 
                                       ('loss_delay', [72., 168.]), 
                                       ('temp_min', [0., 10.])])

#qualitative_parameters = OrderedDict([('year', {'default':2004., 'values':[1998., 2003., 2004.]}),
#                                      ('variety', {'default':1, 'values':[1, 2, 3, 4]})])

qualitative_parameters = OrderedDict([('year', {'default':2012, 'values':[2012]})])

list_param_names = qualitative_parameters.keys() + quantitative_parameters.keys()

generate_parameter_set(quantitative_parameters,
                       qualitative_parameters,
                       parameter_range_file = 'param_range_SA.txt',
                       sample_file = 'septo_morris_input.txt',
                       num_trajectories = 10,
                       num_levels = 10)

nb_cpu = cpu_count()
filename = 'septo_morris_input_full.txt'
param_values = np.loadtxt(filename, delimiter=' ').tolist()
samples = map(lambda x: param_values_to_dict(x, list_param_names), param_values)

# Run disease simulation
if __name__ == '__main__':
    pymap(annual_loop, samples, nb_cpu)