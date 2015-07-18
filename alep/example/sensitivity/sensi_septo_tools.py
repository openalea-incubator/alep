""" Functions used for sensitivity analysis of septoria model
"""
from disease_sensi_morris import *
from septo_decomposed import run_disease, make_canopy
from alinea.alep.disease_outputs import get_recorder
try:
    import cPickle as pickle
except:
    import pickle
from alinea.alep.disease_outputs import *
from alinea.echap.disease.septo_data_treatment import (get_date_threshold,
                                                       get_speed)

### Run and save simulation
def variety_code():
    return {'Mercia':1, 'Rht3':2, 'Tremie12':3, 'Tremie13':4}
    
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
        stored_rec = './septoria/'+variety.lower()+'_'+str(year)+'_'+str(int(i_sample))+'.pckl'
        f_rec = open(stored_rec, 'w')
        pickle.dump(recorder, f_rec)
        f_rec.close()
        del recorder
    except:
        print 'evaluation failed'
        
def save_sensitivity_outputs(year = 2012,
                             variety = 'Tremie12',
                             parameter_range_file = 'septo_param_range.txt',
                             input_file = 'septo_morris_input_full.txt',
                             output_file = 'septo_morris_output_tremie_12.csv'):
    parameter_names = pd.read_csv(param_range_file,
                                  header=None, sep = ' ')[0].values
    list_param_names = ['i_sample', 'year', 'variety'] + parameter_names
    df_in = pd.read_csv(input_file, sep=' ',
                        index_col=0, names=list_param_names)
    vc = variety_code()
    df_in = df_in[(df_in['year']==year) & (df_in['variety']==vc[variety])]
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
