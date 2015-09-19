""" Functions used for sensitivity analysis of septoria model
"""
from disease_sensi_morris import *
from alinea.alep.disease_outputs import (AdelWheatRecorder,
                                        get_synthetic_outputs_by_leaf)
from alinea.alep.simulation_tools.septo_decomposed import annual_loop_septo

### Run and save simulation
class SeptoSensiRecorder(AdelWheatRecorder):
    """ Record simulation output on every leaf of main stems in a dataframe during simulation """
    def __init__(self, group_dus = True, 
                 fungus_name = 'septoria', 
                 increment = 1000):
        super(AdelSeptoRecorder, self).__init__(group_dus = group_dus, 
                                                fungus_name = fungus_name,
                                                increment = increment)
        columns = ['date', 'degree_days', 'num_plant', 'num_leaf_bottom', 'leaf_area', 
                   'leaf_green_area', 'fnl', 'surface_spo', 'surface_spo_on_green', 
                   'surface_empty', 'surface_empty_on_green', 'surface_dead']
        self.data = pandas.DataFrame(data = [[np.nan for col in columns] for i in range(self.increment)], 
                                     columns = columns)
    
    def get_values_single_leaf(self, g, date, degree_days, id_list):
        dict_lf = {}
        dict_lf['date'] = date
        dict_lf['degree_days'] = degree_days
        
        # Update leaf properties
        areas = g.property('area')
        green_areas = g.property('green_area')
        lengths = g.property('length')
        senesced_lengths = g.property('senesced_length')
        fnls = g.property('nff')
        a_label_splitted = self.a_labels[id_list[0]].split('_')
        dict_lf['num_plant'] = int(a_label_splitted[0].split('plant')[1])
        dict_lf['num_leaf_bottom'] = int(a_label_splitted[2].split('metamer')[1])
        dict_lf['leaf_area'] = sum([areas[id] for id in id_list])
        dict_lf['leaf_green_area'] = sum([green_areas[id] for id in id_list])
        dict_lf['fnl'] =  fnls[g.complex_at_scale(id_list[0], 2)]

        # Update properties of dispersal units and lesions
        surface_spo = 0.
        surface_spo_on_green = 0.
        surface_empty = 0.
        surface_empty_on_green = 0.
        surface_dead = 0.
        
        for id in id_list:
            leaf = g.node(id)
            if 'dispersal_units' in leaf.properties():
                for du in leaf.dispersal_units:
                    if du.fungus.name == self.fungus_name:
                        if self.group_dus:
                            nb_dus += du.nb_dispersal_units
                        else:
                            nb_dus += 1
                                
            if 'lesions' in leaf.properties():
                for les in leaf.lesions:
                    if les.fungus.name == self.fungus_name:
                        if self.group_dus:
                            nb_les = les.nb_lesions
                            nb_les_on_green = les.nb_lesions_non_sen
                            ratio_green = float(nb_les_on_green)/nb_les if nb_les>0. else 0.
                            surface_nec_on_green += les.surface_nec *  ratio_green
                            surface_spo_on_green += les.surface_spo * ratio_green
                            surface_empty_on_green += les.surface_empty * ratio_green
                        else:
                            nb_lesions += 1
                            if les.position[0][0]>leaf.senesced_length:
                                nb_lesions_on_green += 1
                                surface_nec_on_green = les.surface_nec_on_green
                                surface_spo_on_green = les.surface_spo
                                surface_empty_on_green = les.surface_empty
                        surface_spo += les.surface_spo
                        surface_empty += les.surface_empty
                        surface_dead += les.surface_dead
        dict_lf['surface_spo'] = surface_spo
        dict_lf['surface_empty'] = surface_empty
        dict_lf['surface_dead'] = surface_dead
        dict_lf['surface_nec_on_green'] = surface_nec_on_green
        dict_lf['surface_spo_on_green'] = surface_spo_on_green
        dict_lf['surface_empty_on_green'] = surface_empty_on_green
        return dict_lf

    def leaf_necrotic_area(self):
        self.data['leaf_necrotic_area'] = self.data['surface_spo'] + \
                                          self.data['surface_empty']

    def leaf_necrotic_area_on_green(self):
        self.data['leaf_necrotic_area_on_green'] = self.data['surface_spo_on_green'] + \
                                                   self.data['surface_empty_on_green']

    def _ratio(self, variable='leaf_necrotic_area', against='leaf_area'):
        r = []
        for ind in self.data.index:
            a = self.data[against][ind]
            if a > 0.:
                r.append(self.data[variable][ind]/a)
            else:
                r.append(0.)
        return r
    
    def severity(self):
        """ Necrotic area of lesions compared to total leaf area """
        if not 'leaf_necrotic_area' in self.data:
            self.leaf_necrotic_area()
        self.data['severity'] = self._ratio(variable='leaf_necrotic_area',
                                            against='leaf_area')
        
    def severity_on_green(self):
        """ Necrotic area of lesions on green compared to green leaf area """
        if not 'leaf_necrotic_area_on_green' in self.data:
            self.leaf_necrotic_area_on_green()
        self.data['severity_on_green'] = self._ratio(variable='leaf_necrotic_area_on_green',
                                                     against='leaf_green_area')

    def get_audpc(self, variable='severity'):
        for pl in set(self.data['num_plant']):
            df_pl =  self.data[self.data['num_plant'] == pl]
            for lf in set(df_pl['num_leaf_top']):
                df_lf = df_pl[df_pl['num_leaf_top'] == lf]
                ind_data_lf = df_lf.index
                if round(df_lf['leaf_green_area'][pandas.notnull(df_lf['leaf_disease_area'])].iloc[-1],10)==0.:
                    data = df_lf[variable][df_lf['leaf_green_area']>0]
                    ddays = df_lf['degree_days'][df_lf['leaf_green_area']>0]
                    data_ref = numpy.ones(len(data))
                    if len(data[data>0])>0:
                        audpc = simps(data[data>0], ddays[data>0])
                        audpc_ref = simps(data_ref[data_ref>0], ddays[data_ref>0])
                        if numpy.isnan(audpc):
                            audpc = trapz(data[data>0], ddays[data>0])
                        if numpy.isnan(audpc_ref):
                            audpc_ref = trapz(data_ref[data_ref>0], ddays[data_ref>0])
                    else:
                        audpc = 0.
                        audpc_ref = 0.
                    self.data.loc[ind_data_lf, 'audpc'] = audpc
                    self.data.loc[ind_data_lf, 'normalized_audpc'] = audpc/audpc_ref if audpc_ref>0. else 0.
                else:
                    self.data.loc[ind_data_lf, 'audpc'] = np.nan
                    self.data.loc[ind_data_lf, 'normalized_audpc'] = np.nan    
    
    def post_treatment(self, variety = None):
        self.data = self.data[~pandas.isnull(self.data['date'])]
        self.add_leaf_numbers()
        self.leaf_necrotic_area()
        self.leaf_necrotic_area_on_green()
        self.severity()
        self.severity_on_green()
        self.get_audpc()
        if variety is not None:
            self.add_variety(variety=variety)
            
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
                        proba_inf=0.5, age_infection=True, Smin=0.02, 
                        density_dus_emitted_ref=1e5, growth_rate=0.0006,
                        degree_days_to_chlorosis=130., 
                        degree_days_to_necrosis=110.,
                        degree_days_to_sporulation=30.,
                        sporulating_fraction=7e-2, **sample)

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