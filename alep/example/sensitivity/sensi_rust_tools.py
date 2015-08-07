""" Functions used for sensitivity analysis of septoria model
"""
from disease_sensi_morris import *
from alinea.alep.disease_outputs import AdelSeptoRecorder
from sensi_septo_tools import variety_code
from alinea.alep.simulation_tools.brown_rust_decomposed import annual_loop_rust
import numpy as np
import pandas as pd
from scipy.integrate import simps

class SensiRustRecorder(AdelSeptoRecorder):
    """ Record simulation output on every leaf of main stems in a dataframe during simulation """
    def __init__(self, group_dus = True, fungus_name = 'brown_rust',
                 increment = 1000):
        super(AdelSeptoRecorder, self).__init__(group_dus = group_dus, 
                                                fungus_name = fungus_name,
                                                increment = increment)
        columns = ['date', 'degree_days', 'num_plant', 
                   'num_leaf_bottom', 'fnl',
                   'leaf_area', 'leaf_disease_area']
        self.data = pandas.DataFrame(data = [[np.nan for col in columns] 
                                    for i in range(1000)], columns = columns)
    
    def get_values_single_leaf(self, g, date, degree_days, id_list):
        dict_lf = {}
        dict_lf['date'] = date
        dict_lf['degree_days'] = degree_days
        
        # Update leaf properties
        fnls = g.property('nff')
        areas = g.property('area')
        a_label_splitted = self.a_labels[id_list[0]].split('_')
        dict_lf['num_plant'] = int(a_label_splitted[0].split('plant')[1])
        dict_lf['num_leaf_bottom'] = int(a_label_splitted[2].split('metamer')[1])
        dict_lf['fnl'] =  fnls[g.complex_at_scale(id_list[0], 2)]
        dict_lf['leaf_area'] = sum([areas[id] for id in id_list])

        # Update properties of dispersal units and lesions
        surface = 0.
        for id in id_list:
            leaf = g.node(id)
            if 'lesions' in leaf.properties():
                for les in leaf.lesions:
                    if les.fungus.name == self.fungus_name:
                        surface += les.surface_alive      
        dict_lf['leaf_disease_area'] = surface
        return dict_lf

    def severity(self):
        if not 'leaf_disease_area' in self.data:
            self.leaf_disease_area()
        self.data['severity'] = [self.data['surface'][ind]/self.data['leaf_area'][ind] if self.data['leaf_area'][ind]>0. else 0. for ind in self.data.index]
    
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
        self.severity()
        self.get_audpc()
        if variety is not None:
            self.add_variety(variety=variety)           

def get_stored_rec(variety, year, i_sample, i_boot):
    return './brown_rust/'+variety.lower()+'_'+ str(year)+ \ 
            '_'+str(int(i_sample))+'_boot'+str(i_boot)+'.csv'
            
def run_brown_rust(sample):
    i_sample = sample.pop('i_sample')
    i_boot = sample.pop('i_boot')
    year = sample.pop('year')
    variety = sample.pop('variety')
        sowing_date = '10-15'
    if year == 2012:
        sowing_date = '10-21'
    elif year == 2013:
        sowing_date = '10-29'
    output_file = get_stored_rec(variety, year, i_sample, i_boot)
    annual_loop_rust(year = year, variety = variety, sowing_date=sowing_date,
                    nplants = 5, output_file = output_file, **kwds)

def get_rust_morris_path(year = 2012, variety = 'Tremie12', nboots = 5):
    return './brown_rust/rust_morris_output_'+variety.lower()+'_'+str(year)+'.csv'
                    
def save_sensitivity_outputs(year = 2012, variety = 'Tremie12',
                             parameter_range_file = 'rust_param_range.txt',
                             input_file = './brown_rust/rust_morris_input_full.txt',
                             nboots = 5):
    parameter_names = pd.read_csv(parameter_range_file, header=None, sep = ' ')[0].values.tolist()
    list_param_names = ['i_sample', 'year', 'variety'] + parameter_names
    for boot in range(nboots):
        input_boot = input_file[:-9]+'_boot'+str(boot)+input_file[-9:]
        df_in = pd.read_csv(input_boot, sep=' ',
                            index_col=0, names=list_param_names)
        vc = variety_code()
        df_in = df_in[(df_in['year']==year) & (df_in['variety']==vc[variety])]
        i_samples = df_in.index
        df_out = pd.DataFrame(columns = ['num_leaf_top'] + list(df_in.columns) + ['normalized_audpc'])
        for i_sample in i_samples:
            stored_rec = get_stored_rec(variety, year, i_sample, boot)
            df_reco = pd.read_csv(stored_rec)
            for lf in np.unique(df_reco['num_leaf_top']):
                df_reco_lf = df_reco[(df_reco['num_leaf_top']==lf)]
                output = {}
                output['num_leaf_top'] = lf
                for col in df_in.columns:
                    output[col] = df_in.loc[i_sample, col]
                output['normalized_audpc'] = df_reco_lf.normalized_audpc.mean()
                df_out = df_out.append(output, ignore_index = True)
        output_file = get_rust_morris_path(year=year, variety=variety, boot=boot)
        df_out.to_csv(output_file)
    
def plot_septo_morris_by_leaf(year = 2012, variety = 'Tremie12',
                             parameter_range_file = './brown_rust/rust_param_range.txt',
                             input_file = './brown_rust/rust_morris_input_full.txt'):
    output_file = get_rust_morris_path(year=year, variety=variety)
    df_out = pd.read_csv(output_file)
    plot_morris_by_leaf(df_out, variable=variable,
                        parameter_range_file=parameter_range_file,
                        input_file=input_file)