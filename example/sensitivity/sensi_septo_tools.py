""" Functions used for sensitivity analysis of septoria model
"""
from disease_sensi_morris import *
from alinea.alep.disease_outputs import (AdelWheatRecorder,
                                         get_synthetic_outputs_by_leaf)
from alinea.alep.simulation_tools.septo_decomposed import annual_loop_septo


### Run and save simulation
def variety_code():
    return {'Mercia': 1, 'Rht3': 2, 'Tremie12': 3, 'Tremie13': 4, 'Custom': 5}


def variety_decode():
    return {v: k for k, v in variety_code().items()}


def get_stored_rec(variety, year, i_sample, i_boot):
    if variety == 'Custom':
        return './septo_wheat/' + variety.lower() + '_' + str(year) + \
               '_' + str(int(i_sample)) + '_boot' + str(int(i_boot)) + '.csv'
    else:
        return './septoria/' + variety.lower() + '_' + str(year) + \
               '_' + str(int(i_sample)) + '_boot' + str(int(i_boot)) + '.csv'


def run_septoria(sample):
    i_sample = sample.pop('i_sample')
    print('------------------------------')
    print('i_sample %d' % i_sample)
    print('------------------------------')
    i_boot = sample.pop('i_boot')
    year = int(sample.pop('year'))
    variety = variety_decode()[sample.pop('variety')]
    sowing_date = '10-15'
    if year == 2012:
        sowing_date = '10-21'
    elif year == 2013:
        sowing_date = '10-29'
    output_file = get_stored_rec(variety, year, i_sample, i_boot)
    annual_loop_septo(year=year, variety=variety, sowing_date=sowing_date,
                      nplants=15, output_file=output_file, **sample)


def run_custom_septoria(sample):
    i_sample = sample.pop('i_sample')
    print('------------------------------')
    print('i_sample %d' % i_sample)
    print('------------------------------')
    i_boot = sample.pop('i_boot')
    year = int(sample.pop('year'))
    variety = variety_decode()[sample.pop('variety')]
    sowing_date = '10-15'
    if year == 2012:
        sowing_date = '10-21'
    elif year == 2013:
        sowing_date = '10-29'
    output_file = get_stored_rec(variety, year, i_sample, i_boot)
    annual_loop_septo(year=year, variety=variety, sowing_date=sowing_date,
                      nplants=15, output_file=output_file,
                      **sample)


def get_septo_morris_path(year=2012, variety='Tremie12', folder='septoria'):
    return './' + folder + '/septo_morris_output_' + variety.lower() + '_' + str(year) + '.csv'


def save_sensitivity_outputs(year=2012, variety='Tremie12',
                             parameter_range_file='septo_param_range.txt',
                             input_file='septo_morris_input_full.txt',
                             folder='septoria',
                             nboots=5):
    parameter_names = pd.read_csv('./' + folder + '/' + parameter_range_file, header=None, sep=' ')[0].values.tolist()
    list_param_names = ['i_sample', 'i_boot', 'year', 'variety'] + parameter_names
    df_out = pd.DataFrame()
    for boot in range(nboots):
        input_boot = './' + folder + '/' + input_file[:-9] + '_boot' + str(boot) + input_file[-9:]
        df_in = pd.read_csv(input_boot, sep=' ',
                            index_col=0, names=list_param_names)
        vc = variety_code()
        df_in = df_in[(df_in['year'] == year) & (df_in['variety'] == vc[variety])]
        i_samples = df_in.index
        df_out_b = pd.DataFrame(
            columns=['num_leaf_top'] + list(df_in.columns) + ['normalized_audpc', 'audpc', 'audpc_500'])
        for i_sample in i_samples:
            stored_rec = get_stored_rec(variety, year, i_sample, boot)
            df_reco = pd.read_csv(stored_rec)
            cols_to_del = np.unique([col for col in df_reco.columns if col.startswith('Unnamed')])  # UGLY HACK
            for col in cols_to_del:
                df_reco = df_reco.drop(col, 1)
            df_s = get_synthetic_outputs_by_leaf(df_reco)
            for lf in np.unique(df_reco['num_leaf_top']):
                df_reco_lf = df_reco[(df_reco['num_leaf_top'] == lf)]
                output = {}
                output['num_leaf_top'] = lf
                for col in df_in.columns:
                    output[col] = df_in.loc[i_sample, col]
                # output['normalized_audpc'] = df_s[df_s['num_leaf_top']==lf].normalized_audpc.astype(float).values[0]
                output['audpc'] = df_s[df_s['num_leaf_top'] == lf].audpc.astype(float).values[0]
                # output['audpc_500'] = df_reco_lf.audpc_500.mean()
                output['max_severity'] = df_s[df_s['num_leaf_top'] == lf].max_severity.astype(float).values[0]
                df_out_b = df_out_b.append(output, ignore_index=True)
        df_out_b['i_boot'] = boot
        df_out = pd.concat([df_out, df_out_b])
    output_file = get_septo_morris_path(year=year, variety=variety, folder=folder)
    df_out.to_csv(output_file)


def plot_septo_morris_by_leaf(year=2012, variety='Tremie12',
                              variable='normalized_audpc',
                              parameter_range_file='septo_param_range.txt',
                              input_file='septo_morris_input.txt',
                              nboots=5, ylims=None, force_rename={}, markers_SA={},
                              folder='septoria'):
    output_file = get_septo_morris_path(year=year, variety=variety, folder=folder)
    df_out = pd.read_csv(output_file)
    plot_morris_by_leaf(df_out, variable=variable,
                        parameter_range_file=parameter_range_file,
                        input_file=input_file, nboots=nboots,
                        ylims=ylims, force_rename=force_rename,
                        markers_SA=markers_SA, folder=folder)


def plot_septo_morris_by_leaf_by_boot(year=2012, variety='Tremie12',
                                      variable='normalized_audpc',
                                      parameter_range_file='./septoria/septo_param_range.txt',
                                      input_file='./septoria/septo_morris_input.txt',
                                      nboots=5, ylims=None):
    output_file = get_septo_morris_path(year=year, variety=variety)
    df_out = pd.read_csv(output_file)
    plot_morris_by_leaf_by_boot(df_out, variable=variable,
                                parameter_range_file=parameter_range_file,
                                input_file=input_file, nboots=nboots, ylims=ylims)


def plot_septo_morris_3_leaves(year=2012, variety='Tremie12',
                               leaves=[10, 5, 1], variable='normalized_audpc',
                               parameter_range_file='./septoria/septo_param_range.txt',
                               input_file='./septoria/septo_morris_input.txt',
                               nboots=5, ylims=None, force_rename={}, axs=None):
    output_file = get_septo_morris_path(year=year, variety=variety)
    df_out = pd.read_csv(output_file)
    plot_morris_3_leaves(df_out, leaves=leaves, variable=variable,
                         parameter_range_file=parameter_range_file,
                         input_file=input_file, nboots=nboots, ylims=ylims,
                         force_rename=force_rename, axs=axs)


def septo_scatter_plot_by_leaf(year=2012, variety='Tremie12',
                               variable='normalized_audpc',
                               parameter='Smax'):
    output_file = get_septo_morris_path(year=year, variety=variety)
    df_out = pd.read_csv(output_file)
    scatter_plot_by_leaf(df_out, variable=variable, parameter=parameter)


def septo_boxplot_by_leaf(year=2012, variety='Tremie12',
                          variable='normalized_audpc',
                          parameter='Smax', ylims=None):
    output_file = get_septo_morris_path(year=year, variety=variety)
    df_out = pd.read_csv(output_file)
    boxplot_by_leaf(df_out, variable=variable,
                    parameter=parameter, ylims=ylims)


def septo_boxplot_3_leaves(year=2012, leaves=[10, 5, 1],
                           variety='Tremie12',
                           variable='normalized_audpc',
                           parameter='Smax', ylims=None):
    output_file = get_septo_morris_path(year=year, variety=variety)
    df_out = pd.read_csv(output_file)
    boxplot_3_leaves(df_out, leaves=leaves, variable=variable,
                     parameter=parameter, ylims=ylims)


def force_rename_SA_wheat():
    return {'tiller_probability': r"$\mathit{Tiller}$",
            'proba_main_nff': r"$\mathit{FLN}$",
            'scale_HS': r"$\mathit{Phyllochron}$",
            'scale_leafDim_length': r"$\mathit{Length}_{leaf}$",
            'scale_leafDim_width': r"$\mathit{Width}_{leaf}$",
            'scale_leafRate': r"$\mathit{Elongation}_{leaf}$",
            'scale_stemDim': r"$\mathit{Length}_{stem}$",
            'scale_stemRate': r"$\mathit{Elongation}_{stem}$",
            'scale_fallingRate': r"$\mathit{Curvature}_{leaf}$",
            'scale_leafSenescence': r"$\mathit{Senescence}_{leaf}$"}


def markers_SA():
    return {'scale_HS': 'o', 'scale_leafSenescence': '^',
            'scale_stemDim': 's', 'scale_stemRate': 'p', 'scale_tillering': '*',
            'scale_leafDim_length': 'd', 'scale_leafDim_width': 'D',
            'scale_leafRate': '<', 'scale_fallingRate': '>',
            'proba_main_nff': 'v'}


def plot_morris_3_leaves_2_years(years=[2011, 2013], variety='Custom',
                                 leaves=[10, 5, 1], variable='audpc',
                                 parameter_range_file='2011/septo_param_range.txt',
                                 input_files=['2011/septo_morris_input.txt',
                                              '2013/septo_morris_input.txt'],
                                 nboots=5, ylims=None, force_rename={},
                                 axs=None, save_fig=True, folder='septo_wheat'):
    fig, axs = plt.subplots(2, 3, figsize=(10, 6))
    for yr, ax, input_file in zip(years, axs, input_files):
        df_out = pd.read_csv(get_septo_morris_path(year=yr, variety='Custom', folder=folder))
        sfx = '- %d' % yr
        plot_morris_3_leaves(df_out, leaves=leaves, variable=variable,
                             parameter_range_file=parameter_range_file,
                             input_file=input_file, nboots=nboots, ylims=ylims,
                             force_rename=force_rename,
                             axs=ax, annotation_suffix=sfx,
                             folder=folder)
    if save_fig:
        fig.savefig('morris_3_leaves_2_years', bbox_inches='tight')


def septo_boxplot_3_leaves_3_params(year=2013, leaves=[10, 5, 1],
                                    variety='Custom',
                                    variable='audpc',
                                    parameters=['scale_HS', 'scale_leafSenescence', 'scale_stemRate'],
                                    ylims=None):
    fig, axs = plt.subplots(3, 3, figsize=(15, 6))
    output_file = get_septo_morris_path(year=year, variety=variety)
    df_out = pd.read_csv(output_file)
    for parameter, ax in zip(parameters, axs):
        boxplot_3_leaves(df_out, leaves=leaves, variable=variable,
                         parameter=parameter, ylims=ylims, axs=ax)


def plot_septo_morris_two_years_two_leaves(years=[2011, 2013], leaves=[1, 2],
                                           variety='Custom',
                                           variable='audpc',
                                           parameter_range_file='septo_param_range.txt',
                                           input_file='septo_morris_input.txt',
                                           nboots=5, ylims=None, force_rename=force_rename_SA_wheat(),
                                           markers_SA=markers_SA(),
                                           folder='septoria'):
    # import matplotlib.pyplot as plt
    plt.rcParams['text.usetex'] = True
    fig, axs = plt.subplots(2, 2, figsize=(20, 15))
    for i_year, year in enumerate(years):
        problem = read_param_file('./' + folder + '/' + str(year) + '/' + parameter_range_file)
        output_file = get_septo_morris_path(year=year, variety=variety, folder=folder)
        df_out = pd.read_csv(output_file)
        axs_ = axs[i_year]
        for i_leaf, leaf in enumerate(leaves):
            ax = axs_[i_leaf]
            if ax == axs[0][-1]:
                lgd = plot_morris_one_leaf(df_out, problem,
                                           leaf=leaf, ax=ax,
                                           variable=variable,
                                           input_file=input_file,
                                           nboots=nboots, ylims=ylims,
                                           force_rename=force_rename,
                                           markers_SA=markers_SA,
                                           add_legend=True,
                                           annotation_suffix=' - ' + str(year),
                                           colored=True,
                                           folder=folder + '/' + str(year),
                                           markersize=20,
                                           return_legend=True)
            else:
                plot_morris_one_leaf(df_out, problem,
                                     leaf=leaf, ax=ax,
                                     variable=variable,
                                     input_file=input_file,
                                     nboots=nboots, ylims=ylims,
                                     force_rename=force_rename,
                                     markers_SA=markers_SA,
                                     add_legend=False,
                                     annotation_suffix=' - ' + str(year),
                                     colored=True,
                                     folder=folder + '/' + str(year),
                                     markersize=20,
                                     return_legend=False)

    plt.savefig('Figure4.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.rcParams['text.usetex'] = False


if __name__ == '__main__':
    plot_septo_morris_two_years_two_leaves()
    plt.show()
