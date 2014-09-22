""" Weather utilities for alep.

Gather methods concerning weather management specific to alep functioning.
"""
from pylab import *
import pandas
import numpy
from alinea.weather.global_weather import *
from alinea.weather.mini_models import leaf_wetness_rapilly
from matplotlib.ticker import FuncFormatter
from datetime import timedelta, date

def is_raining(rain_eval):
    """ Check if it is raining or not
    
    Parameters
    ----------
    rain_eval: TimeControl EvalValue instance
        Weather data divided according to rain occurences
    
    Returns
    -------
    True or False
    """
    if rain_eval:
        if rain_eval.value.rain.sum() > 0:
            return True
        else:
            return False
    else:
        return False

def get_septoria_weather(data_file='meteo00-01.txt', sep = ';'):
    """ Read weather data and adapt it to septoria model (add wetness)
    
    Parameters
    ----------
    data_file: str
        Name of weather data file (in txt)
    
    Returns
    -------
    weather: pandas dataframe
        Dataframe of weather indexed by date and with explicitely named columns
        See `astk.Weather`
    """
    meteo_path = shared_data(alinea.septo3d, 'meteo00-01.txt')
    weather = Weather(data_file=meteo_path)
    weather.check(varnames=['wetness'], models={'wetness':wetness_rapilly})
    
    # weather = Weather(data_file=data_file, sep = sep)
    # weather = add_wetness(weather)
    # weather = add_rain_dispersal_events(weather)
    return weather

def wetness_rapilly(data):
    """ Compute leaf wetness as in Rapilly et Jolivet, 1976 as a 
        function of rain or relative humidity and PAR.
    
    Parameters
    ----------
    data: pandas dataframe 
        Weather data as read by Weather class
    
    Returns
    -------
    sequence of True or False
        True if the leaf is wet, False oherwise
    """
    return np.greater(data[['rain']].values, 0.) + np.less(data[['PPFD']].values,644.) * np.greater_equal(data[['relative_humidity']].values, 85.)

def mean_rain(data):
    """ Create a sequence of rain events with only mean values."""
    pass
    
def add_wetness(weather):
    """ Complete weather data with wetness.
    """
    wet = dict(wetness=[])
    for i_line in range(len(weather.data)):
        wet['wetness'].append(leaf_wetness_rapilly(weather.data.rain[i_line], 
                              weather.data.relative_humidity[i_line],
                              weather.data.PPFD[i_line]))
    wetness = pandas.DataFrame(wet)
    wetness.index = weather.data.rain.index
    weather.data = weather.data.join(wetness)
    return weather
    
def linear_degree_days(data, start_date, base_temp=0, max_temp=25.):
    df = data['temperature_air'].copy()
    df[df<base_temp]=0.
    df[df>max_temp]=0.
    dd = np.zeros(len(df))
    if isinstance(start_date, str):
        start_date = datetime.strptime(start_date, '%Y-%m-%d %H:%M:%S')
    ind_start = len(df.ix[:start_date+timedelta(0,60)])
    seq = pd.date_range(start=df.index[0], end=start_date-timedelta(0,60), freq='H')
    seq = seq.order(ascending=False)
    dd[ind_start:]=np.cumsum((df.ix[start_date+timedelta(0,60):].values-base_temp)/24.)
    dd[:len(seq)]=-np.cumsum((df.ix[seq].values[::-1]-base_temp)/24.)[::-1]
    return dd

def add_julian_days(data, sowing_date="2010-10-15"):
    indexes = np.array([ind.strftime('%Y-%m-%d') for ind in data.index])
    ind_sowing = np.where(indexes==sowing_date)[0][0]
    return [(data.index[t]-data.index[ind_sowing]).days for t in range(len(data.index))]
        
def add_septoria_infection_risk(data, temp_min=10, temp_max=25):
    """ Add True or False if there is a risk of infection for septoria.
    
    Infection is possible if wetness duration > 10h and 10<temp<30 deg C. 
    """
    temp1 = data.temperature_air>=temp_min
    temp2 = data.temperature_air<=temp_max
    infect_cond = data.wetness * temp1 * temp2
    septo_infection_risk = np.zeros(len(data))
    counter = 0.
    for i_line in range(len(data)):
        if infect_cond[i_line]==True:
            counter += 1.
            if counter >= 10.:
                septo_infection_risk[i_line] = True
            else:
                septo_infection_risk[i_line] = False
        else:
            counter = 0.
            septo_infection_risk[i_line] = False   
    return septo_infection_risk
    
def plot_septo_infection_risk(weather, start_date="2010-10-15 12:00:00", adel=None, axis = 'degree_days', ax = None, xlims=None, only_with_event=True):
    def form_tick(x, pos):
        t = date.fromordinal(int(x))
        return t.strftime('%b')+'\n'+str(int(weather.get_variable('degree_days', t)))
    
    if ax == None:
        fig, ax = plt.subplots(figsize=(8,1.8))
    
    if not 'degree_days' in weather.data.columns:
        weather.check(varnames=['degree_days'], models={'degree_days':linear_degree_days}, start_date=start_date, base_temp=0., max_temp=30.)
    if not 'wetness' in weather.data.columns:
        weather.check(varnames=['wetness'], models={'wetness':wetness_rapilly})
    if not 'septo_infection_risk' in weather.data.columns:
        weather.check(varnames=['septo_infection_risk'], models={'septo_infection_risk':add_septoria_infection_risk})
        
    if axis == 'degree_days':
        index = weather.data.degree_days
    elif axis == 'date':
        index = weather.data.index  
    if only_with_event==True:
        if not 'septo_risk_with_event' in weather.data.columns:
            weather.check(varnames=['septo_risk_with_event'], models={'septo_risk_with_event':add_septoria_risk_with_event})
        ax.vlines(index, [0], weather.data.septo_risk_with_event, 'k', alpha=0.1)
    else:
        ax.vlines(index, [0], weather.data.septo_infection_risk, 'k', alpha=0.1)
        
    if adel!=None:
        df = adel.phenT()
        emergence_dates = [numpy.mean(df[df['n']==lf]['tip']) for lf in range(int(max(df['n'])))]
        for lf in range(int(max(df['n']))):
            ax.annotate('', xy=(numpy.mean(df[df['n']==lf]['tip']), 0.7), xycoords='data',
                        xytext=(numpy.mean(df[df['n']==lf]['tip']), 1.), arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color='g'))
    
    ax.set_yticks([])
    ax.set_ylim([0,1])
    ax.set_title(str(weather.data.index[0].year)+'-'+str(weather.data.index[-1].year))
    ax.set_xlabel('Degree days', fontsize=16)
    if axis == 'date':
        formatter = FuncFormatter(form_tick)
        ax.xaxis.set_major_formatter(FuncFormatter(formatter))
    if xlims!=None:
        ax.set_xlim(xlims)       
    plt.tight_layout()
    return ax
    
def add_septoria_efficient_rain(data, viability = 5):
    """ viability (in days) is the period after which inoculum is dead if no infection. """
    df = data.ix[:, ['rain', 'septo_infection_risk']]
    df['first_risk'] = np.append(np.zeros(1), data['septo_infection_risk'].values[1:]-data['septo_infection_risk'].values[:-1])
    df['efficient_rain'] = np.zeros(len(df))
    for date, row in df.iterrows():
        if row['rain']>0 and sum(df.ix[date:date+timedelta(viability,0), 'first_risk'])>0:
            df.set_value(date, 'efficient_rain', 1)
    return df['efficient_rain']
    
def add_septoria_risk_with_event(data, viability = 5):
    """ viability (in days) is the period after which inoculum is dead if no infection. """
    df = data.ix[:, ['rain', 'septo_infection_risk']]
    df['bool_rain'] = (df['rain']>0)*1
    df['septo_risk_with_event'] = np.zeros(len(df))
    for date, row in df.iterrows():
        if row['septo_infection_risk']>0 and sum(df.ix[date-timedelta(days=viability-1, hours=23):date+timedelta(hours=1), 'bool_rain'])>0:
            df.set_value(date, 'septo_risk_with_event', 1)   
    return df['septo_risk_with_event']
    
def add_rain_dispersal_events(weather):
    """ Add a column to indicate the dispersal events at hours with max rain for each dispersal event,
        Add a column to indicate the value of average rain during dispersal event,
        Add a column to indicate the duration of each dispersal event.
    """
    max_rain = 0.
    ind_max = 0.
    rain_counter = 0.
    rain_amount = 0.
    
    dispersal = zeros(len(weather.data))
    rain = zeros(len(weather.data))
    rain_duration = zeros(len(weather.data))
    for i_line in range(len(weather.data)):
        if weather.data.rain[i_line] > 0.:
            rain_counter += 1.
            rain_amount += weather.data.rain[i_line]
            if weather.data.rain[i_line] > max_rain:
                max_rain = weather.data.rain[i_line]
                ind_max = i_line
        else:
            if ind_max > 0.:
                rain[ind_max] = rain_amount/rain_counter
                rain_duration[ind_max] = rain_counter
                dispersal[ind_max] = 1
                max_rain = 0.
                ind_max = 0.
                rain_counter = 0.
                rain_amount = 0.
    
    del weather.data['rain']
    rain = pandas.DataFrame(dict(rain=rain))
    rain.index = weather.data.temperature_air.index
    weather.data = weather.data.join(rain)
    
    rain_duration = pandas.DataFrame(dict(rain_duration=rain_duration))
    rain_duration.index = weather.data.temperature_air.index
    weather.data = weather.data.join(rain_duration)
    
    dispersal_event = pandas.DataFrame(dict(dispersal_event=dispersal))
    dispersal_event.index = weather.data.temperature_air.index
    weather.data = weather.data.join(dispersal_event)
    return weather