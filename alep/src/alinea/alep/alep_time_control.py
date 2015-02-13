""" Strategies of time control used in alep projects. """
import numpy
import pandas
from datetime import datetime
from alinea.astk.TimeControl import *

def evaluation_sequence(delays, eval_time='start'):
    """ retrieve evaluation filter from sequence of delays
    """
    if eval_time == 'start':
        seq = [[True if i == 0 else False for i in range(int(d))] for d in delays]
    elif eval_time == 'end':
        seq = [[True if i == int(d)-1 else False for i in range(int(d))] for d in delays]
    return reduce(lambda x,y: x + y, seq)

class CustomIterWithDelays(IterWithDelays):

    def __init__(self, values = [None], delays = [1], eval_time='start'):
        super(CustomIterWithDelays, self).__init__(values=values, delays=delays)
        self.eval_time = eval_time
        self._evalseq = iter(evaluation_sequence(delays, eval_time))
        self.count = 0.

    def __iter__(self):
        return CustomIterWithDelays(self.values, self.delays, self.eval_time)
        
    def next(self):
        self.ev = self._evalseq.next()
        if self.ev : 
            try: #prevent value exhaustion to stop iterating
                self.val = self._iterable.next()
                self.dt = self._iterdelays.next()
                self.count += 1
            except StopIteration:
                pass
        if self.count == 0.:
            self.val=None
            self.dt = None
        return EvalValue(self.ev, self.val, self.dt)

def add_notation_dates(data, notation_dates_file):
    import re
    if data.index[-1].year != int(re.findall(r'\d+', notation_dates_file)[0]):
        raise NameError("Year of notation_dates_file does not correspond to weather data")
    df = pandas.read_csv(notation_dates_file)
    dates = [datetime.strptime(df['notation_dates'][d], '%d/%m/%Y %H:%M') for d in df.index]
    df2 = pandas.DataFrame(numpy.zeros(len(data)), index = data.index)
    df2.ix[dates, :] = 1
    return df2.values == 1
    
def septoria_filter(seq, weather, degree_days=10., base_temp = 0., 
                    Tmin = 10., Tmax = 30., WDmin = 10., rain_min = 0.2, start_date=None):
    if not 'septo_infection_risk_with_event' in weather.data.columns:
        from alinea.alep.alep_weather import add_septoria_risk_with_event
        weather.check(varnames=['septo_infection_risk_with_event'], 
                      models={'septo_infection_risk_with_event':add_septoria_risk_with_event}, 
                      wetness_duration_min = WDmin, temp_min = Tmin, temp_max = Tmax)
    if not 'septo_degree_days' in weather.data.columns:
        from alinea.alep.alep_weather import linear_degree_days
        weather.check(varnames=['septo_degree_days'], models={'septo_degree_days':linear_degree_days}, start_date=start_date, base_temp=base_temp)
    wetness = weather.data.wetness[seq]
    temperature = weather.data.temperature_air[seq]
    rain = weather.data.rain[seq]
    rain[rain <= rain_min] = 0
    cond_inf = weather.data.septo_infection_risk_with_event[seq]
    ddays = weather.data.degree_days[seq].diff()
    ddays.ix[0] = 0.
    count_wet = 0.
    count_rain = 0.
    count_ddays = 0.
    df = pandas.Series([False for i in range(len(wetness))], index = wetness.index)
    for i, row in df.iteritems():
        # if wetness[i] and Tmin <= temperature[i] < Tmax:
            # if count_wet == 0. and i < df.index[-1]-9 and cond_inf[i+9] == True:
                # df[i] = True
            # count_ddays = 0.
            # count_wet += 1.
            # if count_wet >= WDmin and cond_inf[i] == True:
                # df[i] = True
        # else:
            # count_wet = 0.
        if wetness[i] and Tmin <= temperature[i] < Tmax:
            if count_wet == 0. and i < df.index[-1]-9 and cond_inf[i+9] == True:
                df[i] = True
            count_wet += 1.
        else:
            if i > df.index[0] and count_wet >= WDmin and cond_inf[i-1] == True:
                df[i - 1] = True
            count_wet = 0.
    
    for i, row in df.iteritems():
        if rain[i] > 0:
            if count_rain == 0.:
                df[i] = True
            count_rain += 1
        else:
            if count_rain > 0.:
                df[i] = True
            count_rain = 0.

    for i, row in df.iteritems():
        if row == True:
            count_ddays = 0.
        count_ddays += ddays[i]
        if count_ddays >= degree_days and i > df.index[0]:
            df[i - 1] = True
            count_ddays = 0.
    df[df.index[0]] = True
    return df