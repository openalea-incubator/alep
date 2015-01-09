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
    
def septo_infection_filter(seq, weather, rain_filter, degree_days=20., base_temp = 0., start_date=None):
    if not 'septo_infection_risk_with_event' in weather.data.columns:
        from alinea.alep.alep_weather import add_septoria_risk_with_event
        weather.check(varnames=['septo_infection_risk_with_event'], models={'septo_infection_risk_with_event':add_septoria_risk_with_event})
    if not 'septo_degree_days' in weather.data.columns:
        from alinea.alep.alep_weather import linear_degree_days
        weather.check(varnames=['septo_degree_days'], models={'septo_degree_days':linear_degree_days}, start_date=start_date, base_temp=base_temp)
        
    cond_inf = weather.data['septo_infection_risk_with_event'][seq].values
    ddays = weather.data['septo_degree_days'][seq].values
    ind = []
    count=0
    for i in range(len(cond_inf)):
        if (rain_filter[i]==True and round(weather.data['rain'][seq][i], 14)>0) or ('notation_dates' in weather.data.columns and weather.data['notation_dates'][seq][i] == True):
            cond_inf[i] = True
            count = ddays[i]
        else:
            if i == 0:
                count += ddays[i]
            else:
                count += (ddays[i] - ddays[i-1])
            if count >= degree_days:
                cond_inf[i] = True
                count = (ddays[i] - ddays[i-1])
    return cond_inf==1