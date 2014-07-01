""" Strategies of time control used in alep projects. """
import numpy
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

def septo_infection_filter(seq, weather, rain_filter, degree_days=20, start_date=None):
    if not 'septo_infection_risk' in weather.data.columns:
        from alinea.alep.alep_weather import add_septoria_infection_risk
        weather.check(varnames=['septo_infection_risk'], models={'septo_infection_risk':add_septoria_infection_risk})
    if not 'septo_degree_days' in weather.data.columns:
        from alinea.alep.alep_weather import basic_degree_days
        weather.check(varnames=['septo_degree_days'], models={'septo_degree_days':basic_degree_days}, start_date=start_date, base_temp=-2.)
    cond_inf = weather.data['septo_infection_risk'][seq].values
    ddays = weather.data['septo_degree_days'][seq].values
    ind = []
    count=0
    previous_count =0.
    for i in range(len(cond_inf)):
        if rain_filter[i]==True:
            cond_inf[i] = True
        if cond_inf[i]==True:
            count = 0.
        else:
            count += (ddays[i] - previous_count)
            if count >= degree_days:
                ind.append(i-1)
                count = 0.
            previous_count = ddays[i]
    cond_inf[ind] = True
    return cond_inf==1