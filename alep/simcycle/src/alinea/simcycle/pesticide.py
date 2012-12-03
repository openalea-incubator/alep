"""
This module implements the model of Milne et al (2007) of the efficay of a (mix of) pesticide compounds

"""

#from math import log
from numpy import exp,log
from itertools import groupby
from operator import itemgetter
from numpy import recfromcsv

products_parameters = [{'compound': 'Epoxiconazole', 'Ap': 0.71, 'Kp': 6.0, 'Ae': 0.5, 'type_code': 0, 'Ke': 7.01, 'decay_rate': 0.069, 'dose_max_ha': 125}]

def compounds_from_csv(csvname, delimiter = ';') :
    """ 
    Read a csv of compounds parameters and import them in a dict.
    Expected columns are :
        - 'compound'
        - 'dose_max_ha'
        - 'type_code'
        - 'Ap'
        - 'Kp'
        - 'Ae'
        - 'Ke'
        - 'decay_rate'
    
    """
    
    tab = recfromcsv(csvname, delimiter = delimiter, case_sensitive = True)
    d = [dict(zip(tab.dtype.names, data)) for data in tab]
    return d
 
#products_parameters = compounds_from_csv('products_parameters.csv')

class SimcyclePesticide(Exception): pass       

def _dose_norm(dose, dose_max_ha):
    """ normalise doses(g.m-2) """
    dn = 0
    if dose_max_ha > 0 :
        dn = float(dose * 1e4) / dose_max_ha
    return dn
 
def _dose_response(dose, dose_max_ha, A, K): 
    return A * (1 - exp(-K * _dose_norm(dose,dose_max_ha)))
    
def _dose_decay(decay_rate, initial_dose, days):
    return initial_dose * exp(-decay_rate * days) 

def _dose_equivalent(dose_other, dmax_ref, dmax_other, A_ref, A_other, K_ref, K_other):
    dneq = 0 
    dn_other = _dose_norm(dose_other,dmax_other)
    if A_ref > 0 and A_other > 0 :
        #print((float(A_other)/float(A_ref)) * (1 - exp(- K_other * dn_other)))
        dneq = - (1./K_ref) * log(1 - (float(A_other)/A_ref) * (1 - exp(- K_other * dn_other)))  
    return dneq * dmax_ref / 1e4
    
    
def add_doses(pref,pother, effect='protectant'):
    """ add equivalent dose of other in pref."""
    deq = 0
    if effect == 'protectant':
        deq = _dose_equivalent(pother['dose'], pref['dose_max_ha'], pother['dose_max_ha'], pref['Ap'],pother['Ap'], pref['Kp'],pother['Kp'])
    elif effect == 'eradicant':
        deq = _dose_equivalent(pother['dose'], pref['dose_max_ha'], pother['dose_max_ha'], pref['Ae'],pother['Ae'], pref['Ke'],pother['Ke'])
    else :
        raise SimcyclePesticide('Unknown effect: %s'%effect)
    
    pref['dose'] += deq
    return pref
           
def global_efficacy(doses, compound_parameters = products_parameters):
    """
    compute global efficacy of a mixture of coumpounds
    
    :Parameters:
      - `doses` : A dict of ('compound' : doses) items giving the amount of product 'compound' present on the surface of the leaves (g.m-2)
      - `coumpound_parameters` : A list of dict of parameters. Parameters are :
          - 'compound' (string) - Name of the active substance
          - `dose_max_ha`(float) - Maximum recommended dose of coumpound for a single dose application (g.ha-1)
          - `type_code` (int) - Code for the mode of action of coumpound
          - `Ap` and `Kp` (floats) - Parameters for dose/response curve for protectant effect of compound
          - `Ae` and `Ke` (floats) - Parameters for dose/response curve for eradicant effect of compound
          - `decay_rate` (float, [0,1]) - Decay rate of the active substance over time

  :Returns:
      - `efficacy` : A dict with the following items :
          - `protectant` (float, [0,1]) - Protectant efficacy of the active subsance. Affects the number of successful infections. No infection will occur when `protectant` = 1. 
          - `eradicant` (float, [0,1]) - Eradicant efficacy of the active subsance. Affects the fungal development. Fungal development will be stopped when `eradicant` = 1.
          
    """

    efficacy = {'protectant': None, 'eradicant': None}
    ptable = dict([(p['compound'],p) for p in compound_parameters])
    d = {}
    for k in doses:
        try:
            d[k]=ptable[k]
        except KeyError:
            raise SimcyclePesticide('product %s not found in parameter dict'%k)
        d[k].update(dose=doses[k])
    # Sort products in order to group them by mode of action with 'groupby'
    d = sorted(d.values(),key= itemgetter('type_code')) 
    d = [list(g) for k,g in groupby(d,key= itemgetter('type_code'))]
    
    A = {'eradicant':'Ae','protectant':'Ap'}
    K = {'eradicant':'Ke','protectant':'Kp'}
    for mode in ('eradicant','protectant'):
        # Sort products inside the lists to get reference product (highest A) in first position
        ds = [sorted(g, key=itemgetter(A[mode]), reverse=True) for g in d]
        # Reduce lists with equivalent doses computation (additive effect)
        dc = [reduce(lambda x,y: add_doses(x,y,effect=mode),g) for g in ds]
        # 1 - efficacy computation for reference products with updated doses
        effs = [(1-_dose_response(p['dose'], p['dose_max_ha'],p[A[mode]], p[K[mode]])) for p in dc]
        #print effs
        # Application of the multiplicative effect
        efficacy[mode] = 1 - reduce(lambda x,y: x * y, effs)
                
    return efficacy

def erode_products(erosion_rates, compound_parameters = products_parameters):
    """ 
    Compute loss of activity of compounds due to evolution of resistance of strains
    :Parameters:
      - `erosion_rates` (dict) - A dict of ('compound': erosion_rate) where:
          - `erosion_rate` (float, [0,1]) - Rate of loss of efficacy of compounds due to evolution of resistance of strains
      - `coumpound_parameters` :see 'global_efficiency'
          
    :Returns:
      - `new_pars` : Same dict as `coumpound_parameters` with updated `Ap` and `Ae` linked to erosion of efficacy due to evolution of resistance of strains
      
    """
    
    ptable = dict([(p['compound'],p) for p in compound_parameters])
    new_pars = ptable.copy()
    for k in erosion_rates.keys() :
        if not k in ptable.keys():
            raise SimcyclePesticide('%s not found in compound_parameters'%k)
        
        new_pars[k]['Ap'] *= (1 - erosion_rates[k])
        new_pars[k]['Ae'] *= (1 - erosion_rates[k])
    
    return ptable.values(),
    
def decay(doses, days, compound_parameters = products_parameters):
    """ 
    Update doses after decay over time 
    :Parameters:
      - `coumpound_parameters` : see 'global_efficiency'
      - `doses` : A dict of ('compound' : doses) items giving the amount of product 'compound' present on the surface of the leaves (g.m-2)
      - `days` (float) - Number of days since the product application
          
      :Returns:
      - `current_dose` (float) - Updated mass of product on the surface of the leaves (g.m-2)
    
    """
    ptable = dict([(p['compound'],p) for p in compound_parameters])
    for k in doses.keys():
        try:
            doses[k] = _dose_decay(ptable[k]['decay_rate'], doses[k], days)
        except KeyError:
            raise SimcyclePesticide('product %s not found in parameter dict'%k)
        
    return doses,
