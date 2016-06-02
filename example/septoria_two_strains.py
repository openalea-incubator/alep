""" Example of use of septoria model for two strains of septoria fungus
    with different parameters """
    
from septo_decomposed import *

class InfectionControlTwoStrains:
    def __init__(self, strain_names = ['septoria_strain1', 'septoria_strain2']):
        self.strain_names = strain_names

    def control(self, g, label = 'LeafElement'):
        """ Control if dispersal units can infect at their position """
        labels = g.property('label')
        bids = (v for v,l in labels.iteritems() if l.startswith('blade'))
        dispersal_units = {k:v for k,v in g.property('dispersal_units').iteritems() if len(v)>0.}
        areas = g.property('area')
        lesions = g.property('lesions')
        for blade in bids:
            leaf = [vid for vid in g.components(blade) if labels[vid].startswith(label)]
            leaf_lesions = sum([lesions[lf] for lf in leaf if lf in lesions], []) 
            les_surf = sum([les.surface for les in leaf_lesions])
            leaf_area = sum([areas[lf] for lf in leaf])
            ratio_les_surface = min(1, round(les_surf,3)/round(leaf_area,3)) if round(leaf_area,3)>0. else 0.
            for vid in set(leaf) & set(dispersal_units):
                dus_both_strains = []
                for strain in self.strain_names:
                    dus_to_keep = []
                    dus = []
                    for du in dispersal_units[vid]:
                        if du.is_active:
                            if du.status == 'deposited':
                                dus_to_keep.append(du)
                            elif du.status == 'emitted':
                                dus.append(du)
                                
                    if len(dus)>0 and dus[0].fungus.group_dus == True:
                        total_nb_dus = len(sum([du.position for du in dus],[]))
                    else:
                        total_nb_dus = len(dus)
                    nb_on_lesions = int(total_nb_dus*ratio_les_surface)
                    for du in range(nb_on_lesions):
                        random.shuffle(dus)
                        dus[0].position = dus[0].position[1:]
                        if dus[0].nb_dispersal_units==0.:
                            dus[0].disable()
                            dus = dus[1:]
                    for du in dus:
                        du.set_status(status = 'deposited')
                    dus_both_strains = dus_to_keep + dus
                if sum([DU.nb_dispersal_units for DU in dus_both_strains])>10 and len(dus)>0:
                    import pdb
                    pdb.set_trace()
                dispersal_units[vid] = dus_both_strains

def septo_disease_two_strains(adel, strain_names = ['septoria_strain1', 'septoria_strain2'],
                                sporulating_fraction = [0.5e-4, 0.5e-4], 
                                distri_chlorosis = [None, None], **kwds):
    """ Choose models to assemble the disease model.

    Parameters for the disease cycle are given in kwds in the following format:
        - kwds = {param1:[value_strain1, value_strain_2]}
    """
    domain = adel.domain
    domain_area = adel.domain_area
    convUnit = adel.convUnit
    inoculum = []
    for i_str, strain in enumerate(strain_names):
        if distri_chlorosis[i_str] is not None:
            fungus_strain = variable_septoria(distri_chlorosis = distri_chlorosis[i_str])
            mutable = True
        else:
            fungus_strain = plugin_septoria()
            mutable = False
        kwds_strain = {k:v[i_str] for k,v in kwds.iteritems()}
        fungus_strain.parameters(name = strain, group_dus=True, nb_rings_by_state=1, **kwds_strain)
        inoculum.append(SoilInoculum(fungus_strain, 
                                     sporulating_fraction=sporulating_fraction[i_str],
                                     domain_area=domain_area, mutable = mutable))
    contaminator = PopDropsSoilContamination(domain=domain, domain_area=domain_area)
    growth_controler = PriorityGrowthControl()
    infection_controler = InfectionControlTwoStrains(strain_names = strain_names)
    # infection_controler = BiotrophDUPositionModel()
    emitter = PopDropsEmission(domain=domain, compute_star = False)
    transporter = PopDropsTransport(domain = domain, domain_area = domain_area, dh = 0.01, convUnit = convUnit)
    return inoculum, contaminator, infection_controler, growth_controler, emitter, transporter
                
def run_disease_two_strains(start_date = "2011-10-21 12:00:00", end_date = "2012-07-18 01:00:00", 
                            variety = 'Tremie2', nplants = 1, nsect = 7,
                            disc_level = 5, dir = './adel/tremie_2012_1pl_7sect', 
                            strain_names = ['septoria_strain1', 'septoria_strain2'],
                            sporulating_fraction = [0.5e-4, 0.5e-4], 
                            distri_chlorosis = [None, None], **kwds):
    """ Simulate epidemics, loop over campaign

    Parameters for the disease cycle are given in kwds in the following format:
        - kwds = {param1:[value_strain1, value_strain_2]}
    """
    # Setup simulation
    adel, weather, seq, rain_timing, canopy_timing, septo_timing = setup(start_date = start_date,
                        end_date = end_date, variety = variety, nplants = nplants, nsect = nsect, 
                        disc_level = disc_level, Tmin = 0.)
    it_wheat = 0
    g,TT = adel.load(it_wheat, dir=dir)
    leaf_ids = adel_ids(g)
    
    recorders = [initiate_all_adel_septo_recorders(g, nsect, 
                    date_sequence = [v.index[-1] for v in septo_timing.values],
                    fungus_name = strain) for strain in strain_names]
                        
    # Assemble disease models
    inoculum, contaminator, infection_controler, growth_controler, emitter, transporter = septo_disease_two_strains(adel, strain_names,
                sporulating_fraction, distri_chlorosis, **kwds)
    
    for i, controls in enumerate(zip(canopy_timing, rain_timing, septo_timing)):
        canopy_iter, rain_iter, septo_iter = controls
        
        # Grow wheat canopy
        if canopy_iter:
            it_wheat += 1
            newg,TT = adel.load(it_wheat, dir=dir)
            move_properties(g, newg)
            g = newg
            leaf_ids = adel_ids(g)
        
        # Get weather for date and add it as properties on leaves
        if septo_iter:
            set_properties(g,label = 'LeafElement',
                           temperature_sequence = septo_iter.value.temperature_air,
                           wetness_sequence = septo_iter.value.wetness,
                           dd_sequence = septo_iter.value.degree_days)
        if rain_iter:
            set_properties(g,label = 'LeafElement',
                           rain_intensity = rain_iter.value.rain.mean(),
                           rain_duration = len(rain_iter.value.rain) if rain_iter.value.rain.sum() > 0 else 0.)

        # External contamination
        geom = g.property('geometry')
        if rain_iter and len(geom)>0 and rain_iter.value.rain.mean()>0.:
            for i_str in range(len(strain_names)):
                g = external_contamination(g, inoculum[i_str], contaminator, rain_iter.value)

        # Develop disease (infect for dispersal units and update for lesions)
        if septo_iter:
            infect(g, septo_iter.dt, infection_controler, label='LeafElement')
            update(g, septo_iter.dt, growth_controler, senescence_model=None, label='LeafElement')
            
        # Disperse and wash
        if rain_iter and len(geom)>0 and rain_iter.value.rain.mean()>0.:
            for strain in strain_names:
                g = disperse(g, emitter, transporter, strain, label='LeafElement', weather_data=rain_iter.value)
        
        # Save outputs
        if septo_iter:     
            date = septo_iter.value.index[-1]
            for recorder_strain in recorders:
                for plant in recorder_strain:
                    deads = []
                    for lf, recorder in recorder_strain[plant].iteritems():
                        recorder.update_vids_with_labels(adel_ids = leaf_ids)
                        recorder.record(g, date, degree_days = septo_iter.value.degree_days[-1])
                        if recorder.date_death != None and sum(recorder.data.leaf_green_area)>0.:
                            deads += [s for s in lf.split() if s.isdigit()]
                    if len(deads) > 0:
                        map(lambda x: recorder_strain[plant][x].inactivate(date), map(lambda x: 'F%d' % x, range(1, min(deads)+1)))
                    
    for recorder_strain in recorders:
        for plant in recorder_strain:
            for recorder in recorder_strain[plant].itervalues():
                recorder.get_complete_dataframe()
                recorder.get_normalized_audpc(variable='necrosis_percentage')
                recorder.get_audpc(variable='necrosis_percentage')
    
    return g, recorders