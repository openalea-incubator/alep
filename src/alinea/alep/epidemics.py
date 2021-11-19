"""Basic epidemic"""

class Epidemics(object):

    def __init__(self, fungus, dispersal):
        #todo fungus as a list of fungus
        self.fungus = fungus
        self.dispersal = dispersal

    def update_dus(self, deposits: Dict[Any, list], global_env=None, local_env = None):
        for vid, dus in deposits.items():
            for du in dus:
                if local_env is None:
                    du.update(global_env)
                else:
                    du.update(local_env[vid])
                if du.n_objects <= 0:
                    del du

    def update_lesions(self, lesions: Dict[Any, list], global_env=None, local_env = None):
        for vid, lesion_list in lesions.items():
            for les in lesion_list:
                if local_env is None:
                    les.update(global_env)
                else:
                    les.update(local_env[vid])
                if les.n_objects <= 0:
                    del les
                # TODO : merge dead lesions

    def infection(self, deposits: Dict[Any, list]):
        lesions = {}
        for vid, dus in deposits.items():
            if vid not in lesions:
                lesions[vid] = []
            for du in dus:
                response, mutations = du.infection_response()
                if du.n_objects <= 0:
                    del du
                lesion = du.produce(response, **mutations)
                lesion.setup()
                if lesion.n_objects > 0:
                    found = False
                    for l in lesions[vid]:
                        if lesion.is_like(l):
                            l.merge(lesion)
                            found = True
                            break
                    if not found:
                        lesions[vid].append(lesion)
        return lesions

    def growth(self, lesions, ressources):
        pass

    def step(self, lesions, dispersal_units, dispersal_targets=None, ressources=None, global_env=None, local_env = None, grow=False, disperse=False, infect=False, **kwds):
        loss = 0
        self.update_dus(dispersal_units, global_env, local_env)
        self.update_lesions(lesions, global_env, local_env)
        if disperse:
            deposits, loss =self.dispersal.disperse(self, lesions, targets = dispersal_targets, ** kwds)
            self.update_dus(deposits, global_env, local_env)
            #TODO : merge dispersal_units and deposits
        if infect :
            self.infection(dispersal_units)
        if grow:
            self.growth(lesions, ressources)

        return lesions, deposits, loss