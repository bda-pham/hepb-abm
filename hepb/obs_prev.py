import tables as tb
import random as rd
from hepb.constants import Origin
"""
Observer for prevalence of disease in a certain age group.
"""
from disease.observers.obs_base import Observer

class PrevalenceObserver(Observer):
    def __init__(self, h5file, state_labels, observed_states):
        self.state_labels = state_labels
        self.observed_states = observed_states
        desc = dict((x, tb.UInt32Col(pos=i+1)) for i, x in enumerate(self.state_labels))
        desc['t'] = tb.UInt32Col(pos=0)
        desc['I_uncertain'] = tb.Float32Col()
        desc['Prevalence_all'] = tb.Float32Col()
        super(PrevalenceObserver, self).__init__(h5file, 'prevalence', desc, 'Prevalence Observer')
        
    def create_storage(self, description, title):
        """
        Called when initialised with a h5file opened in 'w'rite mode.
        Adds tables for storage of data (plus local links)
        """

        Observer.create_storage(self, description, title)
        group = self.h5file.get_node('/', 'prevalence')
        self.prev_group = self.h5file.create_table(
            group, 'prev_group',
            {'t': tb.UInt32Col(), 'prev_age': tb.Float32Col(shape=(6)), 'prev_household': tb.Float32Col(shape=(9))},
             'Prevalence Group')
   
    def load_storage(self):
        """
        Called when initialised with a h5file opened in 'r'ead or 'a'ppend mode.
        Creates local links to existing tables.
        """

        Observer.load_storage(self)
        self.prev_group = self.h5file.root.cases.prev_group

    def update(self, t, pop, disease, cases, new_I, **kwargs):
        prev_age = [0] * 6

        age_i_count = [0] * 6
        age_count = [0] * 6
        age_thresholds = [10, 25, 35, 45, 55, 200]

        i_count = 0
        Thai_count = 0

        for ind in pop.I.values():
            age_group = -1
            if ind.origin == Origin.MIGRANT:
                continue
            
            for i, age_threshold in enumerate(age_thresholds):
                if ind.age < age_threshold:
                    age_group = i
                    break
            # age_group = min(7, int(ind.age / 10))
            # household_size = len(pop.groups['household'][ind.groups['household']])
            # for i, house_threshold in enumerate(house_thresholds):
            #     if household_size < house_threshold:
            #         household_group = i
            #         break
            
            age_count[age_group] += 1
            Thai_count += 1
            if ind.state.label in self.observed_states:
                # prob = 1 - min(1, 0.03*len(ind.infections))
                # print(f'age: {ind.age}, infections no: {len(ind.infections)} -- detection prob: {prob}')
                # if rd.random() < prob:
                i_count += 1
                age_i_count[age_group] += 1
                # print(f'-----> DETECTED!')
        for i in range(6):
            if age_count[i] > 0:
                prev_age[i] = age_i_count[i] / age_count[i]

        # print(f"t: {t} -- Thai prevalence: {i_count/Thai_count}")
        self.prev_group.row['t'] = t
        self.prev_group.row['prev_age'] = prev_age
        self.prev_group.row.append()

        self.h5file.flush()
