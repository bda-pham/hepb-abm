import tables as tb
import random as rd
"""
Observer for prevalence of disease in a certain age group.
"""
from disease.observers.obs_base import Observer

class PrevalenceObserver(Observer):
    def __init__(self, h5file, state_labels, age_min=0, age_max=99):
        self.state_labels = state_labels
        self.age_group = range(age_min, age_max+1)
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
        count = {label:0 for label in self.state_labels}
        prev_age = [0] * 6
        prev_house = [0] * 9
        count_I_uncertain = 0
        for age in self.age_group:
            for p in pop.I_by_age[age].values():
                count[p.state.label] += 1
                # if p.state.label == "I":
                #     prob = 1 - min(0.5, 0.01*len(p.infections))
                #     if rd.random() < prob:
                #         count_I_uncertain += 1

        count_all = {label:0 for label in self.state_labels}

        age_i_count = [0] * 6
        age_count = [0] * 6
        age_thresholds = [10, 25, 35, 45, 55, 200]

        house_i_count = [0] * 9
        house_count = [0] * 9
        house_thresholds = [2, 4, 6, 50]

        for ind in pop.I.values():
            count_all[ind.state.label] += 1
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
            
            household_group = min(9, len(pop.groups['household'][ind.groups['household']]))-1
            age_count[age_group] += 1
            house_count[household_group] += 1
            
            if ind.state.label in disease.infectious_states:
                # prob = 1 - min(1, 0.03*len(ind.infections))
                # print(f'age: {ind.age}, infections no: {len(ind.infections)} -- detection prob: {prob}')
                # if rd.random() < prob:
                age_i_count[age_group] += 1
                house_i_count[household_group] += 1
                # print(f'-----> DETECTED!')
        for i in range(6):
            if age_count[i] > 0:
                prev_age[i] = age_i_count[i] / age_count[i]
        for i in range(9):
            if house_count[i] > 0:
                prev_house[i] = house_i_count[i] / house_count[i]

        # print(f't: {t} --- Age prev: {prev_age}')
        # print(f'Household prev: {prev_house}')
        self.row['t'] = t
        self.row['I_uncertain'] = count_I_uncertain
        self.row['Prevalence_all'] = 1.0 * sum([count_all[s] for s in disease.infectious_states]) / len(pop.I)
        for label in self.state_labels:
            self.row[label] = count[label]
        self.row.append()

        self.prev_group.row['t'] = t
        self.prev_group.row['prev_age'] = prev_age
        self.prev_group.row['prev_household'] = prev_house
        self.prev_group.row.append()

        self.h5file.flush()
