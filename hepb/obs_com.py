import tables as tb
import random as rd
from hepb.constants import Origin
"""
Observer for prevalence of disease in a certain age group.
"""
from disease.observers.obs_base import Observer

class ComObserver(Observer):
    def __init__(self, h5file, state_labels, age_min=0, age_max=100, label="community", title='Community Observer'):
        self.state_labels = state_labels
        self.age_group = range(age_min, age_max+1)
        desc = dict((x, tb.UInt32Col(pos=i+2)) for i, x in enumerate(self.state_labels))
        desc['t'] = tb.UInt32Col(pos=0)
        desc['community'] = tb.UInt8Col(pos=1)
        desc['origin'] = tb.UInt8Col(pos=1)
        super(ComObserver, self).__init__(h5file, label, desc, title)
        
    def create_storage(self, description, title):
        """
        Called when initialised with a h5file opened in 'w'rite mode.
        Adds tables for storage of data (plus local links)
        """

        Observer.create_storage(self, description, title)

    def load_storage(self):
        """
        Called when initialised with a h5file opened in 'r'ead or 'a'ppend mode.
        Creates local links to existing tables.
        """

        Observer.load_storage(self)

    def update(self, t, pop, disease, cases, new_I, **kwargs):
        count = {com_id:{origin:{label:0 for label in self.state_labels} for origin in Origin} for com_id in pop.groups['community']}
        for age in self.age_group:
            for ind in pop.I_by_age[age].values():
                count[ind.groups['community']][ind.origin][ind.state.label] += 1


        # print(f't: {t} --- Age prev: {prev_age}')
        # print(f'Household prev: {prev_house}')
        for com in pop.groups['community']:
            for origin in Origin:
                self.row['t'] = t
                self.row['community'] = com
                self.row['origin'] = origin.value
                for label in self.state_labels:
                    self.row[label] = count[com][origin][label]
                self.row.append()
        self.h5file.flush()
