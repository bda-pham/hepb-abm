import tables as tb
import numpy as np

"""
Observer for infections, as well as metrics to calculate DALY.
"""
from disease.observers.obs_base import Observer

class IncidenceObserver(Observer):
    def __init__(self, h5file):
        desc = {}
        desc['t'] = tb.UInt32Col()
        desc['infections'] = tb.UInt32Col()
        #test
        
        super(IncidenceObserver, self).__init__(h5file, 'incidence', desc, 'Incidence Observer')
        
    def update(self, t, pop, disease, cases, new_I, **kwargs):
        self.row['t'] = t
        self.row['infections'] = len(new_I)
        self.row.append()
        self.h5file.flush()
