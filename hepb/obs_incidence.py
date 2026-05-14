import tables as tb
import numpy as np
from hepb.constants import Route

"""
Observer for infections, as well as metrics to calculate DALY.
"""
from disease.observers.obs_base import Observer

class IncidenceObserver(Observer):
    def __init__(self, h5file):
        desc = {}
        desc['t'] = tb.UInt32Col()
        desc['infections'] = tb.UInt32Col()
        desc['route'] = tb.UInt8Col()
        #test
        
        super(IncidenceObserver, self).__init__(h5file, 'incidence', desc, 'Incidence Observer')
        
    def update(self, t, pop, disease, cases, new_I, **kwargs):
        
        routes = {r:0 for r in Route}
        for ind in new_I:
            # print(f't = {t}, route: {ind.route}, age: {ind.age}')
            routes[ind.route] += 1
        for r in Route:
            self.row['t'] = t
            self.row['infections'] = routes[r]
            self.row['route'] = r.value
            self.row.append()
        self.h5file.flush()
