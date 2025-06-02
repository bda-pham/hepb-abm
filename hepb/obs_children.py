import tables as tb
import random as rd
from hepb.constants import Origin
"""
Observer for prevalence of disease in a certain age group.
"""
from hepb.obs_com import ComObserver

class ChildrenObserver(ComObserver):
    def __init__(self, h5file, state_labels, age_min=0, age_max=5):
        super(ChildrenObserver, self).__init__(h5file, state_labels, age_min=age_min, age_max=age_max, label='children', title='Children Observer')
        