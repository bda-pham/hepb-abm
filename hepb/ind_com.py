"""
Base individual class for epidemic simulations, adding state and counters.
"""
from disease.general.ind_epi import IndEpi


class IndCom(IndEpi):
    __slots__ = 'origin', 'dying', 'generation'

    def __init__(self, new_id, age=0, sex=0, bootstrap=False, logging=True):
        super(IndCom, self).__init__(new_id, age, sex, bootstrap, logging)
        self.origin = None
        self.dying = False
        self.generation = 0
