"""
Basic disease states.
"""

from math import exp
import numpy as np

from disease.general.state_base import State, StateTimed

class Susceptible(State):
    """Susceptible state"""

    def __init__(self, order):
        super(Susceptible, self).__init__(order, 'S', 'green')
        self.at_risk = True
        self.infectious = False

    # record infection at infection, not exposure
    # TODO update infection probability for Trachoma
    def test_exposure(self, states, ind, foi, rng):
        if rng.random() < 1.0 - exp(-foi):
            ind.next_state = states['A']
            return "boosting"
        else:
            return False


###############################################################################

# class HepBExposed(StateTimed):
#     """Exposed state"""

#     def __init__(self, order, duration):
#         super(HepBExposed, self).__init__(order, duration, 'E', 'o')
#         self.at_risk = False
#         self.infectious = False

#     def update(self, t, ind, states):
#         if super(HepBExposed, self).update(t, ind, states):
#             ind.next_state = states['I']

class DurationGeneratorGeometric(object):
    def __init__(self, rate, rng):
        self.rate = rate
        self.rng = np.random.RandomState(rng.randint(0, 99999999))

    def get_duration(self):
        return self.rng.geometric(min(1, self.rate))

###############################################################################

class Acute(StateTimed):
    """Infected state"""

    def __init__(self, order, duration, rng):
        super(Acute, self).__init__(order, duration, 'A', 'yellow')
        self.at_risk = False
        self.infectious = True
        self.current = set()
        self.rng = np.random.RandomState(rng.randint(0, 99999999))

    def enter(self, t, ind):
        super(Acute, self).enter(t, ind)
        self.current.add(ind.ID)

    def exit(self, t, ind):
        super(Acute, self).exit(t, ind)
        self.current.remove(ind.ID)

    def update(self, t, ind, states):
        super(Acute, self).update(t, ind, states)
        # calculate chronic chance
        # simply 90% for < 5 and 10% for 5+ for now
        if super(Acute, self).update(t, ind, states):
            chr_prob = 0.1
            if ind.age <= 5:
                chr_prob = 0.9 - 0.7 * ind.age / 5
            elif ind.age < 6:
                chr_prob = 0.2 - 0.1 * (ind.age - 5)
            if self.rng.random() < chr_prob:
                ind.next_state = states['C']
            else:
                ind.next_state = states['R']


###############################################################################

class Chronic(State):
    """Diseased state"""

    def __init__(self, order, death_rate, rng):
        super(Chronic, self).__init__(order, 'C', 'red')
        self.at_risk = False
        self.infectious = True
        self.current = set()
        self.death_rate = death_rate
        self.rng = np.random.RandomState(rng.randint(0, 99999999))

    def update(self, t, ind, states):
        super(Chronic, self).update(t, ind, states)
        if self.rng.random() < self.death_rate:
            ind.dying = True

    def enter(self, t, ind):
        super(Chronic, self).enter(t, ind)
        self.current.add(ind.ID)

    def exit(self, t, ind):
        super(Chronic, self).exit(t, ind)
        self.current.remove(ind.ID)


class Recovered(State):
    """Recovered state"""

    def __init__(self, order):
        super(Recovered, self).__init__(order, 'R', 'blue')
        self.at_risk = False
        self.infectious = False


class Treated(State):
    """Treated state"""

    def __init__(self, order):
        super(Treated, self).__init__(order, 'T', 'blue')
        self.at_risk = False
        self.infectious = False

class Vaccinated(State):
    """Vaccinated state"""

    def __init__(self, order):
        super(Vaccinated, self).__init__(order, 'V', 'blue')
        self.at_risk = False
        self.infectious = False
