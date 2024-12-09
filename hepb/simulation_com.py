"""
.. module:: simulation
.. moduleauthor:: Nic Geard <ngeard@unimelb.edu.au>
"""

import os
import tables as tb
import itertools
from random import Random
from math import exp
from collections import defaultdict

from hepb.pop_hhcom import PopHHCom
from disease.general.sim_epi import SimEpi
from disease.general.ind_epi import IndEpi
from population.pop_gen import allocate_couples, split_age_probs
from population.utils import sample_table, load_probs, \
    load_probs_new, load_age_rates, load_prob_tables, load_prob_list


def _adjust_prob(rate, t_per_year):
    """
    convert from an annual rate to a per-time-period probabiliy
    where time-period = 1/t_per_year

    :param rate: annual rate to convert.
    :type rate: double
    :param t_per_year: number of time periods per year.
    :type t_per_year: int
    """
    tp = 1.0 / t_per_year
    return 1 - pow(1 - rate, tp)

def gen_hh_com_age_structured_pop(pop, pop_size, comm_dist, hh_probs, age_probs_i,
                              cutoffs, rng):
    """
    Generate a population of individuals with age structure and household 
    composition.

    Household composition here is approximated by the number of 
    individuals who are:

    - pre-school age (0--4)
    - school age (5--17)
    - adult (18+)

    This is a bit ad hoc, but serves current purposes, and populations are
    intended to be 'burnt in' to create more natural structures.

    :param pop: The population
    :type pop: PopHH
    :param pop_size: The size of the population to be generated.
    :type pop_size: int
    :param hh_probs: A table of household size probabilities.
    :type hh_probs: list
    :param age_probs_i: A table of age probabilities.
    :type age_probs_i: list
    :param cutoffs: A tuple of values that specify ages cutoffs between pre-school, school and adults.
    :type cutoffs: tuple
    :param rng: The random number generator to use.
    :type rng: :class:`random.Random`
    """

    age_probs = [[x, int(y[0])] for x, y in age_probs_i]
    split_probs = split_age_probs(age_probs, cutoffs)
    split_probs.reverse()

    for comm_percent in comm_dist:
        cur_comm = []
        i = 0
        while i < pop_size * comm_percent:
            # get list of [adults, school age, preschool age]
            hh_type = [int(x) for x in sample_table(hh_probs, rng)]
            #            hh_type = [int(x) for x in sample_uniform(hh_probs, rng)]
            cur_hh = []
            for cur_hh_type, cur_prob in zip(hh_type, split_probs):
                for _ in itertools.repeat(None, cur_hh_type):
                    sex = 0 if rng.random() < 0.5 else 1
                    cur_age = sample_table(cur_prob, rng)
                    cur_ind = pop.add_individual(
                        cur_age, sex, adam=True, logging=pop.logging)
                    if pop.logging:
                        cur_ind.add_log(0, 'f', "Individual (bootstrap)")
                    cur_hh.append(cur_ind)
                    cur_comm.append(cur_ind)
                    i += 1
            hh_id = pop.add_group('household', cur_hh)
            for ind in cur_hh:
                ind.groups['household'] = hh_id

            if pop.logging:
                pop.households[hh_id] = Household(0, adam=True)
                for cur_ind in cur_hh:
                    cur_ind.household = pop.households[hh_id]
                pop.households[hh_id].add_log(
                    0, 'f', "Household (bootstrap)",
                    len(pop.groups['household'][hh_id]))
        comm_id = pop.add_group('community', cur_comm)
        print(f"Added community {comm_id}")
        print(pop.groups.keys())
        pass


class SimEpiCom(SimEpi):
    """
    Basic demographic simulation object.

    Handles updating of births, deaths, aging, immigration, and household structure (couple formation and
    disolution, leaving home, etc.)

    :param p: dictionary of simulation parameters.
    :type p: dict
    :param ind_type: the :class:`.individual.Individual` (sub)class stored by this population.
    :type ind_type: class
    :param create_pop: If `True` (default), create a random population; otherwise, this will need to be done later.
    :type create_pop: bool

    """

    def __init__(self, p, disease, rng, ind_type=IndEpi):
        super(SimEpiCom, self).__init__(p, disease, rng, ind_type)

    def create_population(self):
        """
        Create a population according to specified age and household size
        distributions.
        """
        self._setup_params()
        self.P = PopHHCom(self.ind_type, self.p['logging'])
        gen_hh_com_age_structured_pop(self.P, self.p['pop_size'], self.p['com_dist'], self.hh_comp,
                                  self.age_dist, self.p['age_cutoffs'], self.rng)
        c = 0
        for ind in self.P.I.values():
            if 'community' not in ind.groups:
                c += 1
                pass
        print(f'{c} individuals w no comm!')
        allocate_couples(self.P)

    def _choose_partner(self, ind):
        """
        Choose a partner for i_id, subject to parameter constraints.

        :param ind: the first partner in the couple.
        :type ind: ind_type
        :returns: partner if successful, otherwise None.
        """

        mean_age = ind.age + self.p['partner_age_diff'] \
            if ind.sex == 0 else ind.age - self.p['partner_age_diff']
        tgt_age = 0
        candidates = []
        while tgt_age < self.p['min_partner_age']:
            tgt_age = int(self.rng.gauss(mean_age, self.p['partner_age_sd']))
            tgt_set = self.P.individuals_by_age(tgt_age, tgt_age)
            candidates = [
                x for x in tgt_set
                if not x.partner and x.sex != ind.sex
                and x not in self.P.groups['household'][ind.groups['household']]
                and x.groups['community'] == ind.groups['community']
            ]

        # abort if no eligible partner exists
        return None if not candidates else self.rng.choice(candidates)

    def _choose_household(self, ind):
        """
        Process orphans who result when the last remaining adult guardian
        in their household dies.  If they are above 'adult-age' cutoff, place
        them in a new single household, otherwise, reallocate them to an
        existing family household (with at least one other child).
        """
        cur_hh = ind.groups['household']
        # choose destination household from among family households
        candidates = []
        # for hh in self.P.groups_by_min_size('household', 3):
        for hh in self.P.groups_by_min_size('household', 2):
            com = self.P.groups['household'][hh][0].groups['community']
            if hh != cur_hh and com == ind.groups['community']:
                candidates.append(hh)
        if len(candidates) == 0:
            candidates = self.P.groups_by_min_size('household', 2)
        tgt_hh = cur_hh
        while tgt_hh == cur_hh:
            tgt_hh = self.rng.sample(candidates, 1)[0]
        return tgt_hh


    
    def update_all_demo(self, t):
        """
        Carry out a single update of all demographic aspects population.

        :param t: the current time step.
        :type t: int
        :returns: a tuple containing lists of births, deaths, immigrants and birthdays

        """
        births, deaths, immigrants, birthdays = super(SimEpiCom, self).update_all_demo(t)
        # set comm of immigrants
        for ind in immigrants:
            ind.groups['community'] = 0

        # movement between communities
        for hh in self.P.groups['household'].values():
            cur_com = hh[0].groups['community']
            dist = self.p['mobility_rates'][cur_com+1]
            tar_com = proportional_sample(dist, self.rng)
            if tar_com:
                tar_com -= 1
                if tar_com == -1:
                    # leave the region
                    for ind in hh:
                        self.P.remove_individual(ind)
                    pass
                else:
                    for ind in hh:
                        self.P.remove_individual_from_group('community', ind)
                    print(self.P.groups.keys())
                    print(f"Community: household moving from {cur_com} to {tar_com}")
                    self.P.add_individuals_to_group('community', tar_com, hh)

        return births, deaths, immigrants, birthdays
    
def proportional_sample(distribution, rng):
    cul_prob = 1
    for i, prob in enumerate(distribution):
        if rng.random() < prob / cul_prob:
            return i
        cul_prob -= prob
    return None