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
from hepb.constants import Origin
from disease.general.sim_epi import SimEpi
from hepb.ind_com import IndCom
from population.household import Household
from population.pop_gen import allocate_couples, split_age_probs, duplicate_household
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
                    cur_ind.origin = Origin.THAI
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
        # print(f"Added community {comm_id}")
        assert(all([comm_id == ind.groups['community'] for ind in cur_comm]))
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

    def __init__(self, p, disease, rng, ind_type=IndCom):
        super(SimEpiCom, self).__init__(p, disease, rng, ind_type)
        self.imm_residue = 0
 
    def create_population(self):
        """
        Create a population according to specified age and household size
        distributions.
        """
        self._setup_params()
        self.P = PopHHCom(self.ind_type, self.p['logging'])
        gen_hh_com_age_structured_pop(self.P, self.p['pop_size'], self.p['com_dist'], self.hh_comp,
                                  self.age_dist, self.p['age_cutoffs'], self.rng)

        allocate_couples(self.P)

    def _load_demographic_data(self):
        super(SimEpiCom, self)._load_demographic_data()
        self.mobility_rates = [[_adjust_prob(r, self.p['t_per_year']) for r in rates] for rates in self.p["mobility_rates"]]
        self.origin_mobility = self.p["origin_mobility"].copy()

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

            same_origin_candidates = [
                x for x in candidates
                if x.origin == ind.origin
            ]

        # abort if no eligible partner exists
        if not candidates:
            print(f"No eligible partner! origin: {ind.origin}, site: {ind.groups['community']}")
            return None
        else:
            if same_origin_candidates:
                partner = self.rng.choice(same_origin_candidates)
            else:
                partner = self.rng.choice(candidates)
            # print(f'same origin? {partner.origin == ind.origin}')
            return partner

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

    def _main_loop(self, year_begin, years, verbose=False):
        """
        Run simulation.
        """
        t_begin = int(year_begin * self.p['t_per_year'])
        t_end = int((year_begin + years) * self.p['t_per_year'])

        t_now = (year_begin + self.p['year_now']) * self.p['t_per_year']

        if verbose:
            self.start_time = time.time()
            self.print_column_labels()
            self.print_pop_numbers(t_begin)

        self.disease.update_observers(t_begin, disease=self.disease, pop=self.P,
                                      cases=[],
                                      boosting=[],
                                      introduction=False,
                                      new_I=[], rng=self.rng)

        for t in range(t_begin + 1, t_end + 1):
            if t == t_now:
                self.disease.comunity_access[0] = self.p['new_remote_access']
                self.disease.comunity_access[1] = self.p['new_village_access']
                self.disease.origin_access[Origin.MIGRANT] = self.p['new_migrant_access']
                self.origin_mobility[Origin.MIGRANT] = self.p['new_migrant_mobility']
                self.disease.set_annual_treatment_rate(self.p['new_treat_rate'])
            # update demography (if required)
            if self.p['update_demog']:  # and t%52==0:
                births, deaths, imms, birthdays = self.update_all_demo(t)  # *52)
                firstborns = len([x for x in births if x.birth_order == 1])
                self.disease.firstborn_count += firstborns
                self.disease.subsequent_count += (len(births) - firstborns)
                # oh wow, really don't need to be doing THIS all the time!
                # a) no point unless also reinitialising contact matrix;
                # b) no point at all unless population structure is changing over time
                # self.disease.cmatrix.update_age_classes(
                #        births, deaths, imms, birthdays)
                # if self.p['dyn_rates'] and t % (self.p['cm_update_years'] * self.p['t_per_year']) == 0:
                if self.p['cm_update_years'] > 0 and t % (self.p['cm_update_years'] * self.p['t_per_year']) == 0:
                    self._init_contact_matrix(t)
                self.disease.bd_update(t, births, deaths, imms, self.rng)

            # update disease
            if self.disease.update(t, self.P, self.rng):
                if verbose:
                    self.print_pop_numbers(t)
                break  # update returns true if halting upon fade out

            if verbose:
                self.print_pop_numbers(t)

        if verbose:
            self.print_column_labels()
            self.end_time = time.time()
            print("time:", self.end_time - self.start_time)
    
    def update_all_demo(self, t):
        """
        Carry out a single update of all demographic aspects population.

        :param t: the current time step.
        :type t: int
        :returns: a tuple containing lists of births, deaths, immigrants and birthdays

        """

        emigrants = []
        # movement between communities
        mobility_rates = self.mobility_rates[1:]
        for hh in self.P.groups['household'].values():
            mobility = self.origin_mobility[hh[0].origin]
            if mobility > 0 and self.rng.random() < mobility:
                cur_com = hh[0].groups['community']
                tar_com = proportional_sample(mobility_rates[cur_com], self.rng)
                if tar_com is not None:
                    if tar_com == 0:
                        # leave the region
                        emigrants.extend(hh)
                    else:
                        for ind in hh:
                            self.P.remove_individual_from_group('community', ind)
                        # print(f"Community: household moving from {cur_com} to {tar_com}")
                        self.P.add_individuals_to_group('community', tar_com-1, hh)
        for ind in emigrants:
            self.P.remove_individual(ind)
        
        # age each individual by appropriate number of days
        birthdays = self.P.age_population(364 // self.p['t_per_year'])

        deaths = []
        births = []

        # calculate index for fertility and mortality rates
        # basically: use first entry for burn-in, then one entry every 
        # 'period' years, then use the final entry for any remaining years.
        index = min(max(
            0, (t - (self.p['demo_burn'] * self.p['t_per_year'])) //
               (self.p['t_per_year'])), self.dyn_years) \
            if self.p['dyn_rates'] else 0

        # print cur_t / self.p['t_per_year'], \
        #     index, self.p_adj['growth_rates'][index], \
        #     len(self.P.I)

        cur_inds = list(self.P.I.values())
        for ind in cur_inds:
            death, birth = self._update_individual_demo(t, ind, index)
            if death == "error" and birth == "error":
                return "error", "error", "error", "error"
            if death:
                deaths.append(death)
            if birth:
                births.append(birth)

        # trigger delayed (due to pregnancy births)
        for mother in self.preg_schedule[t]:
            # create new individuals
            new_ind = self.P.birth(t, mother, 0 if self.rng.random() < 0.5 else 1)
            births.append(new_ind)
            # remove mother from pregnancy list
            del self.P.preg_current[mother]
        del self.preg_schedule[t]


        # immigration
        ##############

        imm_count = 0
        imm_tgt = len(self.P.I) * self.p_adj['imm_rates'][index] + self.imm_residue
        source_hh_ids = []
        immigrants = []
        while imm_count < imm_tgt - 1:
            hh_id = self.rng.choice(list(self.P.groups['household'].keys()))
            imm_count += len(self.P.groups['household'][hh_id])
            source_hh_ids.append(hh_id)
        self.imm_residue = imm_tgt - imm_count
        # print(f'Immigration: target {imm_tgt}, actual {imm_count}, residue {self.imm_residue}')
        for hh_id in source_hh_ids:
            new_hh_id = duplicate_household(self.P, t, hh_id)
            immigrants.extend(self.P.groups['household'][new_hh_id])

        for ind in immigrants:
            ind.origin = Origin.MIGRANT
        immigration_rates = self.mobility_rates[0][1:]

        if sum(immigration_rates) > 0:
            total = sum(immigration_rates)
            normalised_rates = [r/total for r in immigration_rates]
            for ind in immigrants:
                comm = proportional_sample(normalised_rates, self.rng)
                self.P.add_individuals_to_group('community', comm, [ind])
        else:
            for ind in immigrants:
                comm = self.rng.choice(range(len(immigration_rates))) 
                self.P.add_individuals_to_group('community', comm, [ind])

        # natural population growth
        ####################

        # print "growth rate (%d) = %f: %d" % (index, self.p_adj['growth_rates'][index], len(self.P.I))

        # calculate number of new individuals to add (whole and fraction)
        new_individuals = (len(self.P.I) * self.p_adj['growth_rates'][index] + len(deaths) - len(immigrants))

        # get whole part of new individuals
        new_now = int(max(0, new_individuals))
        # add fractional part to residue accumulation
        self.growth_residue += (new_individuals - new_now)
        # grab any new 'whole' individuals
        new_residue = int(self.growth_residue)
        new_now += new_residue
        self.growth_residue -= new_residue

        # create the new individuals
        for _ in range(new_now):
            #mother = self._choose_mother(index)
            mother = self._choose_mother(0)
            if mother == "error":
                return "error", "error", "error", "error"
            new_ind = self.P.birth(t, mother, 0 if self.rng.random() < 0.5 else 1)
            if len(new_ind.parents) > 1:
                new_ind.origin = new_ind.parents[0].origin if self.rng.random() < 0.5 else new_ind.parents[1].origin
            births.append(new_ind)

        # set comm of immigrants
        # immigration
        ##############

        # print(f'immigrants {len(immigrants)}, emigrants {len(emigrants)}, births {len(births)}, deaths {len(deaths)}')
        
        print(f't: {t}, communities: { {k:{origin.name:sum(ind.origin==origin for ind in comm) for origin in Origin} for k, comm in self.P.groups["community"].items()}}, total: {len(self.P.I)}')
        return births, deaths, immigrants, birthdays
    
    def _update_individual_demo(self, t, ind, index=0):
        """
        Update individual ind; check for death, couple formation, leaving home
        or divorce, as possible and appropriate.
        """

        death = None
        birth = None

        # DEATH / BIRTH:
        #if self.rng.random() > exp(-self.death_rates[ind.sex][ind.age][index]):
        if ind.dying or self.rng.random() > exp(-self.death_rates[ind.sex][ind.age][0]):

            # currently pregnant women are 'immune' from dying
            if ind in self.P.preg_current:
                return death, birth
            death = ind
            
            # trigger death and reallocate any orphan children
            orphans = self.P.death(t, ind)
            for cur_dep in orphans:
                if cur_dep.age > self.p['leaving_age']:
                    self.P.leave_home(t, cur_dep)
                else:
                    hh = self._choose_household(cur_dep)
                    self.P.allocate_orphan(t, cur_dep, hh)

            
        # COUPLE FORMATION:
        elif self.p['couple_age'] < ind.age < self.p['couple_age_max'] \
                and not ind.partner \
                and self.rng.random() < self.p_adj['couple_probs'][index]:
            partner = self._choose_partner(ind)
            if partner:
                self.P.form_couple(t, ind, partner)

        # LEAVING HOME:
        elif ind.age > self.p['leaving_age'] \
                and ind.with_parents \
                and not ind.partner \
                and self.rng.random() < self.p_adj['leaving_probs'][index]:
            self.P.leave_home(t, ind)

        # DIVORCE:
        elif self.p['divorce_age'] < ind.age < self.p['divorce_age_max'] \
                and ind.partner \
                and self.rng.random() < self.p_adj['divorce_probs'][index]:
            self.P.separate_couple(t, ind)

        # ELSE: individual has a quiet year...
        return death, birth

def proportional_sample(distribution, rng):
    if 1 in distribution:
        return distribution.index(1)
    else:
        cul_prob = 1
        for i, prob in enumerate(distribution):
            if rng.random() < prob / cul_prob:
                return i
            cul_prob -= prob
        return None