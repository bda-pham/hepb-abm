from disease.general.disease_base import DiseaseBase
from hepb.exposure_hepb import ExposureHepB
from hepb.states_hepb import *
from itertools import chain
from statistics import mean
from disease.general.duration import DurationGeneratorFixed

class DiseaseHepB(DiseaseBase):
    def __init__(self, p, cmatrix, rng, fname, mode):
        super(DiseaseHepB, self).__init__(p, None, cmatrix, fname, mode, basic_susceptible='S',
                                          basic_infection='C', disease_states=('A','C'))
        self.exposure = ExposureHepB(p)

        self.days_per_t = 364 / p['t_per_year']
        self.q_nl = p['q_nl']
        self.nl_power = p['nl_power']
        self.acu_mtc_prob = p['acu_mtc_prob']
        self.chr_mtc_prob = p['chr_mtc_prob']
        self.start_t = (p['burn_in'] + p['epi_burn_in']) * p['t_per_year']
        self.max_cover_t = self.start_t + p['start_period'] * p['t_per_year']
        self.max_treat_cover_t = self.start_t + p['start_treat_period'] * p['t_per_year']
        self.vac_cover = p['vac_cover']
        self.pmtct_cover = p['pmtct_cover']
        self.start_ratio = p['start_ratio']
        self.start_treat_ratio = p['start_treat_ratio']
        self.communities = range(len(p['com_dist']))

        self.imm_prev = p['imm_prevalence']
        self.comunity_access = p['healthcare_access'].copy()
        self.origin_access = p['origin_access'].copy()
        self.set_annual_treatment_rate(p['treat_rate'])

        sus = Susceptible(0)
        acu = Acute(1, DurationGeneratorFixed(p['t_per_year']/2), rng)
        chr = Chronic(2, 1-pow(1-p['death_rate']/364, self.days_per_t), rng)
        rec = Recovered(3)
        tre = Treated(4)
        vac = Vaccinated(5)
        self.migrant_comm_id = p['migrant_comm_id']
        
        self.add_states(sus, acu, chr, rec, tre, vac)

    def set_annual_treatment_rate(self, treatment_rate):
        self._treatment_rate = 1-pow(1-treatment_rate/364, self.days_per_t)

    def set_counts(self, P):
        """
        Initialise state counts and I_by_age on basis of state of 
        individuals in population P.
        """
        # reset counts to zero
        for state in self.states.values():
            state.count = 0
        # for label in self.infectious_states:
        #     if self.cmatrix:
        #         self.by_age[label] = {comm: [0] * len(self.cmatrix.age_classes) for comm in self.communities}
        #     else:
        #         self.by_age[label] = [0]

        # set counts on basis of P
        for ind in P.I.values():
            # make sure individual is using this disease's state objects
            if ind.prev_state:
                ind.prev_state = self.states[ind.prev_state.label]
            ind.state = self.states[ind.state.label]
            if ind.next_state:
                ind.next_state = self.states[ind.next_state.label]
            # update state count
            self.states[ind.state.label].count += 1
            # update I_by_age dict

            if ind.state.label in self.infectious_states:
                # self.by_age[ind.state.label][ind.groups['community']][self.cmatrix.age_map[ind.age]] += 1
                ind.state.current.add(ind.ID)

        # update age-sorted lists of individuals
        self.cmatrix.init_age_classes(P)

    def add_states(self, *states):
        """Add a new disease state."""
        for state in states:
            self.states[state.label] = state

        self.infectious_states = [label for label in self.disease_states if self.states[label].infectious]
        # for label in self.infectious_states:
        #     if self.cmatrix:
        #         self.by_age[label] = {comm: [0] * len(self.cmatrix.age_classes) for comm in self.communities}
        #     else:
        #         self.by_age[label] = [0]

    def seed_infection(self, t, P, cases, rng, seed_inds=None):
        """Seed initial infection (set everyone else to susceptible)."""
        for ind in P.I.values():
            ind.next_state = self.states[self.basic_susceptible]
            self.tick(t, ind)
        if not seed_inds:
            main_pop = [ind for ind in P.I.values() if ind.groups['community'] != self.migrant_comm_id]
            seed_inds = rng.sample(main_pop, cases)
        for ind in seed_inds:
            ind.next_state = self.states[self.basic_infection]
            ind.source = -3
            I_in_age, I_out_age = self.tick(t, ind)
            # if I_in_age is not None and I_in_age >= 0:
            #     self.by_age[ind.state.label][ind.groups['community']][self.cmatrix.age_map[I_in_age]] += 1

    def bd_update(self, t, births, deaths, imms, rng):
        """
        Update state counts for births and deaths and immigration.
        
        Currently immigrants are treated as arriving susceptible.
        
        TODO: update this to handle, e.g., transfer of maternal immunity.
        """
        prev = self.states["C"].count / sum([x.count for x in self.states.values()])
        cover_level = min(1, self.start_ratio + (1-self.start_ratio) * (t - self.start_t) / (self.max_cover_t - self.start_t))
        for ind in chain(imms):
            if rng.random() < prev + self.imm_prev:
                ind.next_state = self.states['C']
            else:
                ind.next_state = self.states[self.basic_susceptible]
            self.tick(t, ind)

        universal_pmtct_cover = cover_level * self.pmtct_cover
        universal_vac_cover = cover_level * self.vac_cover
        for ind in births:
            if ind.groups['community'] == self.migrant_comm_id:
                ind.next_state = self.states[self.basic_susceptible]
                self.tick(t, ind)
                continue
            mother = [p for p in ind.parents if p.sex == 0]
            ind.next_state = self.states[self.basic_susceptible]
            community = ind.groups['community']
            modified_pmtct_cover = universal_pmtct_cover * self.comunity_access[community] * self.origin_access[ind.origin]
            if mother:
                mother = mother.pop()
                acu_prob = 0
                if mother.state.label == "C":
                    acu_prob = self.chr_mtc_prob
                elif mother.state.label == "A":
                    acu_prob = self.acu_mtc_prob


                if acu_prob > 0 and rng.random() < acu_prob:
                    if t > self.start_t and modified_pmtct_cover > 0 and rng.random() < modified_pmtct_cover:
                        ind.next_state = self.states["V"]
                    else:
                    #print(f"Mother {mother.state.label} --> baby infected (acute)")
                        ind.next_state = self.states["A"]


            modified_vac_cover = universal_vac_cover * self.comunity_access[community] * self.origin_access[ind.origin]
            if ind.next_state != self.states["V"] and t > self.start_t and modified_vac_cover > 0:
                if ind.next_state == self.states[self.basic_susceptible] and rng.random() < modified_vac_cover:
                    ind.next_state = self.states['V']
                elif ind.next_state == self.states["A"] and rng.random() < modified_vac_cover / 2:
                    ind.next_state = self.states['V']
            self.tick(t, ind)

        for ind in deaths:
            if ind.state:
                ind.state.exit(t, ind)
                # if ind.state.label in self.infectious_states:
                #     self.by_age[ind.state.label][ind.groups['community']][self.cmatrix.age_map[ind.age_at_infection]] -= 1
                ind.state = None
        self.birth_count = len(births)
        self.death_count = len(deaths)

        print(f"t: {t} -- Prevalence: {prev}")

    def update(self, t, P, rng):
        cover_level = min(1, self.start_treat_ratio + (1-self.start_treat_ratio) * (t - self.start_t) / (self.max_treat_cover_t - self.start_t))

        # treatment
        if t > self.start_t:
            for ind_ID in self.states['C'].current:
                ind = P.I[ind_ID]
                universal_treatment_rate = cover_level * self._treatment_rate
                treatment_prob = universal_treatment_rate * self.comunity_access[ind.groups['community']] * self.origin_access[ind.origin]
                if rng.random() < treatment_prob:
                    ind.next_state = self.states["T"]

        super(DiseaseHepB, self).update(t, P, rng)

    def check_exposure(self, t, P, rng):
        """
        Check all at risk individuals for potential exposure to infection.
        """

        # we will return a list of exposures (broken into infectious cases and boosting)
        cases = dict(infection=[], boosting=[])

        # return empty lists if there is currently no force of infection acting
        if sum([self.states[label].count for label in self.infectious_states]) == 0:
            return cases

        # create a set of individuals currently at risk
        # split depending upon whether network or matrix is being used to calculate community exposure
        pop_at_risk = set()
        comm = {}
        if self.cnetwork:
            # population at risk consists of household members and network neighbours of
            # people currently in an infectious state
            for cur_inf_state in self.infectious_states:
                for cur_I in self.states[cur_inf_state].current:
                    pop_at_risk.update([x for x in P.housemates(P.I[cur_I]) if x.state.at_risk])
                    pop_at_risk.update([x for x in self.cnetwork.get_contacts(P, P.I[cur_I]) if x.state.at_risk])
        else:
            # population at risk consists of potentially everybody
            pop_at_risk = [x for x in P.I.values() if x.groups['community'] != self.migrant_comm_id and x.state.at_risk]
            # compute exposure from community:
            # comm is the force of infection arising from infection in the community
            # EC[i] is a vector containing the contact rates between an individual in age group i
            # and individuals in each age group j, weighted by the (approximate) number of people
            # in age group j (as community mixing is density dependent).
            # I_by_age is the number of infected individuals in each age group j
            # total force of infection for each age class i is equal to the products of
            # weighted contact rate and the number of infected individuals, summed over
            # each age class i.
            # change from Nic's model: calculate by_age every timestep
            by_age = {label:{community: [0]*len(self.cmatrix.age_classes) for community in self.communities} for label in self.infectious_states}
            total = {community: [0]*len(self.cmatrix.age_classes) for community in self.communities}
            
            for ind in P.I.values():
                if ind.groups['community'] != self.migrant_comm_id:
                    if ind.state.infectious:
                        by_age[ind.state.label][ind.groups['community']][self.cmatrix.age_map[ind.age]] += 1
                    total[ind.groups['community']][self.cmatrix.age_map[ind.age]] += 1
            
            for label in self.infectious_states:
                comm[label] = {}
                if self.cmatrix:
                    for community in self.communities:
                        comm[label][community] = np.array([np.sum(self.cmatrix.EC[i] * by_age[label][community] / total[community]) for i in range(101)])
                else:
                    # dummy line for non-age based mixing
                    comm[label] = np.array([float(self.states[label].count) / len(P.I) for _ in range(101)])

        # test exposure for each individual at risk
        for ind in pop_at_risk:
            # split depending upon whether network or matrix is being used to calculate community exposure
            if self.cnetwork:
                foi = self.exposure.calc_foi_fast(t, ind, P, self.cnetwork, rng)
            else:
                foi = self.exposure.calc_foi_fast(t, ind, P, comm, rng)
            exposure_type = ind.state.test_exposure(self.states, ind, foi, rng)
            #                    p, s = P.hh_parents_siblings(ind)
            #                    p_I = len([x for x in p if x.state.infectious])
            #                    s_I = len([x for x in s if x.state.infectious])
            if exposure_type == 'infection':
                ind.infections.append(t)
                cases['infection'].append(ind)
            elif exposure_type == 'boosting':
                cases['boosting'].append(ind)
        return cases
    
    def update_ind_states(self, t, P):
        """
        Update disease status of all individuals.
        Returns a list of newly infectious ('symptomatic' individuals)
        """
        new_I = []
        # second loop updates current state
        for ind in P.I.values():
            old_state = ind.state.label
            I_in_age, I_out_age = self.tick(t, ind)
            new_state = ind.state.label

            if old_state not in self.infectious_states and new_state in self.infectious_states:
                new_I.append(ind)
            if not self.cmatrix: continue

            # if I_in_age is not None and I_in_age >= 0 and new_state in self.infectious_states:
            #     self.by_age[new_state][ind.groups['community']][self.cmatrix.age_map[I_in_age]] += 1
            # if I_out_age is not None and I_out_age >= 0 and old_state in self.infectious_states:
            #     self.by_age[old_state][ind.groups['community']][self.cmatrix.age_map[I_out_age]] -= 1

        return new_I