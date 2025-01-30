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

        days_per_t = 364 / p['t_per_year']
        self.q_nl = p['q_nl']
        self.nl_power = p['nl_power']
        self.acu_mtc_prob = p['acu_mtc_prob']
        self.chr_mtc_prob = p['chr_mtc_prob']
        self.start_t = (p['burn_in'] + p['epi_burn_in']) * p['t_per_year']
        self.end_t = self.start_t + p['years'][1] * p['t_per_year']
        self.vac_cover = p['vac_cover']
        self.pmtct_cover = p['pmtct_cover']
        self.start_ratio = p['start_ratio']

        self.imm_prev = p['imm_prevalence']
        self.healthcare_access = p['healthcare_access']

        sus = Susceptible(0)
        acu = Acute(1, DurationGeneratorFixed(p['t_per_year']/2), rng)
        chr = Chronic(2, 1-pow(1-p['treat_rate']/364, days_per_t), 1-pow(1-p['death_rate']/364, days_per_t), rng, 
                      healthcare_access=p['healthcare_access'])
        rec = Recovered(3)
        tre = Treated(4)
        vac = Vaccinated(5)
        
        self.add_states(sus, acu, chr, rec, tre, vac)

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
        #         self.by_age[label] = [0] * len(self.cmatrix.age_classes)
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
                # self.by_age[ind.state.label][self.cmatrix.age_map[ind.age]] += 1
                ind.state.current.add(ind.ID)

        # update age-sorted lists of individuals
        self.cmatrix.init_age_classes(P)

    
    def bd_update(self, t, births, deaths, imms, rng):
        """
        Update state counts for births and deaths and immigration.
        
        Currently immigrants are treated as arriving susceptible.
        
        TODO: update this to handle, e.g., transfer of maternal immunity.
        """

        for ind in chain(imms):
            if rng.random() < self.imm_prev:
                ind.next_state = self.states['C']
            else:
                ind.next_state = self.states[self.basic_susceptible]
            self.tick(t, ind)

        universal_pmtct_cover = self.start_ratio * self.pmtct_cover + (1-self.start_ratio) * self.pmtct_cover * (t - self.start_t) / (self.end_t - self.start_t)
        universal_vac_cover = self.start_ratio*self.vac_cover + (1-self.start_ratio) * self.vac_cover * (t - self.start_t) / (self.end_t - self.start_t)
        for ind in births:
            mother = [p for p in ind.parents if p.sex == 0]
            ind.next_state = self.states[self.basic_susceptible]
            community = ind.groups['community']
            modified_pmtct_cover = universal_pmtct_cover * self.healthcare_access[community]
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


            modified_vac_cover = universal_vac_cover * self.healthcare_access[community]
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
                #     self.by_age[ind.state.label][self.cmatrix.age_map[ind.age_at_infection]] -= 1
#            ind.state = None
        self.birth_count = len(births)
        self.death_count = len(deaths)

