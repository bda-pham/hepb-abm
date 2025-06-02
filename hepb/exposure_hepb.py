from math import sin, pi

from disease.SEIR.exposure_SEIR_add import ExposureSEIRadd


class ExposureHepB(ExposureSEIRadd):
    def __init__(self, p):
        super(ExposureHepB, self).__init__(p) 
        self.acu_ht_prob = p['acu_ht_prob']

    def calc_foi_fast(self, t, ind, pop, comm_precomp, rng):
        """
        Calculates force of infection acting on an individual based on prevalence
        of infection in household and community

        :param t: current timestep
        :param ind: current individual
        :param pop: population
        :param comm_precomp: pre-computed age-specific community FOI
        :return: FOI acting on ind
        """

        # list of hh_members
        hh_members = pop.groups['household'][ind.groups['household']]

        # list of infectious hh members (I_H)
        hh_I = []
        hh_infectivity = 0
        for x in hh_members:
            hh_I.append(x)
        for x in hh_I:
            if x.state.label == "C":
                hh_infectivity += 1
            elif x.state.label == 'A':
                hh_infectivity += 0.7
        N_sub_1 = len(hh_members) - 1

        # proportion/count of household members infectious (I_H / (N_H - 1)^alpha)
        hh_I_prop = hh_infectivity / pow(N_sub_1, self.alpha) if N_sub_1 else 0
        # calculate household contribution
        hh = self.q_h * hh_I_prop
        # calculate community contribution
        # (uses pre-computed values for \sum_j (\eta_{ij} I_j / N_j) )
        comm = self.q * (comm_precomp['A'][ind.groups['community']][ind.age]*self.acu_ht_prob + comm_precomp['C'][ind.groups['community']][ind.age])
        # print(comm_precomp['I'][ind.age])

        # combined exposure is sum of household and community
        combined = hh + comm

        # apply relative adult susceptibility
        if ind.age > 18:
            combined *= self.rel_susc_adult

        # determine if household source and if
#        ind.hh_frac = hh / (hh + comm) if (hh + comm) > 0.0 else 0.0
        ind.hh_frac = hh / combined if combined > 0.0 else 0.0
        # print(f"HH: {hh}, COMM: {comm}, HH_FRAC = {ind.hh_frac}")

        if ind.hh_frac > 0.0 and rng.random() < ind.hh_frac:
            ind.hh_source = 1
            ind.source = rng.sample(hh_I, 1)[0]
        else:
            ind.hh_source = 0
            ind.source = 0

        # compute seasonal forcing term
        sf = (1.0 + (self.sf_amp * sin(t * self.two_pi))) if self.sf_amp > 0.0 else 1.0
        return sf * combined
