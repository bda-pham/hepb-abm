__author__ = 'ngeard'

p = {
    # directories
    'resource_prefix': 'data/',
    'prefix': 'output/',

    # demographic model parameters
    'hh_composition': 'hh_comp.dat',
    'age_distribution': 'age_dist.dat',
    'fertility_parity_probs': 'fertility_age_parity_probs.dat',
    'fertility_age_probs': 'fertility_age_probs.dat',
    'death_rates_m': 'death_rates_male.dat',
    'death_rates_f': 'death_rates_female.dat',

    'couple_prob': 0.06,
    'leaving_prob': 0.005,
    'divorce_prob': 0.001,
    'couple_age': 15,
    'couple_age_max': 60,
    'leaving_age': 18,
    'divorce_age': 18,
    'divorce_age_max': 60,
    'partner_age_diff': -2,
    'partner_age_sd': 2,
    'min_partner_age': 15,
    'birth_gap_mean': 270,
    'birth_gap_sd': 1,

    'pop_size': 90,
    'growth_rate': 0.0,
    'imm_rate': 0,

    'com_dist': [0.2, 0.8],
    'mobility_rates': [[0, 1, 0],
                       [0, 0.98, 0.02],
                       [0, 0, 1]],
    'healthcare_access': [0.2, 1],

    'preg': False,
    'use_parity': False,
    'dyn_rates': True,
    'update_demog': True,

    'growth_rate_file': 'growth_rates.dat',
    'imm_rate_file': 'imm_rates.dat',
    'divorce_prob_file': 'divorce_probs.dat',
    'leaving_prob_file': 'leaving_probs.dat',
    'couple_prob_file': 'couple_probs.dat',

    # contact model parameters
    'cm_gauss': False,
    'cm_smooth': True,
    'sigma_2': 10.0,
    'epsilon': 0.8,
    'phi': 0.7,
    'cm_update_years': 5,

    # disease model parameters. rates are annual
    'infectious_rate': 1.0/5,
    'acu_mtc_prob': 0.7,
    'chr_mtc_prob': 0.95,
    'treat_rate': 0.009,
    'death_rate': 00.000189,
    'imm_prevalence': 0.1,

    # intervention parameters
    'pmtct_cover': 0.5,
    'vac_cover' : 0.95,
    'treat_cover': 0.05,

    'start_ratio': 0.2,

    'q': 0.000,
    'q_nl': 2.7,
    'nl_power': 2,
    'q_h': 0.25,

    'cover': 0.8,

    't_offset': 65,

    'external_exposure_rate': 0,#5e-6,

    'random_seed': True,
    'seed': 0,
    't_per_year': 13,
    'years': [0, 40],
    'burn_in': 100,
    'burn_in_t_per_year': 1,
    'epi_burn_in': 60,

    'halt': False,

    # run parameters
    'num_runs': 1,
    'initial_cases': 1000,
    'output_list': ['all'],
    'save_cp': True,
    'logging': False,
}

p['demo_burn'] = p['burn_in']# + p['epi_burn_in']



