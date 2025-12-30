__author__ = 'ngeard'
from hepb.constants import Origin

p = {
    # directories
    'resource_prefix': 'data/',
    'prefix': 'output/',

    # demographic model parameters
    'hh_composition': 'hh_comp.dat',
    'age_distribution': 'age_dist.dat',
    'fertility_parity_probs': 'fertility_age_parity_probs.dat',
    'fertility_age_probs': 'fertility_age_probs.dat',
    # 'age_fertility_rates_files': {
    #     Origin.THAI: 'fertility_age_rates.dat',
    #     Origin.MIGRANT: 'fertility_age_rates_migrant.dat'},
    # 'total_fertility_rates_file': 'total_fertility_rates.dat',
    # 'total_fertility_rates_migrant_file': 'total_fertility_rates_migrant.dat',
    'death_rates_m': 'death_rates_male.dat',
    'death_rates_f': 'death_rates_female.dat',

    'death_rates_m_migrant': 'death_rates_male_migrant.dat',
    'death_rates_f_migrant': 'death_rates_female_migrant.dat',

    'couple_prob': 0.06, # 0.09
    'leaving_prob': 0.005,
    'divorce_prob': 0.001, # 0.005
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

    'pop_size': 180,
    'migrant_pop_size': 33,
    'migrant_comm_id': 3,
    'growth_rate': 0.0,
    'imm_rate': 0,

    'com_dist': [0.2, 0.2, 0.6],
    'mobility_rates': [[0, 0.5, 0.5, 0],
                       #[0, 0.5, 0.3, 0.2], old row 1
                       [0, 0.65, 0.3, 0.05],
                       [0, 0.3, 0.7, 0],
                       [0, 0.05, 0, 0.95]],
    'healthcare_access': [0.5, 1, 1],
    'new_remote_access': 0.5,
    'new_village_access': 1,
    'thai_access': 1,
    'migrant_access': 0.6,
    'new_migrant_access': 0.6,
    'thai_mobility': 0.02,
    'migrant_mobility': 1,
    'new_migrant_mobility': 1,
    'cross_origin_couple_rate': 0.5,

    'thai_fertility': 1,
    'migrant_fertility': 1,
    'max_migrant_age': 101,
    'generation_threshold': 1,
    'separate_growth': True,
    'separate_death': True,

    'preg': False,
    'use_parity': False,
    'dyn_rates': True,
    'update_demog': True,

    'growth_rate_file': 'growth_rates.dat',
    'migrant_growth_rate_file': 'growth_rates_migrant.dat',
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
    'acu_ht_prob': 0.7,
    'acu_mtc_prob': 0.7,
    'chr_mtc_prob': 0.9,
    'death_rate': 0.000189,
    # 'death_rate': 0,
    'imm_prevalence': 0.03,

    # intervention parameters
    'pmtct_cover': 0.85,
    'vac_cover' : 0.98,
    'treat_rate': 0.075,
    'new_treat_rate': 0.075,

    'start_ratio': 0.80,
    'start_period': 20,

    'start_treat_ratio': 0.0,
    'start_treat_period': 40,

    'q': 0.0000046,
    'q_nl': 2.7,
    'nl_power': 2,
    'q_h': 0.0001,

    'cover': 0.8,

    't_offset': 65,

    'external_exposure_rate': 0,#5e-6,

    'random_seed': False,
    'seed': 0,
    't_per_year': 13,
    'years': [0, 60],
    'year_now': 40,
    'burn_in': 80,
    'burn_in_t_per_year': 1,
    'epi_burn_in': 60,

    'halt': False,

    # run parameters
    'num_runs': 10,
    'initial_cases': 145,
    'output_list': ['all'],
    'save_cp': True,
    'logging': False,
}

p['demo_burn'] = p['burn_in']# + p['epi_burn_in']
p['origin_access'] = {Origin.THAI: p['thai_access'], Origin.MIGRANT: p['migrant_access']}
p['origin_mobility'] = {Origin.THAI: p['thai_mobility'], Origin.MIGRANT: p['migrant_mobility']}
