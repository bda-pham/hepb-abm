"""
A stub for testing/evaluating dynamic demography simulations.
"""
repo_path = 'E:/hepb'
import os, sys
sys.path.append(os.path.join(repo_path, 'simodd-pop'))
sys.path.append(os.path.join(repo_path, 'simodd-dis'))
sys.path.append(os.path.join(repo_path, 'simodd-hepb'))


from population.simulation import Simulation
from disease.general.contact_matrix import ContactMatrix
from population.utils import create_path
from observers.obs_pop import PopulationObserver
# from pop_explore.output_pop import *
from disease.experiments.param_combo import ParamComboIt
from disease.general.ind_epi import IndEpi
from params import p
from hepb.simulation_com import SimEpiCom
from hepb.disease_hepb import DiseaseHepB
from random import Random

def run_single(p, cur_seed):
    create_path(p['prefix'])
    # print("Creating population...")
    iterations = (p['years'][1] + p['epi_burn_in'] + p['burn_in']) * p['burn_in_t_per_year']
    #pop_fname = os.path.join(p['prefix'], 'population.hd5')
        
    rng = Random(cur_seed)
    disease = DiseaseHepB(p, ContactMatrix(), rng, p['prefix']+"/population.hd5", mode='w')
    sim = SimEpiCom(p, disease, rng, ind_type=IndEpi)
    sim.add_observers(PopulationObserver(sim.h5file))

    # print("Running simulation...")
    # print("iter\tyears\tdays\tpeople\thouses\tbirths\tdeaths\tmigrants")
    sim.create_population()
    for i in range(iterations):
        # print "burn in; year", i
        t = i * 364 // p['burn_in_t_per_year']
        b, d, im, b2 = sim.update_all_demo(i)
        sim.update_observers(t, pop=sim.P)
        # print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(i, t // 364, t % 364,
        #     len(sim.P.I), len(sim.P.groups['household']),
        #     len(b), len(d), len(im)))
        if b == "error":
            print("error with birth!")
            break

    # print("Simulation complete!")
    sim.done(complete=True)
    print("Finished, output at "+p["prefix"])
    return sim

if __name__ == '__main__':
    # sweep_params = [
    #     {'name': 'couple_prob', 'values': [0.04]},
    #     {'name': 'leaving_prob', 'values': [0.005]},
    #     {'name': 'divorce_prob', 'values': [0.001]}
    # ]
    sweep_params = [
        # {'name': 'couple_prob', 'values': [0.03]},
        # {'name': 'leaving_prob', 'values': [0.01]},
        # {'name': 'divorce_prob', 'values': [0.01]}
        {'name': 'growth_rate', 'values': [0]}
    ]
    p['t_per_year'] = p['burn_in_t_per_year']
    
    p['num_runs'] = 10

    # generate parameter combinations (converting iterator to list)
    param_combos = list(ParamComboIt(p, sweep_params))


    # just for info, 
    for i, x in enumerate(param_combos):
        cur_sim = run_single(x, i)
