repo_path = 'E:/hepb'
import os, sys
import pandas as pd
sys.path.append(os.path.join(repo_path, 'simodd-pop'))
sys.path.append(os.path.join(repo_path, 'simodd-dis'))
sys.path.append(os.path.join(repo_path, 'simodd-hepb'))
print(sys.path)

from population.utils import create_path

from disease.general.contact_matrix import ContactMatrix
from disease.general.run import go_single
from hepb.simulation_com import SimEpiCom

from disease.experiments.param_combo import ParamComboIt
from hepb.disease_hepb import DiseaseHepB
from hepb.obs_states import StateObserver
from hepb.obs_prev import PrevalenceObserver
from params import p


class DiseaseModel(DiseaseHepB):
    """
    Local version of SEIDS disease, adding observers and vaccines specific 
    to this set of experiments.
    """

    def __init__(self, p, cmatrix, rng, fname, mode='w'):
        super(DiseaseModel, self).__init__(p, cmatrix, rng, fname, mode)

        self.add_observers(StateObserver(self.h5file, self.state_labels()),
                           PrevalenceObserver(self.h5file, self.state_labels()))    # observers track various statistics during sim

# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# - # - MAIN  - # - #
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

if __name__ == '__main__':

    # setup base parameters
    p['pop_prefix'] = p['prefix']
    p['epi_prefix'] = p['prefix']
    p['overwrite'] = True
    p['overwrite_cp'] = True
    p['save_cp'] = False
    # create_path(p['prefix'])
    # print(p)

    # (basic usage) run simulation
    # go_single(p, DiseaseModel, ContactMatrix(), p['seed'], verbose=True)

    # sweep parameters
    sweep_params = [
    ]

    # generate parameter combinations (converting iterator to list)
    param_combos = list(ParamComboIt(p, sweep_params))
    # ethiopia_matrix = pd.read_csv('data/ethiopia-contact-matrix.csv', header=None).to_numpy()

    # just for info, 
    for x in [p]:
        x['q_h'] = x['q'] * 1.8
        # print out prefix and seed (used as output directory) 
        print(x['prefix'], x['seed'])
        # then run simulation
        cmatrix = ContactMatrix()
        go_single(x, DiseaseModel, cmatrix, x['seed'], sim_type=SimEpiCom, verbose=False)

    
