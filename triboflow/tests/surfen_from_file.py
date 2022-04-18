from fireworks import Firework, Workflow
from fireworks import LaunchPad

from triboflow.firetasks.surfen_tools import FT_StartSurfaceEnergyFromFile

encut = 400
k_dens = 6
comp_params_user = {'encut': encut, 'k_dens': k_dens}

common_kwargs = {'functional': 'PBE',
                 'db_file': 'auto',
                 'high_level': 'surfen_test',
                 'comp_params': comp_params_user}

surfen_from_file_inp = [{'filename': '/home/yalcin/scratch/fe310.cif',
                         'miller': (3, 1, 0),
                         'mpid': 'mp-13',
                         'functional': 'PBE',
                         'custom_id': 'fe310_crys',
                         }]

lpad = LaunchPad.auto_load()
for inp in surfen_from_file_inp:
    inp.update(common_kwargs)
    millerstr = ''.join([str(i) for i in inp['miller']])

    FW = Firework(FT_StartSurfaceEnergyFromFile(**inp))
    WF = Workflow.from_Firework(FW)

    lpad.add_wf(WF)
