from fireworks import Firework, Workflow
from fireworks import LaunchPad

from triboflow.firetasks.surfen_tools import FT_SurfEnFromFile

#filename = '/home/fs71411/firaty/al100.cif'
#miller = (1, 0, 0)
mpid = 'mp-149'
functional = 'PBE'

db_file = 'auto'
high_level = True
#custom_id = 'Al_100_cif_CRY_425_12'
comp_params = {'encut': 400, 'k_dens': 6}

millers = [(1,0,0), (1,1,0), (1,1,1), (2,1,0), (2,1,1), (2,2,1)]
lpad = LaunchPad.auto_load()
for miller in millers:
    millerstr = ''.join([str(i) for i in miller])
    filename = f'/home/fs71411/firaty/cifs/Si-({millerstr}).cif'
    custom_id = f'Si_{millerstr}_CRY_400_6'

    FT = FT_SurfEnFromFile(filename=filename,
                       miller=miller,
                       mpid=mpid,
                       functional=functional,
                       db_file=db_file,
                       high_level=high_level,
                       custom_id=custom_id,
                       comp_params=comp_params)

    FW = Firework(FT)
    WF1 = Workflow.from_Firework(FW)

    lpad = LaunchPad.auto_load()
    lpad.add_wf(WF1)
