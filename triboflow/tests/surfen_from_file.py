from fireworks import Firework, Workflow
from fireworks import LaunchPad

from triboflow.firetasks.surfen import FT_SurfEnFromFile

filename = '/home/firat/al100.cif'
miller = (1, 0, 0)
mpid = 'mp-134'
functional = 'PBE'

db_file = 'auto'
high_level = True
custom_id = 'Al_100_cif_from_CRY'
comp_params = {}

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
