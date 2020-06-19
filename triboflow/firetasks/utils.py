""" Utility Firetasks for the triboflow project.

Created on Wed Jun 17 15:59:59 2020
@author: mwo
"""

from fireworks import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize

@explicit_serialize
class FT_PassSpec(FiretaskBase):
    """Update only certain keys in the first level of the spec.
    
    If the key_list contatins only '_all', update the whole spec!
    """
    
    _fw_name = 'Pass Spec'
    required_params=['key_list']
    def run_task(self, fw_spec):
        update = {}
        if self['key_list'] == ['_all']:
            spec = fw_spec
            return FWAction(update_spec = spec)
        else:
            for k in self['key_list']:
                if fw_spec.get(k) is None:
                    raise ValueError('{} can not be passed on to the next'
                                     'FT/FW because it is not in the spec.\n'
                                     'Currently the spec has the keys:\n'
                                     '{}'.format(k, fw_spec.keys()))
                update[k] = fw_spec.get(k)
            return FWAction(update_spec = update)


@explicit_serialize
class FT_PrintSpec(FiretaskBase):
    """Prints the spec of the current Workflow to the screen.
    
    Not only prints the current spec in a pretty way, but also returns a
    FWAction that updates the spec of future to include the current spec.
    """
    
    _fw_name = 'Print Spec'
    
    def run_task(self, fw_spec):
        import pprint
        
        pprint.pprint(fw_spec)
        
        spec = fw_spec
        return FWAction(update_spec = spec)


