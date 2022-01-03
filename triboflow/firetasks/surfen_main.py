from fireworks import explicit_serialize, FiretaskBase, FWAction

from triboflow.workflows.subworkflows import surface_energy_swf


@explicit_serialize
class FT_StartSurfaceEnergy(FiretaskBase):
    _fw_name = 'Starts a subworkflow that calculates surface energies as detour'
    required_params = ['mpid', 'functional', 'sg_params', 'sg_filter']
    optional_params = ['db_file', 'high_level', 'comp_params_user', 'custom_id']

    def run_task(self, fw_spec):
        mpid = self.get('mpid')
        functional = self.get('functional')
        sg_params = self.get('sg_params')
        sg_filter = self.get('sg_filter')

        db_file = self.get('db_file', 'auto')
        high_level = self.get('high_level', True)
        comp_params_user = self.get('comp_params_user', {})
        custom_id = self.get('custom_id', None)

        WF = surface_energy_swf(mpid=mpid,
                                functional=functional,
                                sg_params=sg_params,
                                sg_filter=sg_filter,
                                db_file=db_file,
                                high_level=high_level,
                                comp_params_user=comp_params_user,
                                custom_id=custom_id)

        return FWAction(detours=WF, update_spec=fw_spec)