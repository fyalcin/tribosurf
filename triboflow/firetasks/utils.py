""" Utility Firetasks for the triboflow project.

Created on Wed Jun 17 15:59:59 2020
@author: mwo
"""

from pprint import pprint

from fireworks import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize
from atomate.utils.utils import env_chk

from triboflow.utils.database import Navigator, StructureNavigator
from triboflow.utils.structure_manipulation import interface_name
from triboflow.utils.surfen_tools import move_result


@explicit_serialize
class FT_MoveResults(FiretaskBase):
    """
    Firetask to move the result of a VASP calculation from the Fireworks database to the destination.

    Parameters
    ----------
    tag : str
        Task label to query for in the tasks collection of the Fireworks database.
    fltr : dict
        Filter dictionary that is used to query for the destination in the database.
    coll : str
        Collection to move the results into in the destination database.
    loc : str
        Location in the collection the results should be written to.
    custom_dict : dict, optional
        Custom dictionary to write into the destination. The default is {}.
    db_file : str, optional
        Full path of the db.json. The default is 'auto'.
    high_level : str, optional
        Whether to query the results from the high level database or not. The default is True.

    Returns
    -------
    FWAction that just updates the spec.

    """
    _fw_name = "Move the results from the tag to the entry found by flag."
    required_params = ['tag', 'fltr', 'coll', 'loc']
    optional_params = ['custom_dict', 'db_file', 'high_level']

    def run_task(self, fw_spec):
        tag = self.get('tag')
        fltr = self.get('fltr')
        coll = self.get('coll')
        loc = self.get('loc')

        custom_dict = self.get('custom_dict', {})
        db_file = self.get('db_file', 'auto')
        high_level = self.get('high_level', True)

        move_result(tag, fltr, coll, loc, custom_dict, db_file, high_level)

        return FWAction(update_spec=fw_spec)


@explicit_serialize
class FT_ChooseCompParams(FiretaskBase):
    """Choose the computational parameters for the interface calculations.
    
    Here computational parameters are selected for the interface based on the
    converged parameters of the two constituent materials.
    
    Parameters
    ----------
    loc_1 : list of str
        Location in the spec of the converged parameters of the first material.
    loc_2 : list of str
        Location in the spec of the converged parameters of the second
        material.
    out_loc : list of str
        Location in the spec as specified by a list of keys into which the
        output should be written.
        
    Returns
    -------
    FWAction that modifies the spec for the next firewoks at the out_loc.
    """

    _fw_name = 'Choose Comp Params'
    required_params = ['mp_id_1', 'mp_id_2', 'miller_1', 'miller_2',
                       'functional']
    optional_params = ['db_file', 'high_level_db']

    def run_task(self, fw_spec):
        mp_id_1 = self.get('mp_id_1')
        mp_id_2 = self.get('mp_id_2')
        miller_1 = self.get('miller_1')
        miller_2 = self.get('miller_2')
        functional = self.get('functional')

        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        hl_db = self.get('high_level_db', True)

        nav_structure = StructureNavigator(
            db_file=db_file,
            high_level=hl_db)
        data_1 = nav_structure.get_bulk_from_db(
            mp_id=mp_id_1,
            functional=functional)
        data_2 = nav_structure.get_bulk_from_db(
            mp_id=mp_id_2,
            functional=functional)

        k_dens1 = data_1['comp_parameters']['k_dens']
        k_dens2 = data_2['comp_parameters']['k_dens']
        encut1 = data_1['comp_parameters']['encut']
        encut2 = data_2['comp_parameters']['encut']

        metal_1 = data_1['comp_parameters']['is_metal']
        metal_2 = data_2['comp_parameters']['is_metal']

        k_dens = max(k_dens1, k_dens2)
        encut = max(encut1, encut2)
        metal = any((metal_1, metal_2))

        name = interface_name(mp_id_1, miller_1, mp_id_2, miller_2)

        nav_high = Navigator(
            db_file=db_file,
            high_level=hl_db)
        nav_high.update_data(
            collection=functional + '.interface_data',
            fltr={'name': name},
            new_values={'$set': {'comp_parameters.k_dens': k_dens,
                                 'comp_parameters.encut': encut,
                                 'comp_parameters.is_metal': metal}})
        return


@explicit_serialize
class FT_PrintFromBulkDB(FiretaskBase):
    _fw_name = 'Print bulk data from DB'
    required_params = ['mp_id', 'functional']
    optional_params = ['db_file', 'high_level_db']

    def run_task(self, fw_spec):
        db_file = self.get('db_file', env_chk('>>db_file<<', fw_spec))
        mp_id = self['mp_id']
        functional = self['functional']

        hl_db = self.get('high_level_db', True)

        nav_structure = StructureNavigator(
            db_file=db_file,
            high_level=hl_db)
        bulk_dict = nav_structure.get_bulk_from_db(
            mp_id=mp_id,
            functional=functional)
        print('')
        pprint(bulk_dict)
        print('')


@explicit_serialize
class FT_PassSpec(FiretaskBase):
    """Update only certain keys in the first level of the spec.
    
    If the key_list contatins only '_all', update the whole spec!
    """

    _fw_name = 'Pass Spec'
    required_params = ['key_list']

    def run_task(self, fw_spec):

        update = {}
        if self['key_list'] == ['_all']:
            spec = fw_spec
            return FWAction(update_spec=spec)
        else:
            for k in self['key_list']:
                if fw_spec.get(k) is None:
                    raise ValueError('{} can not be passed on to the next'
                                     'FT/FW because it is not in the spec.\n'
                                     'Currently the spec has the keys:\n'
                                     '{}'.format(k, fw_spec.keys()))
                update[k] = fw_spec.get(k)
            return FWAction(update_spec=update)


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

        return FWAction(update_spec=spec)
