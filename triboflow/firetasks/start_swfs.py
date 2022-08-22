from uuid import uuid4

from atomate.utils.utils import env_chk
from fireworks import explicit_serialize, FiretaskBase, FWAction
from pymatgen.core import Structure
from pymatgen.core.surface import Slab
from pymatgen.core.interface import Interface

from triboflow.utils.database import Navigator, StructureNavigator
from triboflow.utils.structure_manipulation import interface_name
from triboflow.workflows.subworkflows import (adhesion_energy_swf,
                                              surface_energy_swf,
                                              calc_pes_swf, calc_ppes_swf,
                                              converge_swf,
                                              dielectric_constant_swf,
                                              make_and_relax_slab_swf,
                                              charge_analysis_swf)


@explicit_serialize
class FT_StartChargeAnalysisSWF(FiretaskBase):
    """Start an charge redistribution analysis subworkflow.

    Take an interface from the high_level_db and compute
    the charge density redistribution through a subworkflow.

    Parameters
    ----------
    mp_id_1 : str
        MaterialsProject ID number for the first material
    mp_id_2 : str
        MaterialsProject ID number for the second material
    miller_1 : str or [int]
        Miller indices of the first material
    miller_2 : str or [int]
        Miller indices of the second material
    functional : str
        Functional with which the workflow is run. PBE or SCAN.
    db_file : str, optional
        Full path of the db.json file to be used. The default is to use
        env_chk to find the file.
    high_level_db : str or True, optional
        Name of the high_level database to use. Defaults to 'True', in which
        case it is read from the db.json file.
    interface_label : str, optional
        Label that the relaxed interface has in the high-level database. Can
        be either structure@min (default), or structure@max at the moment.
    """
    required_params = ['mp_id_1', 'mp_id_2', 'miller_1', 'miller_2',
                       'functional']
    optional_params = ['db_file', 'high_level_db', 'interface_label']

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
        interface_label = self.get('interface_label', 'relaxed_structure@min')

        nav = Navigator(db_file, high_level=hl_db)

        name = interface_name(mp_id_1, miller_1, mp_id_2, miller_2)

        interface_dict = nav.find_data(collection=functional + '.interface_data',
                                       fltr={'name': name})

        redistribution_was_calculated = interface_dict.get('charge_density_redist')
        comp_params = interface_dict.get('comp_parameters', {})

        if not redistribution_was_calculated:
            interface = Interface.from_dict(interface_dict[interface_label])

            SWF = charge_analysis_swf(interface = interface,
                                      interface_name = name,
                                      functional = functional,
                                      db_file = db_file,
                                      high_level_db = hl_db,
                                      comp_parameters = comp_params)

            return FWAction(detours=SWF)
        else:
            return FWAction(update_spec=fw_spec)

@explicit_serialize
class FT_StartAdhesionSWF(FiretaskBase):
    """Start an adhesion subworkflow.

    Take relaxed top and bottom slabs of an interface, as well as the relaxed
    interface structure (by default the one with the lowest energy) and compute
    the adhesion energy through a subworkflow.

    Parameters
    ----------
    mp_id_1 : str
        MaterialsProject ID number for the first material
    mp_id_2 : str
        MaterialsProject ID number for the second material
    miller_1 : str or [int]
        Miller indices of the first material
    miller_2 : str or [int]
        Miller indices of the second material
    functional : str
        Functional with which the workflow is run. PBE or SCAN.
    db_file : str, optional
        Full path of the db.json file to be used. The default is to use
        env_chk to find the file.
    adhesion_handle: str, optional
        Flag under which the adhesion energy will be saved in the interface_data
        collection of the high_level database.
    high_level_db : str or True, optional
        Name of the high_level database to use. Defaults to 'True', in which
        case it is read from the db.json file.
    """
    required_params = ['mp_id_1', 'mp_id_2', 'miller_1', 'miller_2',
                       'functional']
    optional_params = ['db_file', 'adhesion_handle', 'high_level_db']

    def run_task(self, fw_spec):
        mp_id_1 = self.get('mp_id_1')
        mp_id_2 = self.get('mp_id_2')
        miller_1 = self.get('miller_1')
        miller_2 = self.get('miller_2')
        functional = self.get('functional')
        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        adhesion_handle = self.get('adhesion_handle', 'adhesion_energy@min')
        hl_db = self.get('high_level_db', True)

        nav = Navigator(db_file, high_level=hl_db)

        name = interface_name(mp_id_1, miller_1, mp_id_2, miller_2)

        interface_dict = nav.find_data(collection=functional + '.interface_data',
                                       fltr={'name': name})

        adhesion_was_calculated = interface_dict.get(adhesion_handle)
        comp_params = interface_dict.get('comp_parameters', {})

        if not adhesion_was_calculated:
            top_slab = Slab.from_dict(interface_dict['top_aligned_relaxed'])
            bottom_slab = Slab.from_dict(interface_dict['bottom_aligned_relaxed'])
            interface = Structure.from_dict(interface_dict['relaxed_structure@min'])

            SWF = adhesion_energy_swf(top_slab,
                                      bottom_slab,
                                      interface,
                                      interface_name=name,
                                      functional=functional,
                                      comp_parameters=comp_params)

            return FWAction(detours=SWF)
        else:
            return FWAction(update_spec=fw_spec)


@explicit_serialize
class FT_StartBulkConvoSWF(FiretaskBase):
    """ Starts a convergence subworkflow.

    Starts either an energy cutoff or kpoint density convergence of a material
    with a given MPID and functional through a subworkflow.

    Parameters
    ----------
    conv_type : str
        Either "kpoints" or "encut", depending on what to converge.
    mp_id : str
        MaterialsProject ID number for the material
    functional : str
        Functional with which the workflow is run. PBE or SCAN.
    db_file : str, optional
        Full path of the db.json file to be used. The default is to use
        env_chk to find the file.
    encut_start : float, optional
        Starting encut value for the first run. Defaults to the largest EMIN
        in the POTCAR.
    encut_incr : float, optional
        Increment for the encut during the convergence. Defaults to 25.
    k_dens_start : float, optional
        Starting kpoint density in 1/Angstrom. Defaults to 1.0
    k_dens_incr : float, optional
        Increment for the kpoint convergence. Can be set quite small since
        there is a check in place to see if a new mesh is actually constructed
        for each density. Defaults to 0.1.
    n_converge : int, optional
        Number of calculations that have to be inside the convergence
        threshold for convergence to be reached. Defaults to 3.
    high_level_db : str or True, optional
        Name of the high_level database to use. Defaults to 'True', in which
        case it is read from the db.json file.
    """

    _fw_name = 'Start Encut or Kdensity Convergence'
    required_params = ['conv_type', 'mp_id', 'functional']
    optional_params = ['db_file', 'encut_start', 'encut_incr', 'k_dens_start',
                       'k_dens_incr', 'n_converge', 'high_level_db']

    def run_task(self, fw_spec):

        conv_type = self.get('conv_type')
        mp_id = self.get('mp_id')
        functional = self.get('functional')
        db_file = self.get('db_file')
        n_converge = self.get('n_converge', 3)
        encut_start = self.get('encut_start', None)
        encut_incr = self.get('encut_incr', 25)
        k_dens_start = self.get('k_dens_start', 2.0)
        k_dens_incr = self.get('k_dens_incr', 0.1)
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        hl_db = self.get('high_level_db', True)

        if conv_type not in ['kpoints', 'encut']:
            raise ValueError('"type" input must be either "kpoints" or'
                             '"encut".\nYou have passed {}'.format(conv_type))

        nav_structure = StructureNavigator(
            db_file=db_file,
            high_level=hl_db)
        data = nav_structure.get_bulk_from_db(
            mp_id=mp_id,
            functional=functional)

        if conv_type == 'encut':
            stop_convergence = data.get('encut_info')
        elif conv_type == 'kpoints':
            stop_convergence = data.get('k_dense_info')

        if not stop_convergence:
            structure_dict = data.get('structure_equiVol')
            if not structure_dict:
                structure_dict = data.get('primitive_structure')
                if not structure_dict:
                    structure_dict = data.get('structure_fromMP')
                    if not structure_dict:
                        raise LookupError('No structure found that can be used '
                                          'as input for the convergence swf.')
            structure = Structure.from_dict(structure_dict)
            comp_params = data.get('comp_parameters', {})
            SWF = converge_swf(structure=structure,
                               conv_type=conv_type,
                               flag=mp_id,
                               comp_parameters=comp_params,
                               functional=functional,
                               encut_start=encut_start,
                               encut_incr=encut_incr,
                               k_dens_start=k_dens_start,
                               k_dens_incr=k_dens_incr,
                               n_converge=n_converge,
                               print_help=False)
            return FWAction(detours=SWF, update_spec=fw_spec)
        else:
            return FWAction(update_spec=fw_spec)


@explicit_serialize
class FT_StartSurfaceEnergySWF(FiretaskBase):
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


@explicit_serialize
class FT_StartDielectricSWF(FiretaskBase):
    """ Starts a dielectric subworkflow.

    Starts a subworkflow that calculates and updates the dielectric constant
    for a given bulk material. Will also by default update all slabs with the
    same mpid.

    Parameters
    ----------
    mp_id : str
        MaterialsProject ID number for the material
    functional : str
        Functional with which the workflow is run. PBE or SCAN.
    db_file : str, optional
        Full path of the db.json file to be used. The default is to use
        env_chk to find the file.
    high_level_db : str or True, optional
        Name of the high_level database to use. Defaults to 'True', in which
        case it is read from the db.json file.
    update_bulk : bool, optional
        If the bulk entry for the given mpid should be updated.
        The default is True.
    update_slabs : bool, optional
        If the slab entries matching a given mpid should be updated (all miller
        indices. The default is False.
    """

    _fw_name = 'Start Encut or Kdensity Convergence'
    required_params = ['mp_id', 'functional']
    optional_params = ['db_file', 'high_level_db', 'update_bulk', 'update_slabs']

    def run_task(self, fw_spec):

        mp_id = self.get('mp_id')
        functional = self.get('functional')
        db_file = self.get('db_file')
        update_bulk = self.get('update_bulk', True)
        update_slabs = self.get('update_slabs', True)
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        hl_db = self.get('high_level_db', True)

        nav_structure = StructureNavigator(
            db_file=db_file,
            high_level=hl_db)
        data = nav_structure.get_bulk_from_db(
            mp_id=mp_id,
            functional=functional)

        structure_dict = data.get('structure_equiVol')
        if not structure_dict:
            structure_dict = data.get('primitive_structure')
            if not structure_dict:
                structure_dict = data.get('structure_fromMP')
                if not structure_dict:
                    raise LookupError('No structure found that can be used '
                                      'as input for the convergence swf.')
        structure = Structure.from_dict(structure_dict)
        comp_params = data.get('comp_parameters', {})
        flag = 'test_dielectric_WF_' + str(uuid4())
        SWF = dielectric_constant_swf(structure=structure,
                                      mpid=mp_id,
                                      flag=flag,
                                      comp_parameters=comp_params,
                                      functional='PBE',
                                      db_file=db_file,
                                      hl_db=hl_db,
                                      update_bulk=update_bulk,
                                      update_slabs=update_slabs)
        return FWAction(detours=SWF, update_spec=fw_spec)


@explicit_serialize
class FT_StartPESCalcSWF(FiretaskBase):
    """ Start a PES subworkflow.

    Starts a PES subworkflow using data from the high-level database.
    This is intended to be used to start a PES subworkflow from a main
    workflow.

    Parameters
    ----------
    mp_id_1 : str
        Materials Project database ID for the first material of the interface.
    mp_id_2 : str
        Materials Project database ID for the second material of the interface.
    miller_1 : list of int or str
        Miller index of the first material.
    miller_2 : list of int or str
        Miller index of the second material.
    functional : str
        Which functional to use; has to be 'PBE' or 'SCAN'.
    db_file : str, optional
        Full path to the db.json file that should be used. Defaults to
        '>>db_file<<', to use env_chk.

    Returns
    -------
    FWAction that produces a detour PES subworkflow.
    """
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

        name = interface_name(mp_id_1, miller_1, mp_id_2, miller_2)

        nav_structure = StructureNavigator(
            db_file=db_file,
            high_level=hl_db)
        interface_dict = nav_structure.get_interface_from_db(
            name=name,
            functional=functional
        )
        comp_params = interface_dict['comp_parameters']
        interface = Interface.from_dict(interface_dict['unrelaxed_structure'])
        already_done = interface_dict.get('relaxed_structure@min')

        if not already_done:
            SWF = calc_pes_swf(interface=interface,
                               interface_name=name,
                               functional=functional,
                               comp_parameters=comp_params,
                               output_dir=None)

            return FWAction(detours=SWF, update_spec=fw_spec)

        else:
            return FWAction(update_spec=fw_spec)


@explicit_serialize
class FT_StartPPESWF(FiretaskBase):
    """
    Start a CalcPPES_SWF subworkflow that calculates a PPES.

    The workflow is only added if there are not already relevant results in
    the high-level database.

    Parameters
    ----------
    interface_name : str
        Name of the interface in the high-level database.
    functional : str
        Which functional to use; has to be 'PBE' or 'SCAN'.
    distance_list : list of float, optional
        Modification of the equilibrium distance between the slabs.
        The default is [-0.5, -0.25, 0.0, 0.25, 0.5, 2.5, 3.0, 4.0, 5.0, 7.5].
    out_name : str, optional
        Name for the PPES data in the high-level database. The default is
        'PPES@minimum'.
    structure_name : str, optional
        Name of the structure in the interface entry to the high-level database
        for which the PPES should be calculated. The default is
        'minimum_relaxed'.
    spec : dict, optional
        fw_spec that can be passed to the SWF and will be passed on. The
        default is {}.

    Returns
    -------
    SWF : fireworks.core.firework.Workflow
        Subworkflow to calculate the PPES for a certain interface.

    """

    required_params = ['interface_name', 'functional', 'distance_list']
    optional_params = ['db_file', 'structure_name', 'out_name', 'high_level_db']

    def run_task(self, fw_spec):

        name = self.get('interface_name')
        functional = self.get('functional')
        tag = self.get('tag')

        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)

        structure_name = self.get('structure_name', 'minimum_relaxed')
        out_name = self.get('out_name', 'PPES@minimum')

        d_list = self.get('distance_list')
        hl_db = self.get('high_level_db', True)

        nav_structure = StructureNavigator(
            db_file=db_file,
            high_level=hl_db)
        interface_dict = nav_structure.get_interface_from_db(
            name=name,
            functional=functional)

        calc_PPES = True
        if interface_dict.get('PPES') is not None:
            if interface_dict['PPES'].get(out_name) is not None:
                print('\n A PPES-object with out_name: ' + out_name +
                      '\n has already been created in the interface entry: ' +
                      name + '\n for the ' + functional + ' functional.')
                calc_PPES = False

        if calc_PPES:
            SWF = calc_ppes_swf(interface_name=name,
                                functional=functional,
                                distance_list=d_list,
                                out_name=out_name,
                                structure_name=structure_name,
                                spec=fw_spec)

            return FWAction(additions=SWF, update_spec=fw_spec)
        else:
            return FWAction(update_spec=fw_spec)


@explicit_serialize
class FT_StartSlabRelaxSWF(FiretaskBase):
    """Start a subworkflow as a detour to make a slab and relax it.

    Parameters
    ----------
    mp_id : str
        ID number for structures in the material project.
    miller : list of int or str
        Miller indices for the slab generation. Either single str e.g. '111',
        or a list of int e.g. [1, 1, 1]
    functional : str
        functional that is used for the calculation.
    db_file : str, optional
        Full path to the db.json file which holds the location and access
        credentials to the database. If not given uses env_chk.
    slab_struct_name : str, optional
        Name of the slab to be put in the DB (identified by mp_id and miller).
        Defaults to 'unrelaxed_slab'.
    relax_type : str
        Relaxation type for the get_custom_vasp_relax_settings helper_function.
    bulk_struct_name : str, optional
        Name of the bulk structure in the bulk database (material is
        identified by mp_id, but there might be different structures of the
        same material.) Defaults to 'structure_equiVol'.
    slab_out_name : str, optional
        Name of the slab to be put in the DB (identified by mp_id and miller).
        Defaults to 'relaxed_slab'.

    Returns
    -------
        Starts a new subworkflow as a detour to the current workflow.
    """

    required_params = ['mp_id', 'miller', 'functional']
    optional_params = ['db_file', 'slab_struct_name', 'relax_type',
                       'bulk_struct_name', 'slab_out_name', 'high_level_db']

    def run_task(self, fw_spec):

        mp_id = self.get('mp_id')
        if type(self['miller']) == str:
            miller = [int(k) for k in list(self['miller'])]
        else:
            miller = self['miller']

        functional = self.get('functional')

        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)

        slab_name = self.get('slab_struct_name', 'unrelaxed_slab')
        relax_type = self.get('relax_type', 'slab_pos_relax')
        bulk_name = self.get('bulk_struct_name', 'structure_equiVol')
        slab_out_name = self.get('slab_out_name', 'relaxed_slab')
        hl_db = self.get('high_level_db', True)

        nav_structure = StructureNavigator(
            db_file=db_file,
            high_level=hl_db)

        data = nav_structure.get_bulk_from_db(
            mp_id=mp_id,
            functional=functional)
        bulk_struct = Structure.from_dict(data[bulk_name])

        slab_data = nav_structure.get_slab_from_db(
            mp_id=mp_id,
            functional=functional,
            miller=miller)
        comp_params = slab_data.get('comp_parameters')
        min_thickness = slab_data.get('min_thickness', 10)
        min_vacuum = slab_data.get('min_vacuum', 25)

        WF = make_and_relax_slab_swf(bulk_structure=bulk_struct,
                                     miller_index=miller,
                                     flag=mp_id,
                                     comp_parameters=comp_params,
                                     functional=functional,
                                     min_thickness=min_thickness,
                                     min_vacuum=min_vacuum,
                                     relax_type=relax_type,
                                     slab_struct_name=slab_name,
                                     out_struct_name=slab_out_name,
                                     print_help=False)

        return FWAction(detours=WF)
