from fireworks import FWAction, FiretaskBase, Firework, Workflow
from fireworks.utilities.fw_utilities import explicit_serialize

from triboflow.firetasks.utils import FT_MoveResults
from triboflow.utils.database import Navigator
from triboflow.utils.structure_manipulation import slab_from_file
from triboflow.utils.surfen_tools import get_surfen_inputs_from_mpid, generate_surfen_wfs_from_inputs, \
    write_surface_energies_to_db, put_surfen_inputs_into_db, generate_surfen_entries, get_entry_by_loc, \
    get_surfen_inputs_from_slab


@explicit_serialize
class FT_SurfEnFromFile(FiretaskBase):
    _fw_name = 'Calculates the surface energy of the slab loaded from a file.'
    required_params = ['filename', 'miller', 'mpid', 'functional']
    optional_params = ['db_file', 'high_level', 'custom_id', 'comp_params']

    def run_task(self, fw_spec):
        filename = self.get('filename')
        miller = self.get('miller')
        mpid = self.get('mpid')
        functional = self.get('functional')

        db_file = self.get('db_file', 'auto')
        high_level = self.get('high_level', True)
        custom_id = self.get('custom_id')
        comp_params_user = self.get('comp_params')

        nav = Navigator(db_file, high_level)
        comp_params = nav.find_data('PBE.bulk_data', {'mpid': mpid}).get('comp_parameters')
        comp_params.update(comp_params_user)

        slab, SG = slab_from_file(filename, mpid, functional, miller, db_file, high_level)

        inputs_list = get_surfen_inputs_from_slab(slab, SG, custom_id=custom_id)
        if not isinstance(inputs_list, list):
            inputs_list = [inputs_list]
        inputs_list = generate_surfen_wfs_from_inputs(inputs_list, comp_params)

        fltr = {'mpid': mpid}
        coll = f'{functional}.slab_data.LEO'

        FW1 = Firework(
            FT_PutSurfenInputsIntoDB(inputs_list=inputs_list, sg_params={}, comp_params=comp_params,
                                     fltr=fltr, coll=coll, db_file=db_file, high_level=high_level),
            name=f"Generate surface energy inputs for {mpid} with {functional} and put in DB")

        FW2 = Firework(
            FT_RelaxSurfaceEnergyInputs(inputs_list=inputs_list, fltr=fltr, coll=coll, db_file=db_file,
                                        high_level=high_level),
            name=f"Generate and relax surface energy inputs for {mpid} with {functional}")

        FW3 = Firework(
            FT_WriteSurfaceEnergies(inputs_list=inputs_list, fltr=fltr, coll=coll, db_file=db_file,
                                    high_level=high_level),
            name=f"Calculate the surface energies for {mpid} with {functional} and put into DB")
        WF = Workflow([FW1, FW2, FW3], {FW1: [FW2], FW2: [FW3]})

        return FWAction(detours=WF)


@explicit_serialize
class FT_RelaxSurfaceEnergyInputs(FiretaskBase):
    """ Perform VASP calculations on the input structures in order to calculate
    surface energy.

    Parameters
    ----------
    inputs_list: list
        List with elements containing information about the inputs necessary for surface
        energy calculations for a given material represented by its Materials Project ID.
        Generated by the function generate_surfen_wfs_from_inputs that takes in a
        previous version of the inputs_list generated by get_surfen_inputs_from_slab
        by inputting a Slab and a SlabGenerator corresponding to this Slab object.
    fltr : dict
        Filter to use when looking up results in the database. Generally involves the mpid.
    coll : str
        Collection to query the results from in the database.
    db_file : str, optional
        Full path to the db.json file which holds the location and access
        credentials to the database. The default is 'auto'.
    high_level : str, optional
        Whether to query the results from the high level database or not. The default is True.
    Returns
    -------
        FWAction that detours to a Workflow containing the Fireworks that will perform
        the VASP calculations.
    """
    _fw_name = "Perform various VASP calculations on the structures in order to calculate surface energy."
    required_params = ['inputs_list', 'fltr', 'coll']
    optional_params = ['db_file', 'high_level']

    def run_task(self, fw_spec):
        inputs_list = self.get('inputs_list')
        fltr_high = self.get('fltr')
        coll = self.get('coll')
        db_file = self.get('db_file', 'auto')
        high_level = self.get('high_level', True)

        nav_high = Navigator(db_file, high_level=high_level)
        nav_low = Navigator(db_file, high_level=False)

        WF_list = []
        for slab in inputs_list:
            hkl = slab['slab_params']['hkl']
            for entry in slab.get('inputs'):
                tag = entry.get('tag')
                loc = entry.get('loc')

                # If the calculation is already in the high_level database, we skip this
                # particular slab
                result_high = get_entry_by_loc(nav_high, fltr_high, coll, loc).get('output')
                if result_high:
                    continue

                comp = str(entry.get('struct').composition)
                calc_tag = entry.get('calc_tag')
                WF_calc = Workflow.from_dict(entry.get('WF'))
                WF_move = Workflow([Firework(FT_MoveResults(tag=tag,
                                                            fltr=fltr_high,
                                                            coll=coll,
                                                            loc=loc,
                                                            db_file=db_file,
                                                            high_level=high_level),
                                             name=f"Move {comp}-{hkl}-{calc_tag} results FW.")],
                                   name=f"Move {comp}-{hkl}-{calc_tag} results WF.")

                fltr_low = {'task_label': tag}
                result_low = nav_low.find_data('tasks', fltr_low)

                if result_low:
                    WF_list.append(WF_move)
                else:
                    WF_calc.append_wf(WF_move, WF_calc.leaf_fw_ids)
                    WF_list.append(WF_calc)
        return FWAction(detours=WF_list, update_spec=fw_spec)


@explicit_serialize
class FT_SlabOptOrientation(FiretaskBase):
    """
    Starts a Workflow that will find the lowest energy orientation of a material
    desribed by its MaterialsProject ID given certain constraints.

    Parameters
    ----------
    mpid : str
        ID of the material in the MaterialsProject database.
    functional : TYPE
        Which functional to use; has to be 'PBE' or 'SCAN'.
    db_file : TYPE, optional
        Full path of the db.json. The default is 'auto'.
    high_level : TYPE, optional
        Whether to query the results from the high level database or not.
        The default is True.
    max_index : int, optional
        Maximum miller up to which unique orientations will be searched.
        The default is 2.
    comp_params : dict
        Computational parameters for the VASP simulations. If not set, default
        parameters will be used instead. The default is {}.
    sg_params : dict
        Parameters to be used in the SlabGenerator.
    bvs_method : str, optional
        Filtering method that is used in conjunction with the bond valence sums.
        'threshold' with a 'bvs_param' provided in the kwargs will filter out slabs
        with bond valence sums (1+bvs_param)*bvs_min where bvs_min is the minimum
        bond valence sum of all the slabs.
        'all' will proceed with all the slabs generated regardless of their bond
        valence sum.
        'min_N' with a 'bvs_param' will take the slabs with the bvs_param lowest
        bond valence sums.
        The default is 'threshold'.
    override : bool
        Whether or not to run the workflow even though there exists an
        "opt_orientation" key in the slab entry for this specific material.
        Useful when one wants to have access to the surface energies of different
        terminations in order to generate the Wulff shape for example.

    """

    _fw_name = "Generate necessary structures for surface energy calculation and relax them."
    required_params = ['mpid', 'functional', 'sg_params', 'sg_filter']
    optional_params = ['db_file', 'high_level', 'comp_params', 'override', 'fake']

    def run_task(self, fw_spec):
        mpid = self.get('mpid')
        functional = self.get('functional')
        sg_params = self.get('sg_params')
        sg_filter = self.get('sg_filter')
        print(sg_params)

        db_file = self.get('db_file', 'auto')
        high_level = self.get('high_level', True)
        # max_index = self.get('max_index', 1)
        comp_params_user = self.get('comp_params', {})

        # bvs_method = self.get('bvs_method', 'tolerance')
        # bvs_param = self.get('bvs_param', 0.2)

        override = self.get('override', False)
        fake = self.get('fake', False)

        nav_high = Navigator(db_file, high_level=True)

        comp_params = nav_high.find_data(f'{functional}.bulk_data', {'mpid': mpid})['comp_parameters']
        comp_params.update(comp_params_user)

        slab_dict = nav_high.find_data(f'{functional}.slab_data.LEO', {'mpid': mpid})

        if slab_dict:
            opt_orientation = slab_dict.get('opt_slab')
            if opt_orientation and not override:
                return FWAction(update_spec=fw_spec)

        inputs_list = get_surfen_inputs_from_mpid(mpid,
                                                  functional,
                                                  sg_params,
                                                  sg_filter,
                                                  db_file,
                                                  high_level)

        inputs_list = generate_surfen_wfs_from_inputs(inputs_list, comp_params, fake)

        coll = f'{functional}.slab_data.LEO'
        fltr = {'mpid': mpid}

        FW1 = Firework(
            FT_PutSurfenInputsIntoDB(inputs_list=inputs_list, sg_params=sg_params, comp_params=comp_params,
                                     fltr=fltr, coll=coll, db_file=db_file, high_level=high_level),
            name=f"Generate surface energy inputs for {mpid} with {functional} and put in DB")

        FW2 = Firework(
            FT_RelaxSurfaceEnergyInputs(inputs_list=inputs_list, fltr=fltr, coll=coll, db_file=db_file,
                                        high_level=high_level),
            name=f"Generate and relax surface energy inputs for {mpid} with {functional}")

        FW3 = Firework(
            FT_WriteSurfaceEnergies(inputs_list=inputs_list, fltr=fltr, coll=coll, db_file=db_file,
                                    high_level=high_level),
            name=f"Calculate the surface energies for {mpid} with {functional} and put into DB")

        FW4 = Firework(
            FT_GenerateSurfenEntries(fltr=fltr, coll=coll, db_file=db_file, high_level=high_level),
            name=f"Update the consolidated surface energy list for {mpid} with {functional}")

        WF = Workflow([FW1, FW2, FW3, FW4], {FW1: [FW2], FW2: [FW3], FW3: [FW4]},
                      name=f'Compute surface energies of likely low energy orientations of {mpid} with {functional}')

        return FWAction(detours=WF)


@explicit_serialize
class FT_WriteSurfaceEnergies(FiretaskBase):
    """ Calculates the surface energies and updates the LEO collection entries
    with them.

    Parameters
    ----------
    inputs_list: list
        List with elements containing information about the inputs necessary for surface
        energy calculations for a given material represented by its Materials Project ID.
        Generated by the function generate_surfen_wfs_from_inputs that takes in a
        previous version of the inputs_list generated by get_surfen_inputs_from_slab
        by inputting a Slab and a SlabGenerator corresponding to this Slab object.
    fltr : dict
        Filter to use when looking up results in the database. Generally involves the mpid.
    coll : str
        Collection to query the results from in the database.
    db_file : str, optional
        Full path to the db.json file which holds the location and access
        credentials to the database. The default is 'auto'.
    high_level : str, optional
        Whether to query the results from the high level database or not. The default is True.
    Returns
    -------
        Calculates the surface energies for the slabs given in inputs_list and writes them to
        the database.
    """

    _fw_name = "Calculate the surface energies and write to DB"
    required_params = ['inputs_list', 'fltr', 'coll']
    optional_params = ['db_file', 'high_level']

    def run_task(self, fw_spec):
        inputs_list = self.get('inputs_list')
        fltr = self.get('fltr')
        coll = self.get('coll')
        db_file = self.get('db_file', 'auto')
        high_level = self.get('high_level', True)

        write_surface_energies_to_db(inputs_list, fltr, coll, db_file, high_level)

        return FWAction(update_spec=fw_spec)


@explicit_serialize
class FT_PutSurfenInputsIntoDB(FiretaskBase):
    """
    Firetask that writes all the inputs and parameters for surface energy calculations described in
    the inputs_list into the database to be later updated with the surface energies.

    Parameters
    ----------
    inputs_list: list
        List with elements containing information about the inputs necessary for surface
        energy calculations for a given material represented by its Materials Project ID.
        Generated by the function generate_surfen_wfs_from_inputs that takes in a
        previous version of the inputs_list generated by get_surfen_inputs_from_slab
        by inputting a Slab and a SlabGenerator corresponding to this Slab object.
    sg_params : dict
        Parameters to be used in the SlabGenerator.
    fltr : dict
        Filter to use when looking up results in the database. Generally involves the mpid.
    coll : str
        Collection to query the results from in the database.
    db_file : str, optional
        Full path to the db.json file which holds the location and access
        credentials to the database. The default is 'auto'.
    high_level : str, optional
        Whether to query the results from the high level database or not. The default is True.

    Returns
    -------
    FWAction that updates the spec.

    """
    _fw_name = "Put the structures and computational parameters for the surface energy calculation into the DB"
    required_params = ['inputs_list', 'sg_params', 'comp_params', 'fltr', 'coll']
    optional_params = ['db_file', 'high_level']

    def run_task(self, fw_spec):
        inputs_list = self.get('inputs_list')
        sg_params = self.get('sg_params')
        comp_params = self.get('comp_params')
        fltr = self.get('fltr')
        coll = self.get('coll')

        db_file = self.get('db_file', 'auto')
        high_level = self.get('high_level', True)

        put_surfen_inputs_into_db(inputs_list, sg_params, comp_params, fltr, coll, db_file, high_level)

        return FWAction(update_spec=fw_spec)


@explicit_serialize
class FT_GenerateSurfenEntries(FiretaskBase):
    """ Firetask that goes through all the calculated surface energies in the database and finds the minimum
    of each orientation, creating a new entry in the same document with an easy-to-navigate
    list of surface energies.

    Parameters
    ----------
    fltr : dict
        Filter to use when looking up results in the database. Generally involves the mpid.
    coll : str
        Collection to query the results from in the database.
    db_file : str, optional
        Full path to the db.json file which holds the location and access
        credentials to the database. The default is 'auto'.
    high_level : str, optional
        Whether to query the results from the high level database or not. The default is True.
    Returns
    -------
        FWAction that updates the fw_spec.
    """
    _fw_name = 'Update the database with the lowest surface energy entries with a given input_dict'
    required_params = ['fltr', 'coll']
    optional_params = ['db_file', 'high_level']

    def run_task(self, fw_spec):
        fltr = self.get('fltr')
        coll = self.get('coll')
        db_file = self.get('db_file', 'auto')
        high_level = self.get('high_level', True)

        generate_surfen_entries(fltr, coll, db_file, high_level)

        return FWAction(update_spec=fw_spec)
