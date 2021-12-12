#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 6 22:25:05 2021

Functions that deal with the inputs for surface energy calculations.

    Author: Fırat Yalçın

"""
import numpy as np
from atomate.vasp.fireworks import StaticFW
from atomate.vasp.powerups import add_modify_incar
from fireworks import Firework, Workflow

from triboflow.firetasks.fv import FT_fake_vasp
from triboflow.phys.shaper import Shaper
from triboflow.utils.database import Navigator
from triboflow.utils.structure_manipulation import get_conv_bulk_from_mpid
from triboflow.utils.utils import dict_to_hash
from triboflow.utils.vasp_tools import get_custom_vasp_relax_settings, get_custom_vasp_static_settings
from triboflow.workflows.base import dynamic_relax_swf


def get_by_path(root, items):
    import operator
    from functools import reduce
    """Access a nested object in root by item sequence."""
    return reduce(operator.getitem, items, root)


def generate_input_dict(struct, calc_type, tag):
    """
    Simple function to generate a dictionary that stores the structure
    and the type of calculation to be performed on it.

    Parameters
    ----------
    struct : pymatgen.core.structure.Structure
        Pymatgen Structure object.
    calc_type : str
        Calculation type, either 'relax' or 'static'.
    tag : str
        Short description of the structure in relation to what
        is being calculated. For example, can be 'slab', 'ouc',
        'sto_slab' and so on.

    Returns
    -------
    input_dict : dict
        Dictionary describing the calculation.

    """
    # struct_type = 'slab' if hasattr(struct, 'miller_index') else 'bulk'
    input_dict = {
        'struct': struct,
        # 'struct_type': struct_type,
        'calc_type': calc_type,
        'calc_tag': tag
    }
    return input_dict


def get_surfen_inputs_from_slab(slab, SG=None, tol=0.1, custom_id=None):
    """
    Generates a dictionary containing all the sub structures needed to compute the
    surface energy of the given slab.

    Parameters
    ----------
    slab : pymatgen.core.surface.Slab
        Pymatgen Slab object.
    SG : pymatgen.core.surface.SlabGenerator
        Pymatgen SlabGenerator object. Necessary to generate substructures for
        non-stoichiometric slabs with a specific termination.
    tol : float
        Tolerance value used in the layering of sites used in various
        methods throughout. Value is in Angstroms.
    custom_id : str
        Unique id to identify the database entry linked to the surface energy
        calculation of the slab.

    Returns
    -------
    inputs_dict : dict
        Dictionary containing a summary of the parameters of the input slab,
        the slab itself, and the substructures needed in surface energy calculation.

    """
    # symmetry and stoichiometry of the input slab is identified in order to add
    # the correct calculations to the inputs_dict
    nn_method = 'all'
    id_slab = Shaper.identify_slab(slab)
    sym = id_slab['symmetric']
    sto = id_slab['stoichiometric']
    ## oriented unit cell is used for the reference bulk energies
    ouc = Shaper.get_constrained_ouc(slab)
    ouc_layers = len(Shaper._get_layers(ouc, tol))
    slab_layers = len(Shaper._get_layers(slab, tol))
    slab_thickness = Shaper._get_proj_height(slab, 'slab')
    ouc_input = generate_input_dict(ouc, 'relax', 'ouc')
    millerstr = ''.join([str(i) for i in slab.miller_index])
    inputs_dict = {'struct': slab,
                   'custom_id': custom_id,
                   'inputs': [ouc_input],
                   'slab_params': {'sym': sym,
                                   'sto': sto,
                                   'layer_tol': tol,
                                   'thickness_layers': slab_layers,
                                   'thickness_A': slab_thickness,
                                   'hkl': millerstr,
                                   'bvs': slab.energy,
                                   'area': slab.surface_area}}

    if sym and sto:
        # Easiest case, a simple relaxation followed by a static calculation are
        # needed for the slab and the oriented unit cell
        slab_input = generate_input_dict(slab, 'relax', 'slab_relax')
        inputs_dict['inputs'] += [slab_input]
    if sym and not sto:
        # for non-stoichiometric slabs, we need the cleavage energy for which we need
        # a stoichiometric version of the slab with terminations that are complementary
        # to each other, this is done by SlabGenerator.get_slab() which always creates
        # a stoichiometric slab.
        sto_slab = SG.get_slab(slab.shift, tol)
        sto_slab_layers = len(Shaper._get_layers(sto_slab, tol))
        # Since SlabGenerator creates a larger than than we want, we remove layers
        # and the number of layers removed is an integer multiple of the number of layers
        # in the oriented unit cell to preserve stoichiometry
        layers_to_remove = int(ouc_layers * np.floor((sto_slab_layers - slab_layers) / ouc_layers))
        sto_slab = Shaper._remove_layers(sto_slab, layers_to_remove, tol, method='layers')
        sto_slab_input = generate_input_dict(sto_slab, 'static', 'sto_slab')
        slab_relax_input = generate_input_dict(slab, 'relax', 'slab_relax')
        slab_static_input = generate_input_dict(slab, 'static', 'slab_static')
        inputs_dict['inputs'] += [slab_relax_input, slab_static_input, sto_slab_input]
    if not sym:
        # For asymmetric slabs, we need the periodicity in the layering to figure out
        # if the top and bottom terminations are complementary.
        bbs = Shaper._bonds_by_shift(SG, nn_method, tol)
        bvs, indices = np.unique(list(bbs.values()), return_index=True)
        periodicity = len(bvs)

        # Asymmetric slabs have different surface energies on the top and the bottom,
        # which means we need to relax those regions separately.
        slab_tf = Shaper.fix_regions(slab, tol, fix_type='top_half')
        slab_bf = Shaper.fix_regions(slab, tol, fix_type='bottom_half')

        slab_tf_input = generate_input_dict(slab_tf, 'relax', 'slab_top_fixed_relax')
        slab_bf_input = generate_input_dict(slab_bf, 'relax', 'slab_bot_fixed_relax')

        inputs_dict['inputs'] += [slab_tf_input, slab_bf_input]

        if not sto:
            slab_static_input = generate_input_dict(slab, 'static', 'slab_static')
            inputs_dict['inputs'] += [slab_static_input]

            if slab_layers % periodicity != 0:
                inputs_dict['slab_params'].update({'comp': True})
                # For non-stoichiometric slabs, we need the stoichiometric versions in order
                # to calculate the cleavage energy which is a component in the surface energy
                # calculation.
                all_shifts = SG._calculate_possible_shifts(tol)
                top_shift_index = all_shifts.index(slab.shift)
                bot_shift_index = (top_shift_index - slab_layers) % len(all_shifts)
                sto_slab_top = SG.get_slab(slab.shift, tol)
                sto_slab_bot = SG.get_slab(all_shifts[bot_shift_index], tol)
                sto_slab_layers = len(Shaper._get_layers(sto_slab_top, tol))
                layers_to_remove = int(ouc_layers * np.floor((sto_slab_layers - slab_layers) / ouc_layers))
                sto_slab_top = Shaper._remove_layers(sto_slab_top, layers_to_remove, tol, method='layers')
                sto_slab_bot = Shaper._remove_layers(sto_slab_bot, layers_to_remove, tol, method='layers')
                sto_slab_top_input = generate_input_dict(sto_slab_top, 'static', 'sto_slab_top')
                sto_slab_bot_input = generate_input_dict(sto_slab_bot, 'static', 'sto_slab_bot')
                inputs_dict['inputs'] += [sto_slab_top_input, sto_slab_bot_input]
            else:
                sto_slab = SG.get_slab(slab.shift, tol)
                sto_slab_layers = len(Shaper._get_layers(sto_slab, tol))
                layers_to_remove = int(ouc_layers * np.floor((sto_slab_layers - slab_layers) / ouc_layers))
                sto_slab = Shaper._remove_layers(sto_slab, layers_to_remove, tol, method='layers')
                sto_slab_input = generate_input_dict(sto_slab, 'static', 'sto_slab')
                inputs_dict['inputs'] += [sto_slab_input]
    return inputs_dict


def get_surfen_inputs_from_mpid(mpid,
                                functional,
                                sg_params,
                                sg_filter,
                                db_file='auto',
                                high_level=True,
                                comp_params_user={}):
    """
    Generates an input dictionary that contains all the necessary information about
    the slabs that are consistent with the given parameters for which surface energy
    calculations will be performed for.

    Parameters
    ----------
    mpid : str
        ID of the material in the MaterialsProject database.
    functional : str
        Which functional to use; has to be 'PBE' or 'SCAN'.
    sg_params : dict
        Dictionary containing the parameters to be used in the SlabGenerator
    db_file : str, optional
        Full path of the db.json. The default is 'auto'.
    high_level : str, optional
        Whether to query the results from the high level database or not.
        The default is True.
    max_index : int, optional
        Maximum miller up to which unique orientations will be searched.
        The default is 2.
    comp_params_user : dict
        Computational parameters for the VASP simulations. If not set, default
        parameters will be used instead. The default is {}.
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
    **kwargs : dict
        Should include a key "bvs_param" with a value that is compatible
        with the bvs_method used.

    Raises
    ------
    ValueError
        Bulk structure for the given MP-id should be in the database.

    Returns
    -------
    input_dict : dict
        Dictionary containing all the required input parameters for the
        surface energy calculations of the material described by its MP-id.

    """

    coll = f'{functional}.bulk_data'
    bulk_conv = get_conv_bulk_from_mpid(mpid, coll, db_file='auto', high_level=True)
    candidates = generate_candidate_slabs(bulk_conv, sg_params, sg_filter)
    tol = sg_params.get('tol')
    inputs_list = [get_surfen_inputs_from_slab(c[0], c[1], tol) for c in candidates]
    ## TODO: should I put the mpid info somewhere in inputs_list? What about comp_params? Does it belong here?
    return inputs_list


# def update(d, u):
#     for k, v in u.items():
#         if isinstance(v, collections.abc.Mapping):
#             d[k] = update(d.get(k, {}), v)
#         else:
#             d[k] = v
#     return d

def move_result(tag, fltr, coll, loc, custom_dict={}, db_file='auto', high_level=True):
    """
    Moves the result of a VASP calculation from the Fireworks database to the destination.

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
    None.

    """
    # function to edit DB entries with the location set with loc
    nav = Navigator(db_file)
    calc = nav.find_data('tasks', {'task_label': tag})
    # out is the whole output of the calculation including everything,
    fltr_out = ['input', 'output', 'orig_inputs', 'custodian']
    out = {f: calc[f] for f in fltr_out}
    if custom_dict:
        out.update(custom_dict)

    nav_high = Navigator(db_file='auto', high_level=high_level)

    loc = '.'.join(loc)

    nav_high.update_data(
        collection=coll,
        fltr=fltr,
        new_values={'$set': {loc: out}},
        upsert=True)


def write_surface_energies_to_db(inputs_list, fltr, coll, db_file='auto', high_level=True):
    """
    Calculates the surface energies of the slabs given in the inputs_list and writes
    the surface energies to the database under the relevant entries.

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
    None.

    """
    nav = Navigator(db_file, high_level)

    for slab_dict in inputs_list:
        surf_en = calculate_surface_energy_gen(slab_dict, fltr, coll)
        loc = slab_dict['inputs'][0]['loc'][:3]
        loc = '.'.join(loc)
        nav.update_data(
            collection=coll,
            fltr=fltr,
            new_values={'$set': {loc + '.surface_energy': surf_en}},
            upsert=True)


def generate_surfen_entries(fltr, coll, db_file='auto', high_level=True):
    """
    Goes through all the calculated surface energies in the database and finds the minimum
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
    None.

    """
    nav = Navigator(db_file, high_level)
    data = nav.find_data(coll, fltr)
    surfen_dict = {}
    for hkl, miller_data in data['miller_list'].items():
        tmp_hkl = {}
        rel_str = {}
        for uid, slab_data in miller_data.items():
            surfen = slab_data['surface_energy']
            surfen_arr = [surfen['top'], surfen['bottom']]
            tmp_hkl[uid] = surfen_arr
            rel_str[uid] = {k: v['output']['structure'] for k, v in
                            slab_data['calcs'].items() if k.endswith('relax')}
        min_uid = min(tmp_hkl, key=tmp_hkl.get)
        min_surfen = min(tmp_hkl.values())
        params = ['structure', 'comp_params', 'slab_params']
        surfen_dict[hkl] = {'top': min_surfen[0],
                            'bottom': min_surfen[1],
                            'uid': min_uid}
        params_dict = {k: v for k, v in miller_data[min_uid].items() if k in params}
        relaxed_structure = rel_str[min_uid]
        surfen_dict[hkl].update({'relaxed_structure': relaxed_structure}, **params_dict)

    nav.update_data(collection=coll,
                    fltr=fltr,
                    new_values={'$set': {'surfen_list': surfen_dict}},
                    upsert=True)


def put_surfen_inputs_into_db(inputs_list, sg_params, comp_params, fltr, coll, db_file='auto', high_level=True):
    """
    Writes all the inputs and parameters for surface energy calculations described in
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
    comp_params : dict
        Computational parameters used in the calculations. Saved here for easy troubleshooting
        possible calculation issues.
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
    None.

    """
    nav = Navigator(db_file, high_level)
    for slab in inputs_list:
        inputs = slab.get('inputs')
        slab_params = slab.get('slab_params')
        slab = slab.get('struct')
        for calc in inputs:
            struct = calc['struct']
            loc = calc['loc']
            loc_input = '.'.join(loc)
            nav.update_data(
                collection=coll,
                fltr=fltr,
                new_values={'$set': {loc_input + '.structure': struct.as_dict()}},
                upsert=True)
        loc_slab = '.'.join(loc[:3])
        nav.update_data(
            collection=coll,
            fltr=fltr,
            new_values={'$set': {loc_slab + '.structure': slab.as_dict(),
                                 loc_slab + '.slab_params': slab_params,
                                 loc_slab + '.sg_params': sg_params,
                                 loc_slab + '.comp_params': comp_params
                                 }},
            upsert=True)


def get_vis(struct, comp_params, calc_type):
    """
    Generates a VaspInputSet given the structure, computational parameters and calculation type.

    Parameters
    ----------
    struct : pymatgen.core.structure.Structure or pymatgen.core.surface.Slab
        Main structure object of Pymatgen. Can represent either bulk or surface.
    comp_params : dict
        Dictionary of computational parameters including the optimized cutoff energy
        and kpoint mesh density among other parameters.
    calc_type : str
        Type of calculation, can be 'static' or 'relax'.

    Returns
    -------
    pymatgen.io.vasp.sets.MPRelaxSet or pymatgen.io.vasp.sets.MPStaticSet
        VaspInputSet object for the given parameters.

    """
    # Simple function to generate a vasp input set from a structure, comp_params, and calc_type
    vis = get_custom_vasp_relax_settings if calc_type == 'relax' else \
        get_custom_vasp_static_settings
    struct_type = 'slab' if hasattr(struct, 'miller_index') else 'bulk'
    calc_subtype = 'pos_relax' if calc_type == 'relax' else 'from_scratch'
    return vis(struct, comp_params, f'{struct_type}_{calc_subtype}')


def get_calc_wf(struct, vis, tag, fake=False):
    """
    Generates a Workflow for the given parameters.

    Parameters
    ----------
    struct : pymatgen.core.structure.Structure or pymatgen.core.surface.Slab
        Main structure object of Pymatgen. Can represent either bulk or surface.
    vis : pymatgen.io.vasp.sets.MPRelaxSet or pymatgen.io.vasp.sets.MPStaticSet
        VaspInputSet object for the given structure.
    tag : str
        Unique ID for the calculation to later query for in the tasks.

    Returns
    -------
    WF : fireworks.core.firework.Workflow
        Workflow object for a VASP calculation.

    """
    # Simple function to make a Workflow from a structure, vis, and a tag
    if fake:
        WF = Workflow.from_Firework(Firework(FT_fake_vasp(tag=tag)))
        return WF
    if not hasattr(vis, 'prev_kpoints'):
        WF = dynamic_relax_swf([[struct, vis, tag]], add_static=True)
    else:
        FW = StaticFW(structure=struct, name=tag, vasp_input_set=vis)
        WF = Workflow.from_Firework(FW, name=tag)
        WF = add_modify_incar(WF)
    return WF


def generate_candidate_slabs(bulk_conv,
                             sg_params,
                             sg_filter):
    """
    Generates slabs within certain constraints from a given bulk structure.

    Parameters
    ----------
    bulk_conv : pymatgen.core.structure.Structure
        Conventional standard bulk structure used in the SlabGenerator.
    sg_params : dict
        Dict of parameters used in the SlabGenerator. For info about required
        and optional keys, refer to Shaper.generate_slabs()
    sg_filter : dict
        Dict of parameters used to filter the generated slabs.
        Required keys are:
            method
            <method_prefix>_param
        where currently the only implemented methods are 'bvs_threshold' and
        'bvs_min_N'. 'bvs_threshold' filters out slabs with energy attributes
        within a threshold of the minimum value of the energy across all slabs,
        while 'bvs_min_N' filters out N lowest energy slabs. Energy here refers
        to the bond valence sums of broken bonds when the slabs are generated,
        and can be accessed by slab.energy.

    Returns
    -------
    list
        List of all the (Slab, SlabGenerator) tuples that satisfy the constraints.

    """
    reconstruct = sg_params.get('reconstruct', True)
    slabs_list, SG_dict = Shaper.generate_slabs(bulk_conv, sg_params, reconstruct)
    bvs = [slab.energy for slab in slabs_list]

    # lll = sg_params.get('lll_reduce')
    # prim = sg_params.get('prim')
    # sym = sg_params.get('symmetrize')
    # mns = sg_params.get('mns')
    # formula = slabs_list[0].composition.reduced_formula
    # for index, slab in enumerate(slabs_list):
    #     miller = ''.join([str(m) for m in slab.miller_index])
    #     slab.to('poscar', f'~/{formula}_{miller}_{index}_lll_{lll}_prim_{prim}_sym_{sym}_mns_{mns}.vasp')

    method = sg_filter.get('method')
    if method == 'bvs_threshold':
        bvs_tol = sg_filter.get('bvs_param')
        min_bvs = min(bvs)
        filtered_slabs = [slab for slab in slabs_list if slab.energy / min_bvs - 1 < bvs_tol]
    elif method == 'bvs_min_N':
        N = sg_filter.get('bvs_param')
        if N < len(slabs_list):
            sorted_ind = np.argsort(bvs)
            filtered_slabs = [slabs_list[i] for i in sorted_ind[:N]]
        else:
            filtered_slabs = slabs_list
    else:
        filtered_slabs = slabs_list

    return [(slab, SG_dict[slab.miller_index]) for slab in filtered_slabs]


def set_by_path(root, items, value):
    """Set a value in a nested object in root by item sequence."""
    get_by_path(root, items[:-1])[items[-1]] = value


def get_entry_by_loc(nav, fltr, coll, loc):
    """
    Locates and returns the entry at the given location.

    Parameters
    ----------
    nav : triboflow.utils.database.Navigator
        Navigator object for the database being queried.
    fltr : dict
        Filter dictionary that is used to query for the destination in the database.
    coll : str
        MongoDB collection in which to search for.
    loc : list
        Location in which to search the entry at.

    Returns
    -------
    entry : any
        The database entry queried or None.

    """
    data = nav.find_data(coll, fltr)
    if data:
        try:
            entry = get_by_path(data, loc)
        except KeyError:
            # print(f'result not found at {loc} in {coll} for {fltr}')
            return None
        if entry:
            # print(f'result found at {loc} in {coll} for {fltr}')
            return entry


def calculate_surface_energy_gen(slab_dict, fltr, coll, db_file='auto', high_level=True):
    """
    Calculates the surface energies of top and bottom terminations of the slab given as a dictionary.

    Parameters
    ----------
    slab_dict : dict
        Dictionary containing the information about the slab. Generated by the function
        get_surfen_inputs_from_slab by inputting a Slab and a SlabGenerator corresponding
        to this Slab object.
    fltr : dict
        Filter to use when looking up results in the database. Generally involves the mpid.
    coll : str
        Collection to query the results from in the database.
    db_file : str, optional
        Full path of the db.json. The default is 'auto'.
    high_level : str, optional
        Whether to query the results from the high level database or not. The default is True.

    Returns
    -------
    dict
        Surface energies of the top and bottom terminations of the Slab described by the slab_dict.

    """
    nav = Navigator(db_file, high_level)
    slab_params = slab_dict.get('slab_params')

    sym = slab_params.get('sym')
    sto = slab_params.get('sto')
    comp = slab_params.get('comp')
    area = slab_params.get('area')
    en_dict = {}

    inputs = slab_dict.get('inputs')
    for entry in inputs:
        loc = entry['loc']
        calc_tag = entry['calc_tag']
        result = get_entry_by_loc(nav, fltr, coll, loc)
        output = result['output']
        en_dict[calc_tag] = {'nsites': len(output['structure']['sites']),
                             'energy': output['energy'],
                             'energy_per_atom': output['energy_per_atom']}

    bulk_en = en_dict['ouc']['energy_per_atom']
    if sym and sto:
        slab_en = en_dict['slab_relax']['energy']
        nsites = en_dict['slab_relax']['nsites']
        surf_en_top = slab_en - nsites * bulk_en
        surf_en_bot = surf_en_top

    if sym and not sto:
        slab_relax_en = en_dict['slab_relax']['energy']
        slab_static_en = en_dict['slab_static']['energy']
        sto_slab_en = en_dict['sto_slab']['energy']

        sto_slab_nsites = en_dict['sto_slab']['nsites']

        E_cle = sto_slab_en - sto_slab_nsites * bulk_en
        E_rel = slab_relax_en - slab_static_en

        surf_en_top = E_cle + E_rel
        surf_en_bot = surf_en_top

    if not sym:
        slab_tf_relax_en = en_dict['slab_top_fixed_relax']['energy']
        slab_bf_relax_en = en_dict['slab_bot_fixed_relax']['energy']

        nsites = en_dict['slab_bot_fixed_relax']['nsites']

        if sto:
            surf_en_top = slab_bf_relax_en - nsites * bulk_en
            surf_en_bot = slab_tf_relax_en - nsites * bulk_en
        else:
            slab_static_en = en_dict['slab_static']['energy']
            E_rel_top = slab_bf_relax_en - slab_static_en
            E_rel_bot = slab_tf_relax_en - slab_static_en

            if comp:
                sto_slab_en = en_dict['sto_slab']['energy']
                sto_slab_nsites = en_dict['sto_slab']['nsites']

                E_cle = sto_slab_en - sto_slab_nsites * bulk_en

                surf_en_top = E_cle + E_rel_top
                surf_en_bot = E_cle + E_rel_bot
            else:
                sto_slab_top_en = en_dict['sto_slab_top']['energy']
                sto_slab_bot_en = en_dict['sto_slab_bot']['energy']

                sto_slab_top_nsites = en_dict['sto_slab_top']['nsites']
                sto_slab_bot_nsites = en_dict['sto_slab_bot']['nsites']

                E_cle_top = sto_slab_top_en - sto_slab_top_nsites * bulk_en
                E_cle_bot = sto_slab_bot_en - sto_slab_bot_nsites * bulk_en

                surf_en_top = E_cle_top + E_rel_top
                surf_en_bot = E_cle_bot + E_rel_bot

    surf_en_top *= 16.02176565 / (2 * area)
    surf_en_bot *= 16.02176565 / (2 * area)

    return {'top': surf_en_top, 'bottom': surf_en_bot}


def generate_surfen_wfs_from_inputs(inputs_list, comp_params, fake=False):
    """
    Generates Workflows for the calculations needed for surface energy and
    updates the input dictionary with them.

    Parameters
    ----------
    inputs_list : list
        List containing information about the inputs necessary for surface
        energy calculations for a given material represented by its Materials Project ID.
        Generated by the function get_surfen_inputs_from_slab by inputting a Slab and
        a SlabGenerator corresponding to this Slab object.
    comp_params : dict
        Dictionary containing the computational parameters to be used in the VASP
        calculations.
    Returns
    -------
    inputs_list with each element updated with the Workflows, unique IDs, locations, and tags.

    """
    WFs = []
    for slab in inputs_list:
        frac_coords = slab.get('struct').frac_coords
        inputs = slab.get('inputs')
        slab_params = slab.get('slab_params')
        all_params = {**comp_params, **slab_params, 'coords': frac_coords}
        uid = slab.get('custom_id')
        if not uid:
            uid = dict_to_hash(all_params)
        slab_params.update({'uid': uid})
        for calc in inputs:
            struct = calc.get('struct')
            calc_type = calc.get('calc_type')
            calc_tag = calc.get('calc_tag')
            millerstr = slab_params.get('hkl')
            loc = ['miller_list', millerstr, uid, 'calcs', calc_tag]
            formula = struct.composition.reduced_formula
            tag = f'{formula}_{millerstr}_{calc_tag}_{uid}'
            vis = get_vis(struct, comp_params, calc_type)
            WF = get_calc_wf(struct, vis, tag, fake)
            WFs.append(WF)
            tag_dict = {'WF': WF, 'tag': tag, 'loc': loc}
            calc.update(tag_dict)

    return inputs_list
