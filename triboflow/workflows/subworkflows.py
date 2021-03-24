"""SubWorkflows to be used by the TriboFlow Project.

Created on Wed Jun 17 15:47:39 2020
@author: mwo
"""

from uuid import uuid4
import numpy as np

from fireworks import Workflow, Firework

from atomate.vasp.fireworks import StaticFW

from triboflow.fireworks.common import run_pes_calc_fw, make_pes_fw
from triboflow.firetasks.encut_convergence import FT_EnergyCutoffConvo
from triboflow.firetasks.kpoint_convergence import FT_KpointsConvo
from triboflow.firetasks.structure_manipulation import FT_MakeSlabInDB, \
    FT_StartSlabRelax, FT_GetRelaxedSlab
from triboflow.firetasks.PPES import FT_DoPPESCalcs, FT_FitPPES
from triboflow.firetasks.adhesion import FT_CalcAdhesion
from triboflow.utils.database import Navigator, NavigatorMP
from triboflow.utils.vasp_tools import get_emin, get_custom_vasp_static_settings


def adhesion_energy_swf(top_slab,
                        bottom_slab,
                        interface,
                        interface_name=None,
                        functional='PBE',
                        comp_parameters={}):
    """Create a subworkflow to compute the adhesion energy for an interface.
    
    This workflow takes two matched slabs (their cells must be identical) and
    a relaxed interface structure of those slabs and computes the andhesion
    energy. The two matched slabs must be relaxed as well to get correct
    results.
    Output are saved in a high-level database, but may also be also written
    as files and copied to a chosen location. Note that this copy operation
    is generally dependent on which machine the calculations are executed,
    and not on the machine where the workflow is submitted. Also ssh-keys need
    to be set up for remote_copy to work!
    
    Parameters
    ----------
    top_slab : pymatgen.core.surface.Slab
        Relaxed top slab of the interface.
    bottom_slab : pymatgen.core.surface.Slab
        Relaxed bottom slab of the interface.
    interface : pymatgen.core.surface.Slab
        Relaxed interface structure.
    interface_name : str, optional
        Unique name to find the interface in the database with.
        The default is None, which will lead to an automatic interface_name
        generation which will be printed on screen.
    bottom_mpid : str, optional
        ID of the bulk material of the top slab in the MP database.
        The default is None.
    functional : str, optional
        Which functional to use; has to be 'PBE' or 'SCAN'. The default is 'PBE'
    comp_parameters : dict, optional
        Computational parameters to be passed to the vasp input file generation.
        The default is {}.

    Returns
    -------
    SWF : fireworks.core.firework.Workflow
        A subworkflow intended to compute the adhesion of a certain interface.

    """
    try:
        top_miller = list(top_slab.miller_index)
    except:
        raise AssertionError("You have used {} as an input for <top_slab>.\n"
                             "Please use <class 'pymatgen.core.surface.Slab'>"
                             " instead.".format(type(top_slab)))
        
    try:
        bot_miller = list(bottom_slab.miller_index)
    except:
        raise AssertionError("You have used {} as an input for <bot_slab>.\n"
                             "Please use <class 'pymatgen.core.surface.Slab'>"
                             " instead.".format(type(bottom_slab)))
        
    if not interface_name:
        mt = ''.join(str(s) for s in top_miller)
        mb = ''.join(str(s) for s in bot_miller)
        interface_name = (top_slab.composition.reduced_formula+'_'+mt+'_'+
                          bottom_slab.composition.reduced_formula+'_'+mb+
                          '_AutoGen')
        print('\nYour interface name has been automatically generated to be:'
              '\n {}'.format(interface_name))
        
    
    if comp_parameters == {}:
        print('\nNo computational parameters have been defined!\n'
              'Workflow will run with:\n'
              '   ISPIN = 1\n'
              '   ISMEAR = 0\n'
              '   ENCUT = 520\n'
              '   kpoint density kappa = 5000\n'
              'We recommend to pass a comp_parameters dictionary'
              ' of the form:\n'
              '   {"use_vdw": <True/False>,\n'
              '    "use_spin": <True/False>,\n'
              '    "is_metal": <True/False>,\n'
              '    "encut": <float>,\n'
              '    "k_dens": <int>}\n')
    
    tag = interface_name+'_'+str(uuid4())
    
    vis_top = get_custom_vasp_static_settings(top_slab,
                                              comp_parameters,
                                              'slab_from_scratch')
    vis_bot = get_custom_vasp_static_settings(bottom_slab,
                                              comp_parameters,
                                              'slab_from_scratch')
    vis_interface = get_custom_vasp_static_settings(interface,
                                              comp_parameters,
                                              'slab_from_scratch')
    
    FW_top = StaticFW(structure=top_slab, vasp_input_set=vis_top,
                      name=tag+'top')
    FW_bot = StaticFW(structure=bottom_slab, vasp_input_set=vis_bot,
                      name=tag+'bottom')
    FW_interface = StaticFW(structure=interface, vasp_input_set=vis_interface,
                      name=tag+'interface')
    
    FW_results = Firework(FT_CalcAdhesion(interface_name=interface_name,
                                 functional=functional,
                                 top_label=tag+'top',
                                 bottom_label=tag+'bottom',
                                 interface_label=tag+'interface'))
    WF = Workflow(fireworks=[FW_top, FW_bot, FW_interface, FW_results],
                  links_dict={FW_top: [FW_results],
                              FW_bot: [FW_results],
                              FW_interface: [FW_results]},
                  name='Calculate adhesion SWF for {}'.format(interface_name))
    
    return WF
    
    
def calc_pes_swf(top_slab, bottom_slab,
                 interface_name=None,
                 functional='PBE',
                 comp_parameters={},
                 file_output=False,
                 output_dir=None,
                 remote_copy=False,
                 server=None, 
                 user=None, 
                 port=None):
    """Create a subworkflow to compute the PES for an interface of two slabs.
    
    This workflow takes two matched slabs (their cells must be identical) as
    input and computes the potential energy surface (PES) for the interface.
    Output are saved in a high-level database, but may also be also written
    as files and copied to a chosen location. Note that this copy operation
    is generally dependent on which machine the calculations are executed,
    and not on the machine where the workflow is submitted. Also ssh-keys need
    to be set up for remote_copy to work!
    
    Parameters
    ----------
    top_slab : pymatgen.core.surface.Slab
        Top slab of the interface.
    bottom_slab : pymatgen.core.surface.Slab
        Bottom slab of the interface.
    interface_name : str, optional
        Unique name to find the interface in the database with.
        The default is None, which will lead to an automatic interface_name
        generation which will be printed on screen.
    bottom_mpid : str, optional
        ID of the bulk material of the top slab in the MP database.
        The default is None.
    functional : str, optional
        Which functional to use; has to be 'PBE' or 'SCAN'. The default is 'PBE'
    comp_parameters : dict, optional
        Computational parameters to be passed to the vasp input file generation.
        The default is {}.
    file_output : bool, optional
        Toggles file output. The default is False.
    output_dir : str, optional
        Defines a directory the output is to be copied to. (Do not use a
        trailing / and/or relative location symbols like ~/.)
        The default is None.
    remote_copy : bool, optional
        If true, scp will be used to copy the results to a remote server. Be
        advised that ssh-key certification must be set up between the two
        machines. The default is False.
    server : str, optional
        Fully qualified domain name of the server the output should be copied
        to. The default is None.
    user : str, optional
        The user name on the remote server.
    port : int, optional
        On some machines ssh-key certification is only supported for certain
        ports. A port may be selected here. The default is None.

    Returns
    -------
    SWF : fireworks.core.firework.Workflow
        A subworkflow intended to compute the PES of a certain interface.

    """  
    try:
        top_miller = list(top_slab.miller_index)
    except:
        raise AssertionError("You have used {} as an input for <top_slab>.\n"
                             "Please use <class 'pymatgen.core.surface.Slab'>"
                             " instead.".format(type(top_slab)))
        
    try:
        bot_miller = list(bottom_slab.miller_index)
    except:
        raise AssertionError("You have used {} as an input for <bot_slab>.\n"
                             "Please use <class 'pymatgen.core.surface.Slab'>"
                             " instead.".format(type(bottom_slab)))
        
    if not interface_name:
        mt = ''.join(str(s) for s in top_miller)
        mb = ''.join(str(s) for s in bot_miller)
        interface_name = (top_slab.composition.reduced_formula+'_'+mt+'_'+
                          bottom_slab.composition.reduced_formula+'_'+mb+
                          '_AutoGen')
        print('\nYour interface name has been automatically generated to be:'
              '\n {}'.format(interface_name))
        
    
    if comp_parameters == {}:
        print('\nNo computational parameters have been defined!\n'
              'Workflow will run with:\n'
              '   ISPIN = 1\n'
              '   ISMEAR = 0\n'
              '   ENCUT = 520\n'
              '   kpoint density kappa = 5000\n'
              'We recommend to pass a comp_parameters dictionary'
              ' of the form:\n'
              '   {"use_vdw": <True/False>,\n'
              '    "use_spin": <True/False>,\n'
              '    "is_metal": <True/False>,\n'
              '    "encut": <float>,\n'
              '    "k_dens": <int>}\n')
    
    tag = interface_name+'_'+str(uuid4())
    
    FW_1 = run_pes_calc_fw(top_slab=top_slab, 
                           bottom_slab=bottom_slab, 
                           interface_name=interface_name,
                           functional=functional,
                           comp_parameters=comp_parameters,
                           tag=tag,
                           FW_name='Start PES calcs for '+interface_name)
    
    FW_2 = make_pes_fw(interface_name=interface_name,
                       functional=functional,
                       tag=tag,
                       file_output=file_output,
                       output_dir=output_dir,
                       remote_copy=remote_copy,
                       server=server, 
                       user=user, 
                       port=port,
                       FW_name='Parse PES calcs for '+interface_name)
    
    SWF = Workflow([FW_1, FW_2], {FW_1: [FW_2]},
                   name='Calc PES for '+interface_name+' SWF')
    return SWF
    
def calc_ppes_swf(interface_name, functional, distance_list=[-0.5, -0.25, 0.0, 
                  0.25, 0.5, 2.5, 3.0, 4.0, 5.0, 7.5],
                  out_name='PPES@minimum', structure_name='minimum_relaxed',
                  spec={}):
    """
    Generate a subworkflow that calculates a PPES using static calculations.
    
    For a given interface in the high-level database this subworkflow performs
    static calculations for different distances of the two slabs modeling
    brittle cleavage under mode 1 loading using the rigid separation model.
    The results are saved in a energy vs distance array and saved i the high-
    level database alongside a fit to a UBER curve.

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
    tag = interface_name+'_'+str(uuid4())
    
    FW_1 = Firework(FT_DoPPESCalcs(interface_name=interface_name,
                                   functional=functional,
                                   distance_list = distance_list,
                                   tag=tag,
                                   structure_name = structure_name),
                    spec=spec, name='PPES Calculations for '+interface_name)
    
    FW_2 = Firework(FT_FitPPES(interface_name=interface_name,
                                   functional=functional,
                                   distance_list = distance_list,
                                   out_name = out_name,
                                   tag = tag),
                    spec=spec, name='PPES Fitting for '+interface_name)
    
    SWF = Workflow([FW_1, FW_2], {FW_1: [FW_2]},
                   name = 'Calc PPES for '+interface_name+' SWF')
    return SWF

def make_and_relax_slab_swf(bulk_structure,
                            miller_index,
                            flag,
                            comp_parameters={},
                            functional='PBE',
                            min_thickness=10.0,
                            min_vacuum=15.0,
                            relax_type='slab_pos_relax',
                            slab_struct_name='unrelaxed_slab',
                            out_struct_name='relaxed_slab',
                            spec={},
                            file_output=False,
                            output_dir=None,
                            remote_copy=False,
                            server=None, 
                            user=None, 
                            port=None,
                            print_help=True):
    """
    Make and relax a slab.

    Parameters
    ----------
    bulk_structure : pymatgen.core.structure.Structure
        Bulk structure that is used to construct the slab out of.
    miller_index : list of int or str
        Miller indices of the slab to make.
    flag : str
        An identifier to find the results in the database. It is strongly
        suggested to use the proper Materials-ID from the MaterialsProject
        if it is known for the specific input structure. Otherwise use something
        unique which you can find again.
    comp_parameters : dict, optional
        Computational parameters to be passed to the vasp input file generation.
        The default is {}.
    functional : str, optional
        Which functional to use; has to be 'PBE' or 'SCAN'. The default is 'PBE'.
    min_thickness : float, optional
        Minimal thickness of the unrelaxed slab in Angstrom. The default is 10.0.
    min_vacuum : float, optional
        Minimal thickness of the vacuum layer in Angstrom. The default is 25.0.
    relax_type : str, optional
        Which type of relaxation to run. See get_custom_vasp_relax_settings from
        triboflow.utils.vasp_tools. The default is 'slab_pos_relax'.
    slab_struct_name : str, optional
        Name of the unrelaxed slab in the high-level database.
        The default is 'unrelaxed_slab'.
    out_struct_name : TYPE, optional
        DESCRIPTION. The default is 'relaxed_slab'.
    spec : dict, optional
        DESCRIPTION. The default is {}.
    file_output : bool, optional
        Toggles file output. The default is False.
    output_dir : str, optional
        Defines a directory the output is to be copied to. (Do not use a
        trailing / and/or relative location symbols like ~/.)
        The default is None.
    remote_copy : bool, optional
        If true, scp will be used to copy the results to a remote server. Be
        advised that ssh-key certification must be set up between the two
        machines. The default is False.
    server : str, optional
        Fully qualified domain name of the server the output should be copied
        to. The default is None.
    user : str, optional
        The user name on the remote server.
    port : int, optional
        On some machines ssh-key certification is only supported for certain
        ports. A port may be selected here. The default is None.

    Returns
    -------
    SWF : fireworks.core.firework.Workflow
        A subworkflow intended to make and relax a slab from a bulk structure.

    """
    if type(miller_index) == str:
        miller_str = miller_index
        miller = [int(k) for k in list(miller_index)]
    else:
        miller = miller_index
        miller_str = ''.join(str(s) for s in miller_index)
        
    formula = bulk_structure.composition.reduced_formula
    
    if flag.startswith('mp-') and flag[3:].isdigit():
        nav_mp = NavigatorMP()
        formula_from_flag = nav_mp.get_property_from_mp(
            mp_id=flag,
            properties=['pretty_formula'])
        formula_from_flag = formula_from_flag['pretty_formula']
        if not formula_from_flag == formula:
            raise SystemExit('The chemical formula of your structure ({}) '
                             'does not match the chemical formula of the flag '
                             '(mp-id) you have chosen which corresponds '
                             'to {}.\n'.format(
                                 formula, formula_from_flag))
    
    if comp_parameters == {}:
        print('\nNo computational parameters have been defined!\n'
              'Workflow will run with:\n'
              '   ISPIN = 1\n'
              '   ISMEAR = 0\n'
              '   ENCUT = 520\n'
              '   kpoint density kappa = 5000\n'
              'We recommend to pass a comp_parameters dictionary'
              ' of the form:\n'
              '   {"use_vdw": <True/False>,\n'
              '    "use_spin": <True/False>,\n'
              '    "is_metal": <True/False>,\n'
              '    "encut": <int>,\n'
              '    "k_dens": <int>}\n')
    
    if print_help:
        nav = Navigator()
        db_file = nav.path
        print('Once you workflow has finished you can access the '
              'results from the database using this code:\n\n'
              'import pprint\n'
              'from triboflow.utils.database import GetSlabFromDB\n'
              'results = GetBulkFromDB("{}", "{}", "{}", "{}")\n'
              'pprint.pprint(results)\n'.format(flag, db_file, miller, functional))
    
    tag = formula+miller_str+'_'+str(uuid4())
            
    FTs = []
    
    FTs.append(FT_MakeSlabInDB(bulk_structure=bulk_structure,
                               miller=miller,
                               flag=flag,
                               functional=functional,
                               min_thickness=min_thickness,
                               min_vacuum=min_vacuum))
    
    FTs.append(FT_StartSlabRelax(flag=flag, miller=miller,
                                 functional=functional, tag=tag,
                                 comp_parameters=comp_parameters,
                                 slab_struct_name=slab_struct_name,
                                 relax_type=relax_type))
    
    FW = Firework(FTs, spec=spec,
                  name='Make and relax '+formula+miller_str+' slab')
    
    FW2 = Firework(FT_GetRelaxedSlab(flag=flag,
                                     miller=miller,
                                     functional=functional,
                                     tag=tag,
                                     struct_out_name=out_struct_name,
                                     file_output=file_output,
                                     output_dir=output_dir,
                                     remote_copy=remote_copy,
                                     server=server, 
                                     user=user, 
                                     port=port),
                   spec=spec,
                   name='Put relaxed '+formula+miller_str+' slab in DB')
    
    SWF = Workflow([FW, FW2], {FW: [FW2]},
                   name='Make and relax '+formula+miller_str+' SWF')
    return SWF


def converge_kpoints_swf(structure,
                         flag,
                         comp_parameters={},
                         spec={},
                         functional='PBE',
                         k_dens_start=500,
                         k_dens_incr=50,
                         n_converge=3,
                         db_file=None,
                         file_output=False,
                         output_dir=None,
                         remote_copy=False,
                         server=None, 
                         user=None, 
                         port=None,
                         print_help=True):
    """Subworkflows that converges the the k-mesh density via total energy.
    
    Takes a given structure, computational parameters (which includes the
    convergence criterion in eV/atom) and runs static vasp calculations with
    a denser and denser mesh (larger k_dens parameter) until convergence in the
    total energy is achieved. Output is printed to the screen and saved in the
    high-level triboflow database where it can be queried using the flag set.

    Parameters
    ----------
    structure : pymatgen.core.structure.Structure
        The structure for which to converge the K-pint grids.
    flag : str
        An identifier to find the results in the database. It is strongly
        suggested to use the proper Materials-ID from the MaterialsProject
        if it is known for the specific input structure. Otherwise use something
        unique which you can find again.
    comp_parameters : dict
        Dictionary of computational parameters for the VASP calculations.
        Convergence criterion in eV/atom can be given here under the key:
        'energy_tolerance' and defaults to 0.001 (1meV/atom).
    spec : dict
        Previous fw_spec that will be updated and/or passed on for child
        Fireworks.
    k_dens_start : int, optional
        Starting density value for the first run. Defaults to 500.
    k_dens_incr : int, optional
        Increment for the k-mesh density during the convergence. Defaults to
        50. The increment might actually be larger if otherwise no new mesh
        would be formed!
    n_converge : int, optional
        Number of calculations that have to be inside the convergence
        threshold for convergence to be reached. Defaults to 3.
    db_file : str
        Full path to the db.json file that should be used. Defaults to
        None, in which case env_chk will be used in the FT.
    file_output : bool, optional
        Toggles file output. The default is False.
    output_dir : str, optional
        Defines a directory the output is to be copied to. (Do not use a
        trailing / and/or relative location symbols like ~/.)
        The default is None.
    remote_copy : bool, optional
        If true, scp will be used to copy the results to a remote server. Be
        advised that ssh-key certification must be set up between the two
        machines. The default is False.
    server : str, optional
        Fully qualified domain name of the server the output should be copied
        to. The default is None.
    user : str, optional
        The user name on the remote server.
    port : int, optional
        On some machines ssh-key certification is only supported for certain
        ports. A port may be selected here. The default is None.
    print_help : bool, optional
        Prints a few lines of code that shows how to retrieve the results from
        the database. The default is True.

    Returns
    -------
    WF : fireworks.core.firework.Workflow
        A subworkflow intended to find the converged k_distance for a given
        structure.

    """
    formula = structure.composition.reduced_formula
    name = 'Kpoint Convergence SWF of '+formula
    
    if flag.startswith('mp-') and flag[3:].isdigit():
        formula_from_struct = structure.composition.reduced_formula
        nav_mp = NavigatorMP()
        formula_from_flag = nav_mp.get_property_from_mp(
            mp_id=flag,
            properties=['pretty_formula'])
        formula_from_flag = formula_from_flag['pretty_formula']

        if not formula_from_flag == formula_from_struct:
            raise SystemExit('The chemical formula of your structure ({}) '
                             'does not match the chemical formula of the flag '
                             '(mp-id) you have chosen which corresponds '
                             'to {}.\n'.format(
                                 formula_from_struct, formula_from_flag))
    
    if comp_parameters == {}:
        print('\nNo computational parameters have been defined!\n'
              'Workflow will run with:\n'
              '   ISPIN = 1\n'
              '   ISMEAR = 0\n'
              '   ENCUT = 520\n'
              'We recommend to pass a comp_parameters dictionary'
              ' of the form:\n'
              '   {"use_vdw": <True/False>,\n'
              '    "use_spin": <True/False>,\n'
              '    "is_metal": <True/False>,\n'
              '    "encut": <int>}\n')
    
    if print_help:
        nav = Navigator()
        db_file = nav.path
        print('Once you workflow has finished you can access the '
              'results from the database using this code:\n\n'
              'import pprint\n'
              'from triboflow.utils.database import GetBulkFromDB\n'
              'results = GetBulkFromDB("{}", "{}", "{}")\n'
              'pprint.pprint(results)\n'.format(flag, db_file, functional))
    
    tag = "Kpoints group for {} - {}".format(formula, str(uuid4()))
        
    FT_KptsConvo = FT_KpointsConvo(structure=structure,
                                   comp_params=comp_parameters,
                                   tag=tag,
                                   flag=flag,
                                   functional=functional,
                                   db_file=db_file,
                                   k_dens_incr=k_dens_incr,
                                   k_dens_start=k_dens_start,
                                   n_converge=n_converge,
                                   file_output=file_output,
                                   output_dir=output_dir,
                                   remote_copy=remote_copy,
                                   server=server,
                                   user=user,
                                   port=port)
    
    FW_CE = Firework(FT_KptsConvo, spec=spec,
                     name='Kpoint Convergence')
    WF = Workflow([FW_CE], name=name)

    return WF

def converge_encut_swf(structure, 
                       flag, 
                       comp_parameters={}, 
                       spec={},
                       functional='PBE', 
                       deformations=None, 
                       encut_start=None,
                       encut_incr=25, 
                       n_converge=3, 
                       db_file=None,
                       file_output=False,
                       output_dir=None,
                       remote_copy=False,
                       server=None, 
                       user=None, 
                       port=None,
                       print_help=True):
    """Subworkflows that converges the Encut using a fit to an BM-EOS.
    
    Takes a given structure, computational parameters, and a optional list
    of deformations and uses these deformations to compute an
    Birch-Murnaghan equation of state for higher and higher energy cutoffs.
    Once bulk modulus and equilibrium volume do not change any longer,
    convergence is reached. Output is printed to the screen and saved in the
    high-level triboflow database where it can be queried using the mp_id
    of the material.

    Parameters
    ----------
    structure : pymatgen.core.structure.Structure
        The structure for which to converge the energy cutoff parameter.
    flag : str
        An identifier to find the results in the database. It is strongly
        suggested to use the proper Materials-ID from the MaterialsProject
        if it is known for the specific input structure. Otherwise use something
        unique which you can find again.
    comp_parameters : dict, optional
        Dictionary of computational parameters for the VASP calculations. The
        default is {}.
    spec : dict, optional
        Previous fw_spec that will be updated and/or passed on for child
        Fireworks. The default is {}.
    deformations: list of lists, optional
        List of deformation matrices for the fit to the EOS. Defaults to None,
        which results in 5 volumes from 90% to 110% of the initial volume.
    encut_start : float, optional
        Starting encut value for the first run. Defaults to the largest EMIN
        in the POTCAR.
    encut_incr : float, optional
        Increment for the encut during the convergence. Defaults to 25.
    n_converge : int, optional
        Number of calculations that have to be inside the convergence
        threshold for convergence to be reached. Defaults to 3.
    db_file : str, optional
        Full path to the db.json file that should be used. Defaults to
        None, in which case env_chk will be used in the FT.
    file_output : bool, optional
        Toggles file output. The default is False.
    output_dir : str, optional
        Defines a directory the output is to be copied to. (Do not use a
        trailing / and/or relative location symbols like ~/.)
        The default is None.
    remote_copy : bool, optional
        If true, scp will be used to copy the results to a remote server. Be
        advised that ssh-key certification must be set up between the two
        machines. The default is False.
    server : str, optional
        Fully qualified domain name of the server the output should be copied
        to. The default is None.
    user : str, optional
        The user name on the remote server.
    port : int, optional
        On some machines ssh-key certification is only supported for certain
        ports. A port may be selected here. The default is None.
    print_help : bool, optional
        Prints a few lines of code that shows how to retrieve the results from
        the database. The default is True.

    Returns
    -------
    WF : fireworks.core.firework.Workflow
        A subworkflow intended to find the converged ENCUT for a given
        structure.

    """   
    name = 'Encut Convergence SWF of '+structure.composition.reduced_formula
    
    tag = "BM group: {}".format(str(uuid4()))
        
    if flag.startswith('mp-') and flag[3:].isdigit():
        formula_from_struct = structure.composition.reduced_formula
        nav_mp = NavigatorMP()
        formula_from_flag = nav_mp.get_property_from_mp(
            mp_id=flag, 
            properties=['pretty_formula'])
        formula_from_flag = formula_from_flag['pretty_formula']

        if not formula_from_flag == formula_from_struct:
            raise SystemExit('The chemical formula of your structure ({}) '
                             'does not match the chemical formula of the flag '
                             '(mp-id) you have chosen which corresponds '
                             'to {}.\n'.format(
                                 formula_from_struct, formula_from_flag))
    
    if not encut_start:
        # Get the largest EMIN value of the potcar and round up to the
        # next whole 25.
        vis = get_custom_vasp_static_settings(structure, comp_parameters,
                                              'bulk_from_scratch')
        emin = get_emin(vis.potcar)
        encut_start = int(25 * np.ceil(emin/25))
        
    if comp_parameters == {}:
        print('\nNo computational parameters have been defined!\n'
              'Workflow will run with:\n'
              '   ISPIN = 1\n'
              '   ISMEAR = 0\n'
              '   kpoint density kappa = 5000\n'
              'We recommend to pass a comp_parameters dictionary'
              ' of the form:\n'
              '   {"use_vdw": <True/False>,\n'
              '    "use_spin": <True/False>,\n'
              '    "is_metal": <True/False>,\n'
              '    "k_dens": <int>}\n')
    
    if print_help:
        nav = Navigator()
        db_file = nav.path
        print('Once you workflow has finished you can access the '
              'results from the database using this code:\n\n'
              'import pprint\n'
              'from triboflow.utils.database import GetBulkFromDB\n'
              'results = GetBulkFromDB("{}", "{}", "{}")\n'
              'pprint.pprint(results)\n'.format(flag, db_file, functional))
    
    FT_EncutConvo = FT_EnergyCutoffConvo(structure=structure,
                                         comp_params=comp_parameters,
                                         tag=tag,
                                         flag=flag,
                                         functional=functional,
                                         deformations=deformations,
                                         db_file=db_file,
                                         encut_incr=encut_incr,
                                         encut_start=encut_start, 
                                         n_converge=n_converge,
                                         file_output=file_output,
                                         output_dir=output_dir,
                                         remote_copy=remote_copy,
                                         server=server, 
                                         user=user, 
                                         port=port)        
    
    FW_CE = Firework(FT_EncutConvo, spec=spec,
                     name='Encut Convergence')
    WF = Workflow([FW_CE], name=name)

    return WF