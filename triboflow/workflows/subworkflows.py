"""SubWorkflows to be used by the TriboFlow Project.

Created on Wed Jun 17 15:47:39 2020
@author: mwo
"""

from uuid import uuid4
import numpy as np
from fireworks import Workflow, Firework
from triboflow.fireworks.common import RunPESCalcsFW, MakePESFW
from triboflow.firetasks.encut_convergence import FT_EnergyCutoffConvo
from triboflow.firetasks.kpoint_convergence import FT_KpointsConvo
from triboflow.firetasks.structure_manipulation import FT_MakeSlabInDB, \
    FT_StartSlabRelax, FT_GetRelaxedSlab
from triboflow.firetasks.PPES import FT_DoPPESCalcs, FT_FitPPES
from triboflow.utils.database import GetPropertyFromMP
from triboflow.utils.structure_manipulation import InterfaceName
from triboflow.utils.vasp_tools import GetEmin, GetCustomVaspStaticSettings

def CalcPES_SWF(top_slab, bottom_slab,
                top_mpid = None,
                bottom_mpid = None,
                functional = 'PBE',
                comp_parameters = {},
                file_output = False,
                output_dir = None,
                remote_copy = False,
                server = None, 
                user = None, 
                port = None):
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
    top_mpid : str, optional
        ID of the bulk material of the top slab in the MP database.
        The default is None.
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
        
    if top_mpid and bottom_mpid:
        interface_name = InterfaceName(top_mpid, top_miller,
                                       bottom_mpid, bot_miller)
    else:
        mt = ''.join(str(s) for s in top_miller)
        mb = ''.join(str(s) for s in bot_miller)
        interface_name = (top_slab.formula+'_'+mt+'_'+
                         bottom_slab.formula+'_'+mb+'_NO-MPIDs')
        
    
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
    
    FW_1 = RunPESCalcsFW(top_slab=top_slab, 
                         bottom_slab=bottom_slab, 
                         interface_name=interface_name,
                         functional=functional,
                         comp_parameters=comp_parameters,
                         tag=tag,
                         FW_name='Start PES calcs for '+interface_name)
    
    FW_2 = MakePESFW(interface_name=interface_name,
                     functional=functional,
                     tag=tag,
                     file_output=file_output,
                     output_dir=output_dir,
                     remote_copy = remote_copy,
                     server = server, 
                     user = user, 
                     port = port,
                     FW_name='Parse PES calcs for '+interface_name)
    
    SWF = Workflow([FW_1, FW_2], {FW_1: [FW_2]},
                   name = 'Calc PES for '+interface_name+' SWF')
    return SWF
    
def CalcPPES_SWF(interface_name, functional, distance_list = [-0.5, -0.25, 0.0,
                                0.25, 0.5, 2.5, 3.0, 4.0, 5.0, 7.5],
             out_name = 'PPES@minimum', structure_name = 'minimum_relaxed',
             spec = {}):
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

def MakeAndRelaxSlab_SWF(mp_id, miller, functional,
                       relax_type = 'slab_pos_relax',
                       bulk_struct_name = 'structure_equiVol',
                       slab_struct_name = 'unrelaxed_slab',
                       out_struct_name = 'relaxed_slab',
                       spec = {}):
    """
    Make and relax a slab.

    Parameters
    ----------
    mp_id : str
        Materials project ID tag to identify the bulk and slab entries in the
        high-level database.
    miller : list of int or str
        Miller indices of the slab to make.
    functional : str
        Which functional to use; has to be 'PBE' or 'SCAN'.
    relax_type : str, optional
        Which type of relaxation to run. See helper function:
        GetCustomVaspRelaxSettings. The default is 'slab_pos_relax'.
    bulk_struct_name : str, optional
        Name of the bulk to use as the basis for the slab generation in the 
        high-level database. The default is 'structure_equiVol'.
    slab_struct_name : str, optional
        Name of the unrelaxed slab in the high-level database.
        The default is 'unrelaxed_slab'.
    out_struct_name : TYPE, optional
        DESCRIPTION. The default is 'relaxed_slab'.
    spec : TYPE, optional
        DESCRIPTION. The default is {}.

    Returns
    -------
    SWF : TYPE
        DESCRIPTION.

    """
    if type(miller) == str:
        miller_str = miller
        miller = [int(k) for k in list(miller)]
    else:
        miller = miller
        miller_str = ''.join(str(s) for s in miller)
        
    formula = GetPropertyFromMP(mp_id, 'pretty_formula')
    
    tag = formula+miller_str+'_'+str(uuid4())
            
    FTs = []
    
    FTs.append(FT_MakeSlabInDB(mp_id = mp_id, miller = miller,
                               functional = functional,
                               bulk_struct_name = bulk_struct_name))
    
    FTs.append(FT_StartSlabRelax(mp_id = mp_id, miller = miller,
                                 functional = functional, tag = tag,
                                 slab_struct_name = slab_struct_name,
                                 relax_type = relax_type))
    
    FW = Firework(FTs, spec = spec,
                  name = 'Make and relax '+formula+miller_str+' slab')
    
    FW2 = Firework(FT_GetRelaxedSlab(mp_id = mp_id, miller = miller,
                                 functional = functional, tag = tag,
                                 struct_out_name = out_struct_name),
                   spec = spec,
                   name = 'Put relaxed '+formula+miller_str+' slab in DB')
    
    SWF = Workflow([FW, FW2], {FW: [FW2]},
                   name = 'Make and relax '+formula+miller_str+' SWF')
    return SWF


def ConvergeKpoints_SWF(structure, comp_parameters, spec, mp_id, functional,
                      k_dens_start=500, k_dens_incr=50,
                      n_converge=3, db_file=None):
    """Subworkflows that converges the the kmesh density via total energy.
    
    Takes a given structure, computational parameters (which includes the
    convergence criterion in eV/atom) and runs static vasp calculations with
    a denser and denser mesh (larger k_dens parameter) until convergence in the
    total energy is achieved. Output is printed to the screen and saved in the
    high-level triboflow database where it can be queried using the mp_id
    of the material.

    Parameters
    ----------
    structure : pymatgen.core.structure.Structure
        The structure for which to converge the K-pint grids.
    comp_parameters : dict
        Dictionary of computational parameters for the VASP calculations.
        Convergence criterion in eV/atom can be given here under the key:
        'energy_tolerence' and defaults to 0.001 (1meV/atom).
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

    Returns
    -------
    WF : fireworks.core.firework.Workflow
        A subworkflow intended to find the converged k_distance for a given
        structure.

    """
    formula = structure.composition.reduced_formula
    name = 'Kpoint Convergence SWF of '+formula
    
    tag = "Kpoints group for {} - {}".format(formula, str(uuid4()))
        
    FT_KptsConvo = FT_KpointsConvo(structure = structure,
                                   comp_params = comp_parameters,
                                   tag = tag,
                                   mp_id = mp_id,
                                   functional = functional,
                                   db_file = db_file,
                                   k_dens_incr = k_dens_incr,
                                   k_dens_start = k_dens_start,
                                   n_converge = n_converge)
    
    FW_CE = Firework(FT_KptsConvo, spec=spec,
                     name='Kpoint Convergence')
    WF = Workflow([FW_CE], name=name)
    return WF

def ConvergeEncut_SWF(structure, comp_parameters, spec, mp_id, functional,
                      deformations=None, encut_start=None, encut_incr=25,
                      n_converge=3, db_file=None):
    """Subworkflows that converges the Encut using a fit to an BM-EOS.
    
    Takes a given structure, computational parameters, and a optional list
    of deformations and uses these deformations to compute an
    Birch-Murgnahan equation of state for higher and higher energy cutoffs.
    Once bulk modulus and equilibrium volume do not change any longer,
    convergence is reached. Output is printed to the screen and saved in the
    high-level triboflow database where it can be queried using the mp_id
    of the material.

    Parameters
    ----------
    structure : pymatgen.core.structure.Structure
        The structure for which to converge the energy cutoff parameter.
    comp_parameters : dict
        Dictionary of computational parameters for the VASP calculations.
    spec : dict
        Previous fw_spec that will be updated and/or passed on for child
        Fireworks.
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
    db_file : str
        Full path to the db.json file that should be used. Defaults to
        None, in which case env_chk will be used in the FT.

    Returns
    -------
    WF : fireworks.core.firework.Workflow
        A subworkflow intended to find the converged k_distance for a given
        structure.

    """   
    name = 'Encut Convergence SWF of '+structure.composition.reduced_formula
    
    tag = "BM group: {}".format(str(uuid4()))
        
    
    if not encut_start:
        #Get the largest EMIN value of the potcar and round up to the
        #next whole 25.
        vis = GetCustomVaspStaticSettings(structure, comp_parameters,
                                          'bulk_from_scratch')
        emin = GetEmin(vis.potcar)
        encut_start = int(25 * np.ceil(emin/25))
    
    FT_EncutConvo = FT_EnergyCutoffConvo(structure = structure,
                                         comp_params = comp_parameters,
                                         tag = tag,
                                         mp_id = mp_id,
                                         functional = functional,
                                         db_file = db_file,
                                         encut_incr = encut_incr,
                                         encut_start = encut_start, 
                                         n_converge = n_converge)        
    
    FW_CE = Firework(FT_EncutConvo, spec=spec,
                     name='Encut Convergence')
    WF = Workflow([FW_CE], name=name)
    return WF
