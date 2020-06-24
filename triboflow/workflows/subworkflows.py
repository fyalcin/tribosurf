"""SubWorkflows to be used by the TriboFlow Project.

Created on Wed Jun 17 15:47:39 2020
@author: mwo
"""

from uuid import uuid4
from fireworks import Workflow, Firework
from triboflow.firetasks.encut_convergence import FT_EnergyCutoffConvo
from triboflow.firetasks.kpoint_convergence import FT_KpointsConvo

def ConvergeKpoints_SWF(structure, comp_parameters, spec, mp_id, functional,
                      kappa_start=500, kappa_incr=50,
                      n_converge=3, db_file=None):
    """Subworkflows that converges the the kmesh density via total energy.
    
    Takes a given structure, computational parameters (which includes the
    convergence criterion in eV/atom) and runs static vasp calculations with
    a denser and denser mesh (larger kappa parameter) until convergence in the
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
    kappa_start : int, optional
        Starting density value for the first run. Defaults to 500.
    kappa_incr : int, optional
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
                                   k_dens_incr = kappa_incr,
                                   k_dens_start = kappa_start,
                                   n_converge = n_converge)
    
    FW_CE = Firework(FT_KptsConvo, spec=spec,
                     name='Kpoint Convergence')
    WF = Workflow([FW_CE], name=name)
    return WF

def ConvergeEncut_SWF(structure, comp_parameters, spec, mp_id, functional,
                      deformations=None, encut_start=200, encut_incr=25,
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
        Starting encut value for the first run. Defaults to 200.
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
