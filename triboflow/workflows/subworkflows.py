"""SubWorkflows to be used by the TriboFlow Project.

Created on Wed Jun 17 15:47:39 2020
@author: mwo
"""

from uuid import uuid4
from fireworks import Workflow, Firework
from triboflow.firetasks.encut_convergence import FT_EnergyCutoffConvo, \
    FT_PutEncutInfoInDB

from triboflow.helper_functions import GetCustomVaspRelaxSettings, \
    GetCustomVaspStaticSettings

def ConvergeEncut_SWF(structure, comp_parameters, spec, mp_id, functional,
                      deformations=None, encut_start=200, encut_incr=25,
                      n_converge=3, db_file='>>db_file<<'):
    """Subworkflows that converges the Encut using a fit to an BM-EOS.
    
    Takes a given structure, computational parameters, and a optional list
    of deformations and uses these deformations to compute an
    Birch-Murgnahan equation of state for higher and higher energy cutoffs.
    Once bulk modulus and equilibrium volume do not change any longer,
    convergence is reached.

    Parameters
    ----------
    structure : pymatgen.core.structure.Structure
        The structure for which to converge the K-pint grids.
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
        Number of calculations that have to show the same energy as the last
        one as to signify convergence, Defaults to 3.
    db_file : str
        Full path to the db.json file that should be used. Defaults to
        '>>db_file<<', to use env_chk.

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
                                         encut_start = encut_start)        
    
    FW_CE = Firework(FT_EncutConvo, spec=spec,
                     name='Encut Convergence')
    WF = Workflow([FW_CE], name=name)
    return WF
