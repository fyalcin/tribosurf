#! /.fs/data/wolloch/atomate_test/atomate_env/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 14:33:47 2020

@author: mwo
"""
from fireworks.core.firework import FWAction, FiretaskBase

# =============================================================================
# Custom FireTasks
# =============================================================================
class MakeHeteroStructure(FiretaskBase):
    """
    Matches two slab systems to form a heterogeneous interface.
    
    Args:
        bottom_slab_name (str): Name of the structure of the bottom slab that
                                needs to be aligned (structure itself must be
                                in fw_spec).
        top_slab_name (str):    Name of the structure of the top slab that
                                needs to be aligned (structure itself must be
                                in fw_spec).
        parameters (dict):  Dictionary containing the key (values):
                            interface_distance (float): Distance between the
                                two materials.
                            max_area (float): Maximal cross section
                                area of the matched cell in Angstrom squared.
                                Defaults to 200.,
                            max_mismatch (float): Maximal allowed mismatch
                                between lattice vector length in %.
                                Defaults to 0.01,
                            max_angle_diff (float): Maximal allowed mismatch
                                between lattice vector angle in degrees.
                                Defaults to 1.
                            r1r2_tol (float): Tolerance between factors r1 and
                                r2 which relate to the matched cell areas as:
                                abs(r1*area1-r2*area2)/max_area <= r1r2_tol.
                                Defaults to 0.02
        
    Returns a matched interface structure
    """
    _fw_name = 'Make Hetero Structure'
    required_params = ['bottom_slab_name', 'top_slab_name', 'parameters']
    
    def run_task(self, fw_spec):
        from mpinterfaces.transformations import get_aligned_lattices
        from mpinterfaces.transformations import generate_all_configs
        
        bottom_aligned, top_aligned = get_aligned_lattices(
                fw_spec[self['bottom_slab_name']],
                fw_spec[self['top_slab_name']],
                max_area = self['parameters']['max_area'],
                max_mismatch = self['parameters']['max_mismatch'],
                max_angle_diff = self['parameters']['max_angle_diff'],
                r1r2_tol = self['parameters']['r1r2_tol'])
               
        #TODO: Find out if this is actually useful for the PES to return
        #all the interface structures we need, or if another stacking function
        #should be written with another input of lateral shifts.
        hetero_interfaces = generate_all_configs(top_aligned, bottom_aligned,
                            nlayers_2d=1, nlayers_substrate=1,
                            seperation=self['parameters']['interface_distance'])
        
        hetero_name = self['bottom_slab_name']+self['top_slab_name']
        hetero_interfaces[0].to(fmt='poscar', filename=
              'POSCAR_'+hetero_name+'.vasp')
        
        return FWAction(update_spec={'matched_interface': hetero_interfaces[0]}, 
                        stored_data={'matched_interface': hetero_interfaces[0]})


class MakeSlabFromFormula(FiretaskBase):
    """
    Makes a slab with defined thickness out of a bulk structure loaded
    from the MaterialsProject database.
    
    Args:
        formula (str):      Chemical formula of the structure (e.g. 'Fe2O3').
                            The structure with the lowest formation energy for
                            this formula is used unless MP_ID is specified.
        MP_ID (str):        Materials Project ID of the selected bulk structure
                            (e.g. 'mp-990448').
        miller (list):      List of miller indices to describe the selected
                            surface of the slab (e.g. [1, 1, 0]).
        thickness (float):  Minimal thickness for the slab in Angstrom.
        vacuum (float):     Thickness of the vacuum layer between repeating
                            slabs. Defaults to 25.0 Angstrom
        
    Returns a slab 
    """
    _fw_name = 'Make Slab from Formula'
    required_params = ['formula', 'miller', 'thickness']
    optional_params = ['MP_ID', 'vacuum']
    
    def run_task(self, fw_spec):
        from HelperFunctions import GetLowEnergyStructure
        from mpinterfaces.interface import Interface
        
        if self['MP_ID']:
            MP_ID = self['MP_ID']
        else:
            MP_ID = None
        
        if self['vacuum']:
            vacuum = self['vacuum']
        else:
            vacuum = 25.0
            
        #Get bulk structure from MP database
        bulk_structure, MP_ID_bulk = GetLowEnergyStructure(self['formula'],
                                               MP_ID=MP_ID,
                                               PrintInfo=False)

        
        #Construct the slab out of the bulk_structure
        slab = Interface(bulk_structure, hkl = self['miller'],
                   min_thick = self['thickness'],
                   min_vac = vacuum, primitive = False,
                   from_ase = True)
        
        miller = ''.join(str(e) for e in self['miller'])
        slab_name = self['formula'] + miller
        return FWAction(update_spec={slab_name: slab}, 
                        stored_data={slab_name: slab})

class WriteGeneralKpointsInputs(FiretaskBase):
    """
    Prepares the necessary files (POSCAR, PRECALC, maybe INCAR)
    for the K-Point Grid Generator of the Mueller group at John Hopkins
    http://muellergroup.jhu.edu/K-Points.html

    Args:
        structure (pymatgen structure object): Required input structure for
                                               which the Kpoints mesh should
                                               be created. 
        PrecalcDict (dict): Required parameters for the Grid Generator, should
                            definately include the MINDISTANCE flag.
        IncarDict (dict):   Optional INCAR file as a dictionary to pass
                            information about ISYM, MAGMOM, etc...
        WorkDir (str):      Optional working directory where everything happens
    """
    _fw_name = "Write Generalized Kpoints Inputs"
    required_params = ['structure', 'PrecalcDict']
    optional_params = ['IncarDict', 'WorkDir']

    def run_task(self, fw_spec):
        from HelperFunctions import WriteFileFromDict
        
        if self['WorkDir']:
            workdir = self['WorkDir']
        else:
            workdir = './'
        
        self['structure'].to(fmt='poscar', filename=workdir+'POSCAR')
        WriteFileFromDict(self['PrecalcDict'], workdir+'PRECALC')

        if self['IncarDict']:
            WriteFileFromDict(self['IncarDict'], workdir+'INCAR')
            

class GetKpoints(FiretaskBase):
    """
    Reads a KPOINTS file and pushes it to the spec

    Args:
        workdir (str): Optional working directory. Defaults
                       to './'
    Returns:
        FWAction
    """
    _fw_name = 'Get Kpoints'
    optional_params = ['WorkDir']

    def run_task(self, fw_spec):
        import os
        from pymatgen.io.vasp.inputs import Kpoints
        
        if self['WorkDir']:
            workdir = self['WorkDir']
        else:
            workdir = './'

        #Make a pymatgen Kpoints Object from a KPOINTS file
        KPTS = Kpoints().from_file(workdir+'KPOINTS')

        #Check for 'KPOINTS.log in workdir and read if it is there'
        if os.path.isfile(workdir+'KPOINTS.log'):
            with open(workdir+'KPOINTS.log') as f:
                Kpoints_Info = f.readlines()
            return FWAction(stored_data={'KPOINTS_Info': Kpoints_Info},
                             mod_spec=[{'_push': {'KPOINTS': KPTS}}])
        else:
            return FWAction(mod_spec=[{'_push': {'KPOINTS': KPTS}}])


class CleanUpGeneralizedKpointFiles(FiretaskBase):
    """
    Cleans up all the files needed and created by the 
    K-Point Grid Generator of the Mueller group at John Hopkins
    http://muellergroup.jhu.edu/K-Points.html
    """
    _fw_name = 'Clean Up Generalized Kpoint Files'
    optional_params = ['WorkDir']

    def run_task(self, fw_spec):
        from HelperFunctions import RemoveMatchingFiles
        
        if self['WorkDir']:
            workdir = self['WorkDir']
        else:
            workdir = './'

        RemoveMatchingFiles([workdir+'KPOINTS*',
                            workdir+'POSCAR*',
                            workdir+'INCAR', 
                            workdir+'PRECALC'])

class OLD_MakeGeneralizedKpoints(FiretaskBase):
    """
    This Firestask creates a generalized Kgrid using the 
    K-Point Grid Generator of the Mueller group at John Hopkins
    http://muellergroup.jhu.edu/K-Points.html
    Prerequisits:
                The shellscript "getKPoints" must be in the path
                of the computer and configered correctly
    Inputs:
                A pymatgen structure object needs to be provided in
                the spec as 'structure' as well as a dictionary containing
                PRECALC information in 'PrecalcDict'
                If an IncarDict input is found in the spec, The information
                in the INCAR (ISYM, MAGMOM, etc.) will also be used.
    Outputs:
                A pymatgen KPOINTS object is pushed to the spec, and
                saved alongside some generation data in the database.
    """

    from shutil import which

    _fw_name = "MakeGeneralizedKpoints"

    if which('getKPoints') is None:
        print('"getKPoints" could not be found. Make sure it is in the path'
              ' and configured correctly!\nCheck http://muellergroup.jhu.'
              'edu/K-Points.html for more information.')
        raise SystemExit('Could not find "getKPoints"')

    def run_task(self, fw_spec):
        """
        Prepares the necessary files (POSCAR, PRECALC, and INCAR) and then
        runs the getKPoints script. Kpoints generation data is read and
        stored, and the KPOINTS object is pushed to the spec.
        """
        import subprocess
        from pymatgen.io.vasp.inputs import Kpoints
        from HelperFunctions import WriteFileFromDict, RemoveMatchingFiles

        structure = fw_spec['structure']
        PrecalcDict = fw_spec['PrecalcDict']
        
        #Write the PRECALC and POSCAR, and if applicable, also the INCAR files
        WriteFileFromDict(PrecalcDict, 'PRECALC')
        structure.to(fmt='poscar', filename='POSCAR')
        if 'IncarDict' in fw_spec.keys():
            WriteFileFromDict(fw_spec['IncarDict'], 'INCAR')

        #get the generalized MonkhorstPack mesh using getKPoints
        get_kpoints_file = subprocess.Popen('getKPoints', cwd='./')
        get_kpoints_file.communicate()

        #Make a pymatgen Kpoints Object and rean info about creation process
        KPTS = Kpoints().from_file('KPOINTS')
        with open('KPOINTS.log') as f:
            Kpoints_Info = f.readlines()

        #Clean up unneccessary files
        RemoveMatchingFiles(['KPOINTS*', 'POSCAR*', 'INCAR', 'PRECALC'])
                  
        return FWAction(stored_data={'KPOINTS_Info': Kpoints_Info},
                        mod_spec=[{'_push': {'KPOINTS': KPTS}}])