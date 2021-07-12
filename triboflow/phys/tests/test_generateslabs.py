#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 10:54:09 2021

Test the slab generations.

@author: glosi000

"""

from pymatgen.core.structure import Structure

from triboflow.utils.database import NavigatorMP
from triboflow.phys.solid_state import generate_slabs


thickness=[0, 2, 3, 4, 5, 6, 7, 8]
miller = [1, 1, 1]
miller_str = str(miller[0]) + str(miller[1]) + str(miller[2])

nav_mp = NavigatorMP()
Mg, _ = nav_mp.get_low_energy_structure(
   chem_formula='Mg', 
   mp_id='mp-110')

Cu, _ = nav_mp.get_low_energy_structure(
   chem_formula='Cu', 
   mp_id='mp-30')

Cu_local = Structure.from_file('POSCAR')

a = generate_slabs(structure=Mg, 
                   miller=miller, 
                   thickness=thickness, 
                   vacuum=10, 
                   thick_bulk=2,
                   ext_index=0,
                   in_unit_planes=True)

b = generate_slabs(structure=Cu_local, 
                   miller=miller, 
                   thickness=thickness, 
                   vacuum=10, 
                   thick_bulk=2,
                   ext_index=0,
                   in_unit_planes=True)

c = generate_slabs(structure=Cu_local, 
                   miller=miller, 
                   thickness=thickness, 
                   vacuum=10, 
                   thick_bulk=2,
                   ext_index=0,
                   in_unit_planes=True)

for el, thk in zip(a, thickness):
    print(el.lattice)
    print('')
    el.to('poscar', 'POSCAR_'+str(thk)+'_Mg_'+miller_str+'.vasp')
    
for el, thk in zip(b, thickness):
    print(el.lattice)
    print('')
    el.to('poscar', 'POSCAR_'+str(thk)+'_Cu_'+miller_str+'.vasp')

for el, thk in zip(c, thickness):
    print(el.lattice)
    print('')
    el.to('poscar', 'POSCAR_'+str(thk)+'_Cu_local_'+miller_str+'.vasp')
