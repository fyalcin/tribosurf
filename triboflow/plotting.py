#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 17:10:39 2020

@author: mwo
"""
import numpy as np
import os
from matplotlib import pyplot as plt
from triboflow.utils.database import GetHighLevelDB

db_file='/home/mwo/FireWorks/config/db.json'

def plot_kmesh_convo(mp_id, functional, db_file=db_file):
    tribo_db = GetHighLevelDB(db_file)
    coll = tribo_db[functional+'.bulk_data']
    data = coll.find_one({'mpid': mp_id})
    k_data = data['k_dens_info']
    formula = data['formula']
    y = np.asarray(k_data['energy_list'])
    x = np.asarray(k_data['k_dens_list'])
    tolerance = k_data['Energy_tol_abs']*1000
    n_conv = k_data.get('n_converge', 3)
    y = (y-y[-1])*1000
    plt.grid(True)
    plt.xlim(min(x)-min(x)*0.1, max(x)+max(x)*0.1)
    plt.ylim(-tolerance*5,tolerance*5)
    plt.axvline(x[-n_conv],color='k',ls='--')
    plt.axhline(0,color='k')
    plt.axhline(tolerance, color='k', ls = ':')
    plt.axhline(-tolerance, color='k', ls = ':')
    plt.plot(x,y,'ro--', linewidth=1.5)
    plt.title('Kpoint convergence for {}.{} with MP-ID {}'
              .format(functional, formula, mp_id))
    plt.xlabel(r'kpoint density')
    plt.ylabel(r'$\Delta E$ [meV]')
    plt.tight_layout()
    plt.savefig("Kpoint_convergence.png", dpi=300)
    plt.gcf().clear()
    print('Fig is saved at {}'.format(os.getcwd()))
    return