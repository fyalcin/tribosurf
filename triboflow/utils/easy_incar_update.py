# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 10:37:43 2021

@author: Michael Wolloch

Tool to quickly update INCAR parameters for fizzled fireworks. More comfortable
than the "lpad update_fws" syntax. Put it in the path and run it like:
    python easy_incar_update -i <fw_id_1> <fw_id_2> ... -u '{"ALGO": "Fast"}' -r False 
"""
import argparse
import json

from fireworks import LaunchPad


def GetUserInput():
    """ Read the fw_ids, incar update and rerun toggle from the command line.

    Returns
    -------
    args : Namespace
        Parsed input parameters

    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--fw_id", dest = "fw_id",
                        nargs="+", type = int,
                        help="Space separated list of IDs of the fireworks to be updated")
    parser.add_argument("-u", "--incar_update", dest = "update_dict",
                        type = str,
                        help="Str holding a dictionary with VASP INCAR parameters to set. "
                             "Has to have single quotes outside and double quotes inside.")
    parser.add_argument("-f", "--functional", dest = "functional",
                        type = str,
                        help="Functional")
    parser.add_argument("-r", "--rerun", dest = "rerun",
                        type = str,
                        help="If 'True' will update the spec and then rerun the FWs")
    args = parser.parse_args()
    return args


def update_fws(id_list, update, rerun):
    """
    

    Parameters
    ----------
    id_list : [int]
        List of fw_ids
    update : dict
        Dictionary of VASP INCAR flags and corresponding values.
    rerun : bool
        Selects if the fireworks only get updated or also set to rerun.

    Returns
    -------
    None.

    """
    for i in id_list:
        fw = lp.get_fw_by_id(i)
        try:
            metagga = fw.spec['_tasks'][1]['vasp_input_set_params']['user_incar_settings']['METAGGA']
        except:
            metagga = False
        if metagga:
            #print('FW_ID {} is MetaGGA calc'.format(i))
            for k, v in update.items():
                lp.update_spec([i],
                    {"_tasks.1.vasp_input_set_params.user_incar_settings.{}".format(k): v})
        else:
            #print('FW_ID {} is GGA calc'.format(i))
            for k, v in update.items():
                lp.update_spec([i],
                    {"_tasks.0.vasp_input_set.user_incar_settings.{}".format(k): v})
        if rerun:
            lp.rerun_fw(i)

if __name__ == "__main__":
    lp = LaunchPad.auto_load()
    args = GetUserInput()
    id_list = args.fw_id
    update = json.loads(args.update_dict)
    rerun = ( args.rerun == 'True' )
    update_fws(id_list, update, rerun)
        
    