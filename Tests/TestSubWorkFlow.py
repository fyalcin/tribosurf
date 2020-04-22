#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 12:13:00 2020

@author: mwo
"""

from fireworks import FWAction, FiretaskBase, Firework, Workflow, LaunchPad
from fireworks.utilities.fw_utilities import explicit_serialize
from fireworks.core.rocket_launcher import rapidfire

@explicit_serialize
class FT_PrintSpec(FiretaskBase):
    """Print the spec of the current Workflow to the screen."""
    
    _fw_name = 'Print Spec'
    
    def run_task(self, fw_spec):
        import pprint
        
        pprint.pprint(fw_spec)
        
        spec = fw_spec
        return FWAction(update_spec = spec)
    

@explicit_serialize
class FT_AddNumbers(FiretaskBase):
    """Add two numbers (a, b) from the spec and puts the result c in the spec.
    """
    
    _fw_name = 'Add a b'
    required_params = ['name_a', 'name_b', 'res_name']
    
    def run_task(self, fw_spec):
        
        c = fw_spec[self['name_a']] + fw_spec[self['name_b']]
        
        spec = fw_spec
        spec.update({self['res_name']: c})
        return FWAction(update_spec = spec)

@explicit_serialize
class FT_PutInputInSpec(FiretaskBase):
    """Adds parameters to the spec
    """
    
    _fw_name = 'Add Parameters'
    required_params = ['start_1', 'start_2', 'stop_name']
    
    def run_task(self, fw_spec):
        
        if self['start_1'] < self['start_2']:
            smaller = self['start_1']
            larger = self['start_2']
        elif self['start_2'] < self['start_1']:
            smaller = self['start_2']
            larger = self['start_1']
        else:
            smaller = self['start_1']
            larger = self['start_1']
            
        stop = fw_spec[self['stop_name']]
            
        spec = fw_spec
        spec.update({'smaller': smaller, 
                     'larger': larger,
                     'stop_point': stop})
        return FWAction(update_spec = spec)
    
@explicit_serialize
class FibonacciAdderTask(FiretaskBase):

    _fw_name = "Fibonacci Adder Task"

    def run_task(self, fw_spec):
        smaller = fw_spec['smaller']
        larger = fw_spec['larger']
        stop_point = fw_spec['stop_point']

        m_sum = smaller + larger
        if m_sum < stop_point:
            print('The next Fibonacci number is: {}'.format(m_sum))
            # create a new Fibonacci Adder to add to the workflow
            spec = fw_spec
            spec.update({'smaller': larger, 'larger': m_sum, 'stop_point': stop_point})
            new_fw = Firework(FibonacciAdderTask(), spec=spec)
            return FWAction(stored_data={'next_fibnum': m_sum}, additions=new_fw)

        else:
            print('We have now exceeded our limit; (the next Fibonacci number would have been: {})'.format(m_sum))
            spec = fw_spec
            spec.update({'LastFiboNr': larger})
            return FWAction(update_spec = spec)
        
@explicit_serialize
class FT_StartFibo(FiretaskBase):
    """Start a Fibonacci Adder Workflow if the required inputs are in spec.
    """
    
    _fw_name = 'Fibo Start'
    
    def run_task(self, fw_spec):
        
        spec = fw_spec
        if 'smaller' in spec and 'larger' in spec and 'stop_point' in spec:
            print('')
            print('')
            print(' Trying to start new Workflow!')
            print('')
            print('')
            Fibo_FW = Firework(FibonacciAdderTask)
            WF = Workflow([Fibo_FW])
            return FWAction(detours=WF, update_spec=spec)
        else:
            return FWAction(update_spec = spec)

def DoStuff_FW(inputs):
    FT1 = FT_AddNumbers(name_a='to_add_1', name_b='to_add_2',
                        res_name='sum')
    FT2 = FT_PutInputInSpec(start_1=inputs['Inputs']['a'],
                            start_2=inputs['Inputs']['b'],
                            stop_name='sum')
    FT3 = FT_StartFibo()
    FT4 = FT_PrintSpec()
    return Firework([FT1, FT2, FT3, FT4], name='Initialize')


def PrintFW(name):
    FW=Firework(FT_PrintSpec(), name=name)
    return FW

input_dict = {'Other_Params': {'a': 1,
                               'b': 2,
                               'c': 3},
              'Inputs': {'a': 1,
                         'b': 0},
              'to_add_1': 12,
              'to_add_2': 88
              }

#Fib_FW = Firework(FibonacciAdderTask(), name='fib FW')    

#Print_FW = Firework([FT_PrintSpec()], spec=input_dict, name = 'Print Spec Firework')#
Print_FW = PrintFW(name='Print')
Initialize = Firework(FT_PrintSpec(), spec=input_dict,
                                name= 'Initialize Workflow')

FW2 = DoStuff_FW(input_dict) 
WF = Workflow([Initialize, FW2, Print_FW], {Initialize: [FW2], 
                                            FW2: [Print_FW]},
              name='test')

lpad = LaunchPad.auto_load() # loads this based on the FireWorks configuration
lpad.add_wf(WF)

#rapidfire(lpad)