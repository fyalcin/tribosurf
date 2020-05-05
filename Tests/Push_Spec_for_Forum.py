#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from fireworks import Workflow, LaunchPad, FWAction, FiretaskBase, Firework
from fireworks.core.rocket_launcher import rapidfire
from fireworks.utilities.fw_utilities import explicit_serialize


def GetValueFromNestedDict(dictionary, key_list):
    """Return the value of a nested dict from a list of keys."""
    if len(key_list) > 1:
        if key_list[0] in dictionary:
            return GetValueFromNestedDict(dictionary[key_list[0]],
                                          key_list[1:])
        else:
            return None
    return dictionary.get(key_list[0])

@explicit_serialize
class FT_PrintSpec(FiretaskBase):
    """Prints the spec of the current Workflow to the screen."""
    
    def run_task(self, fw_spec):
        import pprint
        pprint.pprint(fw_spec)

@explicit_serialize
class FT_PushNumber(FiretaskBase):
    """Push a number to the end of an array in the spec using mod_spec _push."""
    required_params = ['number', 'numbers_loc']
    def run_task(self, fw_spec):
        out_str = '->'.join(self['numbers_loc'])
        number = self['number']
        return FWAction(mod_spec = [{'_push': {out_str: number}}])
    
@explicit_serialize
class FT_PushNumberandUpdate(FiretaskBase):
    """Push a number to the end of an array in the spec using mod_spec _push."""
    required_params = ['number', 'numbers_loc']
    def run_task(self, fw_spec):
        out_str = '->'.join(self['numbers_loc'])
        number = self['number']
        return FWAction(mod_spec = [{'_push': {out_str: number}}],
                        update_spec = fw_spec)
    
@explicit_serialize
class FT_PushNumberWithSet(FiretaskBase):
    """Push a number to the end of an array in the spec using mod_spec _set."""
    required_params = ['number', 'numbers_loc']
    def run_task(self, fw_spec):
        numbers = GetValueFromNestedDict(fw_spec, self['numbers_loc'])
        out_str = '->'.join(self['numbers_loc'])
        if numbers is None:
            number = [self['number']]
            return FWAction(mod_spec = [{'_set': {out_str: number}}])
        else:
            numbers.append(self['number'])
            return FWAction(mod_spec = [{'_set': {out_str: numbers}}])
    
def FW_PushNumber(number, number_loc, spec=None):
    FT1 = FT_PushNumber(number=number, numbers_loc=number_loc)
    FT2 = FT_PrintSpec()
    FW = Firework([FT1, FT2], spec=spec, name='Push and Print')
    return FW

def FW_PushNumberAndUpdate(number, number_loc, spec=None):
    FT1 = FT_PushNumberandUpdate(number=number, numbers_loc=number_loc)
    FT2 = FT_PrintSpec()
    FW = Firework([FT1, FT2], spec=spec, name='Push, Update and Print')
    return FW

def FW_PushNumberFourTimes(number1,  number2, number3, number4, number_loc,
                           spec=None):
    FT1 = FT_PushNumber(number=number1, numbers_loc=number_loc)
    FT2 = FT_PushNumber(number=number2, numbers_loc=number_loc)
    FT3 = FT_PushNumber(number=number3, numbers_loc=number_loc)
    FT4 = FT_PushNumber(number=number4, numbers_loc=number_loc)
    FT5 = FT_PrintSpec()
    FW = Firework([FT1, FT5, FT2, FT5, FT3, FT5, FT4, FT5], spec=spec,
                  name='Push 4 times and Print')
    return FW

def FW_PushNumberWithSet(number, number_loc, spec=None):
    FT1 = FT_PushNumberWithSet(number=number, numbers_loc=number_loc)
    FT2 = FT_PrintSpec()
    FW = Firework([FT1, FT2], spec=spec, name='Set and Print')
    return FW

if __name__ == "__main__":
    
    num_loc = ['mynumbers', 'numbers']
    lpad = LaunchPad.auto_load()
    
    wf1 = Workflow([FW_PushNumberFourTimes(1, 2, 3, 4, num_loc)],
                   name='Push in single FW')
    lpad.add_wf(wf1)
    
    FW1 = FW_PushNumber(1, num_loc)
    FW2 = FW_PushNumber(2, num_loc)
    FW3 = FW_PushNumber(3, num_loc)
    FW4 = FW_PushNumber(4, num_loc)
    wf2 = Workflow([FW1, FW2, FW3, FW4], {FW1: [FW2], FW2: [FW3], FW3: [FW4]},
                   name='Push in four separate FWs')
    lpad.add_wf(wf2)
    
    FW1 = FW_PushNumberWithSet(1, num_loc)
    FW2 = FW_PushNumberWithSet(2, num_loc)
    FW3 = FW_PushNumberWithSet(3, num_loc)
    FW4 = FW_PushNumberWithSet(4, num_loc)
    wf3 = Workflow([FW1, FW2, FW3, FW4], {FW1: [FW2], FW2: [FW3], FW3: [FW4]},
                    name='Push using _set in four separate FWs')
    lpad.add_wf(wf3)
    
    FW1 = FW_PushNumberAndUpdate(1, num_loc)
    FW2 = FW_PushNumberAndUpdate(2, num_loc)
    FW3 = FW_PushNumberAndUpdate(3, num_loc)
    FW4 = FW_PushNumberAndUpdate(4, num_loc)
    wf4 = Workflow([FW1, FW2, FW3, FW4], {FW1: [FW2], FW2: [FW3], FW3: [FW4]},
                   name='Push And Update in four separate FWs')
    lpad.add_wf(wf4)
    
    rapidfire(lpad)