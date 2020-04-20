#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 17:04:35 2020

@author: mwo
"""

from fireworks import explicit_serialize, FiretaskBase, FWAction, ScriptTask, LaunchPad
from fireworks import Workflow, Firework
from fireworks.core.rocket_launcher import rapidfire

@explicit_serialize
class OldDetour(FiretaskBase):

    def run_task(self, fw_spec):
        print('Running the INITIAL FW')

        dt1 = Firework(ScriptTask.from_str('echo "this is intermediate job 1"'))
        dt2 = Firework(ScriptTask.from_str('echo "this is intermediate job 2"'))
        dt3 = Firework(ScriptTask.from_str('echo "this is intermediate job 3"'))
        wf = Workflow([dt1, dt2, dt3], links_dict={dt1: [dt2, dt3]})
        return FWAction(detours=wf)

@explicit_serialize
class WorkflowDetourClass(FiretaskBase):
    _fw_name = 'WorkflowDetourTask'
    required_params = []
    optional_params = []
    def run_task(self, fw_spec):
        links = {}
        new_fws_lst = []
        for i in range(3):
            hello_fw = Firework(ScriptTask.from_str("echo hello " + str(i)), name='h'+str(i))
            goodbye_fw = Firework(ScriptTask.from_str("echo goodbye " + str(i)), name='g'+str(i))
            new_fws_lst.append(hello_fw)
            new_fws_lst.append(goodbye_fw)
            links.update({hello_fw : [goodbye_fw]})
        wf = Workflow(new_fws_lst, links, name='sub_wf')
        return FWAction(detours = wf)

if __name__ == '__main__':
    fw0 = Firework(ScriptTask.from_str('echo start'), name = 'start')
    fw1 = Firework(OldDetour(), name = 'split')
    fw2 = Firework(ScriptTask.from_str('echo finish'), name = 'finish')
    wf = Workflow([fw0, fw1, fw2], {fw0 : [fw1], fw1: [fw2]}, name='test workflow')

    lpad = LaunchPad.auto_load() # loads this based on the FireWorks configuration
    lpad.add_wf(wf)
