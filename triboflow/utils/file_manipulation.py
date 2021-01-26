import os, glob
from fireworks import ScriptTask, FileTransferTask

def RemoveMatchingFiles(list_of_patterns):
    """Remove all files matching the patterns in the directory."""
    remove_list=[]
    for pattern in list_of_patterns:
        remove_list.extend(glob.glob(pattern))
    if remove_list is []:
        return
    else:
        for file_to_remove in remove_list:
            os.remove(file_to_remove)
        return
    
def CopyOutputFiles(file_list,
                    output_dir,
                    remote_copy=False,
                    server=None,
                    user=None,
                    port=None):
    """Return a Firetask that copys output files locally or to remote server.
    
    Handles file copy from the work-folder to a chosen output directory which
    might be on a differnet machine. In that case scp will be used via
    a ScriptTaks, as the FileTransferTask was found to be unreliable sometimes
    for remote copies.

    Parameters
    ----------
    file_list : list of str
        Filenames to be copied.
    output_dir : str
        Location the output files are copied to if file_output is selected.
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
    FT : fireworks.user_objects.firetasks.fileio_tasks.FileTransferTask or
         fireworks.user_objects.firetasks.script_task.ScriptTask

    """
    if remote_copy:
        if server and user:
            to_copy = ' '.join(file_list)
            scp_str = 'scp {} {}@{}:{}/.'.format(to_copy, user, server,
                                                 output_dir)
            if port:
                scp_str = 'scp -P {} {} {}@{}:{}/.'.format(port, to_copy,
                                                           user, server,
                                                           output_dir)
            FT = ScriptTask.from_str(scp_str)
        else:
            out_str = ("You have requested remote_copy but "
                       "did not specify a remote server "
                       "and/or username!\n"
                       "No copy will be performed!\n")
            FT = ScriptTask.from_str('echo "{}"'.format(out_str))
            
    else:
        FT = FileTransferTask({'files': file_list,
                               'dest': output_dir,
                               'mode': 'copy'})
    return FT