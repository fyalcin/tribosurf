import os, glob

def RemoveMatchingFiles(list_of_patterns):
    """
    Remove all files matching the patterns (wildcards) in the current
    directory.
    """
    remove_list=[]
    for pattern in list_of_patterns:
        remove_list.extend(glob.glob(pattern))
    if remove_list is []:
        return
    else:
        for file_to_remove in remove_list:
            os.remove(file_to_remove)
        return