import numpy as np

from typing import Union
from pymatgen.core.surface import Slab
from pymatgen.core.interface import Interface

from triboflow.phys.shaper import Shaper


def ext_pressure_to_force_array(
    slab: Union[Slab, Interface], external_pressure: float, tol: float = 0.1
) -> list:
    """
    This function takes a salb object and a pressure and returns a force array.

    Parameters:
    -----------
    slab: pymatgen.core.surface.Slab or pymatgen.core.interface.Interface
        The slab or interface object to which the force array will be applied.
    external_pressure: float
        The pressure to be applied to the structure.
    tol: float
        The tolerance used to determine the top and bottom layers of the slab.

    Returns:
    --------
    force_array: list
        The force array to be applied to the structure.
    """

    try:
        area = slab.surface_area
    except AttributeError:
        area = slab.interface_properties["area"]

    # get a layer dictionary that has the layer fractional z-coordinates as
    # keys, and the layer sites as values:

    layers = Shaper.get_layers(slab, tol)

    # get the sites of the top and bottom layers:
    top_layer = layers[max(layers)]
    bot_layer = layers[min(layers)]

    # get the force on the top and bottom layers:
    force_top = -external_pressure * area / len(top_layer)
    force_bot = external_pressure * area / len(bot_layer)

    # create the force array, making sure that all x and y components are zero.
    # the z component is the force on the top and bottom layers for the right
    # sites, and zero for the rest:

    force_array = []
    for i, site in enumerate(slab.sites):
        if i in top_layer:
            force_array.extend([0.0, 0.0, np.round(force_top, 4)])
        elif i in bot_layer:
            force_array.extend([0.0, 0.0, np.round(force_bot, 4)])
        else:
            force_array.extend([0.0, 0.0, 0.0])

    return force_array


# function for changing the force_array to a string, which groups 0.0s
# together, and removes the brackets and commas. 0.0 values are grouped
# together to reduce the size of the force array string, which is more
# convenient for the user to read and can be understood by VASP.
def force_array_to_string(force_array: Union[list, np.array]) -> str:
    """
    This function takes a force array and returns a string reable by VASP.

    Parameters:
    -----------
    force_array: list
        The force array to be converted to a string.

    Returns:
    --------
    force_string: str
        The force array as a string.
    """

    # first, slice the array in a way that blocks of 0.0 values are grouped
    # together and non-zero values are also grouped together:
    force_array = np.array(force_array)
    split_array = np.split(force_array, np.where(np.diff(force_array))[0] + 1)

    # now loop over the blocks, and if the whole block contains only 0.0s,
    # replace the block by "N*0", where N is the length of the block. If the
    # block contains non-zero values, replace the block by the values in the
    # block, separated by spaces:
    force_string = ""
    for block in split_array:
        if not np.any(block):
            force_string += str(len(block)) + "*0 "
        else:
            for value in block:
                force_string += str(value) + " "

    return force_string.strip()