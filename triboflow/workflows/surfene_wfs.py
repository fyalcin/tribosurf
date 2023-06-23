#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 14:47:15 2021

Workflows for calculating and converging the surface energy of a slab.

The module contains:

    ** SurfEneWF **:
        General class to work on crystalline slabs, workflows are static method.
        It includes the following methods:
            - conv_slabthick_surfene
            - conv_slabthick_alat
            - _check_subwf_params
    
    Functions:
    - check_choice

    Author: Gabriele Losi (glosi000)
    Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna

"""

__author__ = "Gabriele Losi"
__copyright__ = "Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna"
__contact__ = "clelia.righi@unibo.it"
__date__ = "February 9th, 2021"


import numpy as np
from fireworks import Workflow, Firework

from triboflow.utils.utils import create_tags, get_miller_str
from triboflow.utils.database import Navigator
from triboflow.firetasks.core import FT_RelaxStructure, FT_MoveTagResults
from triboflow.firetasks.slabs import FT_GenerateSlabs
from triboflow.firetasks.surfene import FT_SurfaceEnergy
from triboflow.utils.errors import SubWFError


class SurfEneWF:
    """
    Collection of static methods to calculate the surface energy of a slab.
    It also contains methods to converge the surface energy of a slab.

    """

    @staticmethod
    def conv_surface_energy(
        structure,
        mp_id,
        miller,
        functional="PBE",
        spec={},
        db_file=None,
        low_level=None,
        high_level=None,
        relax_type="slab_pos_relax",
        comp_params={},
        thick_min=4,
        thick_max=12,
        thick_incr=2,
        vacuum=10,
        in_unit_planes=True,
        ext_index=0,
        cluster_params={},
        parallelization=None,
        recursion=0,
        override=False,
        add_static=True,
    ):
        """
        Description of the method...

        """

        # Define the first Firework
        # =====================================================================

        # Define the thicknesses to create and simulate the slabs
        if parallelization == "low" or parallelization is None:
            if recursion:
                thickness = [thick_min]
            else:
                thickness = [0, thick_min, thick_max]

        elif parallelization == "high":
            try:
                thickness = np.arange(thick_min, thick_max, thick_incr)
                thickness = thickness.tolist()

                # If the array length is zero, the thick values are not well set
                assert len(thickness) > 0

                thickness = [0] + thickness + [thick_max]

            except:
                SubWFError("Thickness parameters are not defined correctly")

        else:
            raise SubWFError(
                "The value passed as parallelization is not known "
                "Allowed value: None, 'low', 'high'"
            )

        # Create the dictionary key where the unrelaxed slab will be saved
        formula = structure.composition.reduced_formula
        miller_str = get_miller_str(miller)
        slab_entry = [
            ["thickness", "data_" + str(thk), "input"] for thk in thickness
        ]

        # Generate the slabs and store them in the low level database, under
        # the dictionary key given by `slab_entry`
        ft_gen_slabs = FT_GenerateSlabs(
            structure=structure,
            mp_id=mp_id,
            miller=miller,
            collection=functional + ".slab_data",
            db_file=db_file,
            database=low_level,
            thickness=thickness,
            thick_bulk=thick_min,
            vacuum=vacuum,
            ext_index=ext_index,
            in_unit_planes=in_unit_planes,
            entry=slab_entry,
        )

        fw_gen_slabs = Firework(
            [ft_gen_slabs],
            spec=spec,
            name="Generate " + str(len(thickness)) + " slabs",
        )

        # Define the second Firework
        # ==================================================

        # Create the tags to store the calculation of the slabs
        tag_prefix = [
            formula + "_slab_" + miller_str + "_" + str(t) for t in thickness
        ]
        tags = create_tags(tag_prefix)

        # Start the navigator
        nav = Navigator(db_file=db_file, high_level=low_level)

        # Create the Firetasks to relax the structures
        fw_relax = []
        fw_move = []
        for thk, n, t in zip(thickness, slab_entry, tags):
            # Define different minor options to set correctly the firetasks
            res = "bulk_pos_relax" if thk == 0 else relax_type
            struct_kind = "bulk" if "bulk" in res else "slab"
            chkrel, chkmove = check_choice(
                t,
                nav,
                mp_id,
                miller,
                thk,
                functional,
                tol=1e-6,
                override=override,
            )

            ft_1 = FT_RelaxStructure(
                mp_id=mp_id,
                functional=functional,
                collection=functional + ".slab_data",
                entry=n,
                tag=t,
                db_file=db_file,
                database=low_level,
                relax_type=res,
                comp_params=comp_params,
                miller=miller,
                struct_kind=struct_kind,
                check_key=chkrel,
                add_static=add_static,
            )

            ft_2 = FT_MoveTagResults(
                mp_id=mp_id,
                collection_from="tasks",
                collection_to=functional + ".slab_data",
                tag=t,
                db_file=db_file,
                database_from=low_level,
                database_to=low_level,
                miller=miller,
                check_entry=chkmove,
                entry_to=[
                    ["thickness", "data_" + str(thk), "output", "structure"],
                    ["thickness", "data_" + str(thk), "output", "nsites"],
                    ["thickness", "data_" + str(thk), "output", "density"],
                    ["thickness", "data_" + str(thk), "output", "energy"],
                    [
                        "thickness",
                        "data_" + str(thk),
                        "output",
                        "energy_per_atom",
                    ],
                    ["thickness", "data_" + str(thk), "output", "bandgap"],
                    ["thickness", "data_" + str(thk), "output", "forces"],
                    ["thickness", "data_" + str(thk), "output", "stress"],
                    ["thickness", "data_" + str(thk), "output", "task_label"],
                ],
                entry_from=[
                    ["output", "structure"],
                    ["nsites"],
                    ["output", "density"],
                    ["output", "energy"],
                    ["output", "energy_per_atom"],
                    ["output", "bandgap"],
                    ["output", "forces"],
                    ["output", "stress"],
                    ["task_label"],
                ],
                struct_kind=struct_kind,
                override=False,
                cluster_params=cluster_params,
            )

            fw_1 = Firework(
                [ft_1],
                spec=spec,
                name="Relax structure: "
                + formula
                + " "
                + miller_str
                + ", thickness: "
                + str(thk),
            )
            fw_2 = Firework(
                [ft_2],
                spec=spec,
                name="Store structure in DB: "
                + formula
                + " "
                + miller_str
                + ", thickness: "
                + str(thk),
            )

            fw_relax.append(fw_1)
            fw_move.append(fw_2)

        # Define the third Firework(s)
        # ==================================================

        # Set the location of the energies in the high level DB
        entry_surfene = [["thickness", "data_0", "output"]]
        thk_loop = thickness if recursion else thickness[1:]
        for thk in thk_loop:
            entry_surfene.append(["thickness", "data_" + str(thk), "output"])

        ft_surfene = FT_SurfaceEnergy(
            mp_id=mp_id,
            collection=functional + ".slab_data",
            miller=miller,
            entry=entry_surfene,
            db_file=db_file,
            high_level_db=False,
        )

        fw_surfene = Firework(
            [ft_surfene], spec=spec, name="Calculate the Surface Energies"
        )

        # Define and return the Workflow
        # ==================================================

        # Build the workflow list and the links between elements
        wf_list = [fw_gen_slabs]

        # Case of parallelization = None, make serial calculations
        if parallelization is None and not recursion:
            links = {fw_gen_slabs: fw_relax[0]}

            for n, fw in enumerate(zip(fw_relax, fw_move)):
                fwr, fwm = fw
                wf_list.append(fwr)
                wf_list.append(fwm)
                links.update({fwr: fwm})
                if n < len(fw_relax) - 1:
                    links.update({fwm: fw_relax[n + 1]})

            links.update({fw_move[-1]: fw_surfene})

        # Other cases or when recursion is True (one element)
        else:
            links = {fw_gen_slabs: fw_relax}
            for fwr, fwm in zip(fw_relax, fw_move):
                wf_list.append(fwr)
                wf_list.append(fwm)
                links.update({fwr: fwm})
                links.update({fwm: fw_surfene})

        wf_list.append(fw_surfene)

        wf = Workflow(wf_list, links, name="Converge surface energy WF")

        return wf


def check_choice(
    tag, nav, mp_id, miller, thk, functional, tol=1e-5, override=False
):
    """
    Check if the results are already stored in the database. This is necessary
    at the moment to avoid that FT_MoveTagResults goes in conflicts with
    FT_RelaxStructure.

    """

    check_relax = None
    check_move = None

    if override:
        return check_relax, check_move

    data = nav.find_data(
        functional + ".slab_data", {"mpid": mp_id, "miller": miller}
    )

    if data is not None:
        try:
            d = data["thickness"]["data_" + str(thk)]["output"]
            old_tag = d["task_label"]
        except:
            d = None
            old_tag = None

        if d is not None:
            if (
                "energy" in d.keys()
                and "energy_per_atom" in d.keys()
                and "nsites" in d.keys()
            ):
                if old_tag == tag:
                    check_relax = ["thickness", "data_" + str(thk), "output"]
                    check_move = ["thickness", "data_" + str(thk), "output"]
                else:
                    calc = nav.find_data("tasks", {"task_label": old_tag})
                    if calc is not None:
                        if "output" in calc.keys():
                            out = calc["output"]
                            is_eq = (
                                abs(out["energy"] - d["energy"]) < tol
                                and abs(
                                    out["energy_per_atom"]
                                    - d["energy_per_atom"]
                                )
                                < tol
                                and d["nsites"] == calc["nsites"]
                            )
                        else:
                            is_eq = False

                    else:
                        is_eq = False

                    if is_eq:
                        check_relax = [
                            "thickness",
                            "data_" + str(thk),
                            "output",
                        ]
                        check_move = [
                            "thickness",
                            "data_" + str(thk),
                            "output",
                        ]

        else:
            calc = nav.find_data("tasks", {"task_label": tag})
            if calc is not None:
                if "output" in calc.keys():
                    out = calc["output"]
                    if (
                        "energy" in out.keys()
                        and "energy_per_atom" in out.keys()
                        and "nsites" in calc.keys()
                    ):
                        check_relax = [
                            "thickness",
                            "data_" + str(thk),
                            "input",
                        ]

    else:
        calc = nav.find_data("tasks", {"task_label": tag})
        if calc is not None:
            if "output" in calc.keys():
                out = calc["output"]
                if (
                    "energy" in out.keys()
                    and "energy_per_atom" in out.keys()
                    and "nsites" in calc.keys()
                ):
                    check_relax = ["thickness", "data_" + str(thk), "input"]

    return check_relax, check_move
