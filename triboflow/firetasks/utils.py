""" Utility Firetasks for the triboflow project.

Created on Wed Jun 17 15:59:59 2020
@author: mwo
"""

from pprint import pprint

from fireworks import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize
from atomate.utils.utils import env_chk

from triboflow.utils.database import Navigator, StructureNavigator
from triboflow.utils.structure_manipulation import interface_name
from triboflow.utils.utils import move_result


@explicit_serialize
class FT_CopyHomogeneousSlabs(FiretaskBase):
    """
    Firetask to copy aligned slabs from unrelaxed to relaxed.

    Since for homogeneous interfaces there cannot be any strain on the already
    relaxed slabs (from convergence), they do not have to be relaxed again.
    We just use this Firetask to copy them over in the interface_data
    collections.

    Parameters
    ----------
    mpid : str
        Material Project's material identifier ID.
    functional : str
        Functional for the identification of the high_level db.
    miller : str
        Miller indices of the slab
    external_pressure : float
        External pressure in GPa.
    db_file : str, optional
        Full path of the db.json. The default is 'auto'.
    high_level_db : str or bool, optional
        If a string is given, the high-level database will be chosen based on
        that string. If True, the db.json file will be used to determine the
        name of the high_level_db. The default is True.

    """

    _fw_name = "Copy slabs for homogeneous interfaces"
    required_params = ["mp_id", "functional", "miller", "external_pressure"]
    optional_params = ["db_file", "high_level_db"]

    def run_task(self, fw_spec):
        mpid = self.get("mp_id")
        functional = self.get("functional")
        miller = self.get("miller")
        pressure = self.get("external_pressure")
        db_file = self.get("db_file")
        if not db_file:
            db_file = env_chk(">>db_file<<", fw_spec)
        hl_db = self.get("high_level_db", True)

        nav = Navigator(db_file, high_level=hl_db)

        inter_name = interface_name(mpid1=mpid, mpid2=mpid, miller1=miller, miller2=miller)
        interface_data = nav.find_data(
            collection=functional + ".interface_data",
            fltr={"name": inter_name, "pressure": pressure},
        )
        top_slab = interface_data["top_aligned"]
        bot_slab = interface_data["bottom_aligned"]

        nav.update_data(
            collection=functional + ".interface_data",
            fltr={"name": inter_name, "pressure": pressure},
            new_values={
                "$set": {
                    "top_aligned_relaxed": top_slab,
                    "bottom_aligned_relaxed": bot_slab,
                }
            },
        )


@explicit_serialize
class FT_UpdateCompParams(FiretaskBase):
    """
    Firetask to update computational parameters for bulk and/or slabs

    Parameters
    ----------
    mpid : str
        Material Project's material identifier ID.
    functional : str
        Functional for the identification of the high_level db.
    new_params : list
        List of strings that identify the new keys that should be written to
        the computational parameters.
    db_file : str, optional
        Full path of the db.json. The default is 'auto'.
    update_bulk : bool, optional
        If the bulk entry for the given mpid should be updated.
        The default is True.
    update_slabs : bool, optional
        If the slab entries matching a given mpid should be updated (all miller
        indices). The default is False.
    high_level_db : str or bool, optional
        If a string is given, the high-level database will be chosen based on
        that string. If True, the db.json file will be used to determine the
        name of the high_level_db. The default is True.

    """

    _fw_name = "Update computational parameters"
    required_params = ["mpid", "functional", "new_params"]
    optional_params = [
        "db_file",
        "update_bulk",
        "update_slabs",
        "high_level_db",
    ]

    def run_task(self, fw_spec):
        mpid = self.get("mpid")
        functional = self.get("functional")
        new_params = self.get("new_params")
        update_bulk = self.get("update_bulk", True)
        update_slabs = self.get("update_slabs", False)
        db_file = self.get("db_file", "auto")
        hl_db = self.get("high_level_db", True)

        nav_high = StructureNavigator(db_file=db_file, high_level=hl_db)
        # get values from spec:
        new_data = {"$set": {}}
        for param in new_params:
            new_data["$set"][f"comp_parameters.{param}"] = fw_spec.get(param, None)

        if update_bulk:
            nav_high.update_data(
                collection=functional + ".bulk_data",
                fltr={"mpid": mpid},
                new_values=new_data,
                upsert=False,
            )

        if update_slabs:
            nav_high.update_many_data(
                collection=functional + ".slab_data",
                fltr={"mpid": mpid},
                new_values=new_data,
                upsert=False,
            )
        return


@explicit_serialize
class FT_MoveResults(FiretaskBase):
    """
    Firetask to move the result of a VASP calculation from the Fireworks database to the destination.

    Parameters
    ----------
    tag : str
        Task label to query for in the tasks collection of the Fireworks database.
    fltr : dict
        Filter dictionary that is used to query for the destination in the database.
    coll : str
        Collection to move the results into in the destination database.
    loc : str
        Location in the collection the results should be written to.
    custom_dict : dict, optional
        Custom dictionary to write into the destination. The default is {}.
    db_file : str, optional
        Full path of the db.json. The default is 'auto'.
    high_level : str, optional
        Whether to query the results from the high level database or not. The default is True.

    Returns
    -------
    FWAction that just updates the spec.

    """

    _fw_name = "Move the results from the tag to the entry found by flag."
    required_params = ["tag", "fltr", "coll", "loc"]
    optional_params = ["custom_dict", "db_file", "high_level"]

    def run_task(self, fw_spec):
        tag = self.get("tag")
        fltr = self.get("fltr")
        coll = self.get("coll")
        loc = self.get("loc")

        custom_dict = self.get("custom_dict", {})
        db_file = self.get("db_file", "auto")
        high_level = self.get("high_level", True)

        move_result(tag, fltr, coll, loc, custom_dict, db_file, high_level)

        return FWAction(update_spec=fw_spec)


@explicit_serialize
class FT_PrintFromBulkDB(FiretaskBase):
    _fw_name = "Print bulk data from DB"
    required_params = ["mp_id", "functional"]
    optional_params = ["db_file", "high_level_db"]

    def run_task(self, fw_spec):
        db_file = self.get("db_file", env_chk(">>db_file<<", fw_spec))
        mp_id = self["mp_id"]
        functional = self["functional"]

        hl_db = self.get("high_level_db", True)

        nav_structure = StructureNavigator(db_file=db_file, high_level=hl_db)
        bulk_dict = nav_structure.get_bulk_from_db(mp_id=mp_id, functional=functional)
        print("")
        pprint(bulk_dict)
        print("")


@explicit_serialize
class FT_PassSpec(FiretaskBase):
    """Update only certain keys in the first level of the spec.

    If the key_list contains only '_all', update the whole spec!
    """

    _fw_name = "Pass Spec"
    required_params = ["key_list"]

    def run_task(self, fw_spec):
        update = {}
        if self["key_list"] == ["_all"]:
            spec = fw_spec
            return FWAction(update_spec=spec)
        else:
            for k in self["key_list"]:
                if fw_spec.get(k) is None:
                    raise ValueError(
                        "{} can not be passed on to the next"
                        "FT/FW because it is not in the spec.\n"
                        "Currently the spec has the keys:\n"
                        "{}".format(k, fw_spec.keys())
                    )
                update[k] = fw_spec.get(k)
            return FWAction(update_spec=update)


@explicit_serialize
class FT_PrintSpec(FiretaskBase):
    """Prints the spec of the current Workflow to the screen.

    Not only prints the current spec in a pretty way, but also returns a
    FWAction that updates the spec of future to include the current spec.
    """

    _fw_name = "Print Spec"

    def run_task(self, fw_spec):
        import pprint

        pprint.pprint(fw_spec)

        spec = fw_spec

        return FWAction(update_spec=spec)
