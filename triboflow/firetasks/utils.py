""" Utility Firetasks for the triboflow project.

Created on Wed Jun 17 15:59:59 2020
@author: mwo
"""

from fireworks import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize

from hitmen_utils.db_tools import VaspDB


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
    high_level : str or bool, optional
        If a string is given, the high-level database will be chosen based on
        that string. If True, the db.json file will be used to determine the
        name of the high_level database. The default is True.

    """

    _fw_name = "Update computational parameters"
    required_params = ["mpid", "functional", "new_params"]
    optional_params = [
        "db_file",
        "update_bulk",
        "update_slabs",
        "high_level",
    ]

    def run_task(self, fw_spec):
        mpid = self.get("mpid")
        functional = self.get("functional")
        new_params = self.get("new_params")
        update_bulk = self.get("update_bulk", True)
        update_slabs = self.get("update_slabs", False)
        db_file = self.get("db_file", "auto")
        hl_db = self.get("high_level", True)

        nav_high = VaspDB(db_file=db_file, high_level=hl_db)
        # get values from spec:
        new_data = {"$set": {}}
        for param in new_params:
            new_data["$set"][f"comp_parameters.{param}"] = fw_spec.get(
                param, None
            )

        if update_bulk:
            nav_high.update_data(
                collection=functional + ".bulk_data",
                fltr={"mpid": mpid},
                new_values=new_data,
                upsert=False,
            )

        if update_slabs:
            nav_high.update_data(
                collection=functional + ".slab_data",
                fltr={"mpid": mpid},
                new_values=new_data,
                upsert=False,
            )
        return


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
