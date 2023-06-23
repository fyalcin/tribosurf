from fireworks import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize

from triboflow.utils.database import Navigator


@explicit_serialize
class FT_fake_vasp(FiretaskBase):
    _fw_name = "Fake VASP calculation"
    required_params = ["tag"]
    optional_params = []

    def run_task(self, fw_spec):
        tag = self.get("tag")
        nav_low = Navigator("auto")
        result = nav_low.find_data("tasks", {"task_label": tag})
        if result:
            return FWAction(update_spec=fw_spec)

        nav_high = Navigator("auto", "surfen_test")
        sample_data = nav_high.find_data("test_static", {"mpid": "mp-149"})
        sample_out_dict = sample_data["miller_list"]["111"][
            "74e1730eda410beb1ed6e744a7b1815e5fbf196c"
        ]["calcs"]["ouc"]
        sample_out_dict.update({"task_label": tag})

        nav_low.insert_data(collection="tasks", data=sample_out_dict)

        return FWAction(update_spec=fw_spec)
