import json
import os
from hitmen_utils.misc_tools import load_defaults
from pprint import pprint

from Levenshtein import ratio


def check_hetero_wf_inputs(inputs):
    """
    Check if the inputs for the heterostructure workflow are valid.
    """
    checked_inputs = {}
    triboflow_defaults = load_defaults("triboflow")
    surfgen_defaults = load_defaults("surfen")
    for key, value in inputs.items():
        if key.endswith("_1") or key.endswith("_2"):
            check_key = key[:-2]
        else:
            check_key = key
        if check_key in triboflow_defaults.keys():
            checked_inputs[key] = compare_keywords_and_add_defaults(
                input_dict=value,
                default_dict=triboflow_defaults[check_key],
                input_dict_name=key,
            )
        elif check_key in surfgen_defaults.keys():
            checked_inputs[key] = compare_keywords_and_add_defaults(
                input_dict=value,
                default_dict=surfgen_defaults[check_key],
                input_dict_name=key,
            )
        else:
            raise KeyError(
                f"Keyword {key} is not a valid keyword. "
                "Please only use the following keywords: "
                f"{list(triboflow_defaults.keys())+list(surfgen_defaults.keys())}"
            )
    # check if any of the final values are "NO_DEFAULT",
    # in which case raise an error
    for key, value in checked_inputs.items():
        for subkey, subvalue in value.items():
            if subvalue == "NO_DEFAULT":
                raise ValueError(
                    f"Keyword {subkey} in {key} NEEDS to be set by the user.\n"
                    "There is no default value for this keyword. "
                    "Please set it in the input file."
                )

    # extra check if either sg_params["miller"] or sg_params["max_index"]
    # is set (not None) for each material. If not, raise an error.
    for mat in ["_1", "_2"]:
        if (
            checked_inputs["sg_params" + mat]["miller"] is None
            and checked_inputs["sg_params" + mat]["max_index"] is None
        ):
            raise ValueError(
                f"Either sg_params{mat}['miller'] or sg_params{mat}['max_index'] NEEDS to be set by the user.\n"
                "Please set one of them in the input file."
            )
        elif (
            checked_inputs["sg_params" + mat]["miller"]
            and checked_inputs["sg_params" + mat]["max_index"]
        ):
            raise ValueError(
                f"Both sg_params{mat}['miller'] AND sg_params{mat}['max_index'] are set.\n"
                "Please set only one of them in the input file."
            )

    return checked_inputs


def compare_keywords_and_add_defaults(
    input_dict: dict,
    default_dict: dict,
    input_dict_name: str,
    lv_threshold: float = 0.6,
) -> dict:
    """
    Compare the keywords of the input dictionary with the keywords of the
    default dictionary. If a keyword is missing in the input dictionary, add
    the default value to the input dictionary.
    """

    out_dict = input_dict.copy()
    # check if all keywords are valid
    for key in input_dict.keys():
        if key not in default_dict.keys():
            for default_key in default_dict.keys():
                if (
                    ratio(key, default_key, score_cutoff=lv_threshold - 0.01)
                    > lv_threshold
                ):
                    raise KeyError(
                        f"Keyword {key} is not a valid keyword. Did you mean {default_key}?"
                    )
            raise KeyError(
                f"Keyword {key} is not a valid keyword."
                f"Please only use the following keywords in {input_dict_name}:\n"
                f"{pprint(default_dict.keys())}"
            )

    for key, value in default_dict.items():
        if key not in input_dict:
            out_dict[key] = value

    return out_dict
