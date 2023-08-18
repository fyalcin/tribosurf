# Description: This file contains the connection to the Materials Project API
# it replaces parts of the utils.database module in the original code, because
# that module is not compatible with the new API.

import os
from typing import Union
from mp_api.client import MPRester
from pymatgen.core import Structure


class MPConnection:
    def __init__(self, mp_api_key: str = None):
        """
        Class to interact with the Materials Project API.

        Parameters
        ----------
        mp_api_key : str, optional
            Materials Project API key. If None, the key is read from the
            environment variable MP_API_KEY. The default is None.
        """

        if mp_api_key is None:
            self.mp_api_key = os.environ.get(
                "MP_API_KEY", "2ivreLaFVQ12DFP1u9h2ih94PokQdeyE"
            )
        else:
            self.mp_api_key = mp_api_key

    def get_low_energy_structure(
        self, chem_formula: str, mp_id: Union[str, None] = None
    ) -> tuple[Structure, str]:
        """
        Get the lowest energy structure of a material from the MP database.

        Parameters
        ----------
        chem_formula : str
            Chemical formula of the material.
        mp_id : str, optional
            Materials Project ID of the material. If None, the structure with
            the lowest energy for this formula will be returned.
            The default is None.

        Returns
        -------
        structure : pymatgen.core.Structure
            Structure of the material.
        mp_id : str
            Materials Project ID of the material.
        """

        if mp_id is None:
            mp_id = self.get_mpid_from_formula(chem_formula)

        with MPRester(api_key=self.mp_api_key) as mpr:
            structure = mpr.get_structure_by_material_id(mp_id)
        return structure, mp_id

    def get_mpid_from_formula(self, chem_formula: str) -> str:
        """
        Get the Materials Project ID of a material from its chemical formula.

        Parameters
        ----------
        chem_formula : str
            Chemical formula of the material.

        Returns
        -------
        str
            Materials Project ID of the material.
        """
        # if self.connection is None:
        #     self.__mp_connect()

        with MPRester(api_key=self.mp_api_key) as mpr:
            doc = mpr.summary.search(
                formula=chem_formula,
                energy_above_hull=[0.0, 0.0],
                fields=["material_id"],
            )
            if len(doc) == 0:
                raise ValueError(
                    "No material with this formula found in the Materials "
                    "Project database.\n"
                    "Please check the formula and try again."
                )
            elif len(doc) > 1:
                raise ValueError(
                    "More than one material with this formula found in the "
                    "Materials Project database.\n"
                    "Please check the formula and try again."
                )
        return doc[0].material_id

    def get_property_from_mp(self, mp_id: str, properties: list[str]) -> dict:
        """
        Get a property of a material from the Materials Project database.

        Parameters
        ----------
        mp_id : str
            Materials Project ID of the material.
        properties : list[str]
            Properties to get from the database.

        Returns
        -------
        dict
            Dictionary with the requested properties.
        """
        with MPRester(api_key=self.mp_api_key) as mpr:
            available_properties = mpr.summary.available_fields

            if not all([prop in available_properties for prop in properties]):
                raise ValueError(
                    "One or more of the requested properties is not "
                    "available in the SummaryRester of the Materials "
                    "Project database."
                )
            output_dict = mpr.summary.search(material_ids=[mp_id], fields=properties)[
                0
            ].dict()

            output_dict.pop("fields_not_requested")
        return output_dict

    def _print_all_resters(self):
        """
        Print all available resters of the Materials Project API.
        """
        with MPRester(api_key=self.mp_api_key) as mpr:
            for rester in mpr._all_resters:
                print(
                    f"{rester.suffix.replace('/', '_')}: "
                    f"{rester.__class__.__name__}"
                )

    def _print_all_properties(self, rester="summary"):
        """
        Print all available properties of the Materials Project API.
        """
        with MPRester(api_key=self.mp_api_key) as mpr:
            sub_rester = getattr(mpr, rester)
            for prop in sub_rester.available_fields:
                print(prop)


if __name__ == "__main__":
    db = MPConnection()
    mpid = db.get_mpid_from_formula("Fe2O3")

    prop_dict = db.get_property_from_mp(
        mpid,
        [
            "energy_per_atom",
            "formation_energy_per_atom",
            "band_gap",
            "formula_pretty",
            "energy_above_hull",
            "g_vrh",
        ],
    )
    print(prop_dict)

    struct = db.get_low_energy_structure(chem_formula="Fe2O3", mp_id=mpid)
    print(struct)

    # db._print_all_properties("elasticity")
