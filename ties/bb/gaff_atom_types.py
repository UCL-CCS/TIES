import pathlib


def _get_element_map() -> dict[str, str]:
    """
    Read the GAFF atom types dictionary.
    """

    # Get the mapping of atom types to elements
    element_map_filename = (
        pathlib.Path(__file__).parents[1] / "data" / "element_atom_type_map.txt"
    )

    # remove the comments lines with #
    lines = filter(
        lambda line: not line.strip().startswith("#") and not line.strip() == "",
        open(element_map_filename).readlines(),
    )
    # convert into a dictionary

    element_map = {}
    for line in lines:
        element, atom_types = line.split("=")

        for atom_type in atom_types.split():
            element_map[atom_type.strip()] = element.strip()

    return element_map


gaff_element_map = _get_element_map()
