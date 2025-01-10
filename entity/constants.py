from dataclasses import dataclass, field


@dataclass(frozen=True)
class Constants:
    '''
    Dataclass consists of some essential constant data
    '''
    AMINOACID_LIST: tuple = (
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
        "TYR", "VAL", "HIE", "HID", "HIP", "CYX"
    )

    LIGAND: str = "SUS"
    RC_CUTOFF_DISTANCE: int = 12
    STANDARD_ATOMS_MAP =  {"C": 0, "O": 1, "N": 2, "H": 3, "S": 4, "P": 4}

