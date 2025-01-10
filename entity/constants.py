from dataclasses import dataclass


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
    '''
    Value lists for BPS function
    '''
    LAMBDA_VALUES: tuple = (1, -1)
    ZETA_VALUES: tuple = (1., 2., 4., 8.)
    ETA_VALUES: tuple = (0.008, 0.04, 0.2, 1.)
    RS_VALUES: tuple = (2., 4., 6., 8., 10.)

