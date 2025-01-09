from dataclasses import dataclass

@dataclass(frozen=True)
class Constants:
    AMINOACID_LIST : tuple = (
                    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
                    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
                    "TYR", "VAL", "HIE", "HID", "HIP", "CYX"
                    )

    LIGAND: str = "SUS"
    RC_CUTOFF_DISTANCE: int = 12