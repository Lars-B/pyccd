from dataclasses import dataclass
from .base import BaseClade


@dataclass(frozen=True)
class DemeClade(BaseClade):
    """
    Clade with a compartment annotation (deme) that indicates a location or similar

    Attributes:
        deme (string): Indicates a compartment this clade belongs to.
    """
    deme: str
