from dataclasses import dataclass
from .base import BaseClade


@dataclass(frozen=True)
class TransmissionAncestryClade(BaseClade):
    """
    Clade with transmission ancestor, i.e. who infected this clade.

    Attributes:
        transm_ancest (str): Transmission ancestor
    """
    transm_ancest: str
