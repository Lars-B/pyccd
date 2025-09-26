from dataclasses import dataclass
from .base import BaseClade


@dataclass(frozen=True)
class TransmissionBlockClade(BaseClade):
    """
    Clade with a flag indicating whether a transmission block (event) occurred on the edge.

    Attributes:
        has_block (bool): Whether this clade was preceded by a block of transmissions.
    """
    has_block: bool
