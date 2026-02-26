from dataclasses import dataclass
from .base import BaseClade


@dataclass(frozen=True)
class SRangesClade(BaseClade):
    ancestral_range: str
