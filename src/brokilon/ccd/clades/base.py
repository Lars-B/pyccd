from dataclasses import dataclass
from functools import total_ordering


@total_ordering
@dataclass(frozen=True)
class BaseClade:
    """
    Base class representing a clade — a set of taxa/leaves.

    Attributes:
        clade (frozenset): A frozen set of taxa or node labels in this clade.
    """
    clade: frozenset

    def __len__(self) -> int:
        """Return the number of taxa/nodes in the clade."""
        return len(self.clade)

    def __lt__(self, other) -> bool:
        """Compare clade sizes for sorting or ordering."""
        if isinstance(other, BaseClade):
            return len(self.clade) < len(other.clade)
        return NotImplemented

    def __eq__(self, other) -> bool:
        """Clades are equal if their taxon sets are equal."""
        return isinstance(other, BaseClade) and self.clade == other.clade
