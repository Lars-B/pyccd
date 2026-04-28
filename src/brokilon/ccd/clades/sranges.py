from .base import BaseClade
import inspect
from dataclasses import dataclass, field


def get_callsite():
    frame = inspect.currentframe().f_back.f_back
    return f"{frame.f_code.co_filename}:{frame.f_lineno} ({frame.f_code.co_name})"


@dataclass(frozen=True)
class SRangesClade(BaseClade):
    ancestral_range: str
    source: str = field(default_factory=get_callsite, compare=False)
