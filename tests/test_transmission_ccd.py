"""
Tests for the transmission_ccd.py functions
"""
import pytest

from pyccd.transmission_ccd import TypeCCD, get_transmission_maps



def test_valid_enum():
    """
    Testing valid enum creation for TypeCCD
    """
    assert TypeCCD("Blocks") == TypeCCD.BLOCKS
    assert TypeCCD("Ancestry") == TypeCCD.ANCESTRY


def test_invalid_enum():
    """
    Testing invalid enum creation for TypeCCD
    """
    with pytest.raises(ValueError):
        TypeCCD("blocks")


def test_get_transmission_maps_wrong_enum_conversion():
    """
    Testing invalid type_str input for get_transmission_maps()
    """
    with pytest.raises(ValueError) as e:
        get_transmission_maps(trees=[], type_str="Bla")
    assert str(e.value).startswith("Type 'Bla' not recognized.")
