# CCD package for python

Simple implementation of conditional clade distributions in python.

## Currently WIP

- Add CCD0 and CCD2
- WIP transmission tree CCDs from BREATH BEAST2 package


# Installation

## Building from source

This uses the `hatchling` backend for building the package.
Clone the repository and then `cd` into the folder and run
```bash
python -m build
```
which will build the packge wheel file.
You can then use `pip` to install the freshly build package.
```bash
pip install dist/*.whl
```
