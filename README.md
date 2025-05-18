# CCD package for python

Simple implementation of conditional clade distributions in python.

## Currently WIP

- Add CCD0 and CCD2
- WIP transmission tree CCDs from BREATH BEAST2 package

# Usage

This package provide a command line interface to compute the tCCD MAP tree.
The tool is called `transcope` and for more information install the package and run
```bash
transcope --help
```

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

Or alternatively after cloning the repository install directly with the 
following: 

```bash
pip install .
```
