# CCD package for python

Simple implementation of conditional clade distributions in python.

## Currently WIP

- Add CCD0 and CCD2
- WIP transmission tree CCDs from BREATH BEAST2 package

# Usage

## tCCD MAP tree

This package provide a command line interface to compute the tCCD MAP tree.
The tool is called `transcope` and for more information install the package and run
```bash
transcope --help
```

## WIW dates data extraction

A commandline tool to extract date information about who-infected-whom.
After installation you can run
```bash
datesWIW --help
```

This tool will create a CSV file with the following columns
[Infector, Infectee, Start Date of infection, 
End date of infection, type of infection (blockcount), Tree index]

This tool currently assumes that the tip labels in the tree are formated like this:
```
ID+YYYY-MM-DD
2+2007-11-01
```
> [!Caution]
> If this is not the case the tool will not work!

If you have another format you can either change the code
[here](https://github.com/Lars-B/pyccd/blob/main/src/pyccd/wiw_date_data.py#L18-L19) or let me know
which format you have and I can add it.

> [!Important]
> It is also assumed that the float scale is in years, i.e. 1.0 branch length equals 1 year

# Installation

A pre-release can be installed with
```bash
pip install https://github.com/Lars-B/pyccd/releases/download/v1.0.0-alpha/pyccd-0.1.0-py3-none-any.whl
```

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
