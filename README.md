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

> [!Note]
> A detailed example for this tool can be found
> [here](https://gist.github.com/Lars-B/0909ffde2938c6782273719fabf1bb06)


A commandline tool to extract date information about who-infected-whom.
After installation you can get more information by running
```bash
datesWIW --help
```

This tool will create a CSV file with the following columns
[Infector, Infectee, Start Date of infection, blockcount, Tree index]

The tree index is corrected for burnin and reflect the index of a tree in the input file without
burin removed.

To convert the tree distances to dates we need the taxon labels to have date information.
The default assumption is that the tip labels are foramted like this:
```
ID+YYYY-MM-DD
2+2007-11-01
```
> [!Caution]
> Use the flags `--date-sep` and `--date-format` to specify your layout
> Be aware that using a separator that is also present in the date format will not work!

> [!Important]
> It is also currently assumed that the float scale is in years, i.e. 1.0 branch length equals 1
> year

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
