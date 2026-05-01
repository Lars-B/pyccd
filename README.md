# Regular and extended CCDs

This is a python implementation for conditional clade distributions.
This package uses a simple map (dict) based implementation, presented by Larget.

Implemented versions are:

- regular CCD0 and CCD1
- phylogeography CCD1 (assuming fixed leaf labels)
- transmission tree CCD0 and CCD1 using BREATH like trees

Currently in development:
- SRanges implementation for sampled ancestors and fossil ranges

# Usage

> [!Caution]
> If you want to use a specific functionality but can not find documentation please let me know.

## Phylogeography MAP tree

The commandline tool for phylogeography trees is called `geography` and you can use the help to
learn more about the input options.
```bash
geography --help
```

## tCCD MAP tree

This package provide a command line interface to compute the tCCD MAP tree.
The tool is called `transcope` and for more information install the package and run
```bash
transcope --help
```

## Regular CCD MAP

These can be computed with the tool `ccdMAP`, use `--help` to find out how to use it.

## Majority rule consensus tree

Use the tool `mrc` to compute a majority rule consensus tree.
This is available for different verions of extended trees, see help.

## Breath helper infection date extraction

> [!Note]
> A detailed example for this tool can be found
> [here](https://gist.github.com/Lars-B/0909ffde2938c6782273719fabf1bb06)

A commandline tool to extract date information about who-infected-whom.
After installation, you can get more information by running

```bash
breath-helper --help
```

This tool will create a CSV file with the following columns
```csv
Infector, Infectee, Start Date of infection, blockcount, Tree index
```

To convert the tree distances to dates we need the taxon labels to have date information.
The default assumption is that the tip labels are foramted like this:
```
ID+YYYY-MM-DD
2+2007-11-01
```
> [!Caution]
> Use the flags `--date-sep` and `--date-format` to specify your layout
> Be aware that using a separator that is also present in the date format will not work!

If this is not the case the tool will automatically create tip dates.
This is done by calibrating the leaf that is furthest from the root to be the current date.
If this is required, the tool will also produce a `csv` dataframe that contains 
the dates for all taxon labels and all trees in the input.

> [!Important]
> It is also currently assumed that the float scale is in years, i.e. 1.0 branch length equals 1
> year

# Installation

## From GitHub (recommended)

Install the latest version directly from Github:

```bash
pip install brokilon@git+https://github.com/Lars-B/pyccd.git
```

---

## Development install

Clone the repository and install locally:

```
git clone https://github.com/Lars-B/pyccd.git
cd pyccd
pip install -e .
```

---

## Optional: build manually

If you really want to build the wheel yourself:

```
python -m build
pip install dist/*.whl
```

Note that this requires additional packages `build` and `hatchling`.
