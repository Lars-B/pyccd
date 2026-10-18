# Regular and extended CCDs

This is a python implementation for conditional clade distributions.
This package uses a simple map (dict) based implementation, presented by Larget.

Implemented versions are:

- regular CCD0 and CCD1
- phylogeography CCD1 (assuming fixed leaf labels)
- transmission tree CCD0 and CCD1 using BREATH like trees

Currently, under testing:
- SRanges implementation for sampled ancestors and fossil ranges

# Usage

> [!Caution]
> If you want to use a specific functionality but can not find documentation please let me know.

## Phylogeography CCD-MAP tree

The commandline tool for phylogeography trees is called `geography` and you can use the help to
learn more about the input options.
```bash
geography --help
```

## Transmission CCD-MAP tree

This package provide a command line interface to compute the tCCD MAP tree.
The tool is called `transcope` and for more information install the package and run
```bash
transcope --help
```

## Stratigraphic-range CCD-MAP tree

The command-line tool for computing stratigraphic-range CCD-MAP trees is called `fossirange`.
Use the help command to learn more about the available input options.

```bash
fossirange --help
```

### Visualizing the stratigraphic-range CCD graph

The command-line tool `fossigraph` generates an interactive Cytoscape visualization of the stratigraphic-range CCD graph as an HTML file.

Use the help command to learn more about the available input options.

```bash
fossigraph --help
```

## Regular CCD-MAP tree

These can be computed with the tool `ccdMAP`, use `--help` to find out how to use it.

## Majority rule consensus tree

Use the tool `mrc` to compute a majority rule consensus tree.
This is available for different versions of extended trees, see help.

## Breath helper infection date extraction

A commandline tool to extract date information about who-infected-whom.
After installation, you can get more information by running

```bash
breath-helper --help
```

There is a GUI wrapper with more details on this tool available at [Lars-B/brokilon_transmission_app](https://github.com/Lars-B/brokilon_transmission_app).

# Installation

## From GitHub (recommended)

Install the latest version directly from Github:

```bash
pip install brokilon@git+https://github.com/Lars-B/pyccd.git
```

## From a release

Alternatively, download the `.whl` file from the [latest GitHub release](https://github.com/Lars-B/pyccd/releases/latest) and install it using:

```bash
pip install /path/to/brokilon-<version>-py3-none-any.whl
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
