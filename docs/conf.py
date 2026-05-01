# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'brokilon'
copyright = '2025, Lars Berling'
author = 'Lars Berling'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
        'sphinx.ext.autodoc',
        'sphinx.ext.viewcode',
        "sphinx.ext.autosummary",
        "sphinx.ext.coverage"
              ]

autosummary_generate = True

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'private-members': True,
    'show-inheritance': True,
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
html_theme = 'furo'
html_logo = "brokilon.png"
pygments_style = "tango"
pygments_dark_style = "monokai"

import os
import sys

from sphinx.ext.apidoc import main as apidoc_main


def run_apidoc(app):
    src_dir = os.path.abspath("../src/brokilon")
    output_dir = os.path.abspath("./auto_gen_docs")

    apidoc_main([
        "-o", output_dir,
        src_dir,
        "--force",        # overwrite existing files
        "--separate",     # one file per module (nicer)
    ])

def setup(app):
    app.connect("builder-inited", run_apidoc)
