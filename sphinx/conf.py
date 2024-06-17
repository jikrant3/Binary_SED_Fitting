# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('..'))
import src

# -- Project information -----------------------------------------------------

project = 'Binary_SED_Fitting'
copyright = '2024, Vikrant Jadhav'
author = 'Vikrant Jadhav'

# The full version, including alpha/beta/rc tags
release = '3.4.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx_copybutton',
              'sphinx.ext.napoleon',
              'sphinx.ext.autodoc',
              'sphinx.ext.autosectionlabel']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ['_static']
html_logo = "logo.png"
html_theme_options = {
    "sidebar_hide_name": True,
    "navigation_with_keys": True,
}