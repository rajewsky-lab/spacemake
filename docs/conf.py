# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'spacemake'
copyright = '2021-2024, Rajewsky Lab'
author = 'Tamas Ryszard Sztanka-Toth, Marvin Jens, Nikos Karaiskos, Nikolaus Rajewsky'

version = '0.8.0'
release = version

# -- General configuration

extensions = [
    "sphinx_rtd_theme",
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.autosectionlabel',
    'nbsphinx'
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}

intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output
html_theme = "sphinx_rtd_theme"
html_theme_options = {
    'navigation_depth': 3
}

# -- Options for EPUB output
epub_show_urls = 'footnote'

import os
import sys
sys.path.insert(0, os.path.abspath('../'))
