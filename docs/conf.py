# molmass/docs/conf.py

import os
import sys

here = os.path.dirname(__file__)
sys.path.insert(0, os.path.split(here)[0])

import molmass

project = 'molmass'
copyright = '1993-2022, Christoph Gohlke'
author = 'Christoph Gohlke'
version = molmass.__version__

extensions = [
    'sphinx.ext.napoleon',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.doctest',
    # 'sphinxcontrib.spelling',
    # 'sphinx.ext.viewcode',
    # 'sphinx.ext.autosectionlabel',
    # 'numpydoc',
    # 'sphinx_issues',
]

templates_path = ['_templates']

exclude_patterns = []  # type: ignore

html_theme = 'alabaster'

html_static_path = ['_static']
html_css_files = ['custom.css']
html_show_sourcelink = False

autodoc_member_order = 'bysource'  # bysource, groupwise
autodoc_default_flags = ['members']
autodoc_typehints = 'description'
autodoc_type_aliases = {'ArrayLike': 'numpy.ArrayLike'}
autoclass_content = 'class'
autosectionlabel_prefix_document = True
autosummary_generate = True

napoleon_google_docstring = True
napoleon_numpy_docstring = False

html_theme_options = {
    'nosidebar': True,
}


def add_api(app, what, name, obj, options, lines):
    if what == 'module':
        lines.extend(('API', '---'))


def setup(app):
    app.connect('autodoc-process-docstring', add_api)
