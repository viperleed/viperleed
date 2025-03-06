"""Configuration file for the Sphinx documentation builder.

For the full list of built-in configuration values, see the documentation:
https://www.sphinx-doc.org/en/master/usage/configuration.html
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2022-10-17'
__license__ = 'GPLv3+'

import datetime

import sphinx_rtd_theme

import viperleed


# -- Project information ---------------------------------------------- https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'ViPErLEED'
copyright = f'{datetime.date.today().year}, ViPErLEED developers'
author = 'ViPErLEED developers'
release = viperleed.__version__


# -- General configuration -------------------------------------------- https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'myst_nb',                 # For including Jupyter notebooks
    'sphinxarg.ext',           # Auto-generation of argparse CLI docs
    'sphinx_copybutton',       # For a copy button in code blocks
    'sphinx_design',           # For tabs, dropdowns, ...
    'sphinx_rtd_theme',        # HTML theme
    'sphinx.ext.autodoc',      # for API documentation from docstrings
    'sphinx.ext.autosummary',  # for API documentation from docstrings
    'sphinx.ext.extlinks',     
    'sphinx.ext.intersphinx',  # For linking to other docs
    'sphinx.ext.mathjax',      
    'sphinx.ext.napoleon',     # For Numpy-style docstrings
    'sphinx.ext.todo',         # For TODOs
    'sphinx.ext.viewcode',
    'sphinxcontrib.bibtex',    # To use .bib files for bibliography
    'sphinxcontrib.inkscapeconverter',  # For SVG images
    'sphinxcontrib.spelling',  # spell checking for the docs ;)
    ]

source_suffix = '.rst'

numfig = True # enumerate figures

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# Spell-checker options
spelling_warning = True
spelling_word_list_filename = 'spelling_wordlist.txt'  # Ignore these

# For bibliography
bibtex_bibfiles = ['references.bib']
bibtex_default_style = 'unsrt' # use numbers
bibtex_reference_style = 'label' # use numbers in text too


# -- Options for HTML output ----------------------------------------- https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'

html_favicon = '../src/viperleed/guilib/icons/favicon.ico'
html_logo = '_static/viperleed_logo_oneline.svg'
html_static_path = ['_static']

# For tabs
myst_enable_extensions = ['colon_fence']

# RTD theme specific
html_theme_options = {
    'logo_only': True,
    'display_version': True,
    'prev_next_buttons_location': 'both',
    'style_external_links': True,
    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'includehidden': False,
    }

html_css_files = [
    'css/theme_overrides.css',
    ]


# -- Options for jupyter notebook execution ---------------------------
# See https://myst-nb.readthedocs.io/en/latest/index.html
# Use a specific kernel instead of the one stored in the notebook 
# metadata. This allows building the documentation on multiple OSes
nb_kernel_rgx_aliases = {'.+': 'python3'}
nb_execution_mode = "off"  # use cached version for now


# -- Options for LaTeX output -----------------------------------------
# https://www.sphinx-doc.org/en/master/latex.html
latex_engine = 'xelatex'
latex_toplevel_sectioning = 'part'
nb_output_stderr = 'remove-warn' # remove matplotlib user warnings
latex_show_urls = 'inline'
latex_show_pagerefs = False # show page numbers
inkscape_converter_args = ['--export-area-page']
latex_logo = '_static/viperleed_logo_circled.pdf'
latex_elements = {
    'papersize': 'a4paper',
    'passoptionstopackages': r'\PassOptionsToPackage{svgnames}{xcolor}',
    'figure_align': 'tb', # Latex figure (float) alignment
    'preamble': r'''
\usepackage{braket}
\usepackage[overlay,absolute]{textpos}% for header in PDF screen version
\usepackage{everypage}
\usepackage{newunicodechar}
\newunicodechar{α}{$\alpha$}
\newunicodechar{Δ}{$\Delta$}

\textblockorigin{28mm}{16.5mm} % position x,y wrt top-left corner of page
%\setlength{\TPHorizModule}{\pdfpagewidth} % text block width = page width
\setlength{\TPHorizModule}{\textwidth} % text block width = text width
\newlength{\chapterNameLength}%
''' + fr'''
% Modify the size of the ViPErLEED logo. See sphinx-doc/sphinx/issues/11930
\AtBeginDocument{{%
  \renewcommand{{\sphinxlogo}}{{%
    \sphinxincludegraphics[width=0.15\textwidth]{{{latex_logo.rsplit('/', 1)[1]}}}%
    \par%
    }}%
  \renewcommand{{\sphinxcrossref}}[1]{{#1}}%
}}%
''' + r'''
% Fix citations in figure captions. They are broken in sphinx when
% using sphinxcontrib-bibtex (as we do). See original issue #276 at
% https://github.com/mcmtroffaes/sphinxcontrib-bibtex/issues/276
\usepackage{etoolbox}
\AtBeginEnvironment{figure}{\pretocmd{\hyperlink}{\protect}{}{}}
''',
    }
