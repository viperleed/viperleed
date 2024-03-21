# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import datetime
import sphinx_rtd_theme

project = 'ViPErLEED'
copyright = f'{datetime.date.today().year}, ViPErLEED-developers'
author = 'ViPErLEED-developers'
release = '0.11.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
              'sphinx.ext.autodoc',
#              'sphinx.ext.napoleon'        # for Numpy docstrings ?
              'sphinx.ext.extlinks',
              'sphinx.ext.mathjax',
              'sphinx.ext.viewcode',
              'sphinx.ext.intersphinx',
              'sphinx_rtd_theme',
              'sphinxcontrib.bibtex',       # to use .bib files for bibliography
              'sphinxcontrib.inkscapeconverter', # for SVG images
              'sphinxcontrib.spelling',     # spell checking for the docs ;)
              'myst_nb',                    # for including Jupyter notebooks
              ]

source_suffix = '.rst'

numfig = True # enumerate figures

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# Wordlist to ignore for spell checker
spelling_word_list_filename='spelling_wordlist.txt'

# For bibliography
bibtex_bibfiles = ['references.bib']
bibtex_default_style = 'unsrt' # use numbers
bibtex_reference_style = 'label' # use numbers in text too


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
#html_theme = 'sphinx_book_theme'
#html_theme = 'piccolo_theme'

html_favicon = '../guilib/icons/favicon.ico'

html_static_path = ['_static']

# RTD theme specific
html_theme_options = {
    'logo_only': False,
    'display_version': True,
    'prev_next_buttons_location': None,
    'style_external_links': True,
    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'includehidden': False,
}

html_css_files = [
    'css/theme_overrides.css',
]

# -- Options for LaTeX output ------------------------------------------------
# https://www.sphinx-doc.org/en/master/latex.html
latex_engine = 'xelatex'
nb_output_stderr = "remove-warn" # remove matplotlib user warnings
latex_show_urls = 'inline'
latex_show_pagerefs = False # show page numbers
inkscape_converter_args = ['--export-area-page']
latex_elements = {
    'papersize': 'a4paper',
    'passoptionstopackages': r'\PassOptionsToPackage{svgnames}{xcolor}',
    'figure_align': 'H', # Latex figure (float) alignment
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
    '''
    }
