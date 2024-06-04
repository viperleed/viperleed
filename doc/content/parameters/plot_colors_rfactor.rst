.. include:: /substitutions.rst

.. _plot_colors_rfactor:

PLOT_IV
=======

PLOT_IV allows specifying some formatting options for files containing
|IV| plots, like the :ref:`Rfactor_plots.pdf<Rfactorplots>` and
:ref:`THEOBEAMS.pdf<THEOBEAMS>`  files. Possible flags are
``axes``, ``colors``, ``legend``, ``overbar`` and ``perpage``,
explained below.

.. note::
   PLOT_IV was previously called PLOT_RFACTORS. To ensure backwards
   compatibility, PLOT_RFACTORS is still a permissible alias for PLOT_IV.

PLOT_IV axes / borders
----------------------

``PLOT_IV axes`` (or ``borders``) defines which axes are drawn for each panel.

**Default**: ``all``

**Syntax:**

::

   PLOT_IV axes = less     ! draw only left and bottom axes. Equivalent to 'lb'
   PLOT_IV axes = bottom   ! draw only bottom axes. Equivalent to 'b'

**Acceptable values:** ``all`` / ``less`` (or ``lb``) / ``bottom``
(or ``b``) / ``none``

Even if "none" is set, the bottom axis for the bottom-most plots will always
be drawn, and labelled with the energy steps.

PLOT_IV colors
--------------

``PLOT_IV colors`` defines what colors should be used when plotting |IV|
curves. In the :ref:`Rfactor_plots.pdf<Rfactorplots>`  file, the first
value specifies the color of the |IV| curve corresponding to the theoretical
beam, the second value for the experimental beam. Only the first value is used
in :ref:`THEOBEAMS.pdf<THEOBEAMS>`, additional colors will be ignored. If fewer
colors than required are specified, the remaining curves will use matplotlib
defaults.

**Default**: uses default matplotlib colors (tab:blue and tab:orange)

**Syntax:**

::

   PLOT_IV colors = red blue
   PLOT_IV colors = tab:blue tab:orange   ! the default
   PLOT_IV colors = #000000 #1c6d69

**Acceptable values:** HTML hex strings or any
`named color <https://matplotlib.org/3.1.0/gallery/color/named_colors.html>`__
recognized by matplotlib.

Expects string values, separated by space. The values will be passed to
matplotlib as string without further interpretation; if matplotlib cannot
interpret them, colors will fall back on the matplotlib default.

PLOT_IV legend
--------------

``PLOT_IV legend`` defines how often the legend is printed. This affects only
the legend labelling the experimental/theoretical beams. The beam indices and
R factors per beam are always printed.

**Default**: ``all``

**Syntax:**

::

   PLOT_IV legend = first     ! draw legend only on the first plot on each page
   PLOT_IV legend = topright  ! draw legend only in the top-right plot on each page; equivalent to 'tr'
   PLOT_IV legend = none      ! don't draw the legend

**Acceptable values:** ``all`` / ``first`` / ``topright`` (or ``tr``) / ``none``

PLOT_IV font_size
-----------------

``PLOT_IV font_size`` defines the font size used
for the labels and legends in the |IV| plots.

**Default**: 10

**Syntax:**

::

   PLOT_IV font_size = 15

**Acceptable values:** float values > 0

PLOT_IV line_width
------------------

``PLOT_IV line_width`` defines the width of
the lines used to plot the |IV| curves.

**Default**: 1.0

**Syntax:**

::

   PLOT_IV line_width = 2.0   ! use line_width 2.0

**Acceptable values:** float values > 0

PLOT_IV overbar / overline
--------------------------

Replaces the minus sign for negative-valued beam indices by an overline.

**Default**: False

**Syntax:**

::

   PLOT_IV overline = True
   PLOT_IV overbar = True     ! equivalent

**Acceptable values:** T(rue), F(alse), not case sensitive

PLOT_IV perpage / layout
------------------------

``PLOT_IV perpage`` (or ``layout``) defines how many panels (i.e., how
many beams) should be rendered on one page of the
:ref:`Rfactor_plots.pdf<Rfactorplots>` and :ref:`THEOBEAMS.pdf<THEOBEAMS>`
files.

**Default**: 2

**Syntax:**

::

   PLOT_IV perpage = 1     ! separate page per beam
   PLOT_IV perpage = 8     ! 2 columns, 4 rows
   PLOT_IV perpage = 21    ! 3 columns, 7 rows
   PLOT_IV layout = 2 2   ! 2 columns, 2 rows
   PLOT_IV layout = 3 4   ! 3 columns, 4 rows
   PLOT_IV perpage = 3 4   ! 3 columns, 4 rows - equivalent to above

**Acceptable values:** Single positive integer, or tuple of two
positive integers

If two values are given, these will be interpreted as a ``(columns, rows)``
layout instruction, which may or may not work well. If a single integer ``N``
is given, layout is automatically chosen as:

::

   columns = round(sqrt(N/2))
   rows = ceil(N/columns)

The width of the figure is fixed as 7 inch. Height is adapted as needed,
with panels keeping a 2:1 aspect ratio.
