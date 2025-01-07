.. _planegroups:

..
    The |pg caption| replacement is on purpose a single line,
    because replacements cannot span multiple lines in sphinx.

.. |pg caption|  replace:: Overview of the 17 plane groups and possible symmetry operations. For each group, exemplary symmetry-equivalent atoms are shown as filled and empty circles in the left panel. The right panel shows the positions of the mirror, glide, and rotational symmetry elements. Symmetry elements drawn in red are related to those drawn in black via translations by integer multiples of the unit-cell vectors. Blue and black arrows in the left panels show examples of in-plane atomic displacements. Atoms linked by symmetry have symmetry-related displacements. Double arrows in groups pm, cm, and rcm indicate that movement is restricted along the mirror plane. Four-character plane-group designations are also indicated in gray, though ViPErLEED uses short forms in all input and output. Possible subgroups and additional notes are shown next to each group.

=====================
Plane symmetry groups
=====================

The symmetry of crystal surfaces can be classified as one of 17 possible plane
symmetry groups (also sometimes referred to as "wallpaper groups"). For more
information, see, for example, the corresponding
`Wikipedia article <https://en.wikipedia.org/wiki/List_of_planar_symmetry_groups>`_.

These groups are defined by the shape of the unit cell and by the symmetry
operations that leave the surface invariant: rotation axes, mirror planes,
and glide planes. :numref:`fig_plane_groups`, shows an overview of the 17
plane groups, allowed symmetry operations and possible subgroups. There,
an "×" indicates the point used by ViPErLEED as the conventional position
of the origin of the unit cell. See the relevant ViPErLEED publication for
more details :cite:p:`kraushoferViPErLEEDPackageCalculation2025`.

..
    In principle, one would like to do it this way:

    .. only:: not latex

        .. _fig_plane_groups:
        .. figure:: /_static/paper_figures/PlaneGroups_embedded.svg
            :alt: Overview of planegroups and possible symmetry operations.
            :align: center

            |pg caption|

    .. only:: latex

        .. _fig_plane_groups:
        .. figure:: /_static/paper_figures/PlaneGroups_two_columns_embedded.svg
            :alt: Overview of planegroups and possible symmetry operations.
            :align: center
            :height: 780px           <<<<<<<<<<<  TODO: ADJUST!

            |pg caption|

    So that different figures are used for html and latex versions.
    In fact, the current figure does not look great on the PDF. The
    text is very small. This, however, is not yet supported in sphinx,
    see https://github.com/sphinx-doc/sphinx/issues/4242.
    For now, stick to using the same figure for both. Switch after
    the issue is fixed.


.. _fig_plane_groups:

.. figure:: /_static/paper_figures/PlaneGroups_embedded.svg
    :alt: Overview of planegroups and possible symmetry operations.
    :align: center
    :height: 780px

    |pg caption|
