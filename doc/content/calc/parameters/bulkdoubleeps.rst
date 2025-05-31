.. _bulkdoubleeps:

BULKDOUBLING_EPS
================

BULKDOUBLING_EPS defines the convergence criterion used to determine when
doubling of bulk layers is interrupted.

**Default**: BULKDOUBLING_EPS = 0.001

**Syntax**:

::

   BULKDOUBLING_EPS = <eps>

**Accepted values**: A single floating-point number,
0.0001 :math:`\leq` ``<eps>`` < 1. Typical <0.005, although the default value
of 0.001 should be appropriate. Change only in case there are some convergence
issues with the layer doubling, and only after checking whether the issues
arise from an input error in the geometry. Cannot be smaller than 0.0001 due
to Fortran reading it as an F7.4

**What is layer doubling?** During layer doubling, the layers defined as
bulk with :ref:`N_BULK_LAYERS`  will be used to determine
the scattered intensity due to the bulk. This follows a layer doubling
scheme, in which

#. the reflection/transmission matrices of one bulk layer are combined to
   create reflection/transmission matrices for two identical bulk layers.
#. the new reflection/transmission matrices of the 2-layers stack are
   combined to create reflection/transmission matrices for a 4-layers
   stack... and so on.

The doubling process is repeated a maximum number of times equal to
:ref:`BULKDOUBLEITER`, and until the reflection/transmission matrices 
do not change by more that BULKDOUBLING_EPS between two subsequent
doubling iterations.
