.. title:: Workflow YAML Reference
.. _workflow YAML reference:

=======================
Workflow YAML Reference
=======================

Introduction
============

This document is a short reference for the features of atomate Workflows that you can control in YAML files. It aims to express all of the features that can make up a workflow. The benefit of YAML file workflows is that they are easy to understand and share, especially for non-programmers.

For details on the YAML format, refer to the `official YAML specification`_.

.. _official YAML specification: http://www.yaml.org/spec/1.2/spec.html

YAML Files in atomate
=====================

The following illustrates an example of a YAML file that can be used in atwf to run a workflow. Unless there is an existing YAML Workflow for the workflow you are trying to create, you will have to determine which required and optional parameters to set. Every Workflow in atomate is required to have a structure as the first parameter. This is implied in all of the YAML files and does not need to be included.

YAML format for the usual MP bandstructure workflow is given as follows:

.. code-block:: yaml

    fireworks:
    - fw: atomate.vasp.fireworks.core.OptimizeFW
    - fw: atomate.vasp.fireworks.core.StaticFW
      params:
        parents: 0
    - fw: atomate.vasp.fireworks.core.NonSCFUniformFW
      params:
        parents: 1
    - fw: atomate.vasp.fireworks.core.NonSCFLineFW
      params:
        parents: 1
    common_params:
      db_file: db.json
      $vasp_cmd: $HOME/opt/vasp
    name: bandstructure
    metadata:
      tag: testing_workflow

At the top there is often a comment (hashtag) describing the workflow (not shown here).

The `fireworks` key is a list of Fireworks; it is expected that
all such Fireworks have "structure" as the first argument and
other optional arguments following that. Each Firework is specified
via "fw": <explicit path>.

You can pass arguments into the Firework constructor using the special
keyword `params`, which is a dict. Any param starting with a $ will
be expanded using environment variables. If multiple fireworks share
the same `params`, you can use `common_params` to specify a common
set of arguments that are passed to all fireworks. Local params
take precedent over global params.

Another special keyword is `parents`, which provides
the *indices* of the parents of that particular Firework in the
list. The indices start at zero, i.e, the first Firework in your list
has zero. Thus, if you want the second Firework in the list to be a child
of the first Firework, you should specify a parent of 0 for the Firework.
Multiple parents are allowed. This allows you to link the Fireworks into a
logical workflow.

In the above example, we have:
* the first Firework (OptimizeFW) will run before anything else
* the second Firework (StaticFW) will run after the OptimizeFW is complete
* the third and fourth Fireworks (NonSCFUniformFW and NonSCFLineFW) will
run after the StaticFW is complete. Note these two Fireworks can run in parallel.

Next, `name` is used to set the Workflow name (structure formula +
name) which can be helpful in record keeping.

Finally, one can specify a `metadata` key as a YAML dict/hash that will
initialize workflow metadata - this is purely optional and for bookkeeping.

EOS Workflow Example
====================

This example shows what a more complicated workflow can look like using the YAML version of the EOS workflow described in the :ref:`running workflows tutorial`.

In order to use this example, create a file called ``eos.yaml`` with a text editor and enter the following text:

.. code-block:: yaml

    # EOS Workflow
    # An optimization Firework followed by 7 deformed structures based on the optimized structure
    # the deformations are +/- 10% volume of the original cell
    fireworks:
    - fw: atomate.vasp.fireworks.core.OptimizeFW
      user_incar_settings:
        SIGMA: 0.2
        ISMEAR: 1
    - fw: atomate.vasp.fireworks.core.TransmuterFW
      params:
        parents: 0
        transformations:
        - DeformStructureTransformation
        transformation_params:
        - "scaling_matrix": [[0.9655, 0, 0], [0, 0.9655, 0], [0, 0, 0.9655]]
    - fw: atomate.vasp.fireworks.core.TransmuterFW
      params:
        parents: 0
        transformations:
        - DeformStructureTransformation
        transformation_params:
        - "scaling_matrix": [[0.9773, 0, 0], [0, 0.9773, 0], [0, 0, 0.9773]]
    - fw: atomate.vasp.fireworks.core.TransmuterFW
      params:
        parents: 0
        transformations:
        - DeformStructureTransformation
        transformation_params:
        - "scaling_matrix": [[0.9888, 0, 0], [0, 0.9888, 0], [0, 0, 0.9888]]
    - fw: atomate.vasp.fireworks.core.TransmuterFW
      params:
        parents: 0
        transformations:
        - DeformStructureTransformation
        transformation_params:
        - "scaling_matrix": [[1.0000, 0, 0], [0, 1.0000, 0], [0, 0, 1.0000]]
    - fw: atomate.vasp.fireworks.core.TransmuterFW
      params:
        parents: 0
        transformations:
        - DeformStructureTransformation
        transformation_params:
        - "scaling_matrix": [[1.0110, 0, 0], [0, 1.0110, 0], [0, 0, 1.0110]]
    - fw: atomate.vasp.fireworks.core.TransmuterFW
      params:
        parents: 0
        transformations:
        - DeformStructureTransformation
        transformation_params:
        - "scaling_matrix": [[1.0217, 0, 0], [0, 1.0217, 0], [0, 0, 1.0217]]
    - fw: atomate.vasp.fireworks.core.TransmuterFW
      params:
        parents: 0
        transformations:
        - DeformStructureTransformation
        transformation_params:
        - "scaling_matrix": [[1.0323, 0, 0], [0, 1.0323, 0], [0, 0, 1.0323]]
    common_params:
      vasp_cmd: >>vasp_cmd<<
      db_file: >>db_file<<

To add this to your LaunchPad go to the folder containing your ``POSCAR`` (or other structure file) and ``eos.yaml``, run the following command to add the workflow to your LaunchPad:

.. code-block:: bash

    atwf add POSCAR -s eos.yaml

The YAML file format is typically considered easy to read, but it is less practical for more complicated workflows. The Python implementation of the EOS workflow is at :py:mod:`atomate.vasp.workflows.base.bulk_modulus` and it uses the existing deformation workflow to express the same as the above YAML file in less than 20 lines of Python code, including imports. Another advantage of using Python is being able to have more control over Fireworks and create them from Firetasks in the workflow, like the ``FitEOSToDb`` Firetask.
