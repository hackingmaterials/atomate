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