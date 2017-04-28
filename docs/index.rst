.. title:: atomate (Materials Science Workflows)

.. image:: _static/atomate_logo_small.png
   :alt: atomate logo

.. pull-quote:: | "Civilization advances by extending the number of important operations which we can perform without thinking about them."
                |    - Alfred North Whitehead

=======
atomate
=======

Theory and computations are powerful tools for understanding and designing materials, but conventional
software for performing these computations are still difficult to use, understand, and automate.
atomate makes it possible to perform complex materials science computations using
very straightforward statements. Features of atomate include:

* It is built on top of state-of-the-art open-source libraries: **pymatgen**, **custodian**, and **FireWorks**. Building off these libraries means that atomate can not only serve as a simple and friendly introduction to computational materials science, but that it is powerful enough for even the most demanding of theory users that require precise control and massive execution.
* It is easy to get "standard" workflows for a wide variety of desired materials properties - optimized structures, band structures, electronic transport properties, dielectric constants, and much more. Just provide a crystal structure (that's it!) and let atomate set up a complete workflow that provides the property you are interested in. You can do this for a single material, 100 materials, or 100,000 materials.
* One can easily change "standard workflows" - whether that is changing some of the default calculation parameters or recomposing the workflow (adding new calculations, removing steps, etc.) - using a very expressive syntax. One can compose very complex new workflows simply by chaining together pre-built calculation steps.
* A system of "powerups" that let you quickly decorate a bare workflow with useful special properties. Just feed the workflow through the powerup and your workflow will have the feature enabled. A config file allows you to automatically set the powerups you want to apply most often.
* It can build large databases of output properties that you can query, analyze, and share in a systematic way.
* It automatically keeps meticulous records of jobs, their directories, runtime parameters, etc.
* Jobs can be run on a variety of computing systems, queue systems, and architectures.
* atomate uses a standard interface for adding new types of calculations and workflows such that it is possible for users to contribute new features and grow the capabilities of the software over time.

Note: that atomate is currently built to work with the **VASP** electronic structure software, but it is the intention of atomate to support a variety of software.

============
Installation
============

Although atomate itself should be pip-installable, the actual process to get an infrastructure to run VASP at supercomputers might be a bit more involved (not difficult, just more steps). A detailed :ref:`installation guide <installation tutorial>` will walk you through the setup.

If you've installed before and need a quick reference on what to do, the general steps are:

1. Set up pseudopotentials and MongoDB database(s)
#. Create a Python environment
#. Install the atomate Python packages
#. Configure FireWorks
#. Configure pymatgen

An experimental installer to bootstrap your environment with minimal configuration can be found in the `PhasesResearchLab/install-atomate GitHub repo`_. Several preset systems are available including Edison, Cori, and Stampede, but the installer can also be used more generally with some configuration.

.. _PhasesResearchLab/install-atomate GitHub repo: https://github.com/PhasesResearchLab/install-atomate


==============================
Testing the VASP functionality
==============================

In order to use the VASP functionality, make sure you set up PMG_VASP_PSP_DIR variable (see pymatgen docs). Also, make sure you have MongoDB running in order to execute all the tests.

To test the VASP functionality, first run the unit tests in ``atomate.vasp.tests``. These unit tests are designed to run without installing VASP. Some of them start with a VASP workflow but apply the ``use_fake_vasp`` method to replace calling the VASP executable with a "Faker" that verifies basic properties of the inputs and copies pre-stored output files to the current directory, thus simulating the execution of VASP. Anyway this will help make sure your installation is in good shape.

Note that the unit tests in atomate/vasp/tests/test_vasp_workflows.py can be modified to actually run VASP by setting VASP_CMD to a String representing your VASP command. If you need to debug at a later point, this might be something to refer back to. Furthermore, many tests have a DEBUG option that can sometimes help in finding problems. Sometimes you need to toggle DEBUG on/off a couple of times if you are doing this to make sure all the old data is actually cleared between debug runs; the tearDown() and setUp() methods are still a bit finicky.

==========================
Learning to use atomate
==========================

If you are familiar with (i) VASP, (ii) pymatgen, (iii) custodian, and (iv) FireWorks, then most of atomate such be fairly straightforward. For example, the FireTasks implemented in ``atomate/vasp/firetasks`` should look make at least *some* sense, and the Fireworks implemented in ``atomate/vasp/fireworks`` should also seem logical and mostly clear. Workflows are simply chains of Fireworks (technically, DAGs). Normally, they would be implemented in simple Python, i.e. see the FireWorks codebase about how to compose Workflows with Python, but it turns out they are simple enough that one can write them in a simple YAML text file instead of Python code. There is a custom YAML format that is described in the README for the ``atomate/vasp/workflows/base/library`` folder.

In practice, getting prebuilt workflows is easier than this. For this, just look in ``atomate/vasp/workflows/presets``. This folder contains functions where one can simply give a crystal structure and get back an appropriate workflow. Nothing to it!

There are only a couple of new concepts in atomate that you might need to familiarize yourself with, and they are described below.

The "env_chk", e.g. >>db_file<< syntax
======================================

One issue in coding workflows is what to do when different machines require different settings. For example, the path to the VASP executable or the path to a file containing database credentials might be located in different places on different machines. For users wanting to run on multiple machines, such parameters cannot be hard-coded. However, users that are running on a single machine, or those that are testing things out, might prefer to hard-code those parameters.

The ``env_chk`` functionality is a way to support both hard-coding of parameters as well as letting the machine (or more specifically, the FireWorker) set the parameter. Many of the FireTasks in atomate, e.g., ``RunVaspDirect``, state in the docs that they "support ``env_chk``" for a parameter such as ``vasp_cmd``. What this means is that you have two options for creating the FireTask:

Option 1 is to use something like ``my_task = RunVaspDirect(vasp_cmd="vasp")``. This behaves exactly as you would expect in regular Python, i.e., the string literal "vasp" set as the ``vasp_cmd`` parameter.

Option 2 is to use the ``env_chk`` notation which looks like this: ``my_task = RunVaspDirect(vasp_cmd=">>my_vasp_cmd<<")``. If ``env_chk`` parameters like `vasp_cmd`` are enclosed in the ``>><<`` symbols, it is interpreted that the user wants to get the values from the FireWorker's ``env`` value. That is, when executing the workflow, one must use a FireWorker that contains an env that looks like ``{"my_vasp_cmd": "mpirun -n 24 vasp"}``. Here, the ``my_vasp_cmd`` in the dictionary matches the ``>>my_vasp_cmd<<`` string in the env_chk. Thus, when VASP is executed via this FireWorker, it will execute the command ``mpirun -n 24 vasp``. Other FireWorkers, for example located on different computing centers, might execute different VASP commands and can support this by setting a different value of the FireWorker ``env``. The workflow can be kept intact since the workflow is merely pointing to the ``my_vasp_cmd`` env variable and not setting the VASP command explicitly. There are more details about setting the FireWorker env variables in the FireWorks tutorials (in particular the Worker tutorial). The unit tests also use the env_chk feature to find the db configuration file. e.g., see the unit test: ``atomate.vasp.tests.test_vasp_workflows.TestVaspWorkflows#test_single_Vasp_dbinsertion`` and you will have a flavor for how this works. Just remember that if you see something like this ``>>db_file<<``, when running your Workflow your FireWorker will need to set the env like this: ``FWorker(env={"db_file": "path/to/db.json"})`` and you will need to use that FireWorker when launching the jobs.

CalcLocs
========

If you are running multiple VASP jobs that depend on copying the outputs of previous jobs, one issue is how to pass the directory information of previous VASP jobs from Firework to Firework. It is possible to do this manually (as was done in the MPWorks codebase), or using the ``pass_job_info`` keyword built into Fireworks, but the standard way to do this in atomate is *CalcLocs*. Procedurally, all you need to do is add the ```PassCalcLocs`` FireTask to every Firework that contains a VASP job (see ``atomate.vasp.fireworks.core`` for examples). Downstream jobs like ``CopyVaspOutput`` will have a ``calc_loc`` variable that can be set to True, and will automatically get the previous VASP dir parsed from before. Similar with ``VaspToDbTask``. Note that a couple of advantages of this system are:

* It is a general way of passing VASP directories that works with any Firework, and doesn't require you to code the logic of passing VASP directories inside of other functions (e.g., database insertion tasks as was done previously in MPWorks). Thus, the task of reporting and passing the VASP job location is well-separated from the other functions and can just be added in very easily. The only downside is that you have to remember to add in this FireTask.
* The CalcLocs maintains a running dictionary of job type to job location. If you need to grab outputs from multiple jobs (or say, from two jobs back), it is all supported within the framework. Just read the docs, e.g., of ``CopyVaspOutput``.
* Job directories are located across different machines and require ``scp`` or some other complex transfer mechanism are automatically handled by this infrastructure. You don't have to lift a finger! Just tell the parent Firework to pass the calcloc and the child firework to copy the vasp output (which supports the calcloc framework).

Workflow "Powerups"
===================

Workflow powerups are intended to be like function decorators, but for Workflows. For example, let's say you've built a multi-step workflow that computes a band structure. Now, you want to make sure that once a workflow starts running, it is prioritized to finish that particular workflow versus starting other workflows. By passing your workflow through a "powerup", you can get back a decorated workflow that sets the priorities of the Fireworks inside your workflow to endow this behavior (e.g., give all children Fireworks 2X the priority of the root parent). This particular powerup is located in ``atomate.vasp.vasp_powerups.add_priority``. Another powerups allows you to track the status of your jobs (last few lines in output files) in the FireWorks database, for example.

Note that another planned "powerup" is to endow Workflows with duplicate checking, i.e., to make sure the same structure is not run twice. In the past, such duplicate checking logic would be mixed in with the rest of the Workflow (about setting up VASP parameters, running VASP, etc.), and the end result was a very messy workflow code. It was also difficult to turn duplicate checking off and on as desired since all the logic was intermixed. By moving the duplicate checking to a "powerup", one can simply enable duplicate checking by passing the Workflow through the appropriate powerup.

See the ``vasp_powerups.py`` file for examples.

Tuning performance on different machines
========================================

VASP has certain INCAR parameters like NCORE, NPAR, KPAR, etc. that can be tuned
based on your machine. Since the ``ModifyIncar`` firetask supports
``env_chk``, these values can also be set in the fireworker config file
(my_fworker.yaml). E.g.,

.. code-block:: yaml

    env:
      incar_update:
        NCORE: 24

Note that NCORE sets the number of cores that work on a single orbital.
Typically, you want to set this between 1 (higher memory requirements) and
the number of cores per node (lower memory requirements while still
maintaining fast communication times between workers on an a single orbital).
A good starting point might be setting NCORE equal to the square root of
number of cores per node as per the VASP manual. The following information
might come in handy when setting the NCORE parameter on NERSC machines:

* Edison - 24 tasks per node
* Cori - 32 tasks per node
* Matgen - 16 tasks per node

Thus, a good starting point is to set NCORE=4 for Matgen/Edison and NCORE=8 for
Cori. Reduce NCORE if you want to try to increase speed at the risk of having
lower memory available per orbital.

====
Help
====

FAQ:
====

Q:What do I actually need to do to get a job running?
-----------------------------------------------------

:A: First, you need to install and configure atomate (see the :ref:`installation tutorial <installation tutorial>`) for your computing center of interest. Next you need to get some workflows. The easiest way is to throw a pymatgen Structure object into one of the prebuilt workflow functions in ``atomate/vasp/workflows/presets``. Et voilá! You have a workflow object. Next you need to put the workflow into your LaunchPad using the add_wf method in FireWorks. Finally, you need to run the workflow using FireWorks, e.g. using rlaunch, qlaunch or any of the other FireWorks tools.

    Basically, the goal of atomate is to help you get some workflows. e.g., you have a structure and you know you want the dielectric constant - atomate will help you get a workflow to accomplish that. All the details of running workflows, managing them, etc. is handled by FireWorks. Note that there is also an ``mmwf`` script that is intended to help you in putting a Workflow in the LaunchPad, but if you don't really understand what it's doing, it's probably best to ignore this for now.

Q: How do I know what workflows are available?
----------------------------------------------

:A: Browse the library folder in ``atomate/vasp/workflows/base`` for the raw workflows. Browse ``atomate/vasp/workflows/presets`` for preset workflows (just give a Structure, get back a workflow)

Q: I have a workflow that is almost what I want, but I want to tune some settings. How?
---------------------------------------------------------------------------------------

:A: Workflows are composed of Fireworks which are in turn composed of FireTasks. First look at code of the actual Fireworks that your workflow is referring to. Does the Firework contain a parameter for the setting that you want? If so, you can modify the workflow YAML file to set that parameter. If you are sure your Firework does not have the parameter you want, look at the FireTasks inside the Firework. Do those have a parameter for the setting that you want? If yes, the best option is to probably compose the Workflow in Python rather than YAML. It is generally *very* easy to do this. If you don't see the option anywhere, you will need to code it inside the FireTask/Firework.

Q: How do I create a brand new workflow?
----------------------------------------

:A: If you just want to rearrange, add, or delete Fireworks in one of the existing workflows, simply create a new YAML file that contains the sequence of steps you want.

    If the Fireworks that are currently implemented in atomate do not contain the function you want, you will need to write a new Firework (and maybe new FireTasks) and connect them into a workflow. Maybe try referring to how some of the existing workflows are constructed to learn how to do this.

Q: Are there any unit tests to make sure atomate is giving me sensible answers?
-------------------------------------------------------------------------------

:A: We are working on it...

Q: Is there a command line tool?
--------------------------------

:A: The ``atwf`` tool is there but somewhat under development. If you know what you are doing it is probably helpful, if you don't know what you are doing then using this tool probably will not lead to your success in running a workflow.

==============
Citing atomate
==============

We will write and publish a paper on atomate at a later point. For now, you can cite the following two works::

    (1) Jain, A.; Ong, S. P.; Chen, W.; Medasani, B.; Qu, X.; Kocher, M.;
    Brafman, M.; Petretto, G.; Rignanese, G.-M.; Hautier, G.; Gunter, D.;
    Persson, K. A. FireWorks: a dynamic workflow system designed for
    high-throughput applications, Concurr. Comput. Pract. Exp., 2015, 22,
    doi:10.1002/cpe.3505.

    (2) Ong, S. P.; Richards, W. D.; Jain, A.; Hautier, G.; Kocher, M.; Cholia,
    S.; Gunter, D.; Chevrier, V. L.; Persson, K. a.; Ceder, G. Python Materials
    Genomics (pymatgen): A robust, open-source python library for materials
    analysis, Comput. Mater. Sci., 2013, 68, 314–319,
    doi:10.1016/j.commatsci.2012.10.028.


.. _contributing-label:

====================================
Contributing / Contact / Bug Reports
====================================

Want to see something added or changed? There are many ways to make that a reality! Some ways to get involved are:

* Help us improve the documentation - tell us where you got 'stuck' and improve the install process for everyone.
* Let us know if you need support for a queueing system or certain features.
* Point us to areas of the code that are difficult to understand or use.
* Contribute code!

The list of contributors to atomate can be found :doc:`here </contributors>`.

atomate is currently in an alpha release. Although atomate is open source, currently no support is provided for atomate other than for those who are contributing to its development. There is an `atomate Google Group`_ dedicated to discussion and support related to development

.. _atomate Google Group: https://groups.google.com/forum/#!forum/atomate

=========
Changelog
=========

:doc:`Changelog </changelog>`

=======
License
=======

atomate is released under a modified BSD license; the full text can be found :doc:`here</license>`.

===========================
Comprehensive Documentation
===========================

Some comprehensive documentation is listed below (only for the brave!)

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
