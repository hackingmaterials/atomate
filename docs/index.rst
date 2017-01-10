.. title:: atomate (Materials Science Workflows)

.. image:: _static/atomate_logo_small.png
   :alt: atomate logo

.. pull-quote:: | "Civilization advances by extending the number of important operations which we can perform without thinking about them."
                |    - Alfred North Whitehead

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

Although atomate itself should be pip-installable, the actual process to get an infrastructure to run VASP at supercomputers might be a bit more involved (not difficult, just more steps). See "Detailed Installation Notes" at the end of this document for some notes on how to do this.

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

===========================
Detailed Installation Notes
===========================

Here are some notes on how to get atomate up and running in a production system at your supercomputing center. These notes are geared towards the NERSC supercomputing center. You'll need to fill in details and adapt accordingly for other centers.

A. Things you need to do once
=============================

Here are some things you will likely only need to do once (per machine) as an "initial install".

Preliminaries
-------------

1. Make sure you can access to a MongoDB installation from the compute nodes. i.e. you can either start and stop a Mongo server yourself or have credentials to a Mongo server that's always available. Also confirm there are no firewalls from your compute node to your Mongo server. If you are able to get through the FireWorks tutorials on running jobs through a queue, then this step is probably OK. If you are unsure, I recommend actually trying that first before going through all the atomate stuff.
2. Make sure you have access to the VASP executable and pseudopotential files. If you cannot run VASP manually, you cannot do it through this infrastructure. I recommend making sure you know how to run VASP manually on your supercomputer before embarking on this installation.

Set some environment variables
------------------------------

1. Make sure your ``PMG_VASP_PSP_DIR`` environment variable is set to point to your VASP pseudopotential directories (this is a pymatgen thing). Probably you want to put this in your ``.bash_profile`` (or equivalent, e.g., ``.bashrc.ext`` at NERSC) and never have to worry about this again. Otherwise, you will need to do this each and every time.

Install some codes
------------------

1. Load any modules that are needed to do a Python installation.

#. Create a directory in a spot on disk that has relatively fast access from compute nodes. Your Python codes and config files will go here. We will call this place ``<<INSTALL_DIR>>``.

#. It's probably best to make this directory a virtual environment, in case you want to have multiple environments later (for different projects, perhaps for different machines, etc). This will also help in avoiding permissions problems with installing Python codes. So create a virtualenv in the ``<<INSTALL_DIR>>`` using the ``virtualenv`` command. If you know what you are doing, you can probably make things work without virtualenv.

#. Activate your virtualenv, e.g. ``source <<INSTALL_DIR>>/bin/activate``. Now you are ready to install codes.

#. I would suggest making a subdirectory for codes, e.g. ``<<INSTALL_DIR>>/codes`` and then moving to that directory for the remainder.

#. Technically, you just need the atomate code which will contain all the dependencies, and you might be able to get by using the ``pip`` install. What I do is actually install the full source of the atomate code and all of its important dependencies inside ``<<INSTALL_DIR>>/codes``. This includes a ``git clone`` followed by a ``python setup.py develop`` for the following codes:

   - fireworks - ``git clone https://www.github.com/materialsproject/fireworks.git``
   - pymatgen - ``git clone https://www.github.com/materialsproject/pymatgen.git``
   - pymatgen-db - ``git clone https://www.github.com/materialsproject/pymatgen-db.git``
   - custodian - ``git clone https://www.github.com/materialsproject/custodian.git``
   - atomate (from hackingmaterials org) - ``git clone https://www.github.com/hackingmaterials/atomate.git``

#. If all the installation seemed to go smoothly, you are all set! You can try running some unit tests in the code to help confirm things. Note that some of the unit tests in some of the codes will require a MongoDB server.

Configure a bunch of things
---------------------------

In addition to having the code installed, you will need to configure a bunch of settings for running at your computing cluster. This includes setting up your queue adapter and submission script template, providing credentials to your databases, and setting locations of loggers and miscellaneous items.

1. Copy the contents of ``/atomate/vasp/config_example/standard_config`` to ``<<INSTALL_DIR>>/config``. We can work off these files to begin with rather than creating the files from scratch.

There is a lot to configure, so let's tackle the files one by one. We will start simple and get more complex.

Note that all variables enclosed in ``<<>>``, e.g. ``<<HOSTNAME>>``, must be modified by the user.

**my_launchpad.yaml**

As you should know, this file contains the configuration for the FireWorks database (LaunchPad). Make sure to set:

* ``<<HOSTNAME>>`` - the host of your FWS db server
* ``<<PORT_NUM>>`` - the port of your FWS db server
* ``<<DB_NAME>>`` - whatever you want to call your database. If you are not feeling creative, call it ``vasp_calcs``.
* ``<<ADMIN_USERNAME>>`` and ``<<ADMIN_PASSWORD>>`` - the (write) credentials to access your DB. Delete these lines if you do not have password protection in your DB.
* ``<<LOG_DIR>>`` - you can leave this to ``null``. If you want logging, put a directory name str here.
* The other settings, I've left to defaults. Feel free to modify them if you know what you are doing.

You can test whether your connection is running by running ``lpad -l my_launchpad.yaml reset``. This will reset and initialize your FireWorks database. Note that you might see some strange message about ``<<ECHO_STR>>``. We will fix that configuration later - feel free to ignore it for now.

**db.json**

This file contains credentials needed by the pymatgen-db code to insert the results of your VASP calculations. The easiest solution is to use the same database as your FireWorks database, but just use a different collection name. Or, you could use separate databases for FireWorks and VASP results. It is up to you.

For all settings, set to the same as the FireWorks database (``my_launchpad.yaml``) if you're keeping things simple. Or, use the settings for your dedicated database for VASP outputs. Note that since this is a JSON file, you need to use valid JSON conventions. e.g., wrap String values in quotes.

Once you've set up the credentials this file should be good to go.

**FW_config.yaml**

This file contains your global FireWorks settings. Later on (not now), you will set an environment variable called ``FW_CONFIG_FILE`` that points to this file. This file subsequently gives the directory name of where to find the other FWS-related files (my_launchpad.yaml, my_fworker.yaml, and my_qadapter.yaml). Anyway, in terms of setting up this file, set:

* ``<<PATH_TO_CONFIG_DIR>>`` - this is the **full** name of the directory containing the files ``my_launchpad.yaml``, ``my_fworker.yaml``, and ``my_qadapter.yaml``. The easiest way to set this variable is to navigate to ``<<INSTALL_DIR/config>>``, type ``pwd``, and paste the result into this variable.
* ``<<ECHO_TEST>>`` - the simplest thing is to delete this line. If you want, put an identifying string here. Whatever you put will be echoed back whenever you issue a FireWorks command. It is sometimes helpful if you are working with multiple databases and prefer a reminder of which database you are working with.

**my_fworker.yaml**

This file is both simple and complicated. The basic setup is simple. But, setting the ``env`` variable properly requires knowing about the details of the workflows you are going to run. Make sure you understand the ``env_chk`` framework (described elsewhere in the docs) to really know what is going on here.

* ``<<name>>`` - set to any name that describes this Worker. e.g. ``Generic NERSC``.
* ``<<env.db_file>>`` - many of the workflows implemented in atomate use the ``env_chk`` framework to get the path to the tasks database file from here. This allows setting different database files on different systems. Anyway, you want to put the **full** path of ``<<INSTALL_DIR>>/config/db.json``.
* ``<<env.vasp_cmd>>`` - many of the workflows implemented in atomate use the ``env_chk`` framework to get the actual command needed to run VASP because this command differs on different systems and cannot be hard-coded in the workflow itself. So put your full VASP command, e.g. ``mpirun -n 16 vasp`` here.
* ``<<env.scratch_dir>>`` - temporary place where to run VASP calculations using custodian framework. If set to the ``null`` it will simply use the current working directory without using a scratch_dir.

Note that all of these values might depend on the specific system you are running on. The point of the ``my_fworker.yaml`` is precisely to allow for different settings on different
systems. By having a different ``my_fworker.yaml`` file for each intended systems, you can tailor the execution of workflows across systems. This procedure is straightforward but is not covered here. You can set up a second config file, and point your ``FW_CONFIG_FILE`` environment variable to that second config file in order to use different settings (e.g., different ``my_fworker.yaml``).

**my_qadapter.yaml**

This file controls the format of your queue submission script and the commands to submit jobs to the queue (e.g., ``qsub`` versus ``sbatch``). I will not go over how to set this file here. Please refer to the FWS tutorials for that. Note that ``<<CONFIG_DIR>>`` should point to the **full** path of ``<<INSTALL_DIR>>/config``. One further note on this file is that the default uses ``singleshot`` in "reservation" (``-r``) mode. If you want to pack multiple Fireworks into a queue submission you might try turning off reservation mode, and using ``rapidfire`` mode with the appropriate options.

That's it! You've finished basic configuration!

B. Things you need to do each time you log in (or just once if you put it in your .bash_profile)
================================================================================================

In order to run jobs, you must:

1. Load modules for any important libraries (e.g., Python / VASP)
#. Activate your virtualenv (``source <<INSTALL_DIR>>/bin/activate``).
#. set your ``FW_CONFIG_FILE`` env variable to point to ``FW_config.yaml`` (``export FW_CONFIG_FILE=<<INSTALL_DIR>>/config/FW_config.yaml``).

You can put all of these things inside your ``.bash_profile`` or equivalent in order to make them automatic when you log into the cluster. It is up to you.

C. Running some jobs
====================

Ok, you are now ready to test running some jobs! Note that the testing procedure was recently changed and is under development. For now, try getting a workflow from the the ``vasp/workflows/preset`` package in which you can just give a Structure and get back a workflow. Then add that workflow to your LaunchPad and then use FireWorks to launch it in the manner of your desire.

D. Tuning performance on different machines
===========================================

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
FAQS
====

1. **What do I actually need to do to get a job running?**

First, you need to install and configure atomate (see the installation notes above) for your computing center of interest. Next you need to get some workflows. The easiest way is to throw a pymatgen Structure object into one of the prebuilt workflow functions in ``atomate/vasp/workflows/presets``. Et voilá! You have a workflow object. Next you need to put the workflow into your LaunchPad using the add_wf method in FireWorks. Finally, you need to run the workflow using FireWorks, e.g. using rlaunch, qlaunch or any of the other FireWorks tools.
Basically, the goal of atomate is to help you get some workflows. e.g., you have a structure and you know you want the dielectric constant - atomate will help you get a workflow to accomplish that. All the details of running workflows, managing them, etc. is handled by FireWorks. Note that there is also an ``mmwf`` script that is intended to help you in putting a Workflow in the LaunchPad, but if you don't really understand what it's doing, it's probably best to ignore this for now.

2. **How do I know what workflows are available?**

Browse the library folder in ``atomate/vasp/workflows/base`` for the raw workflows. Browse ``atomate/vasp/workflows/presets`` for preset workflows (just give a Structure, get back a workflow)

3. **I have a workflow that is almost what I want, but I want to tune some settings. How?**

Workflows are composed of Fireworks which are in turn composed of FireTasks. First look at code of the actual Fireworks that your workflow is referring to. Does the Firework contain a parameter for the setting that you want? If so, you can modify the workflow YAML file to set that parameter. If you are sure your Firework does not have the parameter you want, look at the FireTasks inside the Firework. Do those have a parameter for the setting that you want? If yes, the best option is to probably compose the Workflow in Python rather than YAML. It is generally *very* easy to do this. If you don't see the option anywhere, you will need to code it inside the FireTask/Firework.

4. **How do I create a brand new workflow?**

If you just want to rearrange, add, or delete Fireworks in one of the existing workflows, simply create a new YAML file that contains the sequence of steps you want.

If the Fireworks that are currently implemented in atomate do not contain the function you want, you will need to write a new Firework (and maybe new FireTasks) and connect them into a workflow. Maybe try referring to how some of the existing workflows are constructed to learn how to do this.

4. **Are there any unit tests to make sure atomate is giving me sensible answers?**

We are working on it...

5. **Is there a command line tool?**

The ``mmwf`` tool is there but somewhat under development. If you know what you are doing it is probably helpful, if you don't know what you are doing then using this tool probably will not lead to your success in running a workflow.

=================
Citing atomate
=================

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

atomate is currently in an alpha release. Although atomate is open source, currently no support is provided for atomate other than for those who are contributing to its development.

=========
Changelog
=========

.. toctree::
   :maxdepth: 1

   changelog

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
