.. title:: Installing atomate
.. _installation tutorial:

==================
Installing atomate
==================

Introduction
============

This guide will get you up and running in an environment for running high-throughput workflows with atomate. atomate is built on pymatgen, custodian, and FireWorks libraries to run materials science workflows. Briefly:

* pymatgen_ is used for creating input and analyzing output of materials science codes
* custodian_ runs VASP and performs error checking/handling and checkpointing
* FireWorks_ enables designing, managing and executing workflows.

Details about how atomate is designed, how these different pieces interact, and how to run and write your own workflows will be covered in later tutorials as well as an academic publication on atomate. For now, these topics will be covered here in enough depth to get you set up and to help you know where to troubleshoot if you are having problems.

It is assumed that you are comfortable with basic Linux shell commands and navigation. If not, `Linux Journey`_ and `Linux Command`_ breifly cover enough to get you started. It will also be helpful if you are familiar with Python, but it is not strictly required for installation.

.. _pymatgen: http://pymatgen.org
.. _custodian: https://materialsproject.github.io/custodian/
.. _FireWorks: https://pythonhosted.org/FireWorks/
.. _Linux Command: http://linuxcommand.org/lc3_learning_the_shell.php
.. _Linux Journey: https://linuxjourney.com/lesson/the-shell


Objectives
==========

* Install and configure atomate on a supercomputing cluster
* Validate the installation with a test workflow

Installation checklist
======================

Completing everything on this checklist should result in a fully functioning environment. Each item will be covered in depth, but this can be used to keep track of the big picture and help reinstall on other systems.

1. Prerequisites_
#. `Create a directory scaffold for atomate`_
#. `Create a Python virtual environment`_
#. `Install Python packages`_
#. `Configure FireWorks`_
#. `Configure pymatgen`_
#. `Run VASP unit tests`_
#. `Run a test workflow`_


Automated installer
===================

The `Phases Research Lab at Penn State`_ has developed an `automated installer`_ to install atomate with minimal user interaction. The installer simply scripts all of the actions given in this installation guide after the user configures their database (see the `Configure FireWorks`_ section). There are a select few preset systems that are handled automatically (include TACC's Stampede, NERSC's Edison and Cori) and otherwise all of the relevant settings can be tweaked in one script and installed. For instructions on the use of the `automated installer`_, see the README. **Disclaimer**: this installer comes with no guarantees or warranty and the authors are not responsible for any problems caused (see the LICENSE). If you run into problems caused by the installer, please open an issue on GitHub.

.. _Phases Research Lab at Penn State: http://www.phases.psu.edu
.. _automated installer: https://github.com/PhasesResearchLab/install-atomate


.. _Prerequisites:

Prerequisites
=============

Before you install, you need to make sure that your "worker" computer (where the simulations will be run, often a computing cluster) that will execute workflows can (i) run the base simulation packages (e.g., VASP, LAMMPs, FEFF, etc) and (ii) connect to a MongoDB database. For (i), make sure you have the appropriate licenses and compilation to run the simulation packages that are needed. For (ii), make sure your computing center doesn't have firewalls that prevent database access. Typically, academic computing clusters as well as systems with a MOM-node style architecture (e.g., NERSC) are OK. High-security government supercomputing centers often require custom treatment and modified execution patterns - some further details are provided later in this installation guide.


VASP
----

Note: you can skip this step if you do not plan to run VASP with atomate.

To get access to VASP on supercomputing resources typically requires that you are added to a user group on the system you work on after your license is verified. You will also need access to the pseudopotentials. For convenience, you might copy these to the same directory you will be installating atomate, but this is not required. Regardless of its location, the directory structure should look like:

::

    pseudopotentials
    ├── POT_GGA_PAW_PBE
    │   ├── POTCAR.Ac.gz
    │   ├── POTCAR.Ac_s.gz
    │   ├── POTCAR.Ag.gz
    │   └── ...
    ├── POT_GGA_PAW_PW91
    │   ├── POTCAR.Ac.gz
    │   ├── POTCAR.Ac_s.gz
    │   ├── POTCAR.Ag.gz
    │   └── ...
    └── POT_LDA_PAW
        ├── POTCAR.Ac.gz
        ├── POTCAR.Ac_s.gz
        ├── POTCAR.Ag.gz
        └── ...


MongoDB
-------

MongoDB_ is a NoSQL database that stores each database entry as a document, which is represented in the JSON format (the formatting is similar to a dictionary in Python). Atomate uses MongoDB to:

* store the workflows that you want to run as well as their state details (through FireWorks - required)
* to parse output files and create database of calculation results (strongly recommended and assumed by most default settings of workflows, but technically optional)

MongoDB must be running and available to accept connections whenever you are running workflows. Thus, it is strongly recommended that you have a server to run MongoDB or (simpler) use a hosting service. Your options are:

* use a commercial service to host your MongoDB instance. These are typically the easiest to use and offer high quality service but require payment for larger databases. mLab_ offers free 500 MB databases with payment required for larger databases or support/backup. The free service is certainly enough to get started for small to medium size projects.
* contact your supercomputing center to see if they offer MongoDB hosting (e.g., NERSC has this)
* self-host a MongoDB server

If you are just starting, we suggest the first (with a free plan) or second option (if available to you). The third option will require you to open up network settings to accept outside connections properly which can sometimes be tricky.

Next, create a new database and set up two new username/password combinations:

- an admin user
- a read-only user

You might choose to have *two* separate databases - one for the workflows and another for the results. We suggest that new users keep both sets of results in a single database and only consider using two databases if they run into specific problems.

 Hang on to your credentials - we will configure FireWorks to connect to them in a later step. Also make sure you know the hostname and port for the MongoDB instance.

.. warning::

    The computers that perform the calculations must have access to your MongoDB server. Some computing resources have firewalls blocking connections. Note that this is not a problem for most computing centers that allow such connections (particularly from MOM-style nodes, e.g. at NERSC, SDSC, etc.), but some of the more security-sensitive centers (e.g., LLNL, PNNL) will run into issues. If you run into connection issues later in this tutorial, some options are:

  * contact your computing center to review their security policy to allow connections from your MongoDB server (best resolution)
  * set up an ssh tunnel to forward connections from allowed machines (the tunnel must be kept alive at all times you are running workflows)
  * use `FireWorks offline mode`_, which is a workable solution but makes the system more difficult to use and limits some features of FireWorks.


.. _MongoDB: https://docs.mongodb.com/manual/
.. _mLab: https://mlab.com
.. _FireWorks offline mode: https://pythonhosted.org/FireWorks/offline_tutorial.html


.. _Create a directory scaffold for atomate:

Create a directory scaffold for atomate
=======================================

Installing atomate includes installation of codes, configuration files, and various binaries and libraries. Thus, it is useful to create a directory structure to hold all these items.

1. Log in to the compute cluster and make sure the Python module you want to use is loaded and added to your rc file (e.g. ``~/.bashrc`` or ``~/.bashrc.ext`` at NERSC)

#. Create a directory in a spot on disk that has relatively fast access from compute nodes *and* that is only accessible by yourself or your collaborators. Your Python codes and config files will go here, including database credentials. We will call this place ``<<INSTALL_DIR>>``. A good name might simply be ``atomate``, but you could also use a project-specific name (e.g., ``atomate-solar``).

#. Now you should scaffold the rest of your ``<<INSTALL_DIR>>`` for the things we are going to do next. Create a directories named ``codes``, ``logs``, and ``config`` so your directory structure looks like:

    ::

        atomate
        ├── codes
        ├── config
        └── logs


.. _Create a Python virtual environment:

Create a Python virtual environment
===================================

We highly recommended that you organize your installation of the atomate and the other Python codes using a virtual environment (e.g. ``virtualenv`` or similar tool such as anaconda). Ultimately, whether you want to use a virtual environment is optional and you don't have to use one if you know what you are doing. Virtual environments allow you to keep an installation of Python and all of the installed packages separate from the installation on the system. Some of the main benefits are:


 * Different Python projects that have conflicting packages can coexist on the same machine.
 * Different versions of Python can exist on the same machine and be managed more easily (e.g. Python 2 and Python 3).
 * You have full rights and control over the environment. If it breaks, you can just delete the folder containing the environment and recreate it. On computing resources, this solves permissions issues with installing and modifying packages.

The easiest way to get a Python virtual environment is to use the ``virtualenv`` tool. Most Python distributions come with ``virtualenv``, but some clusters (e.g., NERSC) have moved towards using Anaconda_, which is a popular distribution of Python designed for scientific computing. If the compute resource you want to access is using Anaconda, you will follow the same general steps, but create your environment with ``conda create``. See the `documentation for the conda command line tool here`_ as well as a `conversion between virtualenv and conda commands <https://conda.io/docs/_downloads/conda-pip-virtualenv-translator.html>`_. To set up your virtual environment:


#. Go to your install directory (``<<INSTALL_DIR>>``) and create a virtual environment there. A good name might be ``atomate_env``. The command to create the environment would be ``virtualenv atomate_env``, which creates a folder ``atomate_env`` in the directory you are in.

#. You can ``ls`` this directory and see that you have the following structure:

    ::

        atomate
        ├──atomate_env/
           ├── bin
           ├── include
           ├── lib
           ├── lib64
           └── pip-selfcheck.json
        ├── codes
        ├── config
        └── logs

#. If you look in the ``bin`` directory, you will see several programs, such as activate, pip, and Python itself. ``lib`` will be where all of your installed packages will be kept, etc. Again, if anything goes wrong in installing Python codes, you can just delete the virtual environment directory (``atomate_env``) and start again.

#. Activate your environment by running ``source <<INSTALL_DIR>>/atomate_env/activate``. This makes it so when you use the command ``python``, the version of ``python`` that you use will be the one in the  ``bin`` directory rather than the system-wide Python. You can read the activation script if you are interested. It's just does a little magic to adjust your path to point towards the ``bin`` and other directories you created.

.. _Anaconda: https://www.continuum.io
.. _documentation for the conda command line tool here: https://conda.io/docs/using/envs.html


.. _Install Python packages:

Install Python packages
=======================

You have successfully set up an environment in which to install atomate! Next, we will download and install all of the atomate-related Python packages. The main tool for installing Python packages is pip and we will use this to install packages (unless you have an Anaconda distribution where again, you'd use conda_). Technically, you could simply use pip to ``pip install atomate`` and pull in atomate and all of the requirements from PyPI_. In that case, you would *not* need to follow the rest of this section's instructions.

However, our recommendation (and the procedure below) is to install directly from GitHub so you can always have the most recent codebase and be able to more easily inspect and modify the code. We'll also do this for the main dependencies of atomate because they often change and evolve together in the source, but not be released to PyPI. Note that this method of installation is required if you will be developing in atomate or any of the other software mentioned here.

1. ``cd`` to your newly created ``codes`` directory.

#. Clone each of the following packages from GitHub using git. You don't have to know the details of how to use git for the installation, but if you are going to be developing code in Python, you should take a look at this `simple git introduction`_. Most Linux distributions include git, so you shouldn't have to install it on the cluster. To download the codes, use the following commands (one command per line)::

        git clone https://www.github.com/materialsproject/fireworks.git
        git clone https://www.github.com/materialsproject/pymatgen.git
        git clone https://www.github.com/atztogo/phonopy.git
        git clone https://www.github.com/materialsvirtuallab/pymatgen-diffusion.git
        git clone https://www.github.com/materialsproject/pymatgen-db.git
        git clone https://www.github.com/materialsproject/custodian.git
        git clone https://www.github.com/hackingmaterials/atomate.git

  Now you should have atomate, custodian, FireWorks, phonopy, pymatgen, pymatgen-db and pymatgen-diffusion folders in your ``codes`` directory.

#. **It is important that you follow the order listed above to install packages**. For example, you must install atomate last so that it uses your custom installation of pymatgen rather than pulling in the requirements as an external library from PyPI. The ``-e`` flag installs as editable. For each of these folders, you ``cd`` into the folders and ``pip install -e .`` (the ``.`` simply means to install from the ``setup.py`` in the current directory) or use the ``conda`` equivalent of this command. There are several clever ways to do this in a one line command as a loop which you can use as an exercise of your shell capabilities [#]_. Note that once installed, if you make changes to the code in these directories, the changes will impact immedately without needing to reinstall. Thus, you can view and modify the code installed in these directories.


.. _conda: https://conda.io/docs/using/pkgs.html
.. _PyPI: https://pypi.python.org/pypi
.. _simple git introduction: http://rogerdudler.github.io/git-guide/

.. _Configure FireWorks:

Configure FireWorks
===================

We've now set up your environment and installed the necessary codes. You're well on your way!

The next step is to configure some the codes for your specific system - e.g., your MongoDB credentials, your computing cluster and its queuing system, etc. - starting with FireWorks. The setup below will be just enough to get your environment bootstrapped. For more details on the installation and specifics of FireWorks, read the `installation guide`_.

.. note:: All of the paths here must be *absolute paths*. For example, the absolute path that refers to ``<<INSTALL_DIR>>`` might be ``/global/homes/u/username/atomate`` (don't use the relative directory ``~/atomate``).

.. warning::

    **Passwords will be stored in plain text!** These files should be stored in a place that is not accessible by unauthorized users. Also, you should make random passwords that are unique only to these databases.

Create the following files in ``<<INSTALL_DIR>>/config``.

db.json
-------

The ``db.json`` file tells atomate where to put the results of parsing calculations from your workflows (i.e., actual property output data on materials). The ``db.json`` file requires you to enter the basic database information as well as what to call the main collection that results are kept in (e.g. ``tasks``) and the authentication information for an admin user and a read only user on the database. Mind that valid JSON requires double quotes around each of the string entries and that all of the entries should be strings except the value of "port", which should be an integer.

.. code-block:: json

    {
        "host": "<<HOSTNAME>>",
        "port": <<PORT>>,
        "database": "<<DB_NAME>>",
        "collection": "tasks",
        "admin_user": "<<ADMIN_USERNAME>>",
        "admin_password": "<<ADMIN_PASSWORD>>",
        "readonly_user": "<<READ_ONLY_PASSWORD>>",
        "readonly_password": "<<READ_ONLY_PASSWORD>>",
        "aliases": {}
    }

my_fworker.yaml
---------------

In FireWorks' distributed `server-worker model`_, each computing resource where you run jobs is a FireWorker (Worker). ``my_fworker.yaml`` controls the environment and settings unique to the cluster, such as the VASP executable. If this is the only cluster you plan on using just one Worker for all of your calculations a minimal setup for the ``my_fworker.yaml`` file is

.. code-block:: yaml

    name: <<WORKER_NAME>>
    category: ''
    query: '{}'
    env:
        db_file: <<INSTALL_DIR>>/config/db.json
        vasp_cmd: <<VASP_CMD>>

Where the <<WORKER_NAME>> is arbitrary and is useful for keeping track of which Worker is running your jobs (an example might be ``Edison`` if you are running on NERSC's Edison resource). ``db_file`` points to the ``db.json`` file that you just configured and contains credentials to connect to the calculation output database. The <<VASP_CMD>> is the command that you would use to run VASP with parallelization (``srun -n 16 vasp``, ``ibrun -n 16 vasp``, ``mpirun -n 16 vasp``, ...). If you don't know which of these to use or which VASP executable is correct, check the documentation for the computing resource you are running on or try to find them interactively by checking the output of ``which srun``, ``which vasp_std``, etc. . If you later want to set up multiple Workers on the same or different machines, you can find information about controlling which Worker can run which job by using the ``name`` field above, or the ``category`` or ``query`` fields that we did not define. For more information on configuring multiple Workers, see the `FireWorks documentation for controlling Workers`_. Such features allow you to use different settings (e.g., different VASP command such as different parallelization amount) for different types of calculations on the same machine or control what jobs are run on various computing centers.

my_launchpad.yaml
-----------------

The ``db.json`` file contained the information to connect to MongoDB for the calculation output database. We must also configure the database for storing and managing workflows within FireWorks using ``my_launchpad.yaml`` as in FireWorks' `server-worker model`_. The LaunchPad is where all of the FireWorks and Workflows are stored. Each Worker can query this database for the status of Fireworks and pull down Fireworks to reserve them in the queue and run them. A ``my_launchpad.yaml`` file with fairly verbose logging (``strm_lvl: INFO``) is below:

.. code-block:: yaml

    host: <<HOSTNAME>>
    port: <<PORT>>
    name: <<DB_NAME>>
    username: <<ADMIN_USERNAME>>
    password: <<ADMIN_PASSWORD>>
    ssl_ca_file: null
    strm_lvl: INFO
    user_indices: []
    wf_user_indices: []

Here's what you'll need to fill out:

* ``<<HOSTNAME>>`` - the host of your MongoDB db server
* ``<<PORT>>`` - the port of your MongoDB db server
* ``<<DB_NAME>>`` - the name of the MongoDB database
* ``<<ADMIN_USERNAME>>`` and ``<<ADMIN_PASSWORD>>`` - the (write) credentials to access your DB. Delete these lines if you do not have password protection in your DB (although you should).

**Note**: If you prefer to use the same database for FireWorks and calculation outputs, these values will largely be duplicated with ``db.json``. If you prefer to use different databases for workflows and calculation outputs, the information here will be different than ``db.json``.

my_qadapter.yaml
----------------

To run your VASP jobs at scale across one or more nodes, you usually submit your jobs through a queue system on the computing resources. FireWorks handles communicating with some of the common queue systems automatically. As usual, only the basic configuration options will be discussed. If you will use atomate as in this tutorial, this basic configuration is sufficient. A minimal ``my_qadapter.yaml`` file for SLURM machines might look like

.. code-block:: yaml

    _fw_name: CommonAdapter
    _fw_q_type: SLURM
    rocket_launch: rlaunch -c <<INSTALL_DIR>>/config singleshot
    nodes: 2
    walltime: 24:00:00
    queue: null
    account: null
    job_name: null
    pre_rocket: null
    post_rocket: null
    logdir: <<INSTALL_DIR>>/logs

The ``_fw_name: CommonAdapter`` means that the queue is one of the built in queue systems and ``_fw_q_type: SLURM`` indicates that the SLURM system will be used. FireWorks supports the following queue systems out of the box:

* PBS/Torque
* SLURM
* SGE
* IBM LoadLeveler

.. note::

  If you aren't sure what queue system the cluster you are setting up uses, consult the documentation for that resource. If the queue system isn't one of these preconfigured ones, consult the `FireWorks documentation for writing queue adapters`_. The FireWorks documentation also has tutorials on setting up your jobs to run on a queue in a way that is more interactive than the minimal details specified here.

``nodes``, ``walltime`` are the default reservations made to the queue as you would expect. ``queue`` refers to the name of the queue you will submit to. Some clusters support this and appropriate values might be ``regular``, ``normal``, ``knl``, etc. as defined by the compute resource you are using. The ``account`` option refers to which account to charge. Again, whether or not you need to set this depends on the resource. ``pre_rocket`` and ``post_rocket`` add lines to before and after you job launches in your queue submission script. One use of this would be to enter directives such as ``#SBATCH -C knl,quad,cache`` to configure SLURM to run on knl nodes.

.. _FireWorks documentation for writing queue adapters: https://pythonhosted.org/FireWorks/qadapter_programming.html?highlight=qadapter


FW_config.yaml
--------------

The ``FW_config.yaml`` file controls different FireWorks settings and also can point to the location of other configuration files. For a more complete reference to the FireWorks parameters you can control see the `FireWorks documentation for modifying the FW config`_. Here you simply need to accomplish telling FireWorks the location of the ``my_launchpad.yaml``, ``my_qadapter.yaml`` and ``my_fworker.yaml`` configuration files.

Create a file called ``FW_config.yaml`` in ``<<INSTALL_DIR>>/config`` with the following contents:

.. code-block:: yaml

    CONFIG_FILE_DIR: <<INSTALL_DIR>>/config

Finishing up
------------

The directory structure of ``<<INSTALL_DIR>>/config`` should now look like

::

    config
    ├── db.json
    ├── FW_config.yaml
    ├── my_fworker.yaml
    ├── my_launchpad.yaml
    └── my_qadapter.yaml

The last thing we need to do to configure FireWorks is add the following line to your RC / bash_profile file to set an environment variable telling FireWorks where to find the ``FW_config.yaml``

.. code-block:: bash

    export FW_CONFIG_FILE=<<INSTALL_DIR>>/config/FW_config.yaml

where ``<<INSTALL_DIR>>`` is your (usual) installation directory.

That's it. You're done configuring FireWorks and most of atomate. You should now perform a check to make sure that you can connect to the database by sourcing your RC file (to set this environment variable) and initializing the database by running the command

.. code-block:: bash

    lpad reset

which should return something like:

.. code-block:: bash

    Are you sure? This will RESET 0 workflows and all data. (Y/N) y
    2015-12-30 18:00:00,000 INFO Performing db tune-up
    2015-12-30 18:00:00,000 INFO LaunchPad was RESET.


.. _installation guide: http://pythonhosted.org/FireWorks/installation.html
.. _server-worker model: https://pythonhosted.org/FireWorks/index.html#centralized-server-and-worker-model
.. _FireWorks documentation for controlling Workers: https://pythonhosted.org/FireWorks/controlworker.html?highlight=category
.. _FireWorks documentation for modifying the FW config: https://pythonhosted.org/FireWorks/config_tutorial.html


.. _Configure pymatgen:

Configure pymatgen
==================

If you are planning to run VASP, the last configuration step is to configure pymatgen to (required) find the pseudopotentials for VASP and (optional) set up your API key from the `Materials Project`_.

1. The pseudopotentials should be in any folder (such as ``<<INSTALL_DIR>>/pps``) as in the `Prerequisites`_.

#. You can get an API key from the `Materials Project`_ by logging in and going to your `Dashboard`_. Enter these into a file called ``.pmgrc.yaml`` in your home folder (i.e., ``~/.pmgrc.yaml``) with the following contents

.. code-block:: yaml

    PMG_VASP_PSP_DIR: <<INSTALL_DIR>>/pps
    PMG_MAPI_KEY: <<YOUR_API_KEY>>

#. If you'd like to use a non-default functional in all of your calculations, you can set the ``DEFAULT_FUNCTIONAL`` key to a functional that is `supported by VASP`_, e.g. ``PS`` to use PBESol.

.. _Materials Project: https://materialsproject.org/dashboard
.. _Dashboard: https://materialsproject.org/dashboard
.. _supported by VASP: https://cms.mpi.univie.ac.at/vasp/vasp/GGA_tag.html


.. _Run VASP unit tests:

Run VASP unit tests
===================

To test atomate and the VASP functionality, make sure you have MongoDB running when executing the tests. You can do this by navigating to ``<<INSTALL_DIR>>/codes/atomate/`` and running ``python setup.py test`` or ``nosetests``.

These unit tests are designed to run without installing VASP. Some of them start with a VASP workflow but apply the ``use_fake_vasp`` method to replace calling the VASP executable with a "Faker" that verifies basic properties of the inputs and copies pre-stored output files to the current directory, thus simulating the execution of VASP. Anyway, this will help make sure your installation is in good shape.

Many tests have a DEBUG option that can sometimes help in finding problems. Sometimes you need to toggle DEBUG on/off a couple of times if you are doing this to make sure all the old data is actually cleared between debug runs.

Note that the unit tests in atomate/vasp/tests/test_vasp_workflows.py can be modified to actually run VASP by setting VASP_CMD to a String representing your VASP command. If you need to debug at a later point, this might be something to refer back to.

If everything looks OK, you are ready to run a test workflow!

.. _Run a test workflow:

Run a test workflow
===================

To make sure that everything is set up correctly and in place, we'll finally run a simple test workflow. In general, two ways to create workflows is using atomate's command line utility ``atwf`` or by creating workflows in Python. More discussion on constructing and running workflows can be found in the :ref:`running workflows tutorial` and details on writing custom workflows can be found in the :ref:`creating workflows`. For now, we will use ``atwf`` to construct a workflow. Ideally you set up a Materials Project API key in the `Configure pymatgen`_ section, otherwise you will need to provide a POSCAR for the structure you want to run. If you have an API key configured, you can run the following to run a structure optimization on Si

.. code-block:: bash

    atwf add -l vasp -s optimize_only.yaml -m mp-149

Alternatively, if you did not set up your API key or want to use a custom POSCAR instead the following command will accomplish the same

.. code-block:: bash

    atwf add -l vasp -s optimize_only.yaml POSCAR

These commands added workflows for running a single structure optimization FireWork to your LaunchPad. You can verify that by using FireWorks' ``lpad`` utility:

.. code-block:: bash

    lpad get_wflows

which should return:

.. code-block:: bash

    [
        {
            "state": "READY",
            "name": "Si--1",
            "created_on": "2015-12-30T18:00:00.000000",
            "states_list": "REA"
        },
    ]

To launch this FireWork and place a reservation in the queue, go to the directory where you would like your calculations to run (e.g. your scratch or work directories) and run the command

.. code-block:: bash

    qlaunch -r rapidfire

**Note**: If you want to run directly rather than through a queue, use ``rlaunch rapidfire`` (go through the FireWorks documentation to understand the details).

If all went well, you can check that the FireWork is in the queue by using the commands for your queue system (e.g. ``squeue`` or ``qstat``) or by checking that the state of the FireWork has changed from ``READY`` to ``RESERVED`` with ``lpad get_wflows``. Once this FireWorks is launched and is completed, you can use pymatgen-db to check that it was entered into your results database by running

.. code-block:: bash

    mgdb query -c <<INSTALL_DIR>>/config/db.json --props task_id formula_pretty output.energy_per_atom

This time, ``<<INSTALL_DIR>>`` can be relative. You should have seen the energy per atom you calculated for Si.

That's it! You've completed the installation tutorial!

See the following pages for more information on the topics we covered here:

* For submitting jobs to the queue in reservation mode see the `FireWorks advanced queue submission tutorial`_
* For using pymatgen-db to query your database see the `pymatgen-db documentation`_
* To see how to run and customize the existing Workflows and FireWorks try the :ref:`running workflows tutorial`

.. _FireWorks advanced queue submission tutorial: https://pythonhosted.org/FireWorks/queue_tutorial_pt2.html
.. _pymatgen-db documentation: https://materialsproject.github.io/pymatgen-db/

Troubleshooting and FAQ:
========================

Q: I can't connect to my LaunchPad database
-------------------------------------------

:A: Make sure the right LaunchPad file is getting selected

  Adding the following line to your ``FW_config.yaml`` will cause the line to be printed every time that configuration is selected

  ::

    ECHO_TEST: Database at <<INSTALL_DIR>>/config/FW_config.yaml is getting selected.

  Then running ``lpad version`` should give the following result if that configuration file is being chosen

  ::

    $ lpad version

    Database at <<INSTALL_DIR>>/config/FW_config.yaml is getting selected.
    FireWorks version: x.y.z
    located in: <<INSTALL_DIR>>/codes/fireworks

  If it's not being found, check that ``echo $FW_CONFIG_FILE`` returns the location of that file (you could use ``cat $FW_CONFIG_FILE`` to check the contents)

:A: Double check all of the configuration settings in ``my_launchpad.yaml``

:A: Have you had success connecting before? Is there a firewall blocking your connection?

:A: You can try following the tutorials of FireWorks which will go through this process in a little more detail.


Q: My job fizzled!
------------------

:A: Check the ``*_structure_optimization.out`` and ``*_structure_optimization.error`` in the launch directory for any errors. Also check the ``FW.json`` to check for a Python traceback.


Q: I made a mistake, how do I cancel my job?
--------------------------------------------

:A: One drawback of using the reservation mode (the ``-r`` in ``qlaunch -r rapidfire``) is that you have to cancel your job in two places: the queue and the LaunchPad. To cancel the job in the queue, use whatever command you usually would (e.g. ``scancel`` or ``qdel``). To cancel or rerun the FireWork, run

    .. code-block:: bash

        lpad defuse_fws -i 1

    or

    .. code-block:: bash

        lpad rerun_fws -i 1

    where `-i 1` means to make perfom the operations on the FireWork at index 1. Run ``lpad -h`` to see all of the options.

The non-reservation mode for qlaunching requires a little less maintenance with certain tradeoffs, which are detailed in the FireWorks documentation.

Q: I honestly tried everything I can to solve my problem. I still need help!
----------------------------------------------------------------------------

:A: There is a Google group for atomate: https://groups.google.com/forum/#!forum/atomate

Footnotes
=========

.. [#] ``for D in */; do cd D && pip install -e . && cd .. ; done``
