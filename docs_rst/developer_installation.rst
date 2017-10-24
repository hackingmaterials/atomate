.. title:: Atomate developer installation
.. _developer installation:

======================================
Installing atomate in development mode
======================================

Introduction
============

This page documents the differences between installing atomate as a user and installing atomate (or other codes) to develop new functionality.
It is only intended as a supplement to the main installation guide and does not stand alone.


Installation
============

.. _codes-develop-mode:

Install Python codes in development mode
----------------------------------------

You can install any code you like in development mode. This would allow you to make changes directly to the source, so you can contribute back changes or new features to atomate (or pymatgen, FireWorks, etc.).

Note that if you previously installed atomate using another tool (e.g., pip or conda), you should uninstall that first before starting this installation. Otherwise, you might have conflicts or unintended behavior resulting from the two different code installations.

The steps for installing pymatgen and atomate in development mode are below.

1. Note: you should have activated your virtual environment or conda environment before proceeding.

#. Create a ``codes`` directory in ``<<INSTALL_DIR>>``

#. ``cd`` to your newly created ``<<INSTALL_DIR>>/codes`` directory.

#. Clone each of the packages you want to install in development mode using git. You don't have to know the details of how to use git for the installation, but if you are going to be developing code in Python, you should take a look at this `simple git introduction <http://rogerdudler.github.io/git-guide/>`_. Most Linux distributions include git, so you shouldn't have to install it on the cluster. To download the codes, use the following commands (one command per line)::

        git clone https://github.com/materialsproject/pymatgen.git
        git clone https://github.com/hackingmaterials/atomate.git

   Now you should have atomate and pymatgen directories in your ``codes`` directory.

#. For each of these folders, you ``cd`` into the folders and ``pip install -e .`` (the ``-e`` flag installs as editable, and the ``.`` simply means to install from the ``setup.py`` in the current directory) or use the ``conda`` equivalent. If you run into problems during this procedure, you may need manually install or load certain dependencies (e.g., use ``module load numpy`` on your system or ``pip install numpy`` or ``conda install numpy``). Note that once installed, if you make changes to the code in these directories, the changes will take effect immediately without needing to reinstall. Thus, you can both view and modify the code installed in these directories.

#. If you want to update these codes later on, execute ``git pull`` followed by ``pip install -e .`` (or again, the ``conda`` equivalent) in the relevant directory. Note that you should update codes in the same order as installation as listed above.


Post-installation
=================

Basic confirmation of installation
----------------------------------

Open up a Python shell using the command ``python``. Confirm that the commands ``import pymatgen`` and ``import atomate`` execute without any issues / errors. Remember that you will need to still be in your virtual environment!

(optional) Run unit tests
-------------------------

If you make changes to atomate, it is a good idea to rerun the unit tests to make sure everything is still working.
The ``db.json`` and ``my_launchpad.yaml`` in the ``<<INSTALL_DIR>>/codes/atomate/atomate/common/test_files`` directory control the database to use for the unit tests. The default is to use a MongoDB running on localhost. You can update these to whatever you like, e.g. a MongoDB instance in the cloud that you use for tests.

.. warning:: Although you can re-use the same Mongo host and port as your production installation for tests, do **not** also use the same database as your production runs! This is why the default configuration uses a database name with ``_unittest`` - so that it won't conflict with any production database. The database and LaunchPad you use in the unit tests **WILL** be reset frequently. **DO NOT USE YOUR PRODUCTION DATABASES FOR TESTING** or you will lose everything!

These unit tests are designed to run without installing VASP or other codes so you can run these on your local machine.
Some of them start with a VASP workflow but apply the ``use_fake_vasp`` method to replace calling the VASP executable with a "Faker" that verifies basic properties of the inputs and copies pre-stored output files to the current directory, thus simulating the execution of VASP.

Note that the unit tests in ``atomate/vasp/tests/test_vasp_workflows.py`` can be modified to actually run VASP by setting VASP_CMD to a String representing your VASP command.
This will just directly execute VASP where you are running the test (i.e., won't submit jobs to a queue).
If you need to debug at a later point, this might be something to refer back to, e.g., by running an interactive job on your supercomputing center (so that you are on a compute node) and trying this procedure.

Many tests have a DEBUG option that can sometimes help in finding problems.
Sometimes you need to toggle DEBUG on/off a couple of times if you are doing this to make sure all the old data is actually cleared between debug runs.

Run the tests by navigating to ``<<INSTALL_DIR>>/codes/atomate/`` and running ``python setup.py test`` or ``nosetests``.