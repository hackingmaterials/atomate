.. title:: Running workflows tutorial
.. _running workflows tutorial:

==========================
Running Workflows Tutorial
==========================

Introduction
============

Once you have a working installation of atomate, you'll want to jump in and start running workflows. There are preset workflows with reasonable settings for many calculations. This tutorial will quickly guide you through customizing and running a preset workflow.


Objectives
==========

* Run a workflow from a YAML file and in Python
* Query for the calculation result in your database
* Visualize the results with matplotlib


Prerequisites
=============

In order for you to complete this tutorial you need

    * A working installation of atomate
    * A text editor

Equation of State Workflows
===========================

An equation of state (EOS) relates state variables to each other. In first-principles calculations at 0K, it is often useful to understand the how the energy of a crystal changes when the crystal volume changes. Here there are instructions for running a workflow to calculate the EOS for fcc Al using both YAML format (text) files and Python programming. Both interfaces can produce equivalent results, so you can use either or both depending on the situation. To calculate the EOS, we need to (i) optimize our ground state structure, (ii) systematically change the crystal volume and calculate the resulting energies, and finally (iii) fit these energies and volumes to one of several EOS functional forms.

.. figure:: _static/elastic_tensor.png
    :alt: Elastic Tensor Workflow
    :scale: 50%

    Atomate elastic tensor workflow. This uses the same pattern as the EOS workflow, but the analysis step (final Firework) calculates the elastic tensor instead of an EOS.

Running an Equation of State Workflow
=====================================

Setup
-----

Make sure you have completed the installation tutorial. Next, create a folder on your HPC resource for this tutorial. It can be located anywhere. You'll keep all of the files for this tutorial there.

If you do not already have important prior results in your atomate environment, start fresh by running::

``lpad reset``

.. warning:: This will reset all of your previous Fireworks and Workflows in your LaunchPad. Do not do this if you have actual results that you want to keep!

If you do not or want to reset your LaunchPad, you can set up a different database to use for tutorials or simply accept mixing your previous results and workflows with the contents of this tutorial in your database.

Create the structure
--------------------

In your text editor, create a file called ``POSCAR`` and enter the following text for a conventional fcc Al structure:

::

    Al4
    1.0
            4.0389294624         0.0000000000         0.0000000000
            0.0000000000         4.0389294624         0.0000000000
            0.0000000000         0.0000000000         4.0389294624
    Al
    4
    Direct
         0.000000000         0.000000000         0.000000000
         0.000000000         0.500000000         0.500000000
         0.500000000         0.000000000         0.500000000
         0.500000000         0.500000000         0.000000000


Note that this does not have to be a POSCAR. You can also supply multiple formats supported by pymatgen, including: Crystallographic Information File (CIF), POSCAR/CONTCAR, CHGCAR, vasprun.xml, CSSR, Netcdf and pymatgen's JSON serialized structures.

There are multiple ways of defining the workflow to execute on the structure. We'll go over 3 options.

Option 1: Use a pre-defined workflow in atomate
-----------------------------------------------

Atomate includes many pre-defined workflows - some defined as YAML files and others as Python functions. Pre-defined YAML files are in ``atomate.vasp.workflows.base.library``. You can add them using the command::

    atwf add [[STRUCTURE_FILE]] -l vasp -s [[NAME_OF_FILE]]

where ``[[STRUCTURE_FILE]]`` is your POSCAR file (or cif, etc.) and ``[[NAME_OF_FILE]]`` is the name of the workflow file. Note that the ``eos`` workflow is not one of the files (an example of a valid file is ``bandstructure.yaml``), so we'll skip this option for this example.

The Python function presets are defined in ``atomate.vasp.workflows.presets``. You can add such workflows using the command::

    atwf add [[STRUCTURE_FILE]] -p [[NAME_OF_PYTHON_FUNCTION]]

An example of a valid Python functions in ``atomate.vasp.workflows.presets`` is ``wf_bulk_modulus``, which is very similar to the workflow we'll actually run. However, this isn't exactly the procedure we'll follow next, so to continue with this example, let's choose one of the other two options for defining the workflow which are more custom than the presets. However, note that you can and are encouraged to use preset workflows where practical and that there exist many such workflows for you to choose from.

Option 2: Create your own workflow file
---------------------------------------

**Define the workflow file**

You can use a text editor to define your own workflow that chains together pre-defined steps in atomate. To get a feeling for this procedure, in your text editor, create a file called ``eos.yaml`` and enter the following text:

.. note:: If your VASP command is anything other than ``vasp_std``, then you'll want to change the last line (common_params). *Don't use any parallelization software (ibrun, srun, etc.)* it's bad behavior and probably won't work.


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
      vasp_cmd: vasp_std
      db_file: >>db_file<<

.. note::
    The YAML file format is typically considered easy to read, but if you want to know more about the YAML format in general you might want to take a look at the `detailed YAML specification`_. If you want to know more specifically about atomate's YAML specification, see the :ref:`Workflow YAML reference`.

.. _detailed YAML specification: http://www.yaml.org/spec/1.2/spec.html


**Add workflow to LaunchPad**

Within the folder containing your ``POSCAR`` (or other structure file) and ``eos.yaml``, run the following command to add the workflow to your LaunchPad:

.. code-block:: bash

    atwf add POSCAR -s eos.yaml

Unless you also want to try making a Python workflow and add it to your LaunchPad, skip ahead to the `Running the workflow`_ section.


Option 3: use Python to generate and add the workflow
-----------------------------------------------------

The YAML version above is more efficient and clear to read and modify than a typical shell script to set up and run these calculations by hand. Even so, this workflow would have been tedious to type out rather than copy-paste. `There must be a better way!`_ Enter Python.

In the installation tutorial, you set up your ``FW_config.yaml``, you indicated the atomate Fireworks can be found at :py:mod:`atomate.vasp.fireworks`. Similarly, atomate preset workflows can be imported from :py:mod:`atomate.vasp.workflows.presets.core`, which thinly wraps the base workflows (:py:mod:`atomate.vasp.workflows.base`) allowing for common settings to be changed with configuration dictionaries. The bulk modulus preset workflow does what the YAML file above does for us. And we can setup the workflow and add it to our LaunchPad ready to run in just a few lines of Python.


.. _There must be a better way!: https://www.youtube.com/watch?v=wf-BqAjZb8M

**Create the workflow script**

In the same directory as the POSCAR, create a Python script name ``eos.py`` with the following contents:

.. note:: If your VASP command is anything other than ``vasp_std``, then you'll want to change the line setting the ``VASP_CMD`` key of the configuration dictionary. *Don't use any parallelization software (ibrun, srun, etc.)* it's bad behavior and probably won't work.

.. code-block:: python

    # Create an EOS from the workflow from the atomate presets
    import numpy as np
    from pymatgen import Structure
    from fireworks import LaunchPad
    from atomate.vasp.workflows.presets.core import wf_bulk_modulus
    from atomate.vasp.powerups import add_modify_incar

    # load structure from file
    struct = Structure.from_file('POSCAR')  # note: many file formats supported, see function docs

    # set up configuration dictionary
    c = {}
    # 7 deformations +/- 10% of the equilibrium volume
    # note that the 1/3 power is so that we scale each direction by (x+1)^(1/3) and the total volume by (x+1)
    c["deformations"] =  [(np.identity(3)*(1+x)**(1.0/3.0)).tolist() for x in np.linspace(-0.1, 0.1, 5)]
    c["VASP_CMD"] = 'vasp_std'

    # create the Workflow
    wf = wf_bulk_modulus(struct, c)

    # now we need to set the correct smearing for the optimization, using the add_modify_incar powerup
    wf = add_modify_incar(wf, {'incar_update': {'SIGMA': 0.2, 'ISMEAR': 1}}, fw_name_constraint='optimization')

    # finally, instatiate the LaunchPad and add the workflow to it
    lpad = LaunchPad.auto_load() # loads this based on the FireWorks configuration
    lpad.add_wf(wf)


**Add workflow to LaunchPad**

If you want to add the workflow to your LaunchPad (e.g., you didn't already go through Option 2 for adding a workflow): from the folder with your ``POSCAR`` and ``eos.py``, run the Python script:

.. code-block:: bash

    python eos.py

.. _Running the workflow:

Running the workflow
--------------------

In both cases, we manually set our ``VASP_CMD`` key to be the plain VASP command for the resource you are using. The reason we did this is because the simulation of EOS for Al is relatively simple, so we can easily run this entire workflow in a couple minutes on a single core without the queue. This is not generally good practice, but we can use it here for demonstration purposes. To run the workflows that you added to the LaunchPad, run the following command. Note the use of ``rlaunch`` rather than ``qlaunch``.

.. code-block:: bash

    rlaunch rapidfire

You should see logging text on the progress of each Firetask in your workflow. The Fireworks have successfully finished launching and running, the results should be added to your database and you can move on.

Analyzing an Equation of State Workflow
=======================================

Querying the results
--------------------

In the Python preset, we get the nice EOS analysis Firework for free. This is not supported in atwf, so we will extract the data for a simple energy vs. volume curve ourselves. With the ``PMGDB_DB_FILE`` varible set in your ``$HOME/.pmgrc.yaml`` file as in the installation instructions, we will be querying the database that the db.json file you created describes.

.. code-block:: bash

    mgdb query --props task_id formula_pretty output.energy_per_atom output.structure.lattice.volume task_label


which will give you an overview of the each Firework you ran. It should look something like

.. code-block:: bash

      task_id  formula_pretty      output.energy_per_atom    output.structure.lattice.volume  task_label
    ---------  ----------------  ------------------------  ---------------------------------  -----------------------------------------------------
            1  Al                                -3.74617                            65.8868  2015-12-30-18-00-00-163825 structure optimization
            2  Al                                -3.69701                            59.2981  2015-12-30-18-00-00-163825 bulk_modulus deformation 0
            3  Al                                -3.73492                            62.5925  2015-12-30-18-00-00-163825 bulk_modulus deformation 1
            4  Al                                -3.74617                            65.8868  2015-12-30-18-00-00-163825 bulk_modulus deformation 2
            5  Al                                -3.73752                            69.1812  2015-12-30-18-00-00-163825 bulk_modulus deformation 3
            6  Al                                -3.71384                            72.4755  2015-12-30-18-00-00-163825 bulk_modulus deformation 4


Now we want to get the results for just our deformations. We add the ``--crit`` option to enable searching based on JSON-formatted criteria. Specifically we just want the deformation. By using the ``--dump`` option and redirection the results to a JSON file, we can load the results in Python for our analysis. We can also simplify the properties are getting, since we are already aware of the other things.

.. code-block:: bash

    mgdb query --crit '{"task_label": {"$regex": "deformation"}}' --props output.energy_per_atom output.structure.lattice.volume --dump > eos-results.json


.. note:: It is important to format your criteria as single quotes on the outside and double quotes on the inside. Double quotes are required for JSON and the single quotes prevent any shell magic that curly braces ('{') usually invoke.


If everything worked, you should have gotten no output, but you should be able to find an ``eos-results.json`` with the following content

.. code-block:: json

    {"output.structure.lattice.volume": 59.29814953343786, "output.energy_per_atom": -3.69701308}
    {"output.structure.lattice.volume": 62.59247826860741, "output.energy_per_atom": -3.7349166475}
    {"output.structure.lattice.volume": 65.88683660010845, "output.energy_per_atom": -3.74616541}
    {"output.structure.lattice.volume": 69.1811925337604, "output.energy_per_atom": -3.73751932}
    {"output.structure.lattice.volume": 72.47551533513209, "output.energy_per_atom": -3.7138362225}


We need to format this file to actual JSON to more easily load the results in Python. Add a "data" property name, extra curly braces and brackets around all of the data and a comma on each line, making the results a list. The file ``eos-results.json`` should look like

.. code-block:: json

    {"data":
    [
    {"output.structure.lattice.volume": 59.29814953343786, "output.energy_per_atom": -3.69701308},
    {"output.structure.lattice.volume": 62.59247826860741, "output.energy_per_atom": -3.7349166475},
    {"output.structure.lattice.volume": 65.88683660010845, "output.energy_per_atom": -3.74616541},
    {"output.structure.lattice.volume": 69.1811925337604, "output.energy_per_atom": -3.73751932},
    {"output.structure.lattice.volume": 72.47551533513209, "output.energy_per_atom": -3.7138362225}
    ]
    }


Analyzing the results
---------------------

Finally, we'll plot the EOS results that we saved in the last section. Simply add the following Python script (``eos-analysis.py``) to your folder and run it

.. code-block:: python

    # eos-analysis.py
    import json
    import matplotlib
    matplotlib.use('Agg') # a little magic for matplotlib to work without a $DISPLAY set
    from matplotlib import pyplot as plt

    # load the results as JSON
    with open('eos-results.json') as f:
        eos_results = json.load(f)

    # get the results into lists of volumes and energies
    volumes = []
    energies = []
    for entry in eos_results['data']:
        volumes.append(entry['output.structure.lattice.volume'])
        energies.append(entry['output.energy_per_atom'])

    # set up the plot, plot the results, and save them to a file
    fig = plt.figure()
    ax = fig.gca()
    ax.plot(volumes, energies, marker='o', linestyle='')
    ax.set_title('Energy vs. Volume for Al')
    ax.set_ylabel('Energy per atom (eV)')
    ax.set_xlabel('Volume (A^3)')
    fig.savefig('eos-energy-volume.png')


If you open the saved figure, ``eos-energy-volume.png``, on your computer you should see the datapoints for your first automated E-V curve plotted!

.. figure:: _static/eos_energy_volume.png
    :alt: Alumninum energy vs. volume

    Energy vs. volume curve for Al created from the EOS volume deformations.

Conclusion
==========

In this tutorial you learned how run a workflow from in a YAML file without writing any code and to do the same in Python. The keys to constructing your own workflows are

We have tried to provide common functionality as preset workflows in Python. Due to some current limitation in the atwf utility, some analysis tasks like the EOS Firework cannot currently be expressed in the YAML, so complete access to full preset workflows can only be achieved in Python.

To see what preset workflows can be run, see the documentation that includes them at :py:mod:`atomate.vasp.workflows.presets`. They can be set up the same way as in this tutorial.

Eventually you may want to create your own workflows that you can use and distribute. The :ref:`creating workflows` article is a guide for writing custom workflows in Python.

