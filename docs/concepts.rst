================
Atomate concepts
================

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